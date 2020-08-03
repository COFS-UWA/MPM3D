#include "Simulations_pcp.h"

#include <cmath>

#include "MaterialModel.h"
#include "Step_T2D_ME_p_Geo.h"

Step_T2D_ME_p_Geo::Step_T2D_ME_p_Geo(const char* _name) :
	Step(_name, "Step_T2D_ME_p_Geo", &solve_substep_T2D_ME_p_geo),
	model(nullptr), damping_ratio(0.0),
	f_ub_ratio_bound(0.0), e_kin_ratio_bound(0.0),
	thread_num(0), not_yet_completed(false)
{

}

Step_T2D_ME_p_Geo::~Step_T2D_ME_p_Geo()
{
	join_all_threads();
}

int Step_T2D_ME_p_Geo::init_calculation()
{
	Model_T2D_ME_p &md = *model;

	if (is_first_step)
	{
		md.init_bcs();
	}

	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element& e = md.elems[e_id];
		e.has_mp = false;
		e.pcls = nullptr;
	}

	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle& pcl = md.pcls[pcl_id];
		pcl.pe = md.find_in_which_element(pcl);
		if (!pcl.pe)
			continue;

		pcl.vol = pcl.m / pcl.density;
		
		Element& e = *pcl.pe;
		e.has_mp = true;
		e.add_pcl(pcl);

		Node& n1 = md.nodes[e.n1];
		n1.has_mp = true;
		Node& n2 = md.nodes[e.n2];
		n2.has_mp = true;
		Node& n3 = md.nodes[e.n3];
		n3.has_mp = true;
	}

	// convergence criteria
	// unbalanced force
	f_ub = 0.0;
	init_f_ub = 0.0;
	f_ub_ratio = 1.0;
	// kinematic energy
	e_kin = 0.0;
	e_kin_max = 0.0;
	e_kin_prev = 0.0;
	e_kin_ratio = 1.0;

	// start calculation
	if (thread_num == 0)
		thread_num = std::thread::hardware_concurrency();

	// init global data
	not_yet_completed.store(true, std::memory_order_relaxed);
	step_barrier.set_thread_num(thread_num);
	cal_barrier.set_thread_num(thread_num);

	// spawn threads
	spawn_all_threads();
	step_barrier.wait_for_others();

	return 0;
}

int Step_T2D_ME_p_Geo::finalize_calculation()
{
	join_all_threads();

	Model_T2D_ME_p& md = *model;
	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle &pcl = md.pcls[pcl_id];
		pcl.vx = 0.0;
		pcl.vy = 0.0;
		pcl.e11 = 0.0;
		pcl.e22 = 0.0;
		pcl.e12 = 0.0;
	}

	return 0;
}

void Step_T2D_ME_p_Geo::spawn_all_threads()
{
	Model_T2D_ME_p& md = *model;
	if (thread_num > 1)
	{
		size_t spawn_thread_num = thread_num - 1;
		cal_threads.resize(spawn_thread_num);
		for (unsigned int th_id = 0; th_id < spawn_thread_num; ++th_id)
			cal_threads[th_id] = std::thread(&Step_T2D_ME_p_Geo::cal_thread_func, this, th_id + 1);
	}
}

void Step_T2D_ME_p_Geo::join_all_threads()
{
	if (not_yet_completed.load(std::memory_order_relaxed))
	{
		not_yet_completed.store(false, std::memory_order_relaxed);
		step_barrier.lift_barrier();
		size_t join_thread_num = thread_num - 1;
		for (unsigned int th_id = 0; th_id < join_thread_num; ++th_id)
			cal_threads[th_id].join();
	}
}

void Step_T2D_ME_p_Geo::cal_thread_func(unsigned int th_id)
{
	step_barrier.wait();

	while (not_yet_completed.load(std::memory_order_relaxed))
	{
		map_pcl_vars_to_nodes_at_elems(th_id);
		cal_barrier.wait();

		update_node_a_and_v(th_id);
		cal_barrier.wait();

		apply_a_and_v_bcs(th_id);
		cal_barrier.wait();

		if (terminate_this_substep.load(std::memory_order_relaxed))
		{
			step_barrier.wait();
			continue;
		}

		cal_de_at_elem(th_id);
		cal_barrier.wait();

		map_de_vol_from_elem_to_node(th_id);
		cal_barrier.wait();

		update_pcl_vars(th_id);
		step_barrier.wait();
	}
}

int solve_substep_T2D_ME_p_geo(void *_self)
{
	typedef Model_T2D_ME_p::Node Node;
	typedef Model_T2D_ME_p::Element Element;
	typedef Model_T2D_ME_p::Particle Particle;

	Step_T2D_ME_p_Geo &self = *(Step_T2D_ME_p_Geo *)(_self);
	Model_T2D_ME_p &md = *self.model;

	self.cur_elem_id.store(0, std::memory_order_relaxed);
	self.step_barrier.lift_barrier();
	self.map_pcl_vars_to_nodes_at_elems(0);
	self.cal_barrier.wait_for_others();

	self.cur_node_id.store(0, std::memory_order_relaxed);
	self.cal_barrier.lift_barrier();
	self.update_node_a_and_v(0);
	self.cal_barrier.wait_for_others();

	self.cur_ax_bc_id.store(0, std::memory_order_relaxed);
	self.cur_ay_bc_id.store(0, std::memory_order_relaxed);
	self.cur_vx_bc_id.store(0, std::memory_order_relaxed);
	self.cur_vy_bc_id.store(0, std::memory_order_relaxed);
	self.cal_barrier.lift_barrier();
	self.apply_a_and_v_bcs(0);
	self.cal_barrier.wait_for_others();

	// unbalanced nodal force and kinetic energy
	self.f_ub = 0.0;
	self.e_kin = 0.0;
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node& n = md.nodes[n_id];
		if (n.has_mp)
		{
			self.f_ub += n.m * n.m * (n.ax * n.ax + n.ay * n.ay);
			self.e_kin += n.m * (n.vx * n.vx + n.vy * n.vy);
		}
	}

	if (self.substep_index == 0) // initial step
	{
		self.init_f_ub = self.f_ub;
		self.e_kin_max = self.e_kin;
	}
	else
	{
		if (self.e_kin_max < self.e_kin)
			self.e_kin_max = self.e_kin;
	}

	self.f_ub_ratio = sqrt(self.f_ub / self.init_f_ub);
	self.e_kin_ratio = sqrt(self.e_kin / self.e_kin_max);

	self.terminate_this_substep.store(false, std::memory_order_relaxed);
	if (self.e_kin < self.e_kin_prev)
	{
		// reset pcl velocity
		for (size_t p_id = 0; p_id < md.pcl_num; ++p_id)
		{
			Particle& pcl = md.pcls[p_id];
			pcl.vx = 0.0;
			pcl.vy = 0.0;
		}
		self.e_kin_prev = 0.0;

		//if (self.e_kin_ratio < self.e_kin_ratio_bound &&
		//	  self.f_ub_ratio < self.f_ub_ratio_bound)
		//		return
		self.terminate_this_substep.store(true, std::memory_order_relaxed);
		self.cal_barrier.lift_barrier();
		self.step_barrier.wait_for_others();
		return 0;
	}
	
	self.e_kin_prev = self.e_kin;

	self.cur_elem_id.store(0, std::memory_order_relaxed);
	self.cal_barrier.lift_barrier();
	self.cal_de_at_elem(0);
	self.cal_barrier.wait_for_others();

	self.cur_node_id.store(0, std::memory_order_relaxed);
	self.cal_barrier.lift_barrier();
	self.map_de_vol_from_elem_to_node(0);
	self.cal_barrier.wait_for_others();

	self.cur_pcl_id.store(0, std::memory_order_relaxed);
	self.cal_barrier.lift_barrier();
	self.update_pcl_vars(0);

	self.step_barrier.wait_for_others();
	return 0;
}

void Step_T2D_ME_p_Geo::map_pcl_vars_to_nodes_at_elems(unsigned int th_id)
{
	Model_T2D_ME_p& md = *static_cast<Model_T2D_ME_p*>(model);

	for (size_t e_id = cur_elem_id.fetch_add(1, std::memory_order_relaxed);
		 e_id < md.elem_num;
		 e_id = cur_elem_id.fetch_add(1, std::memory_order_relaxed))
	{
		Element& e = md.elems[e_id];

		if (!e.has_mp)
			continue;

		e.mi_pcl_vol = 0.0;
		e.s11 = 0.0;
		e.s22 = 0.0;
		e.s12 = 0.0;

		NodeVarAtElem& nv1 = e.node_vars[0];
		NodeVarAtElem& nv2 = e.node_vars[1];
		NodeVarAtElem& nv3 = e.node_vars[2];
		nv1.init();
		nv2.init();
		nv3.init();

		for (Particle* ppcl = e.first_pcl(); e.not_last_pcl(ppcl);
			 ppcl = e.next_pcl(ppcl))
		{
			Particle& pcl = *ppcl;

			// mixed integration
			pcl.vol = pcl.m / pcl.density;
			e.mi_pcl_vol += pcl.vol;
			e.s11 += pcl.s11 * pcl.vol;
			e.s22 += pcl.s22 * pcl.vol;
			e.s12 += pcl.s12 * pcl.vol;

			// map velocity
			double mvx = pcl.m * pcl.vx;
			double mvy = pcl.m * pcl.vy;
			// node 1
			nv1.m += pcl.N1 * pcl.m;
			nv1.vx += pcl.N1 * mvx;
			nv1.vy += pcl.N1 * mvy;
			// node 2
			nv2.m += pcl.N2 * pcl.m;
			nv2.vx += pcl.N2 * mvx;
			nv2.vy += pcl.N2 * mvy;
			// node 3
			nv3.m += pcl.N3 * pcl.m;
			nv3.vx += pcl.N3 * mvx;
			nv3.vy += pcl.N3 * mvy;

			// external load
			if (pcl.has_fx_ext)
			{
				nv1.fx += pcl.N1 * pcl.fx_ext;
				nv2.fx += pcl.N2 * pcl.fx_ext;
				nv3.fx += pcl.N3 * pcl.fx_ext;
			}
			if (pcl.has_fy_ext)
			{
				nv1.fy += pcl.N1 * pcl.fy_ext;
				nv2.fy += pcl.N2 * pcl.fy_ext;
				nv3.fy += pcl.N3 * pcl.fy_ext;
			}
		}

		e.s11 /= e.mi_pcl_vol;
		e.s22 /= e.mi_pcl_vol;
		e.s12 /= e.mi_pcl_vol;
		if (e.mi_pcl_vol > e.area)
			e.mi_pcl_vol = e.area;

		// internal force
		// node 1
		nv1.fx -= (e.dN1_dx * e.s11 + e.dN1_dy * e.s12) * e.mi_pcl_vol;
		nv1.fy -= (e.dN1_dx * e.s12 + e.dN1_dy * e.s22) * e.mi_pcl_vol;
		// node 2
		nv2.fx -= (e.dN2_dx * e.s11 + e.dN2_dy * e.s12) * e.mi_pcl_vol;
		nv2.fy -= (e.dN2_dx * e.s12 + e.dN2_dy * e.s22) * e.mi_pcl_vol;
		// node 3
		nv3.fx -= (e.dN3_dx * e.s11 + e.dN3_dy * e.s12) * e.mi_pcl_vol;
		nv3.fy -= (e.dN3_dx * e.s12 + e.dN3_dy * e.s22) * e.mi_pcl_vol;
	}
}

void Step_T2D_ME_p_Geo::update_node_a_and_v(unsigned int th_id)
{
	Model_T2D_ME_p& md = *static_cast<Model_T2D_ME_p*>(model);

	for (size_t n_id = cur_node_id.fetch_add(1, std::memory_order_relaxed);
		 n_id < md.node_num;
		 n_id = cur_node_id.fetch_add(1, std::memory_order_relaxed))
	{
		Node& n = md.nodes[n_id];

		if (!n.has_mp)
			continue;

		n.m = 0.0;
		n.vx_ori = 0.0;
		n.vy_ori = 0.0;
		n.fx = 0.0;
		n.fy = 0.0;

		for (size_t e_id = 0; e_id < n.n2e_num; ++e_id)
		{
			NodeToElem& n2e = n.n2es[e_id];
			Element& e = md.elems[n2e.e_id];
			if (e.has_mp)
			{
				NodeVarAtElem& nv = e.node_vars[n2e.n_id];
				n.m += nv.m;
				n.vx_ori += nv.vx;
				n.vy_ori += nv.vy;
				n.fx += nv.fx;
				n.fy += nv.fy;
			}
		}

		n.ax = n.fx / n.m;
		n.ay = n.fy / n.m;
		n.vx_ori /= n.m;
		n.vy_ori /= n.m;
		n.vx = n.vx_ori + n.ax * dtime;
		n.vy = n.vy_ori + n.ay * dtime;
		n.dux = n.vx * dtime;
		n.duy = n.vy * dtime;
	}
}

void Step_T2D_ME_p_Geo::apply_a_and_v_bcs(unsigned int th_id)
{
	Model_T2D_ME_p& md = *static_cast<Model_T2D_ME_p*>(model);

	for (size_t a_id = cur_ax_bc_id.fetch_add(1, std::memory_order_relaxed);
		 a_id < md.ax_num;
		 a_id = cur_ax_bc_id.fetch_add(1, std::memory_order_relaxed))
	{
		Node& n = md.nodes[md.axs[a_id].node_id];
		n.ax = md.axs[a_id].a;
		n.vx = n.vx_ori + n.ax * dtime;
		n.dux = n.vx * dtime;
	}

	for (size_t a_id = cur_ay_bc_id.fetch_add(1, std::memory_order_relaxed);
		 a_id < md.ay_num;
		 a_id = cur_ay_bc_id.fetch_add(1, std::memory_order_relaxed))
	{
		Node& n = md.nodes[md.ays[a_id].node_id];
		n.ay = md.ays[a_id].a;
		n.vy = n.vy_ori + n.ay * dtime;
		n.duy = n.vy * dtime;
	}

	for (size_t v_id = cur_vx_bc_id.fetch_add(1, std::memory_order_relaxed);
		 v_id < md.vx_num;
		 v_id = cur_vx_bc_id.fetch_add(1, std::memory_order_relaxed))
	{
		Node& n = md.nodes[md.vxs[v_id].node_id];
		n.ax = 0.0;
		n.vx = md.vxs[v_id].v;
		n.dux = n.vx * dtime;
	}

	for (size_t v_id = cur_vy_bc_id.fetch_add(1, std::memory_order_relaxed);
		 v_id < md.vy_num;
		 v_id = cur_vy_bc_id.fetch_add(1, std::memory_order_relaxed))
	{
		Node& n = md.nodes[md.vys[v_id].node_id];
		n.ay = 0.0;
		n.vy = md.vys[v_id].v;
		n.duy = n.vy * dtime;
	}
}

void Step_T2D_ME_p_Geo::cal_de_at_elem(unsigned int th_id)
{
	Model_T2D_ME_p& md = *static_cast<Model_T2D_ME_p*>(model);

	double de11, de22;
	for (size_t e_id = cur_elem_id.fetch_add(1, std::memory_order_relaxed);
		 e_id < md.elem_num;
		 e_id = cur_elem_id.fetch_add(1, std::memory_order_relaxed))
	{
		Element& e = md.elems[e_id];

		if (!e.has_mp)
			continue;

		Node& n1 = md.nodes[e.n1];
		Node& n2 = md.nodes[e.n2];
		Node& n3 = md.nodes[e.n3];
		de11 = n1.dux * e.dN1_dx + n2.dux * e.dN2_dx + n3.dux * e.dN3_dx;
		de22 = n1.duy * e.dN1_dy + n2.duy * e.dN2_dy + n3.duy * e.dN3_dy;
		e.de_vol_by_3 = (de11 + de22) / 3.0;
		e.dde11 = de11 - e.de_vol_by_3;
		e.dde22 = de22 - e.de_vol_by_3;
		e.de12 = (n1.dux * e.dN1_dy + n2.dux * e.dN2_dy + n3.dux * e.dN3_dy
				+ n1.duy * e.dN1_dx + n2.duy * e.dN2_dx + n3.duy * e.dN3_dx) * 0.5;

		// strain enhancement
		NodeVarAtElem& nv1 = e.node_vars[0];
		NodeVarAtElem& nv2 = e.node_vars[1];
		NodeVarAtElem& nv3 = e.node_vars[2];
		double N_vol;
		for (Particle* ppcl = e.first_pcl(); e.not_last_pcl(ppcl);
			ppcl = e.next_pcl(ppcl))
		{
			Particle& pcl = *ppcl;
			// node 1
			N_vol = pcl.N1 * pcl.vol;
			nv1.de_vol_by_3 += e.de_vol_by_3 * N_vol;
			nv1.se_pcl_vol += N_vol;
			// node 2
			N_vol = pcl.N2 * pcl.vol;
			nv2.de_vol_by_3 += e.de_vol_by_3 * N_vol;
			nv2.se_pcl_vol += N_vol;
			// node 3
			N_vol = pcl.N3 * pcl.vol;
			nv3.de_vol_by_3 += e.de_vol_by_3 * N_vol;
			nv3.se_pcl_vol += N_vol;
		}
	}
}

void Step_T2D_ME_p_Geo::map_de_vol_from_elem_to_node(unsigned int th_id)
{
	Model_T2D_ME_p& md = *static_cast<Model_T2D_ME_p*>(model);

	for (size_t n_id = cur_node_id.fetch_add(1, std::memory_order_relaxed);
		 n_id < md.node_num;
		 n_id = cur_node_id.fetch_add(1, std::memory_order_relaxed))
	{
		Node& n = md.nodes[n_id];

		if (!n.has_mp)
			continue;

		n.de_vol_by_3 = 0.0;
		n.se_pcl_vol = 0.0;
		for (size_t e_id = 0; e_id < n.n2e_num; ++e_id)
		{
			NodeToElem& n2e = n.n2es[e_id];
			Element& e = md.elems[n2e.e_id];
			if (e.has_mp)
			{
				NodeVarAtElem& nv = e.node_vars[n2e.n_id];
				n.de_vol_by_3 += nv.de_vol_by_3;
				n.se_pcl_vol += nv.se_pcl_vol;
			}
		}
		n.de_vol_by_3 /= n.se_pcl_vol;
	}
}

void Step_T2D_ME_p_Geo::update_pcl_vars(unsigned int th_id)
{
	Model_T2D_ME_p& md = *static_cast<Model_T2D_ME_p*>(model);

	double de_vol_by_3, de11, de22, de12;
	size_t pcl_id;
	while ((pcl_id = cur_pcl_id.fetch_add(1, std::memory_order_relaxed)) < md.pcl_num)
	{
		Particle& pcl = md.pcls[pcl_id];
		if (!pcl.pe)
			return;

		Element& e = *pcl.pe;
		Node& n1 = md.nodes[e.n1];
		Node& n2 = md.nodes[e.n2];
		Node& n3 = md.nodes[e.n3];

		// velocity
		pcl.vx += (n1.ax * pcl.N1 + n2.ax * pcl.N2 + n3.ax * pcl.N3) * dtime;
		pcl.vy += (n1.ay * pcl.N1 + n2.ay * pcl.N2 + n3.ay * pcl.N3) * dtime;

		// strain
		de_vol_by_3 = (n1.de_vol_by_3 + n2.de_vol_by_3 + n3.de_vol_by_3) / 3.0;
		de11 = e.dde11 + de_vol_by_3;
		de22 = e.dde22 + de_vol_by_3;
		de12 = e.de12;
		pcl.e11 += de11;
		pcl.e22 += de22;
		pcl.e12 += de12;

		// stress
		// update stress using constitutive model
		double dstrain[6] = { de11, de22, 0.0, de12, 0.0, 0.0 };
		pcl.mm->integrate(dstrain);
		const double* dstress = pcl.mm->get_dstress();
		pcl.s11 += dstress[0];
		pcl.s22 += dstress[1];
		pcl.s12 += dstress[3];
	}
}
