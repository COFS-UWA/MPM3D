#include "Simulations_pcp.h"

#include <cmath>

#include "MaterialModel.h"
#include "Step_T2D_CHM_p_Geo.h"

Step_T2D_CHM_p_Geo::Step_T2D_CHM_p_Geo(const char* _name) :
	Step(_name, "Step_T2D_CHM_p_Geo", &solve_substep_T2D_CHM_p_Geo),
	model(nullptr), damping_ratio(0.0),
	f_ub_ratio_bound(0.0), e_kin_ratio_bound(0.0),
	thread_num(0), not_yet_completed(false)
{

}

Step_T2D_CHM_p_Geo::~Step_T2D_CHM_p_Geo()
{
	join_all_threads();
}

int Step_T2D_CHM_p_Geo::init_calculation()
{
	Model_T2D_CHM_p &md = *model;

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
		Particle &pcl = md.pcls[pcl_id];

		pcl.pe = md.find_in_which_element(pcl);
		if (!pcl.pe)
			continue;
		
		pcl.vol_s = pcl.m_s / pcl.density_s;
		pcl.vol = pcl.vol_s / (1.0 - pcl.n);
		pcl.m_f = pcl.vol * pcl.n * pcl.density_f;

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

int Step_T2D_CHM_p_Geo::finalize_calculation()
{
	join_all_threads();

	Model_T2D_CHM_p &md = *model;
	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle &pcl = md.pcls[pcl_id];
		pcl.vx_s = 0.0;
		pcl.vy_s = 0.0;
		pcl.vx_f = 0.0;
		pcl.vy_f = 0.0;
		pcl.p = 0.0;
	}

	return 0;
}

void Step_T2D_CHM_p_Geo::spawn_all_threads()
{
	Model_T2D_CHM_p& md = *model;
	if (thread_num > 1)
	{
		size_t spawn_thread_num = thread_num - 1;
		cal_threads.resize(spawn_thread_num);
		for (unsigned int th_id = 0; th_id < spawn_thread_num; ++th_id)
			cal_threads[th_id] = std::thread(&Step_T2D_CHM_p_Geo::cal_thread_func, this, th_id + 1);
	}
}

void Step_T2D_CHM_p_Geo::join_all_threads()
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

void Step_T2D_CHM_p_Geo::cal_thread_func(unsigned int th_id)
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

		update_pcl_vars(th_id);
		step_barrier.wait();
	}
}

int solve_substep_T2D_CHM_p_Geo(void *_self)
{
	typedef Model_T2D_CHM_p::Particle Particle;
	typedef Model_T2D_CHM_p::Element Element;
	typedef Model_T2D_CHM_p::Node Node;

	Step_T2D_CHM_p_Geo &self = *(Step_T2D_CHM_p_Geo *)(_self);
	Model_T2D_CHM_p &md = *self.model;

	self.cur_elem_id.store(0, std::memory_order_relaxed);
	self.step_barrier.lift_barrier();
	self.map_pcl_vars_to_nodes_at_elems(0);
	self.cal_barrier.wait_for_others();

	self.cur_node_id.store(0, std::memory_order_relaxed);
	self.cal_barrier.lift_barrier();
	self.update_node_a_and_v(0);
	self.cal_barrier.wait_for_others();

	self.cur_asx_bc_id.store(0, std::memory_order_relaxed);
	self.cur_asy_bc_id.store(0, std::memory_order_relaxed);
	self.cur_vsx_bc_id.store(0, std::memory_order_relaxed);
	self.cur_vsy_bc_id.store(0, std::memory_order_relaxed);
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
			self.f_ub += n.m_s * n.m_s * (n.ax_s * n.ax_s + n.ay_s * n.ay_s);
			self.e_kin += n.m_s * (n.vx_s * n.vx_s + n.vy_s * n.vy_s);
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
			pcl.vx_s = 0.0;
			pcl.vy_s = 0.0;
			pcl.vx_f = 0.0;
			pcl.vy_f = 0.0;
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

	self.cur_pcl_id.store(0, std::memory_order_relaxed);
	self.cal_barrier.lift_barrier();
	self.update_pcl_vars(0);
	self.step_barrier.wait_for_others();

	return 0;
}

void Step_T2D_CHM_p_Geo::map_pcl_vars_to_nodes_at_elems(unsigned int th_id)
{
	Model_T2D_CHM_p& md = *static_cast<Model_T2D_CHM_p*>(model);

	double N_m_s;
	size_t e_id;
	while ((e_id = ATOM_INC(cur_elem_id)) < md.elem_num)
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
			e.mi_pcl_vol += pcl.vol;
			e.s11 += pcl.s11 * pcl.vol;
			e.s22 += pcl.s22 * pcl.vol;
			e.s12 += pcl.s12 * pcl.vol;

			// map velocity
			// node 1
			N_m_s = pcl.N1 * pcl.m_s;
			nv1.m_s += N_m_s;
			nv1.vx_s += N_m_s * pcl.vx_s;
			nv1.vy_s += N_m_s * pcl.vy_s;
			// node 2
			N_m_s = pcl.N2 * pcl.m_s;
			nv2.m_s += N_m_s;
			nv2.vx_s += N_m_s * pcl.vx_s;
			nv2.vy_s += N_m_s * pcl.vy_s;
			// node 3
			N_m_s = pcl.N3 * pcl.m_s;
			nv3.m_s += N_m_s;
			nv3.vx_s += N_m_s * pcl.vx_s;
			nv3.vy_s += N_m_s * pcl.vy_s;

			// external load
			if (pcl.has_fx_s_ext)
			{
				nv1.fx_s += pcl.N1 * pcl.fx_s_ext;
				nv2.fx_s += pcl.N2 * pcl.fx_s_ext;
				nv3.fx_s += pcl.N3 * pcl.fx_s_ext;
			}
			if (pcl.has_fy_s_ext)
			{
				nv1.fy_s += pcl.N1 * pcl.fy_s_ext;
				nv2.fy_s += pcl.N2 * pcl.fy_s_ext;
				nv3.fy_s += pcl.N3 * pcl.fy_s_ext;
			}
		}

		e.s11 /= e.mi_pcl_vol;
		e.s22 /= e.mi_pcl_vol;
		e.s12 /= e.mi_pcl_vol;
		if (e.mi_pcl_vol > e.area)
			e.mi_pcl_vol = e.area;

		// internal force
		// node 1
		nv1.fx_s -= (e.dN1_dx * e.s11 + e.dN1_dy * e.s12) * e.mi_pcl_vol;
		nv1.fy_s -= (e.dN1_dx * e.s12 + e.dN1_dy * e.s22) * e.mi_pcl_vol;
		// node 2
		nv2.fx_s -= (e.dN2_dx * e.s11 + e.dN2_dy * e.s12) * e.mi_pcl_vol;
		nv2.fy_s -= (e.dN2_dx * e.s12 + e.dN2_dy * e.s22) * e.mi_pcl_vol;
		// node 3
		nv3.fx_s -= (e.dN3_dx * e.s11 + e.dN3_dy * e.s12) * e.mi_pcl_vol;
		nv3.fy_s -= (e.dN3_dx * e.s12 + e.dN3_dy * e.s22) * e.mi_pcl_vol;
	}
}

void Step_T2D_CHM_p_Geo::update_node_a_and_v(unsigned int th_id)
{
	Model_T2D_CHM_p& md = *static_cast<Model_T2D_CHM_p*>(model);

	size_t n_id;
	while ((n_id = ATOM_INC(cur_node_id)) < md.node_num)
	{
		Node& n = md.nodes[n_id];

		if (!n.has_mp)
			continue;

		n.m_s = 0.0;
		n.vx_s_ori = 0.0;
		n.vy_s_ori = 0.0;
		n.fx_s = 0.0;
		n.fy_s = 0.0;

		for (size_t e_id = 0; e_id < n.n2e_num; ++e_id)
		{
			NodeToElem& n2e = n.n2es[e_id];
			Element& e = md.elems[n2e.e_id];
			if (e.has_mp)
			{
				NodeVarAtElem& nv = e.node_vars[n2e.n_id];
				n.m_s += nv.m_s;
				n.vx_s_ori += nv.vx_s;
				n.vy_s_ori += nv.vy_s;
				n.fx_s += nv.fx_s;
				n.fy_s += nv.fy_s;
			}
		}

		n.ax_s = n.fx_s / n.m_s;
		n.ay_s = n.fy_s / n.m_s;
		n.vx_s_ori /= n.m_s;
		n.vy_s_ori /= n.m_s;
		n.vx_s = n.vx_s_ori + n.ax_s * dtime;
		n.vy_s = n.vy_s_ori + n.ay_s * dtime;
		n.dux_s = n.vx_s * dtime;
		n.duy_s = n.vy_s * dtime;
	}
}

void Step_T2D_CHM_p_Geo::apply_a_and_v_bcs(unsigned int th_id)
{
	Model_T2D_CHM_p& md = *static_cast<Model_T2D_CHM_p*>(model);

	size_t a_id;
	while ((a_id = ATOM_INC(cur_asx_bc_id)) < md.asx_num)
	{
		Node& n = md.nodes[md.asxs[a_id].node_id];
		n.ax_s = md.asxs[a_id].a;
		n.vx_s = n.vx_s_ori + n.ax_s * dtime;
		n.dux_s = n.vx_s * dtime;
	}
	while ((a_id = ATOM_INC(cur_asy_bc_id)) < md.asy_num)
	{
		Node& n = md.nodes[md.asys[a_id].node_id];
		n.ay_s = md.asys[a_id].a;
		n.vy_s = n.vy_s_ori + n.ay_s * dtime;
		n.duy_s = n.vy_s * dtime;
	}

	size_t v_id;
	while ((v_id = ATOM_INC(cur_vsx_bc_id)) < md.vsx_num)
	{
		Node& n = md.nodes[md.vsxs[v_id].node_id];
		n.ax_s = 0.0;
		n.vx_s = md.vsxs[v_id].v;
		n.dux_s = n.vx_s * dtime;
	}
	while ((v_id = ATOM_INC(cur_vsy_bc_id)) < md.vsy_num)
	{
		Node& n = md.nodes[md.vsys[v_id].node_id];
		n.ay_s = 0.0;
		n.vy_s = md.vsys[v_id].v;
		n.duy_s = n.vy_s * dtime;
	}
}

void Step_T2D_CHM_p_Geo::update_pcl_vars(unsigned int th_id)
{
	Model_T2D_CHM_p& md = *static_cast<Model_T2D_CHM_p*>(model);

	double de_vol_by_3, de11, de22, de12;
	size_t pcl_id;
	while ((pcl_id = ATOM_INC(cur_pcl_id)) < md.pcl_num)
	{
		Particle& pcl = md.pcls[pcl_id];
		if (!pcl.pe)
			return;

		Element& e = *pcl.pe;
		Node& n1 = md.nodes[e.n1];
		Node& n2 = md.nodes[e.n2];
		Node& n3 = md.nodes[e.n3];

		// velocity
		pcl.vx_s += (n1.ax_s * pcl.N1 + n2.ax_s * pcl.N2 + n3.ax_s * pcl.N3) * dtime;
		pcl.vy_s += (n1.ay_s * pcl.N1 + n2.ay_s * pcl.N2 + n3.ay_s * pcl.N3) * dtime;

		// strain
		de11 = n1.dux_s * e.dN1_dx + n2.dux_s * e.dN2_dx + n3.dux_s * e.dN3_dx;
		de22 = n1.duy_s * e.dN1_dy + n2.duy_s * e.dN2_dy + n3.duy_s * e.dN3_dy;
		de12 = (n1.dux_s * e.dN1_dy + n2.dux_s * e.dN2_dy + n3.dux_s * e.dN3_dy
			  + n1.duy_s * e.dN1_dx + n2.duy_s * e.dN2_dx + n3.duy_s * e.dN3_dx) * 0.5;

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

