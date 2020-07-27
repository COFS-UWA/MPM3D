#include "Simulations_pcp.h"

#include <cmath>
#include <iostream>

#include "MaterialModel.h"
#include "Step_T2D_ME_p.h"

Step_T2D_ME_p::Step_T2D_ME_p(const char* _name) :
	Step(_name, "Step_T2D_ME_p", &solve_substep_T2D_ME_p),
	model(nullptr), damping_ratio(0.0),
	thread_num(0), not_yet_completed(false),
	ted_mem(nullptr), pteds(nullptr)
{

}

Step_T2D_ME_p::~Step_T2D_ME_p()
{
	join_thread_and_exit();
	clear_thread_elem_data();
}

int Step_T2D_ME_p::init_calculation()
{
	Model_T2D_ME_p &md = *model;

	if (is_first_step)
	{
		md.init_bcs();
	}

	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle& pcl = md.pcls[pcl_id];
		pcl.pe = md.elems;
		pcl.vol = pcl.m / pcl.density;
		pcl.x_ori = pcl.x;
		pcl.y_ori = pcl.y;
		pcl.ux = 0.0;
		pcl.uy = 0.0;
	}

	// start cal thread
	not_yet_completed = true;
	if (thread_num == 0)
		thread_num = std::thread::hardware_concurrency();
	if (thread_num > 1)
	{
		cal_threads.resize(thread_num-1);
		for (unsigned int th_id = 0; th_id < thread_num-1; ++th_id)
			cal_threads[th_id] = std::thread(&Step_T2D_ME_p::cal_thread_func, this, th_id + 1);
	}
	step_barrier.set_thread_num(thread_num);
	cal_barrier.set_thread_num(thread_num);
	
	alloc_thread_elem_data(md.elem_num);
	step_barrier.wait_for_others();

	return 0;
}

int Step_T2D_ME_p::finalize_calculation()
{
	join_thread_and_exit();
	clear_thread_elem_data();
	return 0;
}

void Step_T2D_ME_p::join_thread_and_exit()
{
	if (not_yet_completed)
	{
		not_yet_completed = false;
		step_barrier.lift_barrier();
		for (unsigned int th_id = 0; th_id < thread_num-1; ++th_id)
			cal_threads[th_id].join();
	}
}

void Step_T2D_ME_p::alloc_thread_elem_data(size_t elem_num)
{
	clear_thread_elem_data();
	ted_mem = new ThreadElemData[thread_num * elem_num];
	pteds = new ThreadElemData*[thread_num];
	for (size_t th_id = 0; th_id < thread_num; ++th_id)
		pteds[th_id] = ted_mem + th_id * elem_num;
}

void Step_T2D_ME_p::clear_thread_elem_data()
{
	if (ted_mem)
	{
		delete[] ted_mem;
		ted_mem = nullptr;
	}
	if (pteds)
	{
		delete[] pteds;
		pteds = nullptr;
	}
}

void Step_T2D_ME_p::cal_thread_func(unsigned int th_id)
{
	step_barrier.wait();

	while (not_yet_completed)
	{
		find_pcls_in_which_elems(th_id);
		cal_barrier.wait();

		map_pcl_vars_to_nodes_at_elems(th_id);
		cal_barrier.wait();

		update_node_a_and_v(th_id);
		cal_barrier.wait();

		apply_a_and_v_bcs(th_id);
		cal_barrier.wait();

		cal_de_at_elem(th_id);
		cal_barrier.wait();
		
		map_de_vol_from_elem_to_node(th_id);
		cal_barrier.wait();

		update_pcl_vars(th_id);
		step_barrier.wait();
	}
}

int solve_substep_T2D_ME_p(void* _self)
{
	typedef Model_T2D_ME_p::Node Node;
	typedef Model_T2D_ME_p::Element Element;
	typedef Model_T2D_ME_p::Particle Particle;

	Step_T2D_ME_p &self = *static_cast<Step_T2D_ME_p *>(_self);
	Model_T2D_ME_p& md = static_cast<Model_T2D_ME_p &>(self.get_model());

	// init calculation globally
	// init nodes
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node& n = md.nodes[n_id];
		n.has_mp = false;
	}

	// init elements
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element& e = md.elems[e_id];
		e.has_mp = false;
	}

	self.cur_pcl_id.store(0);

	self.step_barrier.lift_barrier();
	self.find_pcls_in_which_elems(0);
	self.cal_barrier.wait_for_others();

	self.cur_elem_id.store(0);
	self.cal_barrier.lift_barrier();
	self.map_pcl_vars_to_nodes_at_elems(0);
	self.cal_barrier.wait_for_others();

	self.cur_node_id.store(0);
	self.cal_barrier.lift_barrier();
	self.update_node_a_and_v(0);
	self.cal_barrier.wait_for_others();

	self.cur_ax_bc_id.store(0);
	self.cur_ay_bc_id.store(0);
	self.cur_vx_bc_id.store(0);
	self.cur_vy_bc_id.store(0);
	self.cal_barrier.lift_barrier();
	self.apply_a_and_v_bcs(0);
	self.cal_barrier.wait_for_others();

	self.cur_elem_id.store(0);
	self.cal_barrier.lift_barrier();
	self.cal_de_at_elem(0);
	self.cal_barrier.wait_for_others();

	self.cur_node_id.store(0);
	self.cal_barrier.lift_barrier();
	self.map_de_vol_from_elem_to_node(0);
	self.cal_barrier.wait_for_others();

	self.cur_pcl_id.store(0);
	self.cal_barrier.lift_barrier();
	self.update_pcl_vars(0);

	self.step_barrier.wait_for_others();
	
	return 0;
}

void Step_T2D_ME_p::find_pcls_in_which_elems(unsigned int th_id)
{
	Model_T2D_ME_p& md = *static_cast<Model_T2D_ME_p*>(model);

	ThreadElemData* teds = pteds[th_id];
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
		teds[e_id].init();
	
	size_t pe_id;
	for (size_t pcl_id = cur_pcl_id++; pcl_id < md.pcl_num; pcl_id = cur_pcl_id++)
	{
		Particle &pcl = md.pcls[pcl_id];
		
		if (pcl.pe == nullptr)
			continue;

		pcl.pe = md.find_in_which_element(pcl);
		if (pcl.pe == nullptr)
			continue;
		
		pe_id = pcl.pe->id;
		teds[pe_id].add_pcl(pcl);
		md.elems[pe_id].has_mp = true;
	}
}

void Step_T2D_ME_p::map_pcl_vars_to_nodes_at_elems(unsigned int th_id)
{
	Model_T2D_ME_p& md = *static_cast<Model_T2D_ME_p*>(model);
	
	for (size_t e_id = cur_elem_id++; e_id < md.elem_num; e_id = cur_elem_id++)
	{
		Element &e = md.elems[e_id];

		if (!e.has_mp)
			continue;
		
		Node &n1 = md.nodes[e.n1];
		Node &n2 = md.nodes[e.n2];
		Node &n3 = md.nodes[e.n3];
		NodeVarAtElem& nv1 = e.node_vars[0];
		NodeVarAtElem& nv2 = e.node_vars[1];
		NodeVarAtElem& nv3 = e.node_vars[2];

		n1.has_mp = true;
		n2.has_mp = true;
		n3.has_mp = true;
		
		//e.pcls = nullptr;
		e.pcl_vol = 0.0;
		e.s11 = 0.0;
		e.s22 = 0.0;
		e.s12 = 0.0;
		nv1.init();
		nv2.init();
		nv3.init();
		for (size_t th_id = 0; th_id < thread_num; ++th_id)
		{
			ThreadElemData &ted = pteds[th_id][e_id];
			for (Particle* ppcl = ted.pcls; ppcl; ppcl = ppcl->next_pcl_in_elem)
			{
				// mixed integration
				Particle& pcl = *ppcl;
				pcl.vol = pcl.m / pcl.density;
				e.pcl_vol += pcl.vol;
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
		}

		e.s11 /= e.pcl_vol;
		e.s22 /= e.pcl_vol;
		e.s12 /= e.pcl_vol;
		if (e.pcl_vol > e.area)
			e.pcl_vol = e.area;

		// internal force
		// node 1
		nv1.fx -= (e.dN1_dx * e.s11 + e.dN1_dy * e.s12) * e.pcl_vol;
		nv1.fy -= (e.dN1_dx * e.s12 + e.dN1_dy * e.s22) * e.pcl_vol;
		// node 2
		nv2.fx -= (e.dN2_dx * e.s11 + e.dN2_dy * e.s12) * e.pcl_vol;
		nv2.fy -= (e.dN2_dx * e.s12 + e.dN2_dy * e.s22) * e.pcl_vol;
		// node 3
		nv3.fx -= (e.dN3_dx * e.s11 + e.dN3_dy * e.s12) * e.pcl_vol;
		nv3.fy -= (e.dN3_dx * e.s12 + e.dN3_dy * e.s22) * e.pcl_vol;
	}
}

void Step_T2D_ME_p::update_node_a_and_v(unsigned int th_id)
{
	Model_T2D_ME_p& md = *static_cast<Model_T2D_ME_p*>(model);

	for (size_t n_id = cur_node_id++; n_id < md.node_num; n_id = cur_node_id++)
	{
		Node &n = md.nodes[n_id];

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
			NodeVarAtElem& nv = e.node_vars[n2e.n_id];
			n.m += nv.m;
			n.vx_ori += nv.vx;
			n.vy_ori += nv.vy;
			n.fx += nv.fx;
			n.fy += nv.fy;
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

void Step_T2D_ME_p::apply_a_and_v_bcs(unsigned int th_id)
{
	Model_T2D_ME_p& md = *static_cast<Model_T2D_ME_p*>(model);

	for (size_t a_id = cur_ax_bc_id++; a_id < md.ax_num; a_id = cur_ax_bc_id++)
	{
		Node& n = md.nodes[md.axs[a_id].node_id];
		n.ax = md.axs[a_id].a;
		n.vx = n.vx_ori + n.ax * dtime;
		n.dux = n.vx * dtime;
	}

	for (size_t a_id = cur_ay_bc_id++; a_id < md.ay_num; a_id = cur_ay_bc_id++)
	{
		Node& n = md.nodes[md.ays[a_id].node_id];
		n.ay = md.ays[a_id].a;
		n.vy = n.vy_ori + n.ay * dtime;
		n.duy = n.vy * dtime;
	}

	for (size_t v_id = cur_vx_bc_id++; v_id < md.vx_num; v_id = cur_vx_bc_id++)
	{
		Node& n = md.nodes[md.vxs[v_id].node_id];
		n.ax = 0.0;
		n.vx = md.vxs[v_id].v;
		n.dux = n.vx * dtime;
	}

	for (size_t v_id = cur_vy_bc_id++; v_id < md.vy_num; v_id = cur_vy_bc_id++)
	{
		Node& n = md.nodes[md.vys[v_id].node_id];
		n.ay = 0.0;
		n.vy = md.vys[v_id].v;
		n.duy = n.vy * dtime;
	}
}

void Step_T2D_ME_p::cal_de_at_elem(unsigned int th_id)
{
	Model_T2D_ME_p& md = *static_cast<Model_T2D_ME_p*>(model);

	double de11, de22;
	for (size_t e_id = cur_elem_id++; e_id < md.elem_num; e_id = cur_elem_id++)
	{
		Element &e = md.elems[e_id];
		
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
	}
}

void Step_T2D_ME_p::map_de_vol_from_elem_to_node(unsigned int th_id)
{
	Model_T2D_ME_p& md = *static_cast<Model_T2D_ME_p*>(model);

	for (size_t n_id = cur_node_id++; n_id < md.node_num; n_id = cur_node_id++)
	{
		Node& n = md.nodes[n_id];
		
		if (!n.has_mp)
			continue;
		
		n.de_vol_by_3 = 0.0;
		n.pcl_vol = 0.0;
		for (size_t e_id = 0; e_id < n.n2e_num; ++e_id)
		{
			NodeToElem& n2e = n.n2es[e_id];
			Element& e = md.elems[n2e.e_id];
			n.de_vol_by_3 += e.de_vol_by_3 * e.pcl_vol / 3.0;
			n.pcl_vol += e.pcl_vol / 3.0;
		}
		n.de_vol_by_3 /= n.pcl_vol;
	}
}

void Step_T2D_ME_p::update_pcl_vars(unsigned int th_id)
{
	Model_T2D_ME_p& md = *static_cast<Model_T2D_ME_p *>(model);
	
	double de_vol_by_3, de11, de22, de12;
	for (size_t pcl_id = cur_pcl_id++; pcl_id < md.pcl_num; pcl_id = cur_pcl_id++)
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

		// displacement
		pcl.ux += n1.dux * pcl.N1 + n2.dux * pcl.N2 + n3.dux * pcl.N3;
		pcl.uy += n1.duy * pcl.N1 + n2.duy * pcl.N2 + n3.duy * pcl.N3;

		// update position
		pcl.x = pcl.x_ori + pcl.ux;
		pcl.y = pcl.y_ori + pcl.uy;

		// strain
		//de_vol_by_3 = n1.de_vol_by_3 * pcl.N1 + n2.de_vol_by_3 * pcl.N2 + n3.de_vol_by_3 * pcl.N3;
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

		// density
		pcl.density /= 1.0 + de_vol_by_3 * 3.0;
	}
}
