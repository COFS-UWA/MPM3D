#include "Simulations_pcp.h"

#include <cmath>
#include <iostream>

#include "MaterialModel.h"
#include "Step_T2D_ME_p.h"

Step_T2D_ME_p::Step_T2D_ME_p(const char* _name) :
	Step(_name, "Step_T2D_ME_p", &solve_substep_T2D_ME_p),
	model(nullptr), damping_ratio(0.0),
	thread_num(0), not_yet_completed(false), thread_data_mem(nullptr)
{

}

Step_T2D_ME_p::~Step_T2D_ME_p()
{
	join_all_threads();
	clear_thread_data();
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

	// start calculation
	if (thread_num == 0)
		thread_num = std::thread::hardware_concurrency();

	// init global data
	not_yet_completed.store(true, std::memory_order_relaxed);
	step_barrier.set_thread_num(thread_num);
	cal_barrier.set_thread_num(thread_num);
	// init thread wise data
	alloc_thread_data();
	for (size_t th_id = 0; th_id < thread_num; ++th_id)
	{
		ThreadData& th_data = *pth_datas[th_id];
		th_data.not_combined_yet.test_and_set(std::memory_order_relaxed);
	}

	// spawn threads
	spawn_all_threads();
	
	step_barrier.wait_for_others();
	return 0;
}

int Step_T2D_ME_p::finalize_calculation()
{
	join_all_threads();
	clear_thread_data();
	return 0;
}

void Step_T2D_ME_p::spawn_all_threads()
{
	Model_T2D_ME_p& md = *model;
	if (!md.rigid_circle_is_valid()) // no rigid circle
	{
		cal_substep_func = &solve_substep_T2D_ME_p;
		if (thread_num > 1)
		{
			size_t spawn_thread_num = thread_num - 1;
			cal_threads.resize(spawn_thread_num);
			for (unsigned int th_id = 0; th_id < spawn_thread_num; ++th_id)
				cal_threads[th_id] = std::thread(&Step_T2D_ME_p::cal_thread_func, this, th_id + 1);
		}
	}
	else // has rigid circle
	{
		cal_substep_func = &solve_substep_T2D_ME_p_RigidCircle;
		if (thread_num > 1)
		{
			size_t spawn_thread_num = thread_num - 1;
			cal_threads.resize(spawn_thread_num);
			for (unsigned int th_id = 0; th_id < spawn_thread_num; ++th_id)
				cal_threads[th_id] = std::thread(&Step_T2D_ME_p::cal_thread_func_RigidCircle, this, th_id + 1);
		}
	}
}

void Step_T2D_ME_p::join_all_threads()
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

void Step_T2D_ME_p::alloc_thread_data()
{
	Model_T2D_ME_p& md = *model;
	size_t th_data_pt_size = sizeof(ThreadData **) * thread_num;
	size_t th_elem_data_pt_size = sizeof(ThreadElemData **) * thread_num;
	size_t thread_data_size = sizeof(ThreadData) + sizeof(ThreadElemData) * md.elem_num;
	size_t total_data_size = th_data_pt_size + th_elem_data_pt_size
						   + thread_data_size * thread_num;
	thread_data_mem = new char[total_data_size];

	pth_datas = reinterpret_cast<ThreadData**>(thread_data_mem);
	pth_elem_datas = reinterpret_cast<ThreadElemData **>(thread_data_mem + th_data_pt_size);
	char* data_address = thread_data_mem + th_data_pt_size + th_elem_data_pt_size;
	for (size_t th_id = 0; th_id < thread_num; ++th_id)
	{
		pth_datas[th_id] = reinterpret_cast<ThreadData *>(data_address + thread_data_size * th_id);
		pth_elem_datas[th_id] = reinterpret_cast<ThreadElemData *>((char *)(pth_datas[th_id]) + sizeof(ThreadData));
		
#ifdef _DEBUG
		ThreadElemData *cur_elem_datas = pth_elem_datas[th_id];
		for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
			cur_elem_datas[e_id].elem_id = e_id;
#endif
	}
}

void Step_T2D_ME_p::clear_thread_data()
{
	if (thread_data_mem)
	{
		delete[] thread_data_mem;
		thread_data_mem = nullptr;
	}
}

void Step_T2D_ME_p::cal_thread_func(unsigned int th_id)
{
	step_barrier.wait();

	while (not_yet_completed.load(std::memory_order_relaxed))
	{
		init_cal_vars(th_id);
		cal_barrier.wait();

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
	typedef Step_T2D_ME_p::ThreadData ThreadData;
	typedef Step_T2D_ME_p::ThreadElemData ThreadElemData;

	Step_T2D_ME_p &self = *static_cast<Step_T2D_ME_p *>(_self);
	Model_T2D_ME_p &md = static_cast<Model_T2D_ME_p &>(self.get_model());

	self.cur_node_id.store(0, std::memory_order_relaxed);
	self.cur_elem_id.store(0, std::memory_order_relaxed);
	self.step_barrier.lift_barrier();
	self.init_cal_vars(0);
	self.cal_barrier.wait_for_others();

	self.cur_pcl_id.store(0, std::memory_order_relaxed);
	self.cal_barrier.lift_barrier();
	self.find_pcls_in_which_elems(0);
	self.cal_barrier.wait_for_others();

	self.cur_elem_id.store(0, std::memory_order_relaxed);
	self.cal_barrier.lift_barrier();
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

void Step_T2D_ME_p::init_cal_vars(unsigned int th_id)
{
	Model_T2D_ME_p& md = *static_cast<Model_T2D_ME_p*>(model);

	// init nodes
	size_t n_id;
	while ((n_id = ATOM_INC(cur_node_id)) < md.node_num)
	{
		Node& n = md.nodes[n_id];
		n.has_mp = false;
	}

	// init elements
	size_t e_id;
	while ((e_id = ATOM_INC(cur_elem_id)) < md.elem_num)
	{
		Element& e = md.elems[e_id];
		e.has_mp = false;
	}
}

void Step_T2D_ME_p::find_pcls_in_which_elems(unsigned int th_id)
{
	Model_T2D_ME_p& md = *static_cast<Model_T2D_ME_p*>(model);
	ThreadData& th_data = *pth_datas[th_id];
	ThreadElemData *th_elem_datas = pth_elem_datas[th_id];

	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
		th_elem_datas[e_id].init();
	
	size_t pcl_id, pe_id;
	while ((pcl_id = ATOM_INC(cur_pcl_id)) < md.pcl_num)
	{
		Particle &pcl = md.pcls[pcl_id];
		
		if (pcl.pe == nullptr)
			continue;

		pcl.pe = md.find_in_which_element(pcl);
		if (pcl.pe == nullptr)
			continue;
		
		pe_id = pcl.pe->id;
		th_elem_datas[pe_id].add_pcl(pcl);

		Element& e = md.elems[pe_id];
		e.has_mp = true;
	}

	// combine the pcl list in each element
	size_t stride = 2;
	size_t half_stride = 1;
	size_t oth_id;
	while (half_stride < thread_num)
	{
		if (th_id & (stride - 1))
			break;

		oth_id = th_id + half_stride;

		if (oth_id < thread_num)
		{
			ThreadData& oth_data = *pth_datas[oth_id];
			ThreadElemData* oth_elem_data = pth_elem_datas[oth_id];
			while (oth_data.not_combined_yet.test_and_set(std::memory_order_relaxed));
			// combine the data
			for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
				th_elem_datas[e_id].combine(oth_elem_data[e_id]);
		}

		stride = stride << 1;
		half_stride = half_stride << 1;
	}

	th_data.not_combined_yet.clear(std::memory_order_relaxed);
}

void Step_T2D_ME_p::map_pcl_vars_to_nodes_at_elems(unsigned int th_id)
{
	Model_T2D_ME_p& md = *static_cast<Model_T2D_ME_p*>(model);
	
	size_t e_id;
	while ((e_id = ATOM_INC(cur_elem_id)) < md.elem_num)
	{
		Element &e = md.elems[e_id];

		if (!e.has_mp)
			continue;
		
		e.mi_pcl_vol = 0.0;
		e.s11 = 0.0;
		e.s22 = 0.0;
		e.s12 = 0.0;
		
		Node &n1 = md.nodes[e.n1];
		Node &n2 = md.nodes[e.n2];
		Node &n3 = md.nodes[e.n3];
		n1.has_mp = true;
		n2.has_mp = true;
		n3.has_mp = true;
		
		NodeVarAtElem& nv1 = e.node_vars[0];
		NodeVarAtElem& nv2 = e.node_vars[1];
		NodeVarAtElem& nv3 = e.node_vars[2];
		nv1.init();
		nv2.init();
		nv3.init();

		e.pcls = pth_elem_datas[0][e_id].get_top_pcl();
		for (Particle * ppcl = e.first_pcl(); e.not_last_pcl(ppcl);
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

void Step_T2D_ME_p::update_node_a_and_v(unsigned int th_id)
{
	Model_T2D_ME_p& md = *static_cast<Model_T2D_ME_p*>(model);

	size_t n_id;
	while ((n_id = ATOM_INC(cur_node_id)) < md.node_num)
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

void Step_T2D_ME_p::apply_a_and_v_bcs(unsigned int th_id)
{
	Model_T2D_ME_p& md = *static_cast<Model_T2D_ME_p*>(model);

	size_t a_id;
	while ((a_id = ATOM_INC(cur_ax_bc_id)) < md.ax_num)
	{
		Node& n = md.nodes[md.axs[a_id].node_id];
		n.ax = md.axs[a_id].a;
		n.vx = n.vx_ori + n.ax * dtime;
		n.dux = n.vx * dtime;
	}
	while ((a_id = ATOM_INC(cur_ay_bc_id)) < md.ay_num)
	{
		Node& n = md.nodes[md.ays[a_id].node_id];
		n.ay = md.ays[a_id].a;
		n.vy = n.vy_ori + n.ay * dtime;
		n.duy = n.vy * dtime;
	}

	size_t v_id;
	while ((v_id = ATOM_INC(cur_vx_bc_id)) < md.vx_num)
	{
		Node& n = md.nodes[md.vxs[v_id].node_id];
		n.ax = 0.0;
		n.vx = md.vxs[v_id].v;
		n.dux = n.vx * dtime;
	}
	while ((v_id = ATOM_INC(cur_vy_bc_id)) < md.vy_num)
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
	size_t e_id;
	while ((e_id = ATOM_INC(cur_elem_id)) < md.elem_num)
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

void Step_T2D_ME_p::map_de_vol_from_elem_to_node(unsigned int th_id)
{
	Model_T2D_ME_p& md = *static_cast<Model_T2D_ME_p*>(model);

	size_t n_id;
	while ((n_id = ATOM_INC(cur_node_id)) < md.node_num)
	{
		Node& n = md.nodes[n_id];
		
		if (!n.has_mp)
			continue;

		n.de_vol_by_3 = 0.0;
		n.se_pcl_vol = 0.0;
		for (size_t e_id = 0; e_id < n.n2e_num; ++e_id)
		{
			NodeToElem& n2e = n.n2es[e_id];
			Element &e = md.elems[n2e.e_id];
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

void Step_T2D_ME_p::update_pcl_vars(unsigned int th_id)
{
	Model_T2D_ME_p& md = *static_cast<Model_T2D_ME_p *>(model);
	
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
