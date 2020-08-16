#include "Simulations_pcp.h"

#include <iostream>
#include <assert.h>
#include <cmath>

#include "MaterialModel.h"
#include "Step_T2D_CHM_p.h"

Step_T2D_CHM_p::Step_T2D_CHM_p(const char* _name) :
	Step(_name, "Step_T2D_CHM_p", &solve_substep_T2D_CHM_p),
	model(nullptr), damping_ratio(0.0),
	thread_num(0), not_yet_completed(false), thread_data_mem(nullptr)
{

}

Step_T2D_CHM_p::~Step_T2D_CHM_p()
{
	join_all_threads();
	clear_thread_data();
}

int Step_T2D_CHM_p::init_calculation()
{
	Model_T2D_CHM_p &md = *model;

	if (is_first_step)
	{
		md.init_bcs();
	}

	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle &pcl = md.pcls[pcl_id];
		pcl.pe = md.elems;
		pcl.vol_s = pcl.m_s / pcl.density_s;
		pcl.x_ori = pcl.x;
		pcl.y_ori = pcl.y;
		pcl.ux_s = 0.0;
		pcl.uy_s = 0.0;
		pcl.x_f_ori = pcl.x_f;
		pcl.y_f_ori = pcl.y_f;
		pcl.ux_f = 0.0;
		pcl.uy_f = 0.0;
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

int Step_T2D_CHM_p::finalize_calculation()
{
	join_all_threads();
	clear_thread_data();
	return 0;
}

void Step_T2D_CHM_p::spawn_all_threads()
{
	Model_T2D_CHM_p& md = *model;
	if (thread_num > 1)
	{
		size_t spawn_thread_num = thread_num - 1;
		cal_threads.resize(spawn_thread_num);
		for (unsigned int th_id = 0; th_id < spawn_thread_num; ++th_id)
			cal_threads[th_id] = std::thread(&Step_T2D_CHM_p::cal_thread_func, this, th_id + 1);
	}
}

void Step_T2D_CHM_p::join_all_threads()
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

void Step_T2D_CHM_p::alloc_thread_data()
{
	Model_T2D_CHM_p &md = *model;
	size_t th_data_pt_size = sizeof(ThreadData**) * thread_num;
	size_t th_elem_data_pt_size = sizeof(ThreadElemData**) * thread_num;
	size_t thread_data_size = sizeof(ThreadData) + sizeof(ThreadElemData) * md.elem_num;
	size_t total_data_size = th_data_pt_size + th_elem_data_pt_size
						   + thread_data_size * thread_num;
	thread_data_mem = new char[total_data_size];

	pth_datas = reinterpret_cast<ThreadData**>(thread_data_mem);
	pth_elem_datas = reinterpret_cast<ThreadElemData**>(thread_data_mem + th_data_pt_size);
	char* data_address = thread_data_mem + th_data_pt_size + th_elem_data_pt_size;
	for (size_t th_id = 0; th_id < thread_num; ++th_id)
	{
		pth_datas[th_id] = reinterpret_cast<ThreadData*>(data_address + thread_data_size * th_id);
		pth_elem_datas[th_id] = reinterpret_cast<ThreadElemData*>((char*)(pth_datas[th_id]) + sizeof(ThreadData));

#ifdef _DEBUG
		ThreadElemData* cur_elem_datas = pth_elem_datas[th_id];
		for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
			cur_elem_datas[e_id].elem_id = e_id;
#endif
	}
}

void Step_T2D_CHM_p::clear_thread_data()
{
	if (thread_data_mem)
	{
		delete[] thread_data_mem;
		thread_data_mem = nullptr;
	}
}

void Step_T2D_CHM_p::cal_thread_func(unsigned int th_id)
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

		// rigid circle
		cal_contact_force(th_id);
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

int solve_substep_T2D_CHM_p(void *_self)
{
	typedef Model_T2D_CHM_p::Node Node;
	typedef Model_T2D_CHM_p::Element Element;
	typedef Model_T2D_CHM_p::Particle Particle;
	typedef Step_T2D_CHM_p::ThreadData ThreadData;
	typedef Step_T2D_CHM_p::ThreadElemData ThreadElemData;

	Step_T2D_CHM_p &self = *(Step_T2D_CHM_p *)(_self);
	Model_T2D_CHM_p &md = *self.model;

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

	// rigid circle
	self.cur_elem_id.store(0, std::memory_order_relaxed);
	self.cal_barrier.lift_barrier();
	self.cal_contact_force(0);
	self.cal_barrier.wait_for_others();
	// accumulate contact force on rigid circle
	RigidCircle& rc = md.rigid_circle;
	rc.reset_rf();
	for (size_t th_id = 0; th_id < self.thread_num; ++th_id)
	{
		ThreadData& th_data = *(self.pth_datas[th_id]);
		rc.add_rcf(th_data.rcf);
	}
	rc.update_motion(self.dtime);

	self.cur_node_id.store(0, std::memory_order_relaxed);
	self.cal_barrier.lift_barrier();
	self.update_node_a_and_v(0);
	self.cal_barrier.wait_for_others();

	self.cur_asx_bc_id.store(0, std::memory_order_relaxed);
	self.cur_asy_bc_id.store(0, std::memory_order_relaxed);
	self.cur_vsx_bc_id.store(0, std::memory_order_relaxed);
	self.cur_vsy_bc_id.store(0, std::memory_order_relaxed);
	self.cur_afx_bc_id.store(0, std::memory_order_relaxed);
	self.cur_afy_bc_id.store(0, std::memory_order_relaxed);
	self.cur_vfx_bc_id.store(0, std::memory_order_relaxed);
	self.cur_vfy_bc_id.store(0, std::memory_order_relaxed);
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

void Step_T2D_CHM_p::init_cal_vars(unsigned int th_id)
{
	Model_T2D_CHM_p& md = *static_cast<Model_T2D_CHM_p*>(model);

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

void Step_T2D_CHM_p::find_pcls_in_which_elems(unsigned int th_id)
{
	Model_T2D_CHM_p& md = *static_cast<Model_T2D_CHM_p*>(model);
	ThreadData& th_data = *pth_datas[th_id];
	ThreadElemData* th_elem_datas = pth_elem_datas[th_id];

	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
		th_elem_datas[e_id].init();

	size_t pe_id, pcl_id;
	while ((pcl_id = ATOM_INC(cur_pcl_id)) < md.pcl_num)
	{
		Particle& pcl = md.pcls[pcl_id];

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

void Step_T2D_CHM_p::map_pcl_vars_to_nodes_at_elems(unsigned int th_id)
{
	Model_T2D_CHM_p& md = *static_cast<Model_T2D_CHM_p*>(model);

	size_t e_id;
	double N_m_s, N_m_f;
	double n2_miu_div_k;
	double n2_miu_div_k_vrx_vol;
	double n2_miu_div_k_vry_vol;
	double fx_drag, fy_drag;
	double n_prod_p, one_min_n_prod_p;
	while ((e_id = ATOM_INC(cur_elem_id)) < md.elem_num)
	{
		Element& e = md.elems[e_id];

		if (!e.has_mp)
			continue;

		e.mi_pcl_vol = 0.0;
		e.mi_pcl_n = 0.0;
		e.s11 = 0.0;
		e.s22 = 0.0;
		e.s12 = 0.0;
		e.p = 0.0;

		Node& n1 = md.nodes[e.n1];
		Node& n2 = md.nodes[e.n2];
		Node& n3 = md.nodes[e.n3];
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
		for (Particle* ppcl = e.first_pcl(); e.not_last_pcl(ppcl);
			 ppcl = e.next_pcl(ppcl))
		{
			Particle& pcl = *ppcl;

			pcl.vol = pcl.vol_s / (1.0 - pcl.n);
			pcl.m_f = pcl.n * pcl.density_f * pcl.vol;
			double n2_miu_div_k = pcl.n * pcl.n * md.miu / md.k;
			double n2_miu_div_k_vrx_vol = n2_miu_div_k * (pcl.vx_f - pcl.vx_s) * pcl.vol;
			double n2_miu_div_k_vry_vol = n2_miu_div_k * (pcl.vy_f - pcl.vy_s) * pcl.vol;
			// mixed integration
			e.mi_pcl_vol += pcl.vol;
			e.mi_pcl_n += pcl.vol_s;
			e.s11 += pcl.s11 * pcl.vol;
			e.s22 += pcl.s22 * pcl.vol;
			e.s12 += pcl.s12 * pcl.vol;
			e.p += pcl.p * pcl.vol;

			// node 1
			// mixture phase
			N_m_s = pcl.N1 * pcl.m_s;
			nv1.m_s += N_m_s;
			nv1.vx_s += N_m_s * pcl.vx_s;
			nv1.vy_s += N_m_s * pcl.vy_s;
			// fluid phase
			N_m_f = pcl.N1 * pcl.m_f;
			nv1.m_f += N_m_f;
			nv1.vx_f += N_m_f * pcl.vx_f;
			nv1.vy_f += N_m_f * pcl.vy_f;
			// solid - fluid interaction
			fx_drag = pcl.N1 * n2_miu_div_k_vrx_vol;
			fy_drag = pcl.N1 * n2_miu_div_k_vry_vol;
			nv1.fx_s += fx_drag;
			nv1.fy_s += fy_drag;
			nv1.fx_f -= fx_drag;
			nv1.fy_f -= fy_drag;

			// node 2
			// mixture phase
			N_m_s = pcl.N2 * pcl.m_s;
			nv2.m_s += N_m_s;
			nv2.vx_s += N_m_s * pcl.vx_s;
			nv2.vy_s += N_m_s * pcl.vy_s;
			// fluid phase
			N_m_f = pcl.N2 * pcl.m_f;
			nv2.m_f += N_m_f;
			nv2.vx_f += N_m_f * pcl.vx_f;
			nv2.vy_f += N_m_f * pcl.vy_f;
			// solid - fluid interaction
			fx_drag = pcl.N2 * n2_miu_div_k_vrx_vol;
			fy_drag = pcl.N2 * n2_miu_div_k_vry_vol;
			nv2.fx_s += fx_drag;
			nv2.fy_s += fy_drag;
			nv2.fx_f -= fx_drag;
			nv2.fy_f -= fy_drag;

			// node 3
			// mixture phase
			N_m_s = pcl.N3 * pcl.m_s;
			nv3.m_s += N_m_s;
			nv3.vx_s += N_m_s * pcl.vx_s;
			nv3.vy_s += N_m_s * pcl.vy_s;
			// fluid phase
			N_m_f = pcl.N3 * pcl.m_f;
			nv3.m_f += N_m_f;
			nv3.vx_f += N_m_f * pcl.vx_f;
			nv3.vy_f += N_m_f * pcl.vy_f;
			// solid - fluid interaction
			fx_drag += pcl.N3 * n2_miu_div_k_vrx_vol;
			fy_drag += pcl.N3 * n2_miu_div_k_vry_vol;
			nv3.fx_s += fx_drag;
			nv3.fy_s += fy_drag;
			nv3.fx_f -= fx_drag;
			nv3.fy_f -= fy_drag;

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
			if (pcl.has_fx_f_ext)
			{
				nv1.fx_f += pcl.N1 * pcl.fx_f_ext;
				nv2.fx_f += pcl.N2 * pcl.fx_f_ext;
				nv3.fx_f += pcl.N3 * pcl.fx_f_ext;
			}
			if (pcl.has_fy_f_ext)
			{
				nv1.fy_f += pcl.N1 * pcl.fy_f_ext;
				nv2.fy_f += pcl.N2 * pcl.fy_f_ext;
				nv3.fy_f += pcl.N3 * pcl.fy_f_ext;
			}
		}

		e.mi_pcl_n = 1.0 - e.mi_pcl_n / e.mi_pcl_vol; // 1.0 - Vs / V
		e.s11 /= e.mi_pcl_vol;
		e.s22 /= e.mi_pcl_vol;
		e.s12 /= e.mi_pcl_vol;
		e.p /= e.mi_pcl_vol;
		if (e.mi_pcl_vol > e.area)
			e.mi_pcl_vol = e.area;

		// internal force
		n_prod_p = e.mi_pcl_n * e.p;
		one_min_n_prod_p = e.p - n_prod_p;
		// node 1
		nv1.fx_s -= (e.dN1_dx * (e.s11 - one_min_n_prod_p) + e.dN1_dy * e.s12) * e.mi_pcl_vol;
		nv1.fy_s -= (e.dN1_dx * e.s12 + e.dN1_dy * (e.s22 - one_min_n_prod_p)) * e.mi_pcl_vol;
		nv1.fx_f -= (e.dN1_dx * -n_prod_p) * e.mi_pcl_vol;
		nv1.fy_f -= (e.dN1_dy * -n_prod_p) * e.mi_pcl_vol;
		// node 2
		nv2.fx_s -= (e.dN2_dx * (e.s11 - one_min_n_prod_p) + e.dN2_dy * e.s12) * e.mi_pcl_vol;
		nv2.fy_s -= (e.dN2_dx * e.s12 + e.dN2_dy * (e.s22 - one_min_n_prod_p)) * e.mi_pcl_vol;
		nv2.fx_f -= (e.dN2_dx * -n_prod_p) * e.mi_pcl_vol;
		nv2.fy_f -= (e.dN2_dy * -n_prod_p) * e.mi_pcl_vol;
		// node 3
		nv3.fx_s -= (e.dN3_dx * (e.s11 - one_min_n_prod_p) + e.dN3_dy * e.s12) * e.mi_pcl_vol;
		nv3.fy_s -= (e.dN3_dx * e.s12 + e.dN3_dy * (e.s22 - one_min_n_prod_p)) * e.mi_pcl_vol;
		nv3.fx_f -= (e.dN3_dx * -n_prod_p) * e.mi_pcl_vol;
		nv3.fy_f -= (e.dN3_dy * -n_prod_p) * e.mi_pcl_vol;
	}
}

void Step_T2D_CHM_p::cal_contact_force(unsigned int th_id)
{
	Model_T2D_CHM_p& md = *static_cast<Model_T2D_CHM_p*>(model);
	RigidCircle& rc = md.rigid_circle;
	RigidCircleForce& rcf = pth_datas[th_id]->rcf;

	double rc_x = rc.get_x();
	double rc_y = rc.get_y();
	rcf.reset_rf();
	double dist, norm_x, norm_y;
	double fs_cont, fsx_cont, fsy_cont;
	double ff_cont, ffx_cont, ffy_cont;
	size_t e_id;
	while ((e_id = ATOM_INC(cur_elem_id)) < md.elem_num)
	{
		Element& e = md.elems[e_id];

		if (!e.has_mp /* || not contact circle */)
			continue;

		NodeVarAtElem& nv1 = e.node_vars[0];
		NodeVarAtElem& nv2 = e.node_vars[1];
		NodeVarAtElem& nv3 = e.node_vars[2];

		for (Particle* ppcl = e.first_pcl(); e.not_last_pcl(ppcl);
			 ppcl = e.next_pcl(ppcl))
		{
			Particle& pcl = *ppcl;

			// solid particle
			if (rc.detect_collision_with_point(pcl.x, pcl.y, pcl.vol, dist, norm_x, norm_y))
			{
				fs_cont = md.Ks_cont * dist;
				fsx_cont = fs_cont * norm_x;
				fsy_cont = fs_cont * norm_y;
				// add reaction force to rigid object
				rcf.add_rf(pcl.x, pcl.y, -fsx_cont, -fsy_cont, rc_x, rc_y);
				// node 1
				nv1.fx_s += pcl.N1 * fsx_cont;
				nv1.fy_s += pcl.N1 * fsy_cont;
				// node 2
				nv2.fx_s += pcl.N2 * fsx_cont;
				nv2.fy_s += pcl.N2 * fsy_cont;
				// node 3
				nv3.fx_s += pcl.N3 * fsx_cont;
				nv3.fy_s += pcl.N3 * fsy_cont;
			}

			// fluid particle
			if (rc.detect_collision_with_point(pcl.x_f, pcl.y_f, pcl.vol, dist, norm_x, norm_y))
			{
				ff_cont = md.Kf_cont * dist;
				ffx_cont = ff_cont * norm_x;
				ffy_cont = ff_cont * norm_y;
				// add reaction force to rigid object
				rcf.add_rf(pcl.x_f, pcl.y_f, -ffx_cont, -ffy_cont, rc_x, rc_y);
				// node 1
				nv1.fx_f += pcl.N1 * ffx_cont;
				nv1.fy_f += pcl.N1 * ffy_cont;
				// node 2
				nv2.fx_f += pcl.N2 * ffx_cont;
				nv2.fy_f += pcl.N2 * ffy_cont;
				// node 3
				nv3.fx_f += pcl.N3 * ffx_cont;
				nv3.fy_f += pcl.N3 * ffy_cont;
			}
		}
	}
}

void Step_T2D_CHM_p::update_node_a_and_v(unsigned int th_id)
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
		n.m_f = 0.0;
		n.vx_f_ori = 0.0;
		n.vy_f_ori = 0.0;
		n.fx_f = 0.0;
		n.fy_f = 0.0;

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
				n.m_f += nv.m_f;
				n.vx_f_ori += nv.vx_f;
				n.vy_f_ori += nv.vy_f;
				n.fx_f += nv.fx_f;
				n.fy_f += nv.fy_f;
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

		n.ax_f = n.fx_f / n.m_f;
		n.ay_f = n.fy_f / n.m_f;
		n.vx_f_ori /= n.m_f;
		n.vy_f_ori /= n.m_f;
		n.vx_f = n.vx_f_ori + n.ax_f * dtime;
		n.vy_f = n.vy_f_ori + n.ay_f * dtime;
		n.dux_f = n.vx_f * dtime;
		n.duy_f = n.vy_f * dtime;
	}
}

void Step_T2D_CHM_p::apply_a_and_v_bcs(unsigned int th_id)
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
	while ((a_id = ATOM_INC(cur_afx_bc_id)) < md.afx_num)
	{
		Node& n = md.nodes[md.afxs[a_id].node_id];
		n.ax_f = md.afxs[a_id].a;
		n.vx_f = n.vx_f_ori + n.ax_f * dtime;
		n.dux_f = n.vx_f * dtime;
	}
	while ((a_id = ATOM_INC(cur_afy_bc_id)) < md.afy_num)
	{
		Node& n = md.nodes[md.afys[a_id].node_id];
		n.ay_f = md.afys[a_id].a;
		n.vy_f = n.vy_f_ori + n.ay_f * dtime;
		n.duy_f = n.vy_f * dtime;
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
	while ((v_id = ATOM_INC(cur_vfx_bc_id)) < md.vfx_num)
	{
		Node& n = md.nodes[md.vfxs[v_id].node_id];
		n.ax_f = 0.0;
		n.vx_f = md.vfxs[v_id].v;
		n.dux_f = n.vx_f * dtime;
	}
	while ((v_id = ATOM_INC(cur_vfy_bc_id)) < md.vfy_num)
	{
		Node& n = md.nodes[md.vfys[v_id].node_id];
		n.ay_f = 0.0;
		n.vy_f = md.vfys[v_id].v;
		n.duy_f = n.vy_f * dtime;
	}
}

void Step_T2D_CHM_p::cal_de_at_elem(unsigned int th_id)
{
	Model_T2D_CHM_p& md = *static_cast<Model_T2D_CHM_p*>(model);

	double de11, de22, de_mean;
	size_t e_id;
	while ((e_id = ATOM_INC(cur_elem_id)) < md.elem_num)
	{
		Element& e = md.elems[e_id];

		if (!e.has_mp)
			continue;

		Node& n1 = md.nodes[e.n1];
		Node& n2 = md.nodes[e.n2];
		Node& n3 = md.nodes[e.n3];
		de11 = n1.dux_s * e.dN1_dx + n2.dux_s * e.dN2_dx + n3.dux_s * e.dN3_dx;
		de22 = n1.duy_s * e.dN1_dy + n2.duy_s * e.dN2_dy + n3.duy_s * e.dN3_dy;
		e.de12 = (n1.dux_s * e.dN1_dy + n2.dux_s * e.dN2_dy + n3.dux_s * e.dN3_dy
				+ n1.duy_s * e.dN1_dx + n2.duy_s * e.dN2_dx + n3.duy_s * e.dN3_dx) * 0.5;
		// volumetric strain of solid phase
		e.de_vol_s = de11 + de22;
		de_mean = e.de_vol_s / 3.0;
		e.dde11 = de11 - de_mean;
		e.dde22 = de22 - de_mean;
		// "volumetric strain" of fluid phase, take compression as positive
		e.de_vol_f = (1.0 - e.mi_pcl_n) / e.mi_pcl_n * -e.de_vol_s
				   - (n1.dux_f * e.dN1_dx + n2.dux_f * e.dN2_dx + n3.dux_f * e.dN3_dx
					+ n1.duy_f * e.dN1_dy + n2.duy_f * e.dN2_dy + n3.duy_f * e.dN3_dy);

		// strain enhancement
		NodeVarAtElem& nv1 = e.node_vars[0];
		NodeVarAtElem& nv2 = e.node_vars[1];
		NodeVarAtElem& nv3 = e.node_vars[2];
		double vol_de_vol_s, vol_de_vol_f;
		for (Particle* ppcl = e.first_pcl(); e.not_last_pcl(ppcl);
			 ppcl = e.next_pcl(ppcl))
		{
			Particle& pcl = *ppcl;
			vol_de_vol_s = e.de_vol_s * pcl.vol;
			vol_de_vol_f = e.de_vol_f * pcl.vol;
			// node 1
			nv1.de_vol_s += pcl.N1 * vol_de_vol_s;
			nv1.de_vol_f += pcl.N1 * vol_de_vol_f;
			nv1.se_pcl_vol += pcl.N1 * pcl.vol;
			// node 2
			nv2.de_vol_s += pcl.N2 * vol_de_vol_s;
			nv2.de_vol_f += pcl.N2 * vol_de_vol_f;
			nv2.se_pcl_vol += pcl.N2 * pcl.vol;
			// node 3
			nv3.de_vol_s += pcl.N3 * vol_de_vol_s;
			nv3.de_vol_f += pcl.N3 * vol_de_vol_f;
			nv3.se_pcl_vol += pcl.N3 * pcl.vol;
		}
	}
}

void Step_T2D_CHM_p::map_de_vol_from_elem_to_node(unsigned int th_id)
{
	Model_T2D_CHM_p& md = *static_cast<Model_T2D_CHM_p*>(model);

	size_t n_id;
	while ((n_id = ATOM_INC(cur_node_id)) < md.node_num)
	{
		Node& n = md.nodes[n_id];

		if (!n.has_mp)
			continue;

		n.de_vol_s = 0.0;
		n.de_vol_f = 0.0;
		n.se_pcl_vol = 0.0;
		for (size_t e_id = 0; e_id < n.n2e_num; ++e_id)
		{
			NodeToElem& n2e = n.n2es[e_id];
			Element& e = md.elems[n2e.e_id];
			if (e.has_mp)
			{
				NodeVarAtElem& nv = e.node_vars[n2e.n_id];
				n.de_vol_s += nv.de_vol_s;
				n.de_vol_f += nv.de_vol_f;
				n.se_pcl_vol += nv.se_pcl_vol;
			}
		}
		n.de_vol_s /= n.se_pcl_vol;
		n.de_vol_f /= n.se_pcl_vol;
	}
}

void Step_T2D_CHM_p::update_pcl_vars(unsigned int th_id)
{
	Model_T2D_CHM_p& md = *static_cast<Model_T2D_CHM_p*>(model);

	double de11, de22, de12;
	double de_vol_s, de_vol_f;
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
		pcl.vx_f += (n1.ax_f * pcl.N1 + n2.ax_f * pcl.N2 + n3.ax_f * pcl.N3) * dtime;
		pcl.vy_f += (n1.ay_f * pcl.N1 + n2.ay_f * pcl.N2 + n3.ay_f * pcl.N3) * dtime;

		// displacement
		pcl.ux_s += n1.dux_s * pcl.N1 + n2.dux_s * pcl.N2 + n3.dux_s * pcl.N3;
		pcl.uy_s += n1.duy_s * pcl.N1 + n2.duy_s * pcl.N2 + n3.duy_s * pcl.N3;
		pcl.ux_f += n1.dux_f * pcl.N1 + n2.dux_f * pcl.N2 + n3.dux_f * pcl.N3;
		pcl.uy_f += n1.duy_f * pcl.N1 + n2.duy_f * pcl.N2 + n3.duy_f * pcl.N3;

		// update position
		pcl.x = pcl.x_ori + pcl.ux_s;
		pcl.y = pcl.y_ori + pcl.uy_s;
		pcl.x_f = pcl.x_f_ori + pcl.ux_f;
		pcl.y_f = pcl.y_f_ori + pcl.uy_f;

		// strain enhancement appraoch
		de_vol_s = (n1.de_vol_s + n2.de_vol_s + n3.de_vol_s) / 3.0;
		//de_vol_s = n1.de_vol_s * pcl.N1 + n2.de_vol_s * pcl.N2 + n3.de_vol_s * pcl.N3;
		// porosity
		pcl.n = (de_vol_s + pcl.n) / (1.0 + de_vol_s);

		// strain
		de_vol_s /= 3.0;
		de11 = e.dde11 + de_vol_s;
		de22 = e.dde22 + de_vol_s;
		de12 = e.de12;
		// strain
		pcl.e11 += de11;
		pcl.e22 += de22;
		pcl.e12 += de12;

		// update stress using constitutive model
		double dstrain[6] = { de11, de22, 0.0, de12, 0.0, 0.0 };
		pcl.mm->integrate(dstrain);
		const double* dstress = pcl.mm->get_dstress();
		pcl.s11 += dstress[0];
		pcl.s22 += dstress[1];
		pcl.s12 += dstress[3];

		// pore pressure
		de_vol_f = (n1.de_vol_f + n2.de_vol_f + n3.de_vol_f) / 3.0;
		//de_vol_f = n1.de_vol_f * pcl.N1 + n2.de_vol_f * pcl.N2 + n3.de_vol_f * pcl.N3;
		pcl.p += md.Kf * de_vol_f;
		// fluid density
		pcl.density_f /= 1.0 - de_vol_f;
	}
}
