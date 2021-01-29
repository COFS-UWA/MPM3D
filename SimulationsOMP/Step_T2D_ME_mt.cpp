#include "SimulationsOMP_pcp.h"

#include <fstream>
#include <omp.h>

#include "Step_T2D_ME_mt.h"

#define one_third (1.0/3.0)
#define N_min (1.0e-8)
#define Block_Low(th_id, th_num, data_num) ((th_id)*(data_num)/(th_num))

#ifdef _DEBUG
static std::fstream res_file_t2d_me_mt;
#endif

Step_T2D_ME_mt::Step_T2D_ME_mt(const char* _name) : 
	Step_OMP(_name, "Step_T2D_ME_mt", &substep_func_omp_T2D_ME_mt) {}

Step_T2D_ME_mt::~Step_T2D_ME_mt() {}

int Step_T2D_ME_mt::init_calculation()
{
#ifdef _DEBUG
	res_file_t2d_me_mt.open("me_mt_res.txt", std::ios::out | std::ios::binary);
#endif

	Model_T2D_ME_mt &md = *(Model_T2D_ME_mt *)model;

	omp_set_num_threads(thread_num);

	pcl_m = md.pcl_m;
	pcl_bf = md.pcl_bf;
	pcl_t = md.pcl_t;
	pcl_pos = md.pcl_pos;
	pcl_vol = md.pcl_vol;
	pcl_mat_model = md.pcl_mat_model;

	Model_T2D_ME_mt::SortedPclVarArrays& md_spva0 = md.sorted_pcl_var_arrays[0];
	SortedPclVarArrays& spva0 = sorted_pcl_var_arrays[0];
	spva0.pcl_index = md_spva0.pcl_index;
	spva0.pcl_density = md_spva0.pcl_density;
	spva0.pcl_v = md_spva0.pcl_v;
	spva0.pcl_disp = md_spva0.pcl_disp;
	spva0.pcl_N = md_spva0.pcl_N;
	spva0.pcl_stress = md_spva0.pcl_stress;
	spva0.pcl_strain = md_spva0.pcl_strain;
	spva0.pcl_estrain = md_spva0.pcl_estrain;
	spva0.pcl_pstrain = md_spva0.pcl_pstrain;

	Model_T2D_ME_mt::SortedPclVarArrays& md_spva1 = md.sorted_pcl_var_arrays[1];
	SortedPclVarArrays& spva1 = sorted_pcl_var_arrays[1];
	spva1.pcl_index = md_spva1.pcl_index;
	spva1.pcl_density = md_spva1.pcl_density;
	spva1.pcl_v = md_spva1.pcl_v;
	spva1.pcl_disp = md_spva1.pcl_disp;
	spva1.pcl_N = md_spva1.pcl_N;
	spva1.pcl_stress = md_spva1.pcl_stress;
	spva1.pcl_strain = md_spva1.pcl_strain;
	spva1.pcl_estrain = md_spva1.pcl_estrain;
	spva1.pcl_pstrain = md_spva1.pcl_pstrain;

	elem_num = md.elem_num;
	node_num = md.node_num;

	elem_node_id = md.elem_node_id;
	elem_dN_ab = md.elem_dN_ab;
	elem_dN_c = md.elem_dN_c;
	elem_area = md.elem_area;

	elem_pcl_m = md.elem_pcl_m;
	elem_density = md.elem_density;
	elem_de = md.elem_de;
	elem_m_de_vol = md.elem_m_de_vol;

	elem_node_vm = md.elem_node_vm;
	elem_node_force = md.elem_node_force;

	node_a = md.node_a;
	node_v = md.node_v;
	node_has_vbc = md.node_has_vbc;
	node_am = md.node_am;
	node_de_vol = md.node_de_vol;
	
	cf_tmp.reset();
	//contact_substep_id = md.contact_substep_id;
	//prev_contact_pos = md.prev_contact_pos;
	//prev_contact_tan_force = md.prev_contact_tan_force;
	//memset(contact_substep_id, 0xFF, sizeof(size_t) * md.ori_pcl_num);
	if (md.has_rigid_rect())
	{
		prr = &md.get_rigid_rect();
		prr->reset_f_contact();
	}

	thread_datas = (ThreadData*)thread_mem.alloc(sizeof(ThreadData) * thread_num);

	char* cur_mem = (char*)cal_mem.alloc(
		sizeof(size_t) * (md.pcl_num * 4 + 4)
		+ sizeof(size_t) * (md.elem_num * 13 + 4)
		+ Cache_Alignment
		+ sizeof(size_t) * thread_num * 0x100 * 2);
	pcl_in_elems[0] = ((size_t*)cur_mem) + 1;
	cur_mem += sizeof(size_t) * (md.pcl_num + 2);
	pcl_in_elems[1] = ((size_t*)cur_mem) + 1;
	cur_mem += sizeof(size_t) * (md.pcl_num + 2);
	prev_pcl_ids[0] = (size_t*)cur_mem;
	cur_mem += sizeof(size_t) * md.pcl_num;
	prev_pcl_ids[1] = (size_t*)cur_mem;
	cur_mem += sizeof(size_t) * md.pcl_num;
	valid_elem_id = (size_t*)cur_mem;
	cur_mem += sizeof(size_t) * md.elem_num;
	node_has_elems[0] = ((size_t*)cur_mem) + 1;
	cur_mem += sizeof(size_t) * (md.elem_num * 3 + 2);
	node_has_elems[1] = ((size_t*)cur_mem) + 1;
	cur_mem += sizeof(size_t) * (md.elem_num * 3 + 2);
	node_elem_pairs[0] = (size_t*)cur_mem;
	cur_mem += sizeof(size_t) * md.elem_num * 3;
	node_elem_pairs[1] = (size_t*)cur_mem;
	cur_mem += sizeof(size_t) * md.elem_num * 3;
	cur_mem = cache_aligned(cur_mem);
	elem_count_bin = (size_t*)cur_mem;
	cur_mem += sizeof(size_t) * thread_num * 0x100;
	elem_sum_bin = (size_t*)cur_mem;

	pcl_in_elems[0][-1] = SIZE_MAX;
	pcl_in_elems[1][-1] = SIZE_MAX;
	node_has_elems[0][-1] = SIZE_MAX;
	node_has_elems[1][-1] = SIZE_MAX;

	prev_valid_pcl_num = md.pcl_num;
	valid_pcl_num = 0;
#pragma omp parallel
	{
		size_t my_th_id = size_t(omp_get_thread_num());

		ThreadData& thd = thread_datas[my_th_id];
		new (&thd) ThreadData;
		thd.sorted_pcl_var_id = 1;
		thd.sorted_pcl_in_elem_id = 0;
		//PclVar_T3D_ME_mt& pv_getter = thd.pcl_var_getter;
		//pv_getter.pmodel = &md;

		size_t p_id, ori_p_id, e_id;
		size_t p_id0 = Block_Low(my_th_id, thread_num, prev_valid_pcl_num);
		size_t p_id1 = Block_Low(my_th_id + 1, thread_num, prev_valid_pcl_num);
		size_t pcl_in_mesh_num = 0;
		size_t* pcl_in_elem0 = pcl_in_elems[thd.sorted_pcl_in_elem_id];
		size_t* prev_pcl_id0 = prev_pcl_ids[thd.sorted_pcl_in_elem_id];
		for (p_id = p_id0; p_id < p_id1; ++p_id)
		{
			ori_p_id = spva0.pcl_index[p_id];
			Position& p_p = pcl_pos[ori_p_id];
			Displacement& p_d = spva0.pcl_disp[p_id];
			p_p.x += p_d.ux;
			p_p.y += p_d.uy;
			p_d.ux = 0.0;
			p_d.uy = 0.0;
			ShapeFunc& p_N = spva0.pcl_N[p_id];
			e_id = md.find_pcl_in_which_elem(p_p.x, p_p.y, p_N);
			if (e_id == SIZE_MAX)
				e_id = md.find_pcl_in_which_elem_tol(p_p.x, p_p.y, p_N);
			pcl_in_elem0[p_id] = e_id;
			prev_pcl_id0[p_id] = p_id;
			if (e_id != SIZE_MAX)
				++pcl_in_mesh_num;
		}

#pragma omp critical
		valid_pcl_num += pcl_in_mesh_num;
	}

	pcl_in_elems[0][prev_valid_pcl_num] = SIZE_MAX;
	pcl_in_elems[1][prev_valid_pcl_num] = SIZE_MAX;
	valid_elem_num = 0;
	return 0;
}

int Step_T2D_ME_mt::finalize_calculation()
{
	Model_T2D_ME_mt &md = *(Model_T2D_ME_mt *)model;
	md.pcl_num = valid_pcl_num;
	for (size_t t_id = 0; t_id < thread_num; ++t_id)
		thread_datas[t_id].~ThreadData();
	return 0;
}

int substep_func_omp_T2D_ME_mt(
	void* _self,
	size_t my_th_id,
	double dt,
	double cur_time,
	size_t substp_id)
{
	typedef Model_T2D_ME_mt::SortedPclVarArrays SortedPclVarArrays;
	typedef Model_T2D_ME_mt::Force Force;
	typedef Model_T2D_ME_mt::Acceleration Acceleration;
	typedef Model_T2D_ME_mt::Velocity Velocity;
	typedef Model_T2D_ME_mt::Displacement Displacement;
	typedef Model_T2D_ME_mt::Position Position;
	typedef Model_T2D_ME_mt::Stress Stress;
	typedef Model_T2D_ME_mt::Strain Strain;
	typedef Model_T2D_ME_mt::StrainInc StrainInc;
	typedef Model_T2D_ME_mt::ElemNodeIndex ElemNodeIndex;
	typedef Model_T2D_ME_mt::ShapeFunc ShapeFunc;
	typedef Model_T2D_ME_mt::ShapeFuncAB ShapeFuncAB;
	typedef Model_T2D_ME_mt::ShapeFuncC ShapeFuncC;
	typedef Model_T2D_ME_mt::ElemNodeVM ElemNodeVM;
	typedef Model_T2D_ME_mt::NodeHasVBC NodeHasVBC;

	Step_T2D_ME_mt& self = *(Step_T2D_ME_mt*)(_self);
	Model_T2D_ME_mt& md = *(Model_T2D_ME_mt*)(self.model);

//	SortedPclVarArrays &spcv0 = self.sorted_pcl_var_array[self.pcl_sorted_var_id];
//	size_t* pcl_index0 = pscv0.pcl_index;
//	double *pcl_density0 = pscv0.pcl_density;
//	PclDisp *pcl_disp0 = pscv0.pcl_disp;
//	PclV *pcl_v0 = pscv0.pcl_v;
//	PclShapeFunc* pcl_N0 = pscv0.pcl_N;
//	PclStress* pcl_stress0 = pscv0.pcl_stress;
//	PclStrain *pcl_strain0 = pscv0.pcl_strain;
//	PclStrain* pcl_estrain0 = pscv0.pcl_estrain;
//	PclStrain* pcl_pstrain0 = pscv0.pcl_pstrain;
//
//	PclSortedVarArray& pscv1 = self.pcl_sorted_var_array[self.pcl_sorted_var_id ^ 1];
//	size_t* pcl_index1 = pscv1.pcl_index;
//	double* pcl_density1 = pscv1.pcl_density;
//	PclDisp* pcl_disp1 = pscv1.pcl_disp;
//	PclV *pcl_v1 = pscv1.pcl_v;
//	PclShapeFunc *pcl_N1 = pscv1.pcl_N;
//	PclStress *pcl_stress1 = pscv1.pcl_stress;
//	PclStrain* pcl_strain1 = pscv1.pcl_strain;
//	PclStrain* pcl_estrain1 = pscv1.pcl_estrain;
//	PclStrain* pcl_pstrain1 = pscv1.pcl_pstrain;
//
//	size_t sort_var_id = self.radix_sort_var_id;
//	size_t* new_to_prev_pcl_map = self.new_to_prev_pcl_maps[sort_var_id];
//	size_t* new_to_prev_pcl_map_tmp = self.new_to_prev_pcl_maps[sort_var_id ^ 1];
//	size_t* pcl_in_elem_array = self.pcl_in_elem_arrays[sort_var_id];
//	size_t* pcl_in_elem_array_tmp = self.pcl_in_elem_arrays[sort_var_id ^ 1];
//
//	size_t p_id0 = self.pcl_range[my_th_id].id;
//	size_t p_id1 = self.pcl_range[my_th_id + 1].id;
//	size_t p_id, e_id;
//	size_t prev_pcl_id, ori_pcl_id;
//	double p_m, p_vol, p_N_m;
//	double one_third_bfx, one_third_bfy;
//	for (p_id = p_id0; p_id < p_id1; ++p_id)
//	{
//		prev_pcl_id = new_to_prev_pcl_map[p_id];
//
//		ori_pcl_id = pcl_index1[prev_pcl_id];
//		pcl_index0[p_id] = ori_pcl_id;
//
//		e_id = pcl_in_elem_array[p_id];
//
//		// map pcl mass
//		p_m = self.pcl_m[ori_pcl_id];
//		self.elem_pcl_m[e_id] += p_m;
//
//		// map pcl volume
//		p_vol = p_m / pcl_density1[prev_pcl_id];
//		self.pcl_vol[p_id] = p_vol;
//		self.elem_pcl_vol[e_id] += p_vol;
//
//		// map stress
//		PclStress& p_s1 = pcl_stress1[prev_pcl_id];
//		PclStress& p_s0 = pcl_stress0[p_id];
//		p_s0.s11 = p_s1.s11;
//		p_s0.s22 = p_s1.s22;
//		p_s0.s12 = p_s1.s12;
//		ElemStress & e_s = self.elem_stress[e_id];
//		e_s.s11 += p_s0.s11 * p_vol;
//		e_s.s22 += p_s0.s22 * p_vol;
//		e_s.s12 += p_s0.s12 * p_vol;
//
//		// map velocity
//		PclShapeFunc& p_N1 = pcl_N1[prev_pcl_id];
//		PclShapeFunc& p_N0 = pcl_N0[p_id];
//		p_N0.N1 = p_N1.N1;
//		p_N0.N2 = p_N1.N2;
//		p_N0.N3 = p_N1.N3;
//		PclV& p_v1 = pcl_v1[prev_pcl_id];
//		PclV& p_v0 = pcl_v0[p_id];
//		p_v0.vx = p_v1.vx;
//		p_v0.vy = p_v1.vy;
//		ElemNodeVM& en_vm1 = self.elem_node_vm[e_id * 3];
//		p_N_m = p_N0.N1 * p_m;
//		en_vm1.vm += p_N_m;
//		en_vm1.vmx += p_N_m * p_v0.vx;
//		en_vm1.vmy += p_N_m * p_v0.vy;
//		ElemNodeVM& en_vm2 = self.elem_node_vm[e_id * 3 + 1];
//		p_N_m = p_N0.N2 * p_m;
//		en_vm2.vm += p_N_m;
//		en_vm2.vmx += p_N_m * p_v0.vx;
//		en_vm2.vmy += p_N_m * p_v0.vy;
//		ElemNodeVM& en_vm3 = self.elem_node_vm[e_id * 3 + 2];
//		p_N_m = p_N0.N3 * p_m;
//		en_vm3.vm += p_N_m;
//		en_vm3.vmx += p_N_m * p_v0.vx;
//		en_vm3.vmy += p_N_m * p_v0.vy;
//
//		// External load
//		PclBodyForce& p_bf = self.pcl_bf[ori_pcl_id];
//		one_third_bfx = one_third * p_bf.bfx;
//		one_third_bfy = one_third * p_bf.bfy;
//		PclTraction& p_t = self.pcl_t[ori_pcl_id];
//		ElemNodeForce& en_f1 = self.elem_node_force[e_id * 3];
//		en_f1.fx += one_third_bfx + p_N0.N1 * p_t.tx;
//		en_f1.fy += one_third_bfy + p_N0.N1 * p_t.ty;
//		ElemNodeForce& en_f2 = self.elem_node_force[e_id * 3 + 1];
//		en_f2.fx += one_third_bfx + p_N0.N2 * p_t.tx;
//		en_f2.fy += one_third_bfy + p_N0.N2 * p_t.ty;
//		ElemNodeForce& en_f3 = self.elem_node_force[e_id * 3 + 2];
//		en_f3.fx += one_third_bfx + p_N0.N3 * p_t.tx;
//		en_f3.fy += one_third_bfy + p_N0.N3 * p_t.ty;
//	}
//
//	if (md.has_rigid_rect())
//	{
//		RigidRect& rr = md.get_rigid_rect();
//		RigidRectForce rr_force;
//		rr_force.reset_f_contact();
//		self.apply_rigid_rect_avg(
//			my_th_id, dt,
//			pcl_in_elem_array, pscv0,
//			rr_force
//			);
//#pragma omp critical
//		{
//			rr.combine(rr_force);
//		}
//	}
//#pragma omp barrier
//
//#pragma omp master
//	{
//		if (md.has_rigid_rect())
//		{
//			RigidRect& rr = md.get_rigid_rect();
//			rr.update_motion(dt);
//		}
//	}
//
//	size_t e_id0 = self.elem_range[my_th_id];
//	size_t e_id1 = self.elem_range[my_th_id + 1];
//	double e_pcl_vol;
//	for (e_id = e_id0; e_id < e_id1; ++e_id)
//	{
//		if (self.elem_pcl_vol[e_id] != 0.0)
//		{
//			e_pcl_vol = self.elem_pcl_vol[e_id];
//			self.elem_density[e_id] = self.elem_pcl_m[e_id] / e_pcl_vol;
//
//			ElemStress& e_s = self.elem_stress[e_id];
//			e_s.s11 /= e_pcl_vol;
//			e_s.s22 /= e_pcl_vol;
//			e_s.s12 /= e_pcl_vol;
//			if (e_pcl_vol > self.elem_area[e_id])
//				e_pcl_vol = self.elem_area[e_id];
//			ElemShapeFuncAB& e_sf = self.elem_sf_ab[e_id];
//			ElemNodeForce& en_f1 = self.elem_node_force[e_id * 3];
//			en_f1.fx -= (e_sf.dN1_dx * e_s.s11 + e_sf.dN1_dy * e_s.s12) * e_pcl_vol;
//			en_f1.fy -= (e_sf.dN1_dx * e_s.s12 + e_sf.dN1_dy * e_s.s22) * e_pcl_vol;
//			ElemNodeForce& en_f2 = self.elem_node_force[e_id * 3 + 1];
//			en_f2.fx -= (e_sf.dN2_dx * e_s.s11 + e_sf.dN2_dy * e_s.s12) * e_pcl_vol;
//			en_f2.fy -= (e_sf.dN2_dx * e_s.s12 + e_sf.dN2_dy * e_s.s22) * e_pcl_vol;
//			ElemNodeForce& en_f3 = self.elem_node_force[e_id * 3 + 2];
//			en_f3.fx -= (e_sf.dN3_dx * e_s.s11 + e_sf.dN3_dy * e_s.s12) * e_pcl_vol;
//			en_f3.fy -= (e_sf.dN3_dx * e_s.s12 + e_sf.dN3_dy * e_s.s22) * e_pcl_vol;
//		}
//	}
//#pragma omp barrier
//
//	// update node variables
//	size_t n_id0 = self.node_range[my_th_id];
//	size_t n_id1 = self.node_range[my_th_id + 1];
//	size_t ne_id = self.node_elem_range[my_th_id];
//	size_t n_id, ne_id1, node_var_id, bc_mask;
//	double n_am, n_fx, n_fy, n_vm, n_vmx, n_vmy;
//	for (n_id = n_id0; n_id < n_id1; ++n_id)
//	{
//		n_am = 0.0;
//		n_fx = 0.0;
//		n_fy = 0.0;
//		n_vm = 0.0;
//		n_vmx = 0.0;
//		n_vmy = 0.0;
//		ne_id1 = self.node_elem_list[n_id];
//		for (; ne_id < ne_id1; ++ne_id)
//		{
//			n_am += self.elem_pcl_m[self.elem_id_array[ne_id]];
//			node_var_id = self.node_elem_id_array[ne_id];
//			ElemNodeForce& nf = self.elem_node_force[node_var_id];
//			n_fx += nf.fx;
//			n_fy += nf.fy;
//			ElemNodeVM& nvm = self.elem_node_vm[node_var_id];
//			n_vm += nvm.vm;
//			n_vmx += nvm.vmx;
//			n_vmy += nvm.vmy;
//		}
//		NodeA &node_a = self.node_a[n_id];
//		if (n_am != 0.0)
//		{
//			n_am *= one_third;
//			self.node_am[n_id] = n_am;
//			node_a.ax = n_fx / n_am;
//			node_a.ay = n_fy / n_am;
//		}
//		if (n_vm != 0.0)
//		{
//			NodeV& node_v = self.node_v[n_id];
//			node_v.vx = n_vmx / n_vm + node_a.ax * dt;
//			node_v.vy = n_vmy / n_vm + node_a.ay * dt;
//			NodeHasVBC& node_has_vbc = self.node_has_vbc[n_id];
//			bc_mask = size_t(node_has_vbc.has_vx_bc) + SIZE_MAX;
//			node_a.ax_ui &= bc_mask;
//			node_v.vx_ui &= bc_mask;
//			bc_mask = size_t(node_has_vbc.has_vy_bc) + SIZE_MAX;
//			node_a.ay_ui &= bc_mask;
//			node_v.vy_ui &= bc_mask;
//		}
//	}
//#pragma omp barrier
//
//	// cal element strain and strain enhancement
//	double e_de_vol;
//	for (e_id = e_id0; e_id < e_id1; ++e_id)
//	{
//		if (self.elem_pcl_vol[e_id] != 0.0)
//		{
//			ElemNodeIndex& e_n_id = self.elem_node_id[e_id];
//			NodeV& n_v1 = self.node_v[e_n_id.n1];
//			NodeV& n_v2 = self.node_v[e_n_id.n2];
//			NodeV& n_v3 = self.node_v[e_n_id.n3];
//			ElemShapeFuncAB& e_sf = self.elem_sf_ab[e_id];
//			ElemStrainInc& e_de = self.elem_de[e_id];
//			e_de.de11 = (e_sf.dN1_dx * n_v1.vx + e_sf.dN2_dx * n_v2.vx + e_sf.dN3_dx * n_v3.vx) * dt;
//			e_de.de22 = (e_sf.dN1_dy * n_v1.vy + e_sf.dN2_dy * n_v2.vy + e_sf.dN3_dy * n_v3.vy) * dt;
//			e_de.de12 = (e_sf.dN1_dx * n_v1.vy + e_sf.dN2_dx * n_v2.vy + e_sf.dN3_dx * n_v3.vy
//					   + e_sf.dN1_dy * n_v1.vx + e_sf.dN2_dy * n_v2.vx + e_sf.dN3_dy * n_v3.vx) * dt * 0.5;
//			e_de_vol = e_de.de11 + e_de.de22;
//			self.elem_m_de_vol[e_id] = self.elem_pcl_m[e_id] * e_de_vol;
//			e_de_vol *= one_third;
//			e_de.de11 -= e_de_vol;
//			e_de.de22 -= e_de_vol;
//		}
//	}
//#pragma omp barrier
//
//	double n_am_de_vol;
//	ne_id = self.node_elem_range[my_th_id];
//	for (n_id = n_id0; n_id < n_id1; ++n_id)
//	{
//		if (self.node_am[n_id] != 0.0)
//		{
//			ne_id1 = self.node_elem_list[n_id];
//			n_am_de_vol = 0.0;
//			for (; ne_id < ne_id1; ++ne_id)
//				n_am_de_vol += self.elem_m_de_vol[self.elem_id_array[ne_id]];
//			self.node_de_vol[n_id] = n_am_de_vol * one_third / self.node_am[n_id];
//		}
//		else
//		{
//			ne_id = self.node_elem_list[n_id];
//			self.node_de_vol[n_id] = 0.0;
//		}
//	}
//#pragma omp barrier
//
//	for (e_id = e_id0; e_id < e_id1; ++e_id)
//	{
//		if (self.elem_pcl_vol[e_id] != 0.0)
//		{
//			ElemNodeIndex& e_n_id = self.elem_node_id[e_id];
//			e_de_vol = one_third * 
//				 (self.node_de_vol[e_n_id.n1]
//				+ self.node_de_vol[e_n_id.n2]
//				+ self.node_de_vol[e_n_id.n3]);
//			self.elem_density[e_id] /= (1.0 + e_de_vol);
//			e_de_vol *= one_third;
//			ElemStrainInc& e_de = self.elem_de[e_id];
//			e_de.de11 += e_de_vol;
//			e_de.de22 += e_de_vol;
//		}
//	}
//
//#pragma omp master
//	{
//		self.new_pcl_num = 0;
//	}
//#pragma omp barrier
//
//	// update particle variables
//	double pcl_x, pcl_y;
//	double dstrain[6] = { 0.0 };
//	size_t pcl_in_mesh_num = 0;
//	for (p_id = p_id0; p_id < p_id1; ++p_id)
//	{
//		e_id = pcl_in_elem_array[p_id];
//
//		// update density
//		pcl_density0[p_id] = self.elem_density[e_id];
//
//		prev_pcl_id = new_to_prev_pcl_map[p_id];
//
//		// update stress
//		ori_pcl_id = pcl_index0[p_id];
//		ElemStrainInc& e_de = self.elem_de[e_id];
//		PclStrain& pcl_e0 = pcl_strain0[p_id];
//		PclStrain& pcl_e1 = pcl_strain1[prev_pcl_id];
//		pcl_e0.e11 = pcl_e1.e11 + e_de.de11;
//		pcl_e0.e22 = pcl_e1.e22 + e_de.de22;
//		pcl_e0.e12 = pcl_e1.e12 + e_de.de12;
//		dstrain[0] = e_de.de11;
//		dstrain[1] = e_de.de22;
//		dstrain[3] = e_de.de12;
//		MatModel::MaterialModel &pcl_mm = *self.pcl_mat_model[ori_pcl_id];
//		int32_t mm_res = pcl_mm.integrate(dstrain);
//		const double *dstress = pcl_mm.get_dstress();
//		PclStress& p_s0 = pcl_stress0[p_id];
//		p_s0.s11 += dstress[0];
//		p_s0.s22 += dstress[1];
//		p_s0.s12 += dstress[3];
//		const double* dee = pcl_mm.get_dstrain_e();
//		PclStrain &pcl_ee0 = pcl_estrain0[p_id];
//		PclStrain& pcl_ee1 = pcl_estrain1[prev_pcl_id];
//		pcl_ee0.e11 = pcl_ee1.e11 + dee[0];
//		pcl_ee0.e22 = pcl_ee1.e22 + dee[1];
//		pcl_ee0.e12 = pcl_ee1.e12 + dee[3];
//		const double* dpe = pcl_mm.get_dstrain_p();
//		PclStrain& pcl_pe0 = pcl_pstrain0[p_id];
//		PclStrain& pcl_pe1 = pcl_pstrain1[prev_pcl_id];
//		pcl_pe0.e11 = pcl_pe1.e11 + dpe[0];
//		pcl_pe0.e22 = pcl_pe1.e22 + dpe[1];
//		pcl_pe0.e12 = pcl_pe1.e12 + dpe[3];
//
//		ElemNodeIndex& e_n_id = self.elem_node_id[e_id];
//		PclShapeFunc& p_N = pcl_N0[p_id];
//		
//		// update velocity
//		NodeA& n_a1 = self.node_a[e_n_id.n1];
//		NodeA& n_a2 = self.node_a[e_n_id.n2];
//		NodeA& n_a3 = self.node_a[e_n_id.n3];
//		PclV& p_v0 = pcl_v0[p_id];
//		p_v0.vx += (p_N.N1 * n_a1.ax + p_N.N2 * n_a2.ax + p_N.N3 * n_a3.ax) * dt;
//		p_v0.vy += (p_N.N1 * n_a1.ay + p_N.N2 * n_a2.ay + p_N.N3 * n_a3.ay) * dt;
//		
//		// update displacement
//		NodeV& n_v1 = self.node_v[e_n_id.n1];
//		NodeV& n_v2 = self.node_v[e_n_id.n2];
//		NodeV& n_v3 = self.node_v[e_n_id.n3];
//		PclDisp& p_d1 = pcl_disp1[new_to_prev_pcl_map[p_id]];
//		PclDisp& p_d0 = pcl_disp0[p_id];
//		p_d0.ux = p_d1.ux + (p_N.N1 * n_v1.vx + p_N.N2 * n_v2.vx + p_N.N3 * n_v3.vx) * dt;
//		p_d0.uy = p_d1.uy + (p_N.N1 * n_v1.vy + p_N.N2 * n_v2.vy + p_N.N3 * n_v3.vy) * dt;
//		
//		// update location (in which element)
//		PclPos& p_p = self.pcl_pos[ori_pcl_id];
//		pcl_x = p_p.x + p_d0.ux;
//		pcl_y = p_p.y + p_d0.uy;
//		if (!md.is_in_element(pcl_x, pcl_y, e_id, p_N))
//			e_id = md.find_pcl_in_which_elem(pcl_x, pcl_y, p_N);
//		if (e_id != self.elem_num)
//		{
//			if (p_N.N1 < N_min)
//				p_N.N1 = N_min;
//			if (p_N.N2 < N_min)
//				p_N.N2 = N_min;
//			if (p_N.N3 < N_min)
//				p_N.N3 = N_min;
//			++pcl_in_mesh_num;
//		}
//		new_to_prev_pcl_map[p_id] = p_id;
//		pcl_in_elem_array[p_id] = e_id;
//	}
//
//#pragma omp critical
//	{
//		self.new_pcl_num += pcl_in_mesh_num;
//	}
//#pragma omp barrier
//
//	// sort particle variables
//	size_t* my_cbin = self.elem_count_bin + my_th_id * 0x100;
//	size_t* my_sbin = self.elem_sum_bin + my_th_id * 0x100;
//	size_t data_digit, bin_id, th_id;
//	size_t* other_cbin;
//	for (size_t digit_disp = 0, elem_num_tmp = self.elem_num;
//		 elem_num_tmp; digit_disp += 8, elem_num_tmp >>= 8)
//	{
//		memset(my_cbin, 0, 0x100 * sizeof(size_t));
//
//		for (p_id = p_id0; p_id < p_id1; ++p_id)
//		{
//			data_digit = (pcl_in_elem_array[p_id] >> digit_disp) & 0xFF;
//			++my_cbin[data_digit];
//		}
//
//		my_sbin[0] = my_cbin[0];
//		for (bin_id = 1; bin_id < 0x100; ++bin_id)
//		{
//			my_cbin[bin_id] += my_cbin[bin_id - 1];
//			my_sbin[bin_id] = my_cbin[bin_id];
//		}
//#pragma omp barrier
//
//		for (th_id = 0; th_id < my_th_id; ++th_id)
//		{
//			other_cbin = self.elem_count_bin + th_id * 0x100;
//			for (bin_id = 0; bin_id < 0x100; ++bin_id)
//				my_sbin[bin_id] += other_cbin[bin_id];
//		}
//		for (th_id = my_th_id + 1; th_id < self.thread_num; ++th_id)
//		{
//			other_cbin = self.elem_count_bin + th_id * 0x100;
//			for (bin_id = 1; bin_id < 0x100; ++bin_id)
//				my_sbin[bin_id] += other_cbin[bin_id - 1];
//		}
//
//		for (p_id = p_id1; p_id-- > p_id0;)
//		{
//			data_digit = (pcl_in_elem_array[p_id] >> digit_disp) & 0xFF;
//			pcl_in_elem_array_tmp[--my_sbin[data_digit]] = pcl_in_elem_array[p_id];
//			new_to_prev_pcl_map_tmp[my_sbin[data_digit]] = new_to_prev_pcl_map[p_id];
//		}
//
//		new_to_prev_pcl_map_tmp = new_to_prev_pcl_map;
//		pcl_in_elem_array_tmp = pcl_in_elem_array;
//		sort_var_id ^= 1;
//		new_to_prev_pcl_map = self.new_to_prev_pcl_maps[sort_var_id];
//		pcl_in_elem_array = self.pcl_in_elem_arrays[sort_var_id];
//#pragma omp barrier
//	}
//	
//	// reset element variables
//	e_id0 = self.elem_range[my_th_id];
//	e_id1 = self.elem_range[my_th_id + 1];
//	memset(self.elem_pcl_m + e_id0, 0, (e_id1 - e_id0) * sizeof(double));
//	memset(self.elem_pcl_vol + e_id0, 0, (e_id1 - e_id0) * sizeof(double));
//	memset(self.elem_stress + e_id0, 0, (e_id1 - e_id0) * sizeof(ElemStress));
//	memset(self.elem_node_vm + e_id0 * 3, 0, (e_id1 - e_id0) * 3 * sizeof(ElemNodeVM));
//	memset(self.elem_node_force + e_id0 * 3, 0, (e_id1 - e_id0) * 3 * sizeof(ElemNodeForce));
//	memset(self.elem_m_de_vol + e_id0 * 3, 0, (e_id1 - e_id0) * 3 * sizeof(double));
//
//	if (my_th_id == 0)
//	{
//		self.pcl_range[self.thread_num].id = self.new_pcl_num;
//		self.radix_sort_var_id = sort_var_id;
//		self.pcl_sorted_var_id ^= 1;
//		self.pcl_num = self.new_pcl_num;
//		
//		if (md.has_rigid_rect())
//		{
//			RigidRect &rr = md.get_rigid_rect();
//			self.rr_fx_cont = rr.get_fx_contact();
//			self.rr_fy_cont = rr.get_fy_contact();
//			self.rr_m_cont = rr.get_m_contact();
//			rr.reset_f_contact();
//		}
//
//		if (self.new_pcl_num)
//			self.continue_calculation();
//		else
//			self.exit_calculation();
//	}
//	else
//	{
//		p_id = Block_Low(my_th_id, self.thread_num, self.new_pcl_num);
//		if (p_id < self.new_pcl_num)
//		{
//			e_id = pcl_in_elem_array[p_id];
//			while (e_id == pcl_in_elem_array[p_id + 1])
//				++p_id;
//			self.pcl_range[my_th_id].id = p_id + 1;
//		}
//		else
//		{
//			self.pcl_range[my_th_id].id = self.new_pcl_num;
//		}
//	}
//
//#pragma omp barrier
//	return 0;
//}
//
//int Step_T2D_ME_mt::apply_rigid_rect_avg(
//	size_t my_th_id,
//	double dt,
//	size_t *pcl_in_elem_array,
//	PclSortedVarArray &cur_pscv,
//	RigidRectForce &rr_force
//	)
//{
//	double p_x, p_y;
//	double dist, norm_x, norm_y;
//	double f_cont, fx_cont, fy_cont;
//	size_t e_id, pcl_ori_id;
//	size_t p_id0 = pcl_range[my_th_id].id;
//	size_t p_id1 = pcl_range[my_th_id + 1].id;
//	Model_T2D_ME_mt& md = *(Model_T2D_ME_mt*)model;
//	RigidRect& rr = md.get_rigid_rect();
//	const Point2D &rr_centre = rr.get_centre();
//	for (size_t p_id = p_id0; p_id < p_id1; ++p_id)
//	{
//		pcl_ori_id = cur_pscv.pcl_index[p_id];
//		PclPos& p_p = pcl_pos[pcl_ori_id];
//		PclDisp& p_d = cur_pscv.pcl_disp[p_id];
//		p_x = p_p.x + p_d.ux;
//		p_y = p_p.y + p_d.uy;
//		if (rr.detect_collision_with_point(
//			p_x, p_y, pcl_vol[p_id],
//			dist, norm_x, norm_y))
//		{
//			f_cont = K_cont * dist;
//			fx_cont = f_cont * norm_x;
//			fy_cont = f_cont * norm_y;
//			// apply contact force to rigid body
//			rr_force.add_f_contact(
//				p_x, p_y,
//				-fx_cont, -fy_cont,
//				rr_centre.x, rr_centre.y
//				);
//			// apply contact force to mesh
//			PclShapeFunc &p_N = cur_pscv.pcl_N[p_id];
//			e_id = pcl_in_elem_array[p_id];
//			ElemNodeForce& en_f1 = elem_node_force[e_id * 3];
//			en_f1.fx += p_N.N1 * fx_cont;
//			en_f1.fy += p_N.N1 * fy_cont;
//			ElemNodeForce& en_f2 = elem_node_force[e_id * 3 + 1];
//			en_f2.fx += p_N.N2 * fx_cont;
//			en_f2.fy += p_N.N2 * fy_cont;
//			ElemNodeForce& en_f3 = elem_node_force[e_id * 3 + 2];
//			en_f3.fx += p_N.N3 * fx_cont;
//			en_f3.fy += p_N.N3 * fy_cont;
//		}
//	}

	return 0;
}
