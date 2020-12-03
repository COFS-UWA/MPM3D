#include "SimulationsOMP_pcp.h"

#include <fstream>

#include <omp.h>

#include "Step_T2D_CHM_mt.h"

#define one_third (1.0/3.0)
#define N_min (1.0e-10)
#define Block_Low(th_id, th_num, data_num) ((th_id)*(data_num)/(th_num))

#ifdef _DEBUG
	static std::fstream res_file_t2d_me_mt;
#endif // _DEBUG

Step_T2D_CHM_mt::Step_T2D_CHM_mt(const char* _name) : 
	Step_OMP(_name, "Step_T2D_CHM_mt", &substep_func_omp_T2D_CHM_mt) {}

Step_T2D_CHM_mt::~Step_T2D_CHM_mt() {}

int Step_T2D_CHM_mt::init_calculation()
{
#ifdef _DEBUG
	res_file_t2d_me_mt.open("t2d_chm_mt_res.txt", std::ios::out | std::ios::binary);
#endif

	Model_T2D_CHM_mt &md = *(Model_T2D_CHM_mt *)model;

	omp_set_num_threads(thread_num);

//	const size_t th_num_1 = thread_num + 1;
//	char *mem_range = (char *)task_range_mem.alloc((sizeof(PclRange) + sizeof(size_t) * 3) * th_num_1);
//	pcl_range = (PclRange *)mem_range;
//	mem_range += sizeof(PclRange) * th_num_1;
//	elem_range = (size_t *)mem_range;
//	mem_range += sizeof(size_t) * th_num_1;
//	node_elem_range = (size_t *)mem_range;
//	mem_range += sizeof(size_t) * th_num_1;
//	node_range = (size_t*)mem_range;
//
//	pcl_range[0].id = 0;
//	divide_task_to_thread(thread_num, md.elem_num, elem_range);
//	divide_task_to_thread(thread_num, md.node_num, node_range);
//	node_elem_range[0] = 0;
//	for (size_t th_id = 0; th_id < thread_num; ++th_id)
//		node_elem_range[th_id + 1] = md.node_elem_list[node_range[th_id + 1] - 1];
//
//	radix_sort_var_id = 0;
//	new_to_prev_pcl_maps[0] = (size_t *)radix_sort_var_mem.alloc(sizeof(size_t) * (md.pcl_num * 4 + 2) + Cache_Alignment * 3);
//	new_to_prev_pcl_maps[1] = cache_aligned(new_to_prev_pcl_maps[0] + md.pcl_num);
//	pcl_in_elem_arrays[0] = cache_aligned(new_to_prev_pcl_maps[1] + md.pcl_num);
//	pcl_in_elem_arrays[1] = cache_aligned(pcl_in_elem_arrays[0] + md.pcl_num + 1);
//	pcl_in_elem_arrays[0][md.pcl_num] = md.elem_num;
//	pcl_in_elem_arrays[1][md.pcl_num] = md.elem_num;
//
//	elem_count_bin = (size_t *)elem_bin_mem.alloc(sizeof(size_t) * thread_num * 0x100 * 2);
//	elem_sum_bin = elem_count_bin + thread_num * 0x100;
//
//	pcl_num = md.pcl_num;
//	elem_num = md.elem_num;
//	node_num = md.node_num;
//
//	pcl_m = md.pcl_m;
//	pcl_bf = md.pcl_bf;
//	pcl_t = md.pcl_t;
//	pcl_pos = md.pcl_pos;
//	pcl_vol = md.pcl_vol;
//	pcl_mat_model = md.pcl_mat_model;
//
//	pcl_sorted_var_id = 1;
//
//	PclSortedVarArray &md_pscv0 = md.pcl_sorted_var_array[0];
//	PclSortedVarArray &pscv0 = pcl_sorted_var_array[0];
//	pscv0.pcl_index = md_pscv0.pcl_index;
//	pscv0.pcl_density = md_pscv0.pcl_density;
//	pscv0.pcl_disp = md_pscv0.pcl_disp;
//	pscv0.pcl_v = md_pscv0.pcl_v;
//	pscv0.pcl_N = md_pscv0.pcl_N;
//	pscv0.pcl_stress = md_pscv0.pcl_stress;
//
//	PclSortedVarArray& md_pscv1 = md.pcl_sorted_var_array[1];
//	PclSortedVarArray& pscv1 = pcl_sorted_var_array[1];
//	pscv1.pcl_index = md_pscv1.pcl_index;
//	pscv1.pcl_density = md_pscv1.pcl_density;
//	pscv1.pcl_disp = md_pscv1.pcl_disp;
//	pscv1.pcl_v = md_pscv1.pcl_v;
//	pscv1.pcl_N = md_pscv1.pcl_N;
//	pscv1.pcl_stress = md_pscv1.pcl_stress;
//
//	elem_node_id = md.elem_node_id;
//	elem_area = md.elem_area;
//	elem_sf_ab = md.elem_sf_ab;
//	elem_sf_c = md.elem_sf_c;
//
//	elem_density = md.elem_density;
//	elem_pcl_m = md.elem_pcl_m;
//	elem_pcl_vol = md.elem_pcl_vol;
//	elem_de = md.elem_de;
//	elem_stress = md.elem_stress;
//	elem_m_de_vol = md.elem_m_de_vol;
//
//	elem_node_vm = md.elem_node_vm;
//	elem_node_force = md.elem_node_force;
//
//	elem_id_array = md.elem_id_array;
//	node_elem_id_array = md.node_elem_id_array;
//	node_elem_list = md.node_elem_list;
//	node_a = md.node_a;
//	node_v = md.node_v;
//	node_has_vbc = md.node_has_vbc;
//	node_am = md.node_am;
//	node_de_vol = md.node_de_vol;
//	
//	if (md.has_rigid_rect())
//	{
//		RigidRect &rr = md.get_rigid_rect();
//		rr_fx_cont = 0.0;
//		rr_fy_cont = 0.0;
//		rr_m_cont = 0.0;
//		rr.reset_f_contact();
//	}
//
//	for (size_t th_id = 1; th_id < th_num_1; ++th_id)
//		pcl_range[th_id].id = Block_Low(th_id, thread_num, pcl_num);
//	PclSortedVarArray &psva = md.pcl_sorted_var_array[0];
//	PclSortedVarArray& psva1 = md.pcl_sorted_var_array[1];
//	size_t next_pcl_num = 0;
//#pragma omp parallel
//	{
//		size_t my_th_id = size_t(omp_get_thread_num());
//		
//		// init elem vars
//		size_t e_id0 = elem_range[my_th_id];
//		size_t e_id1 = elem_range[my_th_id + 1];
//		memset(elem_pcl_m + e_id0, 0, (e_id1 - e_id0) * sizeof(double));
//		memset(elem_pcl_vol + e_id0, 0, (e_id1 - e_id0) * sizeof(double));
//		memset(elem_stress + e_id0, 0, (e_id1 - e_id0) * sizeof(ElemStress));
//		memset(elem_node_vm + e_id0 * 3, 0, (e_id1 - e_id0) * 3 * sizeof(ElemNodeVM));
//		memset(elem_node_force + e_id0 * 3, 0, (e_id1 - e_id0) * 3 * sizeof(ElemNodeForce));
//		memset(elem_m_de_vol + e_id0 * 3, 0, (e_id1 - e_id0) * 3 * sizeof(double));
//
//		size_t sort_var_id = radix_sort_var_id;
//		size_t* new_to_prev_pcl_map = new_to_prev_pcl_maps[sort_var_id];
//		size_t* new_to_prev_pcl_map_tmp = new_to_prev_pcl_maps[sort_var_id ^ 1];
//		size_t* pcl_in_elem_array = pcl_in_elem_arrays[sort_var_id];
//		size_t* pcl_in_elem_array_tmp = pcl_in_elem_arrays[sort_var_id ^ 1];
//
//		size_t p_id0 = pcl_range[my_th_id].id;
//		size_t p_id1 = pcl_range[my_th_id + 1].id;
//		size_t pcl_in_mesh_num = 0;
//		size_t e_id, p_id;
//		for (p_id = p_id0; p_id < p_id1; ++p_id)
//		{
//			PclDisp& p_d = psva.pcl_disp[p_id];
//			p_d.ux = 0.0;
//			p_d.uy = 0.0;
//			PclDisp& p_d1 = psva1.pcl_disp[p_id];
//			p_d1.ux = 0.0;
//			p_d1.uy = 0.0;
//
//			PclPos& p_p = md.pcl_pos[p_id];
//			PclShapeFunc& p_N = psva.pcl_N[p_id];
//			e_id = md.find_pcl_in_which_elem(p_p.x, p_p.y, p_N);
//			if (e_id != elem_num)
//			{
//				if (p_N.N1 < N_min)
//					p_N.N1 = N_min;
//				if (p_N.N2 < N_min)
//					p_N.N2 = N_min;
//				if (p_N.N3 < N_min)
//					p_N.N3 = N_min;
//				++pcl_in_mesh_num;
//			}
//			new_to_prev_pcl_map[p_id] = p_id;
//			pcl_in_elem_array[p_id] = e_id;
//		}
//
//#pragma omp critical
//		next_pcl_num += pcl_in_mesh_num;

//#pragma omp barrier
//		}
//
//		if (my_th_id == 0) // master thread
//		{
//			pcl_range[thread_num].id = next_pcl_num;
//			radix_sort_var_id = sort_var_id;
//			pcl_num = next_pcl_num;
//		}
//		else // non master thread
//		{
//			p_id = Block_Low(my_th_id, thread_num, next_pcl_num);
//			if (p_id < next_pcl_num)
//			{
//				e_id = pcl_in_elem_array[p_id];
//				while (e_id == pcl_in_elem_array[p_id + 1])
//					++p_id;
//				pcl_range[my_th_id].id = p_id + 1;
//			}
//			else
//			{
//				pcl_range[my_th_id].id = next_pcl_num;
//			}
//		}
//	}
	
	return 0;
}

int Step_T2D_CHM_mt::finalize_calculation()
{
	Model_T2D_CHM_mt &md = *(Model_T2D_CHM_mt *)model;
	md.pcl_num = pcl_num;
	return 0;
}

int substep_func_omp_T2D_CHM_mt(
	void* _self,
	size_t my_th_id,
	double dt,
	double cur_time,
	size_t substp_id
	)
{
	typedef Model_T2D_CHM_mt::ShapeFunc ShapeFunc;
	typedef Model_T2D_CHM_mt::DShapeFuncAB DShapeFuncAB;
	typedef Model_T2D_CHM_mt::DShapeFuncC DShapeFuncC;
	typedef Model_T2D_CHM_mt::Force Force;
	typedef Model_T2D_CHM_mt::Position Position;
	typedef Model_T2D_CHM_mt::Displacement Displacement;
	typedef Model_T2D_CHM_mt::Velocity Velocity;
	typedef Model_T2D_CHM_mt::Acceleration Acceleration;
	typedef Model_T2D_CHM_mt::Stress Stress;
	typedef Model_T2D_CHM_mt::Strain Strain;
	typedef Model_T2D_CHM_mt::StrainInc StrainInc;
	typedef Model_T2D_CHM_mt::SortedPclVarArrays SortedPclVarArrays;
	typedef Model_T2D_CHM_mt::ElemNodeIndex ElemNodeIndex;
	typedef Model_T2D_CHM_mt::ElemNodeVM ElemNodeVM;
	typedef Model_T2D_CHM_mt::NodeHasVBC NodeHasVBC;
	typedef Step_T2D_CHM_mt::ThreadData ThreadData;

	Step_T2D_CHM_mt& self = *(Step_T2D_CHM_mt*)(_self);
	
	if (self.valid_pcl_num == 0)
	{
#pragma omp master
		self.abort_calculation();

#pragma omp barrier
		return 0;
	}
	
	Model_T2D_CHM_mt& md = *(Model_T2D_CHM_mt*)(self.model);

	const double* const pcl_m_s = self.pcl_m_s;
	const double* pcl_density_s = self.pcl_density_s;
	const Force* const pcl_bf = self.pcl_bf;
	const Force* const pcl_t = self.pcl_t;
	const Position* const pcl_pos = self.pcl_pos;
	double* const pcl_vol = self.pcl_vol;
	MatModel::MaterialModel** const pcl_mat_model = self.pcl_mat_model;
	
//	const size_t thread_num = self.thread_num;
//	ThreadData &thd = self.thread_datas[my_th_id];
//	SortedPclVarArrays &spva0 = self.sorted_pcl_var_arrays[thd.sorted_pcl_var_id];
//	thd.sorted_pcl_var_id ^= 1;
//	SortedPclVarArrays& spva1 = self.sorted_pcl_var_arrays[thd.sorted_pcl_var_id];
//	size_t* pcl_index0 = spva0.pcl_index;
//	double* pcl_n0 = spva0.pcl_n;
//	double* pcl_density_f0 = spva0.pcl_density_f;
//	Velocity* pcl_v0 = spva0.pcl_v; 
//	Displacement* pcl_disp0 = spva0.pcl_disp;
//	Stress* pcl_stress0 = spva0.pcl_stress;
//	Strain* pcl_strain0 = spva0.pcl_strain;
//	Strain* pcl_estrain0 = spva0.pcl_estrain;
//	Strain* pcl_pstrain0 = spva0.pcl_pstrain;
//	ShapeFunc* pcl_N0 = spva0.pcl_N;
//
//	size_t* pcl_index1 = spva1.pcl_index;
//	double* pcl_n1 = spva1.pcl_n;
//	double* pcl_density_f1 = spva1.pcl_density_f;
//	Velocity* pcl_v1 = spva1.pcl_v;
//	Displacement* pcl_disp1 = spva1.pcl_disp;
//	Stress* pcl_stress1 = spva1.pcl_stress;
//	Strain* pcl_strain1 = spva1.pcl_strain;
//	Strain* pcl_estrain1 = spva1.pcl_estrain;
//	Strain* pcl_pstrain1 = spva1.pcl_pstrain;
//	ShapeFunc* pcl_N1 = spva1.pcl_N;
//
//	ElemNodeIndex* elem_node_id;
//	DShapeFuncAB* elem_N_ab;
//	DShapeFuncC* elem_N_c;
//	double* elem_area;
//
//	double* elem_density_f; // elem_num
//	double* elem_pcl_n; // elem_num
//	double* elem_pcl_m_s; // elem_num
//	double* elem_pcl_m_f; // elem_num
//	StrainInc* elem_de; // elem_num
//	double* elem_p; // elem_num
//	double* elem_m_de_vol_s; // elem_num
//	double* elem_m_de_vol_f; // elem_num
//
//	// element-node data
//	ElemNodeVM* elem_node_vm_s; // elem_num * 3
//	ElemNodeVM* elem_node_vm_f; // elem_num * 3
//	Force* elem_node_force_s; // elem_num * 3
//	Force* elem_node_force_f; // elem_num * 3
//
//	Acceleration* node_a_s; // node_num
//	Acceleration* node_a_f; // node_num
//	Velocity* node_v_s; // node_num
//	Velocity* node_v_f; // node_num
//	NodeHasVBC* node_has_vbc_s;
//	NodeHasVBC* node_has_vbc_f; // node_num
//	double* node_am_s; // node_num
//	double* node_am_f; // node_num
//	double* node_de_vol_s; // node_num
//	double* node_de_vol_f; // node_num
//	
//	union
//	{
//		struct
//		{
//			size_t* prev_pcl_id0;
//			size_t* prev_pcl_id1;
//			size_t* pcl_in_elem0;
//			size_t* pcl_in_elem1;
//			size_t* node_has_elem0;
//			size_t* node_has_elem1;
//			size_t* node_elem_pair0;
//			size_t* node_elem_pair1;
//		};
//		struct
//		{
//			size_t prev_pcl_id_ui0;
//			size_t prev_pcl_id_ui1;
//			size_t pcl_in_elem_ui0;
//			size_t pcl_in_elem_ui1;
//			size_t node_has_elem_ui0;
//			size_t node_has_elem_ui1;
//			size_t node_elem_pair_ui0;
//			size_t node_elem_pair_ui1;
//		};
//	};
//
//	prev_pcl_id0 = self.prev_pcl_ids[thd.sorted_pcl_in_elem_id];
//	prev_pcl_id1 = self.prev_pcl_ids[thd.sorted_pcl_in_elem_id ^ 1];
//	pcl_in_elem0 = self.pcl_in_elems[thd.sorted_pcl_in_elem_id];
//	pcl_in_elem1 = self.pcl_in_elems[thd.sorted_pcl_in_elem_id ^ 1];
//	node_has_elem0 = self.node_has_elems[0];
//	node_has_elem1 = self.node_has_elems[1];
//	node_elem_pair0 = self.node_elem_pairs[0];
//	node_elem_pair1 = self.node_elem_pairs[1];
//
//	size_t p_id, bin_id, th_id, pos_id;
//	size_t digit_disp, elem_num_tmp, * other_cbin;
//	size_t p_id0 = Block_Low(my_th_id, thread_num, self.valid_pcl_num);
//	size_t p_id1 = Block_Low(my_th_id + 1, thread_num, self.valid_pcl_num);
//	size_t* const elem_count_bin = self.elem_count_bin;
//	size_t* const my_cbin = elem_count_bin + my_th_id * 0x100;
//	size_t* const elem_sum_bin = self.elem_sum_bin;
//	size_t* const my_sbin = elem_sum_bin + my_th_id * 0x100;
//#define data_digit(num, disp) (((num) >> (disp)) & 0xFF)
//#define swap(a, b) \
//		(a) = (a) ^ (b); \
//		(b) = (a) ^ (b); \
//		(a) = (a) ^ (b)
//	for (digit_disp = 0, elem_num_tmp = self.elem_num;
//		elem_num_tmp; digit_disp += 8, elem_num_tmp >>= 8)
//	{
//		memset(my_cbin, 0, 0x100 * sizeof(size_t));
//
//		for (p_id = p_id0; p_id < p_id1; ++p_id)
//			++my_cbin[data_digit(pcl_in_elem0[p_id], digit_disp)];
//
//		my_sbin[0] = my_cbin[0];
//		for (bin_id = 1; bin_id < 0x100; ++bin_id)
//		{
//			my_cbin[bin_id] += my_cbin[bin_id - 1];
//			my_sbin[bin_id] = my_cbin[bin_id];
//		}
//
//#pragma omp barrier
//
//		for (th_id = 0; th_id < my_th_id; ++th_id)
//		{
//			other_cbin = elem_count_bin + th_id * 0x100;
//			for (bin_id = 0; bin_id < 0x100; ++bin_id)
//				my_sbin[bin_id] += other_cbin[bin_id];
//		}
//		for (th_id = my_th_id + 1; th_id < thread_num; ++th_id)
//		{
//			other_cbin = elem_count_bin + th_id * 0x100;
//			for (bin_id = 1; bin_id < 0x100; ++bin_id)
//				my_sbin[bin_id] += other_cbin[bin_id - 1];
//		}
//
//		for (p_id = p_id1; p_id-- > p_id0;)
//		{
//			pos_id = --my_sbin[data_digit(pcl_in_elem0[p_id], digit_disp)];
//			pcl_in_elem1[pos_id] = pcl_in_elem0[p_id];
//			prev_pcl_id1[pos_id] = prev_pcl_id0[p_id];
//		}
//
//		swap(pcl_in_elem_ui0, pcl_in_elem_ui1);
//		swap(prev_pcl_id_ui0, prev_pcl_id_ui1);
//		thd.sorted_pcl_in_elem_id ^= 1;
//#pragma omp barrier
//	}
//
//	// update p_id0, p_id1
//	size_t e_id;
//	e_id = pcl_in_elem0[p_id0];
//	while (p_id0 != SIZE_MAX && e_id == pcl_in_elem0[--p_id0]);
//	++p_id0;
//	e_id = pcl_in_elem0[p_id1];
//	while (p_id1 != SIZE_MAX && e_id == pcl_in_elem0[--p_id1]);
//	++p_id1;
//	
//	size_t ori_p_id, prev_p_id, ne_id;
//	double p_m, p_vol, p_N_m;
//	double one_third_bfx_s, one_third_bfy_s;
//	double one_third_bfx_f, one_third_bfy_f;
//	double e_p_m = 0.0;
//	double e_p_vol = 0.0;
//	double e_s11 = 0.0;
//	double e_s22 = 0.0;
//	double e_s12 = 0.0;
//	double en1_vm = 0.0;
//	double en1_vmx_s = 0.0;
//	double en1_vmy_s = 0.0;
//	double en1_vmx_f = 0.0;
//	double en1_vmy_f = 0.0;
//	double en2_vm = 0.0;
//	double en2_vmx_s = 0.0;
//	double en2_vmy_s = 0.0;
//	double en2_vmx_f = 0.0;
//	double en2_vmy_f = 0.0;
//	double en3_vm = 0.0;
//	double en3_vmx_s = 0.0;
//	double en3_vmy_s = 0.0;
//	double en3_vmx_f = 0.0;
//	double en3_vmy_f = 0.0;
//	double en1_fx_s = 0.0;
//	double en1_fy_s = 0.0;
//	double en1_fx_f = 0.0;
//	double en1_fy_f = 0.0;
//	double en2_fx_s = 0.0;
//	double en2_fy_s = 0.0;
//	double en2_fx_f = 0.0;
//	double en2_fy_f = 0.0;
//	double en3_fx_s = 0.0;
//	double en3_fy_s = 0.0;
//	double en3_fx_f = 0.0;
//	double en3_fy_f = 0.0;
//	e_id = pcl_in_elem0[p_id0];
//	size_t* const my_valid_elem_id = self.valid_elem_id + e_id;
//	size_t* const my_node_has_elem = node_has_elem1 + e_id * 3;
//	size_t* const my_node_elem_pair = node_elem_pair1 + e_id * 3;
//	size_t my_valid_elem_num = 0;
//	for (p_id = p_id0; p_id < p_id1; ++p_id)
//	{
//		prev_p_id = new_to_prev_pcl_map[p_id];
//
//		ori_p_id = pcl_index1[prev_p_id];
//		pcl_index0[p_id] = ori_p_id;
//
//		e_id = pcl_in_elem0[p_id];
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
//		Stress& p_s1 = pcl_stress1[prev_p_id];
//		Stress& p_s0 = pcl_stress0[p_id];
//		p_s0.s11 = p_s1.s11;
//		p_s0.s22 = p_s1.s22;
//		p_s0.s12 = p_s1.s12;
//		e_s11 += p_s0.s11 * p_vol;
//		e_s22 += p_s0.s22 * p_vol;
//		e_s12 += p_s0.s12 * p_vol;
//
//		// map pore pressure
//		pcl_p0[p_id] = pcl_p1[prev_pcl_id];
//		self.elem_p[p_id] += pcl_p0[p_id] * p_vol;
//
//		// map velocity
//		ShapeFunc& p_N1 = pcl_N1[prev_pcl_id];
//		ShapeFunc& p_N0 = pcl_N0[p_id];
//		p_N0.N1 = p_N1.N1;
//		p_N0.N2 = p_N1.N2;
//		p_N0.N3 = p_N1.N3;
//		PclV& p_v_s1 = pcl_v_s1[prev_pcl_id];
//		PclV& p_v_s0 = pcl_v_s0[p_id];
//		p_v_s0.vx = p_v_s1.vx;
//		p_v_s0.vy = p_v_s1.vy;
//		ElemNodeVM& en_vm_s1 = self.elem_node_vm_s[e_id * 3];
//		p_N_m = p_N0.N1 * p_m;
//		en_vm_s1.vm += p_N_m;
//		en_vm_s1.vmx += p_N_m * p_v0.vx;
//		en_vm_s1.vmy += p_N_m * p_v0.vy;
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
//		Force& p_bf = self.pcl_bf_s[ori_p_id];
//		one_third_bfx_s = one_third * p_bf.fx;
//		one_third_bfy_s = one_third * p_bf.fy;
//
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
//	
//	
//	}
//
//	if (md.has_rigid_circle())
//	{
//		RigidCircle& rc = md.get_rigid_circle();
//		Force2D rc_force;
//		rc_force.reset();
//		self.apply_rigid_circle(
//			p_id0,
//			p_id1,
//			pcl_in_elem0,
//			spva0,
//			rc_force,
//			substp_id,
//			thd
//			);
//
//#pragma omp critical
//		self.cf_tmp.combine(rc_force);
//	}
//#pragma omp barrier
//
//#pragma omp master
//	{
//		if (md.has_rigid_circle())
//		{
//			RigidCircle rc = md.get_rigid_circle();
//			rc.set_rcf(self.cf_tmp);
//			rc.update_motion(dt);
//			self.cf_tmp.reset();
//		}
//	}

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
//	for (n_id = n_id0; n_id < n_id1; ++n_id)
//	{
//		double n_am = 0.0;
//		double n_fx = 0.0;
//		double n_fy = 0.0;
//		double n_vm = 0.0;
//		double n_vmx = 0.0;
//		double n_vmy = 0.0;
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
//			ElemShapeFuncAB& e_sf = self.elem_sf_ab[e_id];
//			ElemStrainInc& e_de = self.elem_de[e_id];
//			NodeV& n_v1 = self.node_v[e_n_id.n1];
//			NodeV& n_v2 = self.node_v[e_n_id.n2];
//			NodeV& n_v3 = self.node_v[e_n_id.n3];
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
//		// update stress
//		ori_pcl_id = pcl_index0[p_id];
//		ElemStrainInc& e_de = self.elem_de[e_id];
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
#pragma omp barrier
	return 0;
}

int Step_T2D_CHM_mt::apply_rigid_circle(
	size_t p_id0,
	size_t p_id1,
	size_t* pcl_in_elem,
	SortedPclVarArrays& cur_spva,
	Force2D& rc_cf,
	size_t substp_id,
	ThreadData& thd
	) noexcept
{
	double p_x, p_y;
	double dist, norm_x, norm_y;
	double f_cont, fx_cont, fy_cont;
	size_t e_id, pcl_ori_id;
	Model_T2D_CHM_mt& md = *(Model_T2D_CHM_mt*)model;
	RigidCircle& rc = md.get_rigid_circle();
	const Point2D &rc_centre = rc.get_centre();
	for (size_t p_id = p_id0; p_id < p_id1; ++p_id)
	{
		pcl_ori_id = cur_spva.pcl_index[p_id];
		Position& p_p = pcl_pos[pcl_ori_id];
		Displacement& p_d = cur_spva.pcl_disp[p_id];
		p_x = p_p.x + p_d.ux;
		p_y = p_p.y + p_d.uy;
		if (rc.detect_collision_with_point(
			p_x, p_y, pcl_vol[p_id],
			dist, norm_x, norm_y))
		{
			f_cont = Kn_cont * dist;
			fx_cont = f_cont * norm_x;
			fy_cont = f_cont * norm_y;
			// apply contact force to rigid body
			rc_cf.add_force(
				p_x, p_y,
				-fx_cont, -fy_cont,
				rc_centre.x, rc_centre.y
				);
			// apply contact force to mesh
			ShapeFunc &p_N = cur_spva.pcl_N[p_id];
			e_id = pcl_in_elem[p_id];
			Force& en_f1 = elem_node_force_s[e_id * 3];
			en_f1.fx += p_N.N1 * fx_cont;
			en_f1.fy += p_N.N1 * fy_cont;
			Force& en_f2 = elem_node_force_s[e_id * 3 + 1];
			en_f2.fx += p_N.N2 * fx_cont;
			en_f2.fy += p_N.N2 * fy_cont;
			Force& en_f3 = elem_node_force_s[e_id * 3 + 2];
			en_f3.fx += p_N.N3 * fx_cont;
			en_f3.fy += p_N.N3 * fy_cont;
		}
	}

	return 0;
}
