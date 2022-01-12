#include "SimulationsOMP_pcp.h"

#include <omp.h>

#include "Step_T3D_CHM_ud_mt_subiter.h"

#include <iostream>
#include <fstream>

std::fstream t3d_chm_ud_mt_subit_db_file;

#define one_fourth (0.25)
#define one_third (1.0/3.0)
#define N_min (1.0e-10)
#define Block_Low(th_id, th_num, data_num) ((th_id)*(data_num)/(th_num))

Step_T3D_CHM_ud_mt_subiter::Step_T3D_CHM_ud_mt_subiter(const char* _name) : 
	Step_OMP(_name, "Step_T3D_CHM_ud_mt_subiter",
		&substep_func_omp_T3D_CHM_ud_mt_subiter),
	max_subiter_num(20), mass_factor(1.0), converge_e_kin_ratio(0.001) {}

Step_T3D_CHM_ud_mt_subiter::~Step_T3D_CHM_ud_mt_subiter() {}

static constexpr double u_div_u_cav_pow_cut_off = 1.0e10;

int Step_T3D_CHM_ud_mt_subiter::init_calculation()
{
	Model_T3D_CHM_mt &md = *(Model_T3D_CHM_mt *)model;

	omp_set_num_threads(thread_num);

	t3d_chm_ud_mt_subit_db_file.open("t3d_chm_ud_mt_subiter.csv", std::ios::binary | std::ios::out);

	pcl_m_s = md.pcl_m_s;
	pcl_density_s = md.pcl_density_s;
	pcl_vol_s = md.pcl_vol_s;
	pcl_bf_s = md.pcl_bf_s;
	pcl_bf_f = md.pcl_bf_f;
	pcl_t = md.pcl_t;
	pcl_pos = md.pcl_pos;
	pcl_mat_model = md.pcl_mat_model;
	pcl_mat_model_copy_offset = md.pcl_mat_model_copy_offset;
	if (md.pcl_mat_model_total_size)
	{
		mat_model_copy_mem.alloc(md.pcl_mat_model_total_size);
		mat_model_copy = (char *)mat_model_copy_mem.aligned_address();
	}

	pcl_vol = md.pcl_vol;
	pcl_prev_stress = md.pcl_prev_stress;
	pcl_destrain = md.pcl_destrain;
	pcl_dpstrain = md.pcl_dpstrain;

	Model_T3D_CHM_mt::SortedPclVarArrays& md_spva0
		= md.sorted_pcl_var_arrays[0];
	SortedPclVarArrays& spva0 = sorted_pcl_var_arrays[0];
	spva0.pcl_index = md_spva0.pcl_index;
	spva0.pcl_n = md_spva0.pcl_n;
	spva0.pcl_density_f = md_spva0.pcl_density_f;
	spva0.pcl_v_s = md_spva0.pcl_v_s;
	spva0.pcl_v_f = md_spva0.pcl_v_f;
	spva0.pcl_u_s = md_spva0.pcl_u_s;
	spva0.pcl_u_f = md_spva0.pcl_u_f;
	spva0.pcl_stress = md_spva0.pcl_stress;
	spva0.pcl_p = md_spva0.pcl_p;
	spva0.pcl_strain = md_spva0.pcl_strain;
	spva0.pcl_estrain = md_spva0.pcl_estrain;
	spva0.pcl_pstrain = md_spva0.pcl_pstrain;
	spva0.pcl_N = md_spva0.pcl_N;

	Model_T3D_CHM_mt::SortedPclVarArrays& md_spva1
		= md.sorted_pcl_var_arrays[1];
	SortedPclVarArrays& spva1 = sorted_pcl_var_arrays[1];
	spva1.pcl_index = md_spva1.pcl_index;
	spva1.pcl_n = md_spva1.pcl_n;
	spva1.pcl_density_f = md_spva1.pcl_density_f;
	spva1.pcl_v_s = md_spva1.pcl_v_s;
	spva1.pcl_v_f = md_spva1.pcl_v_f;
	spva1.pcl_u_s = md_spva1.pcl_u_s;
	spva1.pcl_u_f = md_spva1.pcl_u_f;
	spva1.pcl_stress = md_spva1.pcl_stress;
	spva1.pcl_p = md_spva1.pcl_p;
	spva1.pcl_strain = md_spva1.pcl_strain;
	spva1.pcl_estrain = md_spva1.pcl_estrain;
	spva1.pcl_pstrain = md_spva1.pcl_pstrain;
	spva1.pcl_N = md_spva1.pcl_N;

	elem_num = md.elem_num;
	node_num = md.node_num;

	elem_node_id = md.elem_node_id;
	elem_N_abc = md.elem_N_abc;
	elem_N_d = md.elem_N_d;
	elem_vol = md.elem_vol;

	elem_pcl_int_vol = md.elem_pcl_int_vol;
	elem_density_f = md.elem_density_f;
	elem_pcl_n = md.elem_pcl_n;
	elem_pcl_m = md.elem_pcl_m_s;
	elem_de = md.elem_de;
	elem_p = md.elem_p;
	elem_m_de_vol = md.elem_m_de_vol_s;
	elem_m_pde_vol = md.elem_m_de_vol_f;

	// element-node data
	elem_node_vm = md.elem_node_vm_s;
	elem_node_f_ext = md.elem_node_force_s;
	elem_node_f_int = md.elem_node_force_f;

	node_am = md.node_am_s;
	node_a = md.node_a_s;
	node_v = md.node_v_s;
	node_du = md.node_du_s;
	node_vn = md.node_vn_s;
	node_pv = md.node_pv_s;
	node_pdu = md.node_pdu_s;
	node_de_vol = md.node_de_vol_s;
	node_pde_vol = md.node_de_vol_f;

	node_has_vbc = md.node_has_vbc_s;
	node_vbc_vec = md.node_vbc_vec_s;

	Kf = md.Kf;
	m_cav = md.m_cav;
	u_cav = md.u_cav;
	u_cav0 = md.u_cav0;
	Kf_min_ratio = md.Kf_min_ratio;
	u_div_u_cav_cut_off = pow(u_div_u_cav_pow_cut_off, 1.0 / m_cav);

	thread_datas = (ThreadData*)thread_mem.alloc(sizeof(ThreadData) * thread_num);

	char* cur_mem = (char*)cal_mem.alloc(
		  sizeof(size_t) * (md.pcl_num * 4 + 4)
		+ sizeof(size_t) * (md.elem_num * 17 + 4)
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
	valid_elem_ids = (size_t*)cur_mem;
	cur_mem += sizeof(size_t) * md.elem_num;
	node_has_elems[0] = ((size_t*)cur_mem) + 1;
	cur_mem += sizeof(size_t) * (md.elem_num * 4 + 2);
	node_has_elems[1] = ((size_t*)cur_mem) + 1;
	cur_mem += sizeof(size_t) * (md.elem_num * 4 + 2);
	node_elem_pairs[0] = (size_t*)cur_mem;
	cur_mem += sizeof(size_t) * md.elem_num * 4;
	node_elem_pairs[1] = (size_t*)cur_mem;
	cur_mem += sizeof(size_t) * md.elem_num * 4;
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
		thd.max_pcl_vol = 0.0;

		size_t p_id, ori_p_id, e_id;
		size_t p_id0 = Block_Low(my_th_id, thread_num, prev_valid_pcl_num);
		size_t p_id1 = Block_Low(my_th_id + 1, thread_num, prev_valid_pcl_num);
		size_t pcl_in_mesh_num = 0;
		size_t* pcl_in_elem0 = pcl_in_elems[0];
		size_t* prev_pcl_id0 = prev_pcl_ids[0];
		for (p_id = p_id0; p_id < p_id1; ++p_id)
		{
			ori_p_id = spva0.pcl_index[p_id];
			const double p_vol = pcl_m_s[ori_p_id]
				/ (pcl_density_s[ori_p_id] * (1.0 - spva0.pcl_n[p_id]));
			if (thd.max_pcl_vol < p_vol)
				thd.max_pcl_vol = p_vol;
			Position& p_p = pcl_pos[ori_p_id];
			Displacement& p_u_s = spva0.pcl_u_s[p_id];
			p_p.x += p_u_s.ux;
			p_p.y += p_u_s.uy;
			p_p.z += p_u_s.uz;
			p_u_s.ux = 0.0;
			p_u_s.uy = 0.0;
			p_u_s.uz = 0.0;
			ShapeFunc& p_N = spva0.pcl_N[p_id];
			e_id = md.find_pcl_in_which_elem(p_p.x, p_p.y, p_p.z, p_N);
			if (e_id == SIZE_MAX)
				e_id = md.find_pcl_in_which_elem_tol(p_p.x, p_p.y, p_p.z, p_N);
			pcl_in_elem0[p_id] = e_id;
			prev_pcl_id0[p_id] = p_id;
			if (e_id != SIZE_MAX)
				++pcl_in_mesh_num;
		}

#pragma omp critical
		valid_pcl_num += pcl_in_mesh_num;
	}

	pcm_s = md.pcm_s;
	memset(md.contact_substep_ids_s, 0xFF, sizeof(size_t)* md.ori_pcl_num);
	pcm_f = md.pcm_f;
	memset(md.contact_substep_ids_f, 0xFF, sizeof(size_t)* md.ori_pcl_num);
	
	if (md.has_rigid_cylinder())
	{
		prcy = &md.get_rigid_cylinder();
		prcy->reset_cont_force();
	}
	if (md.has_t3d_rigid_mesh())
	{
		prm = &md.get_t3d_rigid_mesh();
		prm->reset_cont_force();
		// set max dist for efficiency
		double max_pcl_radius = thread_datas[0].max_pcl_vol;
		for (size_t th_id = 1; th_id < thread_num; ++th_id)
		{
			if (max_pcl_radius < thread_datas[th_id].max_pcl_vol)
				max_pcl_radius = thread_datas[th_id].max_pcl_vol;
		}
		max_pcl_radius = 0.5 * pow(max_pcl_radius, one_third) * 4.0;
		prm->init_max_dist(max_pcl_radius);
	}

	pcl_in_elems[0][prev_valid_pcl_num] = SIZE_MAX;
	pcl_in_elems[1][prev_valid_pcl_num] = SIZE_MAX;
	valid_elem_num = 0;

	// subiteration
	subiter_index = 0;
	cur_e_kin = 0.0;
	prev_e_kin = 0.0;
	max_e_kin = -1.0;

	return 0;
}

int Step_T3D_CHM_ud_mt_subiter::finalize_calculation()
{
	Model_T3D_CHM_mt &md = *(Model_T3D_CHM_mt *)model;
	md.pcl_num = prev_valid_pcl_num;
	for (size_t t_id = 0; t_id < thread_num; ++t_id)
		thread_datas[t_id].~ThreadData();
	return 0;
}

int substep_func_omp_T3D_CHM_ud_mt_subiter(
	void* _self,
	size_t my_th_id,
	double dt,
	double cur_time,
	size_t substp_id)
{
	typedef Model_T3D_CHM_mt::ShapeFunc ShapeFunc;
	typedef Model_T3D_CHM_mt::DShapeFuncABC DShapeFuncABC;
	typedef Model_T3D_CHM_mt::DShapeFuncD DShapeFuncD;
	typedef Model_T3D_CHM_mt::Force Force;
	typedef Model_T3D_CHM_mt::Position Position;
	typedef Model_T3D_CHM_mt::Displacement Displacement;
	typedef Model_T3D_CHM_mt::Velocity Velocity;
	typedef Model_T3D_CHM_mt::Acceleration Acceleration;
	typedef Model_T3D_CHM_mt::Stress Stress;
	typedef Model_T3D_CHM_mt::Strain Strain;
	typedef Model_T3D_CHM_mt::StrainInc StrainInc;
	typedef Model_T3D_CHM_mt::SortedPclVarArrays SortedPclVarArrays;
	typedef Model_T3D_CHM_mt::ElemNodeIndex ElemNodeIndex;
	typedef Model_T3D_CHM_mt::ElemNodeVM ElemNodeVM;
	typedef Model_T3D_CHM_mt::NodeHasVBC NodeHasVBC;
	typedef Model_T3D_CHM_mt::NodeVBCVec NodeVBCVec;
	typedef Step_T3D_CHM_ud_mt_subiter::ThreadData ThreadData;

	Step_T3D_CHM_ud_mt_subiter& self = *(Step_T3D_CHM_ud_mt_subiter*)(_self);

	if (self.valid_pcl_num == 0)
	{
#pragma omp master
		self.abort_calculation();

#pragma omp barrier
		return 0;
	}

	Model_T3D_CHM_mt& md = *(Model_T3D_CHM_mt*)(self.model);

	const double* const pcl_m_s = self.pcl_m_s;
	const double* const pcl_density_s = self.pcl_density_s;
	const double* const pcl_vol_s = self.pcl_vol_s;
	const Force* const pcl_bf_s = self.pcl_bf_s;
	const Force* const pcl_bf_f = self.pcl_bf_f;
	const Force* const pcl_t = self.pcl_t;
	const Position* const pcl_pos = self.pcl_pos;
	MatModel::MaterialModel** const pcl_mat_model = self.pcl_mat_model;
	const size_t* const pcl_mat_model_copy_offset = self.pcl_mat_model_copy_offset;

	double* const pcl_vol = self.pcl_vol;
	Stress *const pcl_prev_stress = self.pcl_prev_stress;
	StrainInc *const pcl_destrain = self.pcl_destrain;
	StrainInc *const pcl_dpstrain = self.pcl_dpstrain;

	const size_t thread_num = self.thread_num;
	ThreadData& thd = self.thread_datas[my_th_id];
	SortedPclVarArrays& spva0 = self.sorted_pcl_var_arrays[thd.sorted_pcl_var_id];
	thd.sorted_pcl_var_id ^= 1;
	SortedPclVarArrays& spva1 = self.sorted_pcl_var_arrays[thd.sorted_pcl_var_id];
	// 0
	size_t* const pcl_index0 = spva0.pcl_index;
	double* const pcl_n0 = spva0.pcl_n;
	double* const pcl_density_f0 = spva0.pcl_density_f;
	Displacement* const pcl_u_s0 = spva0.pcl_u_s;
	Displacement* const pcl_u_f0 = spva0.pcl_u_f;
	Velocity* const pcl_v_s0 = spva0.pcl_v_s;
	Velocity* const pcl_v_f0 = spva0.pcl_v_f;
	Stress* const pcl_stress0 = spva0.pcl_stress;
	double* const pcl_p0 = spva0.pcl_p;
	Strain* const pcl_strain0 = spva0.pcl_strain;
	Strain* const pcl_estrain0 = spva0.pcl_estrain;
	Strain* const pcl_pstrain0 = spva0.pcl_pstrain;
	ShapeFunc* const pcl_N0 = spva0.pcl_N;
	// 1
	size_t* const pcl_index1 = spva1.pcl_index;
	double* const pcl_n1 = spva1.pcl_n;
	double* const pcl_density_f1 = spva1.pcl_density_f;
	Displacement* const pcl_u_s1 = spva1.pcl_u_s;
	Displacement* const pcl_u_f1 = spva1.pcl_u_f;
	Velocity* const pcl_v_s1 = spva1.pcl_v_s;
	Velocity* const pcl_v_f1 = spva1.pcl_v_f;
	Stress* const pcl_stress1 = spva1.pcl_stress;
	double* const pcl_p1 = spva1.pcl_p;
	Strain* const pcl_strain1 = spva1.pcl_strain;
	Strain* const pcl_estrain1 = spva1.pcl_estrain;
	Strain* const pcl_pstrain1 = spva1.pcl_pstrain;
	ShapeFunc* const pcl_N1 = spva1.pcl_N;

	const ElemNodeIndex* const elem_node_id = self.elem_node_id;
	const DShapeFuncABC* const elem_N_abc = self.elem_N_abc;
	const DShapeFuncD* const elem_N_d = self.elem_N_d;
	const double* const elem_vol = self.elem_vol;

	double* const elem_density_f = self.elem_density_f;
	double* const elem_pcl_n = self.elem_pcl_n;
	double* const elem_pcl_m = self.elem_pcl_m;
	double *elem_pcl_int_vol = self.elem_pcl_int_vol;
	StrainInc* const elem_de = self.elem_de;
	double* const elem_p = self.elem_p;
	double* const elem_m_de_vol = self.elem_m_de_vol;

	ElemNodeVM* const elem_node_vm = self.elem_node_vm;
	Force* const elem_node_f_ext = self.elem_node_f_ext;
	Force* const elem_node_f_int = self.elem_node_f_int;

	NodeHasVBC* const node_has_vbc = self.node_has_vbc;
	NodeVBCVec* const node_vbc_vec = self.node_vbc_vec;
	
	double* const node_am = self.node_am;
	Acceleration* const node_a = self.node_a;
	Velocity* const node_vn = self.node_vn;
	Velocity* const node_v = self.node_v;
	Displacement* const node_du = self.node_du;
	Velocity* const node_pv = self.node_pv;
	Displacement *const node_pdu = self.node_pdu;
	double* const node_de_vol = self.node_de_vol;
	double* const node_pde_vol = self.node_pde_vol;

	union
	{
		struct
		{
			size_t* prev_pcl_id0;
			size_t* prev_pcl_id1;
			size_t* pcl_in_elem0;
			size_t* pcl_in_elem1;
			size_t* node_has_elem0;
			size_t* node_has_elem1;
			size_t* node_elem_pair0;
			size_t* node_elem_pair1;
		};
		struct
		{
			size_t prev_pcl_id_ui0;
			size_t prev_pcl_id_ui1;
			size_t pcl_in_elem_ui0;
			size_t pcl_in_elem_ui1;
			size_t node_has_elem_ui0;
			size_t node_has_elem_ui1;
			size_t node_elem_pair_ui0;
			size_t node_elem_pair_ui1;
		};
	};

	prev_pcl_id0 = self.prev_pcl_ids[thd.sorted_pcl_in_elem_id];
	prev_pcl_id1 = self.prev_pcl_ids[thd.sorted_pcl_in_elem_id ^ 1];
	pcl_in_elem0 = self.pcl_in_elems[thd.sorted_pcl_in_elem_id];
	pcl_in_elem1 = self.pcl_in_elems[thd.sorted_pcl_in_elem_id ^ 1];
	node_has_elem0 = self.node_has_elems[0];
	node_has_elem1 = self.node_has_elems[1];
	node_elem_pair0 = self.node_elem_pairs[0];
	node_elem_pair1 = self.node_elem_pairs[1];

#pragma omp master
	{
		// init rigid body
		self.cf_tmp.reset();

		// init subiteration
		self.subiter_index = 0;
		self.prev_e_kin = 0.0;
		self.max_e_kin = -1.0;
	}

	size_t p_id, bin_id, th_id, pos_id;
	size_t digit_disp, elem_num_tmp, * other_cbin;
	size_t p_id0 = Block_Low(my_th_id, thread_num, self.prev_valid_pcl_num);
	size_t p_id1 = Block_Low(my_th_id + 1, thread_num, self.prev_valid_pcl_num);
	size_t* const elem_count_bin = self.elem_count_bin;
	size_t* const my_cbin = elem_count_bin + my_th_id * 0x100;
	size_t* const elem_sum_bin = self.elem_sum_bin;
	size_t* const my_sbin = elem_sum_bin + my_th_id * 0x100;
#define data_digit(num, disp) (((num) >> (disp)) & 0xFF)
#define swap(a, b) \
		(a) = (a) ^ (b); \
		(b) = (a) ^ (b); \
		(a) = (a) ^ (b)
	for (digit_disp = 0, elem_num_tmp = self.elem_num;
		elem_num_tmp; digit_disp += 8, elem_num_tmp >>= 8)
	{
		memset(my_cbin, 0, 0x100 * sizeof(size_t));

		for (p_id = p_id0; p_id < p_id1; ++p_id)
		{
			++my_cbin[data_digit(pcl_in_elem0[p_id], digit_disp)];
			assert(pcl_in_elem0[p_id] < self.elem_num ||
				pcl_in_elem0[p_id] == SIZE_MAX);
		}

		my_sbin[0] = my_cbin[0];
		for (bin_id = 1; bin_id < 0x100; ++bin_id)
		{
			my_cbin[bin_id] += my_cbin[bin_id - 1];
			my_sbin[bin_id] = my_cbin[bin_id];
		}

#pragma omp barrier

		for (th_id = 0; th_id < my_th_id; ++th_id)
		{
			other_cbin = elem_count_bin + th_id * 0x100;
			for (bin_id = 0; bin_id < 0x100; ++bin_id)
				my_sbin[bin_id] += other_cbin[bin_id];
		}
		for (th_id = my_th_id + 1; th_id < thread_num; ++th_id)
		{
			other_cbin = elem_count_bin + th_id * 0x100;
			for (bin_id = 1; bin_id < 0x100; ++bin_id)
				my_sbin[bin_id] += other_cbin[bin_id - 1];
		}

		for (p_id = p_id1; p_id-- > p_id0;)
		{
			pos_id = --my_sbin[data_digit(pcl_in_elem0[p_id], digit_disp)];
			pcl_in_elem1[pos_id] = pcl_in_elem0[p_id];
			prev_pcl_id1[pos_id] = prev_pcl_id0[p_id];
			assert((pcl_in_elem0[p_id] < self.elem_num ||
				pcl_in_elem0[p_id] == SIZE_MAX) &&
				(prev_pcl_id0[p_id] < self.prev_valid_pcl_num));
		}

		swap(pcl_in_elem_ui0, pcl_in_elem_ui1);
		swap(prev_pcl_id_ui0, prev_pcl_id_ui1);
		thd.sorted_pcl_in_elem_id ^= 1;
#pragma omp barrier
	}

	// update p_id0, p_id1
	size_t e_id;
	p_id0 = Block_Low(my_th_id, thread_num, self.valid_pcl_num);
	e_id = pcl_in_elem0[p_id0];
	while (p_id0 != SIZE_MAX && e_id == pcl_in_elem0[--p_id0]);
	++p_id0;
	assert(p_id0 <= self.valid_pcl_num);
	p_id1 = Block_Low(my_th_id + 1, thread_num, self.valid_pcl_num);
	e_id = pcl_in_elem0[p_id1];
	while (p_id1 != SIZE_MAX && e_id == pcl_in_elem0[--p_id1]);
	++p_id1;
	assert(p_id1 <= self.valid_pcl_num);

	size_t ori_p_id, prev_p_id, ne_id;
	double p_n, p_m_s, p_m_f, p_vol, p_vol_f, p_N_m;
	double one_fourth_bfx_s;
	double one_fourth_bfy_s;
	double one_fourth_bfz_s;
	double en1_vm_s = 0.0;
	double en1_vmx_s = 0.0;
	double en1_vmy_s = 0.0;
	double en1_vmz_s = 0.0;
	double en2_vm_s = 0.0;
	double en2_vmx_s = 0.0;
	double en2_vmy_s = 0.0;
	double en2_vmz_s = 0.0;
	double en3_vm_s = 0.0;
	double en3_vmx_s = 0.0;
	double en3_vmy_s = 0.0;
	double en3_vmz_s = 0.0;
	double en4_vm_s = 0.0;
	double en4_vmx_s = 0.0;
	double en4_vmy_s = 0.0;
	double en4_vmz_s = 0.0;
	double e_p_m_s = 0.0;
	double e_p_m_f = 0.0;
	double e_n = 0.0;
	double e_p_vol_f = 0.0;
	double e_p_vol = 0.0;
	double e_s11 = 0.0;
	double e_s22 = 0.0;
	double e_s33 = 0.0;
	double e_s12 = 0.0;
	double e_s23 = 0.0;
	double e_s31 = 0.0;
	double e_p = 0.0;
	double en1_fx_ext = 0.0;
	double en1_fy_ext = 0.0;
	double en1_fz_ext = 0.0;
	double en2_fx_ext = 0.0;
	double en2_fy_ext = 0.0;
	double en2_fz_ext = 0.0;
	double en3_fx_ext = 0.0;
	double en3_fy_ext = 0.0;
	double en3_fz_ext = 0.0;
	double en4_fx_ext = 0.0;
	double en4_fy_ext = 0.0;
	double en4_fz_ext = 0.0;
	size_t my_valid_elem_num = 0;
	e_id = pcl_in_elem0[p_id0];
	assert(e_id < self.elem_num || e_id == SIZE_MAX);
	size_t* const my_valid_elem_ids = self.valid_elem_ids + e_id;
	size_t* const my_node_has_elem = node_has_elem1 + e_id * 4;
	size_t* const my_node_elem_pair = node_elem_pair1 + e_id * 4;
	for (p_id = p_id0; p_id < p_id1; ++p_id)
	{
		prev_p_id = prev_pcl_id0[p_id];
		assert(prev_p_id < self.prev_valid_pcl_num);

		// ori_p_id
		ori_p_id = pcl_index1[prev_p_id];
		assert(ori_p_id < md.ori_pcl_num);
		pcl_index0[p_id] = ori_p_id;

		// map pcl mass and volume
		// m_s
		p_m_s = pcl_m_s[ori_p_id];
		e_p_m_s += p_m_s;
		// e_n
		e_n += pcl_vol_s[ori_p_id];
		// vol
		p_n = pcl_n1[prev_p_id];
		p_vol = pcl_vol_s[ori_p_id] / (1.0 - p_n);
		pcl_vol[p_id] = p_vol;
		e_p_vol += p_vol;
		// vol_f
		p_vol_f = p_n * p_vol;
		e_p_vol_f += p_vol_f;
		// m_f
		p_m_f = pcl_density_f1[prev_p_id] * p_vol_f;
		e_p_m_f += p_m_f;

		// map stress
		Stress& p_ps = pcl_prev_stress[p_id];
		Stress& p_s1 = pcl_stress1[prev_p_id];
		p_ps.s11 = p_s1.s11;
		p_ps.s22 = p_s1.s22;
		p_ps.s33 = p_s1.s33;
		p_ps.s12 = p_s1.s12;
		p_ps.s23 = p_s1.s23;
		p_ps.s31 = p_s1.s31;
		e_s11 += p_ps.s11 * p_vol;
		e_s22 += p_ps.s22 * p_vol;
		e_s33 += p_ps.s33 * p_vol;
		e_s12 += p_ps.s12 * p_vol;
		e_s23 += p_ps.s23 * p_vol;
		e_s31 += p_ps.s31 * p_vol;

		// map pore pressure
		e_p += pcl_p1[prev_p_id] * p_vol;

		// map velocity
		ShapeFunc& p_N1 = pcl_N1[prev_p_id];
		ShapeFunc& p_N0 = pcl_N0[p_id];
		p_N0.N1 = p_N1.N1;
		p_N0.N2 = p_N1.N2;
		p_N0.N3 = p_N1.N3;
		p_N0.N4 = p_N1.N4;
		// solid velocity
		Velocity& p_v_s1 = pcl_v_s1[prev_p_id];
		Velocity& p_v_s0 = pcl_v_s0[p_id];
		p_v_s0.vx = p_v_s1.vx;
		p_v_s0.vy = p_v_s1.vy;
		p_v_s0.vz = p_v_s1.vz;
		p_N_m = (p_N0.N1 > N_tol ? p_N0.N1 : N_tol) * (p_m_s + p_m_f);
		en1_vm_s += p_N_m;
		en1_vmx_s += p_N_m * p_v_s0.vx;
		en1_vmy_s += p_N_m * p_v_s0.vy;
		en1_vmz_s += p_N_m * p_v_s0.vz;
		p_N_m = (p_N0.N2 > N_tol ? p_N0.N2 : N_tol) * (p_m_s + p_m_f);
		en2_vm_s += p_N_m;
		en2_vmx_s += p_N_m * p_v_s0.vx;
		en2_vmy_s += p_N_m * p_v_s0.vy;
		en2_vmz_s += p_N_m * p_v_s0.vz;
		p_N_m = (p_N0.N3 > N_tol ? p_N0.N3 : N_tol) * (p_m_s + p_m_f);
		en3_vm_s += p_N_m;
		en3_vmx_s += p_N_m * p_v_s0.vx;
		en3_vmy_s += p_N_m * p_v_s0.vy;
		en3_vmz_s += p_N_m * p_v_s0.vz;
		p_N_m = (p_N0.N4 > N_tol ? p_N0.N4 : N_tol) * (p_m_s + p_m_f);
		en4_vm_s += p_N_m;
		en4_vmx_s += p_N_m * p_v_s0.vx;
		en4_vmy_s += p_N_m * p_v_s0.vy;
		en4_vmz_s += p_N_m * p_v_s0.vz;
	
		// solid external load
		const Force& p_bf_s = pcl_bf_s[ori_p_id];
		const Force& p_bf_f = pcl_bf_f[ori_p_id];
		one_fourth_bfx_s = one_fourth * (p_bf_s.fx + p_bf_f.fx);
		one_fourth_bfy_s = one_fourth * (p_bf_s.fy + p_bf_f.fy);
		one_fourth_bfz_s = one_fourth * (p_bf_s.fz + p_bf_f.fz);
		const Force& p_t = pcl_t[ori_p_id];
		en1_fx_ext += one_fourth_bfx_s + p_N0.N1 * p_t.fx;
		en1_fy_ext += one_fourth_bfy_s + p_N0.N1 * p_t.fy;
		en1_fz_ext += one_fourth_bfz_s + p_N0.N1 * p_t.fz;
		en2_fx_ext += one_fourth_bfx_s + p_N0.N2 * p_t.fx;
		en2_fy_ext += one_fourth_bfy_s + p_N0.N2 * p_t.fy;
		en2_fz_ext += one_fourth_bfz_s + p_N0.N2 * p_t.fz;
		en3_fx_ext += one_fourth_bfx_s + p_N0.N3 * p_t.fx;
		en3_fy_ext += one_fourth_bfy_s + p_N0.N3 * p_t.fy;
		en3_fz_ext += one_fourth_bfz_s + p_N0.N3 * p_t.fz;
		en4_fx_ext += one_fourth_bfx_s + p_N0.N4 * p_t.fx;
		en4_fy_ext += one_fourth_bfy_s + p_N0.N4 * p_t.fy;
		en4_fz_ext += one_fourth_bfz_s + p_N0.N4 * p_t.fz;

		// displacement (for contact)
		Displacement& p_u_s1 = pcl_u_s1[prev_p_id];
		Displacement& p_u_s0 = pcl_u_s0[p_id];
		p_u_s0.ux = p_u_s1.ux;
		p_u_s0.uy = p_u_s1.uy;
		p_u_s0.uz = p_u_s1.uz;

		if (e_id != pcl_in_elem0[p_id + 1])
		{
			// v_s
			ElemNodeVM& en1_v = elem_node_vm[e_id * 4];
			en1_v.vm = en1_vm_s;
			en1_v.vmx = en1_vmx_s;
			en1_v.vmy = en1_vmy_s;
			en1_v.vmz = en1_vmz_s;
			ElemNodeVM& en2_v = elem_node_vm[e_id * 4 + 1];
			en2_v.vm = en2_vm_s;
			en2_v.vmx = en2_vmx_s;
			en2_v.vmy = en2_vmy_s;
			en2_v.vmz = en2_vmz_s;
			ElemNodeVM& en3_v = elem_node_vm[e_id * 4 + 2];
			en3_v.vm = en3_vm_s;
			en3_v.vmx = en3_vmx_s;
			en3_v.vmy = en3_vmy_s;
			en3_v.vmz = en3_vmz_s;
			ElemNodeVM& en4_v = elem_node_vm[e_id * 4 + 3];
			en4_v.vm = en4_vm_s;
			en4_v.vmx = en4_vmx_s;
			en4_v.vmy = en4_vmy_s;
			en4_v.vmz = en4_vmz_s;

			elem_pcl_m[e_id] = e_p_m_s + e_p_m_f;
			e_n = 1.0 - e_n / e_p_vol;
			elem_pcl_n[e_id] = e_n;
			elem_density_f[e_id] = e_p_m_f / e_p_vol_f;
			elem_pcl_int_vol[e_id] = e_p_vol;

			e_s11 /= e_p_vol;
			e_s22 /= e_p_vol;
			e_s33 /= e_p_vol;
			e_s12 /= e_p_vol;
			e_s23 /= e_p_vol;
			e_s31 /= e_p_vol;
			e_p /= e_p_vol;
			elem_p[e_id] = e_p;
			if (e_p_vol > elem_vol[e_id])
				e_p_vol = elem_vol[e_id];

			const DShapeFuncABC& e_dN = elem_N_abc[e_id];
			// node 1
			Force& en1_f_int = elem_node_f_int[e_id * 4];
			en1_f_int.fx = (e_dN.dN1_dx * (e_s11 - e_p) + e_dN.dN1_dy * e_s12 + e_dN.dN1_dz * e_s31) * e_p_vol;
			en1_f_int.fy = (e_dN.dN1_dx * e_s12 + e_dN.dN1_dy * (e_s22 - e_p) + e_dN.dN1_dz * e_s23) * e_p_vol;
			en1_f_int.fz = (e_dN.dN1_dx * e_s31 + e_dN.dN1_dy * e_s23 + e_dN.dN1_dz * (e_s33 - e_p)) * e_p_vol;;
			// node 2
			Force& en2_f_int = elem_node_f_int[e_id * 4 + 1];
			en2_f_int.fx = (e_dN.dN2_dx * (e_s11 - e_p) + e_dN.dN2_dy * e_s12 + e_dN.dN2_dz * e_s31) * e_p_vol;
			en2_f_int.fy = (e_dN.dN2_dx * e_s12 + e_dN.dN2_dy * (e_s22 - e_p) + e_dN.dN2_dz * e_s23) * e_p_vol;
			en2_f_int.fz = (e_dN.dN2_dx * e_s31 + e_dN.dN2_dy * e_s23 + e_dN.dN2_dz * (e_s33 - e_p)) * e_p_vol;
			// node 3
			Force& en3_f_int = elem_node_f_int[e_id * 4 + 2];
			en3_f_int.fx = (e_dN.dN3_dx * (e_s11 - e_p) + e_dN.dN3_dy * e_s12 + e_dN.dN3_dz * e_s31) * e_p_vol;
			en3_f_int.fy = (e_dN.dN3_dx * e_s12 + e_dN.dN3_dy * (e_s22 - e_p) + e_dN.dN3_dz * e_s23) * e_p_vol;
			en3_f_int.fz = (e_dN.dN3_dx * e_s31 + e_dN.dN3_dy * e_s23 + e_dN.dN3_dz * (e_s33 - e_p)) * e_p_vol;
			// node 4
			Force& en4_f_int = elem_node_f_int[e_id * 4 + 3];
			en4_f_int.fx = (e_dN.dN4_dx * (e_s11 - e_p) + e_dN.dN4_dy * e_s12 + e_dN.dN4_dz * e_s31) * e_p_vol;
			en4_f_int.fy = (e_dN.dN4_dx * e_s12 + e_dN.dN4_dy * (e_s22 - e_p) + e_dN.dN4_dz * e_s23) * e_p_vol;
			en4_f_int.fz = (e_dN.dN4_dx * e_s31 + e_dN.dN4_dy * e_s23 + e_dN.dN4_dz * (e_s33 - e_p)) * e_p_vol;

			// external nodal force
			// node 1
			Force& en1_f_ext = elem_node_f_ext[e_id * 4];
			en1_f_ext.fx = en1_fx_ext;
			en1_f_ext.fy = en1_fy_ext;
			en1_f_ext.fz = en1_fz_ext;
			// node 2
			Force& en2_f_ext = elem_node_f_ext[e_id * 4 + 1];
			en2_f_ext.fx = en2_fx_ext;
			en2_f_ext.fy = en2_fy_ext;
			en2_f_ext.fz = en2_fz_ext;
			// node 3
			Force& en3_f_ext = elem_node_f_ext[e_id * 4 + 2];
			en3_f_ext.fx = en3_fx_ext;
			en3_f_ext.fy = en3_fy_ext;
			en3_f_ext.fz = en3_fz_ext;
			// node 4
			Force& en4_f_ext = elem_node_f_ext[e_id * 4 + 3];
			en4_f_ext.fx = en4_fx_ext;
			en4_f_ext.fy = en4_fy_ext;
			en4_f_ext.fz = en4_fz_ext;

			ne_id = my_valid_elem_num * 4;
			my_valid_elem_ids[my_valid_elem_num++] = e_id;

			const ElemNodeIndex& eni = elem_node_id[e_id];
			my_node_has_elem[ne_id] = eni.n1;
			my_node_elem_pair[ne_id] = e_id * 4;
			my_node_has_elem[++ne_id] = eni.n2;
			my_node_elem_pair[ne_id] = e_id * 4 + 1;
			my_node_has_elem[++ne_id] = eni.n3;
			my_node_elem_pair[ne_id] = e_id * 4 + 2;
			my_node_has_elem[++ne_id] = eni.n4;
			my_node_elem_pair[ne_id] = e_id * 4 + 3;

			e_id = pcl_in_elem0[p_id + 1];
			assert(e_id < self.elem_num || e_id == SIZE_MAX);

			e_p_m_s = 0.0;
			e_p_m_f = 0.0;
			e_n = 0.0;
			e_p_vol_f = 0.0;
			e_p_vol = 0.0;
			e_s11 = 0.0;
			e_s22 = 0.0;
			e_s33 = 0.0;
			e_s12 = 0.0;
			e_s23 = 0.0;
			e_s31 = 0.0;
			e_p = 0.0;
			en1_vm_s = 0.0;
			en1_vmx_s = 0.0;
			en1_vmy_s = 0.0;
			en1_vmz_s = 0.0;
			en2_vm_s = 0.0;
			en2_vmx_s = 0.0;
			en2_vmy_s = 0.0;
			en2_vmz_s = 0.0;
			en3_vm_s = 0.0;
			en3_vmx_s = 0.0;
			en3_vmy_s = 0.0;
			en3_vmz_s = 0.0;
			en4_vm_s = 0.0;
			en4_vmx_s = 0.0;
			en4_vmy_s = 0.0;
			en4_vmz_s = 0.0;
			en1_fx_ext = 0.0;
			en1_fy_ext = 0.0;
			en1_fz_ext = 0.0;
			en2_fx_ext = 0.0;
			en2_fy_ext = 0.0;
			en2_fz_ext = 0.0;
			en3_fx_ext = 0.0;
			en3_fy_ext = 0.0;
			en3_fz_ext = 0.0;
			en4_fx_ext = 0.0;
			en4_fy_ext = 0.0;
			en4_fz_ext = 0.0;
		}
	}

#pragma omp critical
	self.valid_elem_num += my_valid_elem_num;

	Force3D rc_force;

	if (md.has_rigid_cylinder())
	{
		rc_force.reset();
		self.apply_rigid_cylinder(
			p_id0, p_id1,
			pcl_in_elem0,
			spva0, rc_force,
			substp_id, thd);

#pragma omp critical
		self.cf_tmp.combine(rc_force);
	}

	if (md.has_t3d_rigid_mesh())
	{
		rc_force.reset();
		self.apply_t3d_rigid_mesh(
			p_id0, p_id1,
			pcl_in_elem0,
			spva0, rc_force,
			substp_id, thd);

#pragma omp critical
		self.cf_tmp.combine(rc_force);
	}

#pragma omp barrier
	// sort node-elem pair according to node id
	size_t ve_id;
	memset(my_cbin, 0, 0x100 * sizeof(size_t));
	for (ve_id = 0; ve_id < my_valid_elem_num; ++ve_id)
	{
		++my_cbin[data_digit(my_node_has_elem[ve_id * 4], 0)];
		++my_cbin[data_digit(my_node_has_elem[ve_id * 4 + 1], 0)];
		++my_cbin[data_digit(my_node_has_elem[ve_id * 4 + 2], 0)];
		++my_cbin[data_digit(my_node_has_elem[ve_id * 4 + 3], 0)];
	}

	my_sbin[0] = my_cbin[0];
	for (bin_id = 1; bin_id < 0x100; ++bin_id)
	{
		my_cbin[bin_id] += my_cbin[bin_id - 1];
		my_sbin[bin_id] = my_cbin[bin_id];
	}

#pragma omp barrier

	for (th_id = 0; th_id < my_th_id; ++th_id)
	{
		other_cbin = elem_count_bin + th_id * 0x100;
		for (bin_id = 0; bin_id < 0x100; ++bin_id)
			my_sbin[bin_id] += other_cbin[bin_id];
	}
	for (th_id = my_th_id + 1; th_id < thread_num; ++th_id)
	{
		other_cbin = elem_count_bin + th_id * 0x100;
		for (bin_id = 1; bin_id < 0x100; ++bin_id)
			my_sbin[bin_id] += other_cbin[bin_id - 1];
	}

	for (ve_id = my_valid_elem_num; ve_id-- > 0;)
	{
		pos_id = --my_sbin[data_digit(my_node_has_elem[ve_id * 4], 0)];
		node_has_elem0[pos_id] = my_node_has_elem[ve_id * 4];
		node_elem_pair0[pos_id] = my_node_elem_pair[ve_id * 4];
		pos_id = --my_sbin[data_digit(my_node_has_elem[ve_id * 4 + 1], 0)];
		node_has_elem0[pos_id] = my_node_has_elem[ve_id * 4 + 1];
		node_elem_pair0[pos_id] = my_node_elem_pair[ve_id * 4 + 1];
		pos_id = --my_sbin[data_digit(my_node_has_elem[ve_id * 4 + 2], 0)];
		node_has_elem0[pos_id] = my_node_has_elem[ve_id * 4 + 2];
		node_elem_pair0[pos_id] = my_node_elem_pair[ve_id * 4 + 2];
		pos_id = --my_sbin[data_digit(my_node_has_elem[ve_id * 4 + 3], 0)];
		node_has_elem0[pos_id] = my_node_has_elem[ve_id * 4 + 3];
		node_elem_pair0[pos_id] = my_node_elem_pair[ve_id * 4 + 3];
	}

#pragma omp barrier

#pragma omp master
	{
		node_has_elem0[self.valid_elem_num * 4] = SIZE_MAX;
		node_has_elem1[self.valid_elem_num * 4] = SIZE_MAX;
	}

	size_t ve_id0 = Block_Low(my_th_id, thread_num, self.valid_elem_num * 4);
	size_t ve_id1 = Block_Low(my_th_id + 1, thread_num, self.valid_elem_num * 4);
	size_t node_num_tmp = self.node_num >> 8;
	for (digit_disp = 8; node_num_tmp; digit_disp += 8, node_num_tmp >>= 8)
	{
		memset(my_cbin, 0, sizeof(size_t) * 0x100);

		for (ve_id = ve_id0; ve_id < ve_id1; ++ve_id)
		{
			++my_cbin[data_digit(node_has_elem0[ve_id], digit_disp)];
			assert(node_has_elem0[ve_id] < self.node_num);
		}

		my_sbin[0] = my_cbin[0];
		for (bin_id = 1; bin_id < 0x100; ++bin_id)
		{
			my_cbin[bin_id] += my_cbin[bin_id - 1];
			my_sbin[bin_id] = my_cbin[bin_id];
		}

#pragma omp barrier

		for (th_id = 0; th_id < my_th_id; ++th_id)
		{
			other_cbin = elem_count_bin + th_id * 0x100;
			for (bin_id = 0; bin_id < 0x100; ++bin_id)
				my_sbin[bin_id] += other_cbin[bin_id];
		}
		for (th_id = my_th_id + 1; th_id < thread_num; ++th_id)
		{
			other_cbin = elem_count_bin + th_id * 0x100;
			for (bin_id = 1; bin_id < 0x100; ++bin_id)
				my_sbin[bin_id] += other_cbin[bin_id - 1];
		}

		for (ve_id = ve_id1; ve_id-- > ve_id0;)
		{
			pos_id = --my_sbin[data_digit(node_has_elem0[ve_id], digit_disp)];
			node_has_elem1[pos_id] = node_has_elem0[ve_id];
			node_elem_pair1[pos_id] = node_elem_pair0[ve_id];
			assert(node_has_elem0[ve_id] < self.node_num);
			assert(node_elem_pair0[ve_id] < self.elem_num * 4);
		}

		swap(node_has_elem_ui0, node_has_elem_ui1);
		swap(node_elem_pair_ui0, node_elem_pair_ui1);
#pragma omp barrier
	}

	// modify ne_id0, ne_id1
	size_t n_id;
	n_id = node_has_elem0[ve_id0];
	while (ve_id0 != SIZE_MAX && n_id == node_has_elem0[--ve_id0]);
	++ve_id0;
	assert(ve_id0 <= self.valid_elem_num * 4);
	n_id = node_has_elem0[ve_id1];
	while (ve_id1 != SIZE_MAX && n_id == node_has_elem0[--ve_id1]);
	++ve_id1;
	assert(ve_id1 <= self.valid_elem_num * 4);

	// update node variables
	size_t bc_mask;
	double n_am = 0.0;
	double n_fx = 0.0;
	double n_fy = 0.0;
	double n_fz = 0.0;
	double n_vm = 0.0;
	double n_vmx = 0.0;
	double n_vmy = 0.0;
	double n_vmz = 0.0;
	n_id = node_has_elem0[ve_id0];
	assert(n_id < self.node_num);
	double vbc_len;
	for (ve_id = ve_id0; ve_id < ve_id1; ++ve_id)
	{
		ne_id = node_elem_pair0[ve_id];
		assert(ne_id < self.elem_num * 4);

		e_id = ne_id / 4;
		n_am += elem_pcl_m[e_id];
		const Force& nf_ext = elem_node_f_ext[ne_id];
		const Force& nf_int = elem_node_f_int[ne_id];
		n_fx += nf_ext.fx - nf_int.fx;
		n_fy += nf_ext.fy - nf_int.fy;
		n_fz += nf_ext.fz - nf_int.fz;

		ElemNodeVM& nvm = elem_node_vm[ne_id];
		n_vm += nvm.vm;
		n_vmx += nvm.vmx;
		n_vmy += nvm.vmy;
		n_vmz += nvm.vmz;

		if (n_id != node_has_elem0[ve_id + 1])
		{
			// solid
			n_am *= one_fourth;
			node_am[n_id] = n_am;
			Acceleration& n_a = node_a[n_id];
			n_a.ax = n_fx / n_am;
			n_a.ay = n_fy / n_am;
			n_a.az = n_fz / n_am;
			Velocity& n_vn = node_vn[n_id];
			n_vn.vx = n_vmx / n_vm;
			n_vn.vy = n_vmy / n_vm;
			n_vn.vz = n_vmz / n_vm;
			Velocity& n_v = node_v[n_id];
			n_v.vx = n_vn.vx + n_a.ax * dt;
			n_v.vy = n_vn.vy + n_a.ay * dt;
			n_v.vz = n_vn.vz + n_a.az * dt;
			NodeVBCVec& n_vbc_v = node_vbc_vec[n_id];
			vbc_len = n_a.ax * n_vbc_v.x + n_a.ay * n_vbc_v.y + n_a.az * n_vbc_v.z;
			n_a.ax -= vbc_len * n_vbc_v.x;
			n_a.ay -= vbc_len * n_vbc_v.y;
			n_a.az -= vbc_len * n_vbc_v.z;
			vbc_len = n_v.vx * n_vbc_v.x + n_v.vy * n_vbc_v.y + n_v.vz * n_vbc_v.z;
			n_v.vx -= vbc_len * n_vbc_v.x;
			n_v.vy -= vbc_len * n_vbc_v.y;
			n_v.vz -= vbc_len * n_vbc_v.z;
			NodeHasVBC& n_has_vbc = node_has_vbc[n_id];
			bc_mask = SIZE_MAX + size_t(n_has_vbc.has_vx_bc);
			n_a.iax &= bc_mask;
			n_v.ivx &= bc_mask;
			bc_mask = SIZE_MAX + size_t(n_has_vbc.has_vy_bc);
			n_a.iay &= bc_mask;
			n_v.ivy &= bc_mask;
			bc_mask = SIZE_MAX + size_t(n_has_vbc.has_vz_bc);
			n_a.iaz &= bc_mask;
			n_v.ivz &= bc_mask;
			Displacement& n_du = node_du[n_id];
			n_du.ux = n_v.vx * dt;
			n_du.uy = n_v.vy * dt;
			n_du.uz = n_v.vz * dt;
			
			Velocity& n_pv = node_pv[n_id];
			n_pv.vx = 0.0;
			n_pv.vy = 0.0;
			n_pv.vz = 0.0;
			Displacement &n_pdu = node_pdu[n_id];
			n_pdu.ux = 0.0;
			n_pdu.uy = 0.0;
			n_pdu.uz = 0.0;

			n_id = node_has_elem0[ve_id + 1];
			assert(n_id < self.node_num || n_id == SIZE_MAX);

			n_am = 0.0;
			n_fx = 0.0;
			n_fy = 0.0;
			n_fz = 0.0;
			n_vm = 0.0;
			n_vmx = 0.0;
			n_vmy = 0.0;
			n_vmz = 0.0;
		}
	}

#pragma omp master
	{
#ifdef _DEBUG
		self.prev_valid_pcl_num_tmp = self.prev_valid_pcl_num;
#endif
		self.prev_valid_pcl_num = self.valid_pcl_num;
		self.valid_pcl_num = 0;
	}

#pragma omp barrier
	// cal element strain and "enhancement"
	double e_de_vol;
	for (ve_id = 0; ve_id < my_valid_elem_num; ++ve_id)
	{
		e_id = my_valid_elem_ids[ve_id];
		assert(e_id < self.elem_num);

		const ElemNodeIndex& eni = elem_node_id[e_id];
		const Displacement &n1_du = node_du[eni.n1];
		const Displacement& n2_du = node_du[eni.n2];
		const Displacement& n3_du = node_du[eni.n3];
		const Displacement& n4_du = node_du[eni.n4];
		const DShapeFuncABC& e_dN = elem_N_abc[e_id];
		StrainInc& e_de = elem_de[e_id];
		e_de.de11 = e_dN.dN1_dx * n1_du.ux + e_dN.dN2_dx * n2_du.ux + e_dN.dN3_dx * n3_du.ux + e_dN.dN4_dx * n4_du.ux;
		e_de.de22 = e_dN.dN1_dy * n1_du.uy + e_dN.dN2_dy * n2_du.uy + e_dN.dN3_dy * n3_du.uy + e_dN.dN4_dy * n4_du.uy;
		e_de.de33 = e_dN.dN1_dz * n1_du.uz + e_dN.dN2_dz * n2_du.uz + e_dN.dN3_dz * n3_du.uz + e_dN.dN4_dz * n4_du.uz;
		e_de.de12 = (e_dN.dN1_dx * n1_du.uy + e_dN.dN2_dx * n2_du.uy + e_dN.dN3_dx * n3_du.uy + e_dN.dN4_dx * n4_du.uy
				   + e_dN.dN1_dy * n1_du.ux + e_dN.dN2_dy * n2_du.ux + e_dN.dN3_dy * n3_du.ux + e_dN.dN4_dy * n4_du.ux) * 0.5;
		e_de.de23 = (e_dN.dN1_dy * n1_du.uz + e_dN.dN2_dy * n2_du.uz + e_dN.dN3_dy * n3_du.uz + e_dN.dN4_dy * n4_du.uz
				   + e_dN.dN1_dz * n1_du.uy + e_dN.dN2_dz * n2_du.uy + e_dN.dN3_dz * n3_du.uy + e_dN.dN4_dz * n4_du.uy) * 0.5;
		e_de.de31 = (e_dN.dN1_dz * n1_du.ux + e_dN.dN2_dz * n2_du.ux + e_dN.dN3_dz * n3_du.ux + e_dN.dN4_dz * n4_du.ux
				   + e_dN.dN1_dx * n1_du.uz + e_dN.dN2_dx * n2_du.uz + e_dN.dN3_dx * n3_du.uz + e_dN.dN4_dx * n4_du.uz) * 0.5;
		e_de_vol = e_de.de11 + e_de.de22 + e_de.de33;
		elem_m_de_vol[e_id] = elem_pcl_m[e_id] * e_de_vol;		
		e_de_vol *= one_third;
		e_de.de11 -= e_de_vol;
		e_de.de22 -= e_de_vol;
		e_de.de33 -= e_de_vol;
	}

#pragma omp barrier
	double n_am_de_vol = 0.0;
	n_id = node_has_elem0[ve_id0];
	assert(n_id < self.node_num);
	for (ve_id = ve_id0; ve_id < ve_id1; ++ve_id)
	{
		e_id = node_elem_pair0[ve_id] / 4;
		assert(e_id < self.elem_num);
		n_am_de_vol += elem_m_de_vol[e_id];
		if (n_id != node_has_elem0[ve_id + 1])
		{
			node_de_vol[n_id] = n_am_de_vol * one_fourth / node_am[n_id];
			n_id = node_has_elem0[ve_id + 1];
			assert(n_id < self.node_num || n_id == SIZE_MAX);
			n_am_de_vol = 0.0;
		}
	}

#pragma omp barrier

	const Displacement *pn1_du, *pn2_du, *pn3_du, *pn4_du;
	StrainInc *pe_de;
	const double* estrain, * pstrain, * dstress;
	e_id = SIZE_MAX;
	for (p_id = p_id0; p_id < p_id1; ++p_id)
	{
		if (e_id != pcl_in_elem0[p_id])
		{
			e_id = pcl_in_elem0[p_id];
			assert(e_id < self.elem_num);

			const ElemNodeIndex& eni = elem_node_id[e_id];
			
			pn1_du = node_du + eni.n1;
			pn2_du = node_du + eni.n2;
			pn3_du = node_du + eni.n3;
			pn4_du = node_du + eni.n4;

			e_de_vol = (node_de_vol[eni.n1] + node_de_vol[eni.n2]
					  + node_de_vol[eni.n3] + node_de_vol[eni.n4]) * one_fourth;
			const double e_de_vol_f = -e_de_vol / elem_pcl_n[e_id];

			elem_pcl_n[e_id] = (e_de_vol + elem_pcl_n[e_id]) / (1.0 + e_de_vol);
			//elem_pcl_int_vol[e_id] *= (1.0 + e_de_vol);

			elem_density_f[e_id] /= (1.0 - e_de_vol_f);
			double Kf_ratio = 1.0;
			// cavitation
			if (self.m_cav != 0.0 && elem_p[e_id] < self.u_cav0)
			{
				const double tmp1 = fabs((elem_p[e_id] - self.u_cav0) / (self.u_cav - self.u_cav0));
				if (tmp1 < self.u_div_u_cav_cut_off)
					Kf_ratio = self.Kf_min_ratio + (1.0 - self.Kf_min_ratio) / (1.0 + pow(tmp1, self.m_cav));
				else
					Kf_ratio = self.Kf_min_ratio + (1.0 - self.Kf_min_ratio) / (1.0 + u_div_u_cav_pow_cut_off);
			}
			elem_p[e_id] += Kf_ratio * self.Kf * e_de_vol_f;
			
			//if (self.substep_index % 1000 == 999 && e_id == 149)
			//	t3d_chm_ud_mt_subit_db_file << self.substep_index << ", "
			//	<< self.current_time << ", " << elem_p[e_id] << ", "
			//	<< Kf_ratio << ", " << Kf_ratio * self.Kf << ",\n";

			pe_de = elem_de + e_id;
			e_de_vol *= one_third;
			pe_de->de11 += e_de_vol;
			pe_de->de22 += e_de_vol;
			pe_de->de33 += e_de_vol;
		}

		ShapeFunc& p_N = pcl_N0[p_id];
		Displacement& p_u_s0 = pcl_u_s0[p_id];
		p_u_s0.ux += p_N.N1 * pn1_du->ux + p_N.N2 * pn2_du->ux + p_N.N3 * pn3_du->ux + p_N.N4 * pn4_du->ux;
		p_u_s0.uy += p_N.N1 * pn1_du->uy + p_N.N2 * pn2_du->uy + p_N.N3 * pn3_du->uy + p_N.N4 * pn4_du->uy;
		p_u_s0.uz += p_N.N1 * pn1_du->uz + p_N.N2 * pn2_du->uz + p_N.N3 * pn3_du->uz + p_N.N4 * pn4_du->uz;
		
		ori_p_id = pcl_index0[p_id];

		// update stress
		MatModel::MaterialModel& pcl_mm = *pcl_mat_model[ori_p_id];
		char* pcl_mm_copy = self.mat_model_copy + pcl_mat_model_copy_offset[ori_p_id];
		pcl_mm.store_to(pcl_mm_copy);
		pcl_mm.integrate(pe_de->de);
		dstress = pcl_mm.get_dstress();
		Stress& p_ps = pcl_prev_stress[p_id];
		Stress& p_s0 = pcl_stress0[p_id];
		p_s0.s11 = p_ps.s11 + dstress[0];
		p_s0.s22 = p_ps.s22 + dstress[1];
		p_s0.s33 = p_ps.s33 + dstress[2];
		p_s0.s12 = p_ps.s12 + dstress[3];
		p_s0.s23 = p_ps.s23 + dstress[4];
		p_s0.s31 = p_ps.s31 + dstress[5];

		estrain = pcl_mm.get_dstrain_e();
		StrainInc& p_dee = pcl_destrain[p_id];
		p_dee.de11 = estrain[0];
		p_dee.de22 = estrain[1];
		p_dee.de33 = estrain[2];
		p_dee.de12 = estrain[3];
		p_dee.de23 = estrain[4];
		p_dee.de31 = estrain[5];

		pstrain = pcl_mm.get_dstrain_p();
		StrainInc& p_dpe = pcl_dpstrain[p_id];
		p_dpe.de11 = pstrain[0];
		p_dpe.de22 = pstrain[1];
		p_dpe.de33 = pstrain[2];
		p_dpe.de12 = pstrain[3];
		p_dpe.de23 = pstrain[4];
		p_dpe.de31 = pstrain[5];
	}

	// start subiteration
	int subiteration_res = 0;
	size_t subiter_index = 0;
	while (subiteration_res == 0 &&
		subiter_index < self.max_subiter_num) // controlled
	{
		subiteration_res = self.subiteration(my_th_id,
			p_id0, p_id1, ve_id0, ve_id1,
			my_valid_elem_ids, my_valid_elem_num,
			spva0, pcl_in_elem0,
			node_has_elem0, node_elem_pair0);
		++subiter_index;
#pragma omp master
		self.subiter_index = subiter_index;
	}

#pragma omp barrier
	const Acceleration* pn1_a, * pn2_a, * pn3_a, * pn4_a;
	//const Velocity* pn1_v, * pn2_v, * pn3_v, * pn4_v;
	double e_density_f, p_x, p_y, p_z;
	size_t p_e_id, pcl_in_mesh_num = 0;
	e_id = SIZE_MAX;
	thd.sorted_pcl_in_elem_id ^= 1;
	for (p_id = p_id0; p_id < p_id1; ++p_id)
	{
		if (e_id != pcl_in_elem0[p_id])
		{
			e_id = pcl_in_elem0[p_id];
			assert(e_id < self.elem_num);

			const ElemNodeIndex& eni = elem_node_id[e_id];
			pn1_a = node_a + eni.n1;
			pn2_a = node_a + eni.n2;
			pn3_a = node_a + eni.n3;
			pn4_a = node_a + eni.n4;
			//pn1_v = node_v + eni.n1;
			//pn2_v = node_v + eni.n2;
			//pn3_v = node_v + eni.n3;
			//pn4_v = node_v + eni.n4;

			e_n = elem_pcl_n[e_id];
			e_density_f = elem_density_f[e_id];
			e_p = elem_p[e_id];
			pe_de = elem_de + e_id;
		}

		// update velocity
		ShapeFunc& p_N = pcl_N0[p_id];
		Velocity& p_v_s0 = pcl_v_s0[p_id];
		p_v_s0.vx += (p_N.N1 * pn1_a->ax + p_N.N2 * pn2_a->ax + p_N.N3 * pn3_a->ax + p_N.N4 * pn4_a->ax) * dt;
		p_v_s0.vy += (p_N.N1 * pn1_a->ay + p_N.N2 * pn2_a->ay + p_N.N3 * pn3_a->ay + p_N.N4 * pn4_a->ay) * dt;
		p_v_s0.vz += (p_N.N1 * pn1_a->az + p_N.N2 * pn2_a->az + p_N.N3 * pn3_a->az + p_N.N4 * pn4_a->az) * dt;
		Velocity& p_v_f0 = pcl_v_f0[p_id];
		p_v_f0.vx = p_v_s0.vx;
		p_v_f0.vy = p_v_s0.vy;
		p_v_f0.vz = p_v_s0.vz;

		// update displacement
		Displacement& p_u_s0 = pcl_u_s0[p_id];
		//p_u_s0.ux += (p_N.N1 * pn1_v->vx + p_N.N2 * pn2_v->vx + p_N.N3 * pn3_v->vx + p_N.N4 * pn4_v->vx) * dt;
		//p_u_s0.uy += (p_N.N1 * pn1_v->vy + p_N.N2 * pn2_v->vy + p_N.N3 * pn3_v->vy + p_N.N4 * pn4_v->vy) * dt;
		//p_u_s0.uz += (p_N.N1 * pn1_v->vz + p_N.N2 * pn2_v->vz + p_N.N3 * pn3_v->vz + p_N.N4 * pn4_v->vz) * dt;
		Displacement& p_u_f0 = pcl_u_f0[p_id];
		p_u_f0.ux = p_u_s0.ux;
		p_u_f0.uy = p_u_s0.uy;
		p_u_f0.uz = p_u_s0.uz;

		// update location (in which element)
		ori_p_id = pcl_index0[p_id];
		assert(ori_p_id < md.ori_pcl_num);
		const Position& p_p = pcl_pos[ori_p_id];
		p_x = p_p.x + p_u_s0.ux;
		p_y = p_p.y + p_u_s0.uy;
		p_z = p_p.z + p_u_s0.uz;
		p_e_id = e_id;
		if (!md.is_in_element(p_x, p_y, p_z, e_id, p_N))
		{
			p_e_id = md.find_pcl_in_which_elem(p_x, p_y, p_z, p_N);
			if (p_e_id == SIZE_MAX)
			{
				if (md.is_in_element_tol(p_x, p_y, p_z, e_id, p_N))
					p_e_id = e_id;
				else
					p_e_id = md.find_pcl_in_which_elem_tol(p_x, p_y, p_z, p_N);
			}
		}
		if (p_e_id != SIZE_MAX) // in mesh
			++pcl_in_mesh_num;
		pcl_in_elem1[p_id] = p_e_id;
		prev_pcl_id1[p_id] = p_id;
		assert(p_e_id < self.elem_num || p_e_id == SIZE_MAX);

		// update n
		pcl_n0[p_id] = e_n;
		// update density
		pcl_density_f0[p_id] = e_density_f;
		// update pore pressure
		pcl_p0[p_id] = e_p;

		prev_p_id = prev_pcl_id0[p_id];
#ifdef _DEBUG
		assert(prev_p_id < self.prev_valid_pcl_num_tmp);
#endif
		Strain& p_e1 = pcl_strain1[prev_p_id];
		Strain& p_e0 = pcl_strain0[p_id];
		p_e0.e11 = p_e1.e11 + pe_de->de11;
		p_e0.e22 = p_e1.e22 + pe_de->de22;
		p_e0.e33 = p_e1.e33 + pe_de->de33;
		p_e0.e12 = p_e1.e12 + pe_de->de12;
		p_e0.e23 = p_e1.e23 + pe_de->de23;
		p_e0.e31 = p_e1.e31 + pe_de->de31;

		const StrainInc& dee = pcl_destrain[p_id];
		const Strain& p_ee1 = pcl_estrain1[prev_p_id];
		Strain& p_ee0 = pcl_estrain0[p_id];
		p_ee0.e11 = p_ee1.e11 + dee.de11;
		p_ee0.e22 = p_ee1.e22 + dee.de22;
		p_ee0.e33 = p_ee1.e33 + dee.de33;
		p_ee0.e12 = p_ee1.e12 + dee.de12;
		p_ee0.e23 = p_ee1.e23 + dee.de23;
		p_ee0.e31 = p_ee1.e31 + dee.de31;

		const StrainInc& dpe = pcl_dpstrain[p_id];
		const Strain& p_pe1 = pcl_pstrain1[prev_p_id];
		Strain& p_pe0 = pcl_pstrain0[p_id];
		p_pe0.e11 = p_pe1.e11 + dpe.de11;
		p_pe0.e22 = p_pe1.e22 + dpe.de22;
		p_pe0.e33 = p_pe1.e33 + dpe.de33;
		p_pe0.e12 = p_pe1.e12 + dpe.de12;
		p_pe0.e23 = p_pe1.e23 + dpe.de23;
		p_pe0.e31 = p_pe1.e31 + dpe.de31;
	}

#pragma omp critical
	self.valid_pcl_num += pcl_in_mesh_num;

#pragma omp master
	{
		// update rigid body motion
		if (md.has_rigid_cylinder())
		{
			RigidCylinder& rcy = md.get_rigid_cylinder();
			rcy.set_cont_force(self.cf_tmp);
			rcy.update_motion(dt);
			//self.cf_tmp.reset();
		}

		if (md.has_t3d_rigid_mesh())
		{
			RigidObjectByT3DMesh& rb = md.get_t3d_rigid_mesh();
			rb.set_cont_force(self.cf_tmp);
			rb.update_motion(dt);
			//self.cf_tmp.reset();
		}

		pcl_in_elem0[self.prev_valid_pcl_num] = SIZE_MAX;
		pcl_in_elem1[self.prev_valid_pcl_num] = SIZE_MAX;
		self.valid_elem_num = 0;
		self.continue_calculation();
	}

#pragma omp barrier
	return 0;
}
