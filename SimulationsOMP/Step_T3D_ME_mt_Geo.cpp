#include "SimulationsOMP_pcp.h"

#include <fstream>

#include <assert.h>
#include <chrono>
#include <omp.h>

#include "Step_T3D_ME_mt_Geo.h"

#ifdef _DEBUG
static std::fstream res_file_T3D_ME_mt_Geo;
#endif

#define one_third (1.0/3.0)
#define one_fourth (0.25)
#define MAX_PCL_RADIUX_SCALE (5.0)
#define Block_Low(th_id, th_num, data_num) ((th_id)*(data_num)/(th_num))

Step_T3D_ME_mt_Geo::Step_T3D_ME_mt_Geo(const char* _name) :
	Step_OMP(_name, "Step_T3D_ME_mt_Geo", &substep_func_omp_T3D_ME_mt_Geo),
	f_ub_ratio_limit(0.001), e_kin_ratio_limit(0.001) {}

Step_T3D_ME_mt_Geo::~Step_T3D_ME_mt_Geo() {}

int Step_T3D_ME_mt_Geo::init_calculation()
{
#ifdef _DEBUG
	res_file_T3D_ME_mt_Geo.open("t3d_stp_mt_geo.csv", std::ios::out | std::ios::binary);
#endif

	Model_T3D_ME_mt &md = *(Model_T3D_ME_mt *)model;

	omp_set_num_threads(thread_num);

	pcl_m = md.pcl_m;
	pcl_bf = md.pcl_bf;
	pcl_t = md.pcl_t;
	pcl_pos = md.pcl_pos;
	pcl_vol = md.pcl_vol;
	pcl_mat_model = md.pcl_mat_model;

	Model_T3D_ME_mt::SortedPclVarArrays& md_spva0 = md.sorted_pcl_var_arrays[0];
	SortedPclVarArrays& spva0 = sorted_pcl_var_arrays[0];
	spva0.pcl_index = md_spva0.pcl_index;
	spva0.pcl_density = md_spva0.pcl_density;
	spva0.pcl_disp = md_spva0.pcl_disp;
	spva0.pcl_v = md_spva0.pcl_v;
	spva0.pcl_stress = md_spva0.pcl_stress;
	spva0.pcl_strain = md_spva0.pcl_strain;
	spva0.pcl_estrain = md_spva0.pcl_estrain;
	spva0.pcl_pstrain = md_spva0.pcl_pstrain;
	spva0.pcl_N = md_spva0.pcl_N;

	Model_T3D_ME_mt::SortedPclVarArrays& md_spva1 = md.sorted_pcl_var_arrays[1];
	SortedPclVarArrays& spva1 = sorted_pcl_var_arrays[1];
	spva1.pcl_index = md_spva1.pcl_index;
	spva1.pcl_density = md_spva1.pcl_density;
	spva1.pcl_disp = md_spva1.pcl_disp;
	spva1.pcl_v = md_spva1.pcl_v;
	spva1.pcl_stress = md_spva1.pcl_stress;
	spva1.pcl_strain = md_spva1.pcl_strain;
	spva1.pcl_estrain = md_spva1.pcl_estrain;
	spva1.pcl_pstrain = md_spva1.pcl_pstrain;
	spva1.pcl_N = md_spva1.pcl_N;

	elem_num = md.elem_num;
	node_num = md.node_num;

	elem_node_id = md.elem_node_id;
	elem_dN_abc = md.elem_dN_abc;
	elem_dN_d = md.elem_dN_d;
	elem_vol = md.elem_vol;

	elem_density = md.elem_density;
	elem_pcl_m = md.elem_pcl_m;
	elem_de = md.elem_de;
	elem_m_de_vol = md.elem_m_de_vol;

	elem_node_vm = md.elem_node_vm;
	elem_node_force = md.elem_node_force;

	node_a = md.node_a;
	node_v = md.node_v;
	node_has_vbc = md.node_has_vbc;
	node_vbc_vec = md.node_vbc_vec;
	node_am = md.node_am;
	node_de_vol = md.node_de_vol;

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
	valid_elem_id = (size_t*)cur_mem;
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
			p_p.z += p_d.uz;
			p_d.ux = 0.0;
			p_d.uy = 0.0;
			p_d.uz = 0.0;
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
	
	pcl_in_elems[0][prev_valid_pcl_num] = SIZE_MAX;
	pcl_in_elems[1][prev_valid_pcl_num] = SIZE_MAX;
	valid_elem_num = 0;
	
	init_f_ub_is_init = false;
	e_kin_max_is_init = false;
	e_kin_prev = 0.0;
	return 0;
}

int Step_T3D_ME_mt_Geo::finalize_calculation()
{
	Model_T3D_ME_mt& md = *(Model_T3D_ME_mt*)model;
	md.pcl_num = valid_pcl_num;
	for (size_t t_id = 0; t_id < thread_num; ++t_id)
		thread_datas[t_id].~ThreadData();
	return 0;
}

int substep_func_omp_T3D_ME_mt_Geo(
	void* _self,
	size_t my_th_id,
	double dt,
	double cur_time,
	size_t substp_id)
{
	typedef Model_T3D_ME_mt::Force Force;
	typedef Model_T3D_ME_mt::Position Position;
	typedef Model_T3D_ME_mt::Displacement Displacement;
	typedef Model_T3D_ME_mt::Velocity Velocity;
	typedef Model_T3D_ME_mt::Acceleration Acceleration;
	typedef Model_T3D_ME_mt::ShapeFunc ShapeFunc;
	typedef Model_T3D_ME_mt::DShapeFuncABC DShapeFuncABC;
	typedef Model_T3D_ME_mt::DShapeFuncD DShapeFuncD;
	typedef Model_T3D_ME_mt::Stress Stress;
	typedef Model_T3D_ME_mt::Strain Strain;
	typedef Model_T3D_ME_mt::StrainInc StrainInc;
	typedef Model_T3D_ME_mt::ElemNodeIndex ElemNodeIndex;
	typedef Model_T3D_ME_mt::ElemNodeVM ElemNodeVM;
	typedef Model_T3D_ME_mt::NodeHasVBC NodeHasVBC;
	typedef Model_T3D_ME_mt::NodeVBCVec NodeVBCVec;
	typedef Model_T3D_ME_mt::SortedPclVarArrays SortedPclVarArrays;
	typedef Step_T3D_ME_mt_Geo::ThreadData ThreadData;

	Step_T3D_ME_mt_Geo& self = *(Step_T3D_ME_mt_Geo*)(_self);

	int efef = self.substep_index;
	if (self.substep_index == 10792)
	{
		int eefef = 0;
	}

	if (self.valid_pcl_num == 0)
	{
#pragma omp master
		self.abort_calculation();

#pragma omp barrier
		return 0;
	}

	Model_T3D_ME_mt& md = *(Model_T3D_ME_mt*)(self.model);

	const double* const pcl_m = self.pcl_m;
	const Force* const pcl_bf = self.pcl_bf;
	const Force* const pcl_t = self.pcl_t;
	const Position* const pcl_pos = self.pcl_pos;
	double* const pcl_vol = self.pcl_vol;
	MatModel::MaterialModel** const pcl_mat_model = self.pcl_mat_model;

	const size_t thread_num = self.thread_num;
	ThreadData& thd = self.thread_datas[my_th_id];
	SortedPclVarArrays& spva0 = self.sorted_pcl_var_arrays[thd.sorted_pcl_var_id];
	thd.sorted_pcl_var_id ^= 1;
	SortedPclVarArrays& spva1 = self.sorted_pcl_var_arrays[thd.sorted_pcl_var_id];

	size_t* const pcl_index0 = spva0.pcl_index;
	double* const pcl_density0 = spva0.pcl_density;
	Displacement* const pcl_disp0 = spva0.pcl_disp;
	Velocity* const pcl_v0 = spva0.pcl_v;
	Stress* const pcl_stress0 = spva0.pcl_stress;
	Strain* const pcl_strain0 = spva0.pcl_strain;
	Strain* const pcl_estrain0 = spva0.pcl_estrain;
	Strain* const pcl_pstrain0 = spva0.pcl_pstrain;
	ShapeFunc* const pcl_N0 = spva0.pcl_N;

	size_t* const pcl_index1 = spva1.pcl_index;
	double* const pcl_density1 = spva1.pcl_density;
	Displacement* const pcl_disp1 = spva1.pcl_disp;
	Velocity* const pcl_v1 = spva1.pcl_v;
	Stress* const pcl_stress1 = spva1.pcl_stress;
	Strain* const pcl_strain1 = spva1.pcl_strain;
	Strain* const pcl_estrain1 = spva1.pcl_estrain;
	Strain* const pcl_pstrain1 = spva1.pcl_pstrain;
	ShapeFunc* const pcl_N1 = spva1.pcl_N;

	ElemNodeIndex* const elem_node_id = self.elem_node_id;
	DShapeFuncABC* const elem_dN_abc = self.elem_dN_abc;
	DShapeFuncD* const elem_dN_d = self.elem_dN_d;
	double* const elem_vol = self.elem_vol;

	double* const elem_density = self.elem_density;
	double* const elem_pcl_m = self.elem_pcl_m;
	StrainInc* const elem_de = self.elem_de;
	double* const elem_m_de_vol = self.elem_m_de_vol;

	ElemNodeVM* const elem_node_vm = self.elem_node_vm;
	Force* const elem_node_force = self.elem_node_force;

	Acceleration* const node_a = self.node_a;
	Velocity* const node_v = self.node_v;
	NodeHasVBC* const node_has_vbc = self.node_has_vbc;
	NodeVBCVec* const node_vbc_vec = self.node_vbc_vec;
	double* const node_am = self.node_am;
	double* const node_de_vol = self.node_de_vol;

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
		self.f_ub = 0.0;
		self.e_kin = 0.0;
	}

	// sort particle variables
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
	double p_m, p_vol, p_N_m;
	double one_fourth_bfx, one_fourth_bfy, one_fourth_bfz;
	double e_p_m = 0.0;
	double e_p_vol = 0.0;
	double e_s11 = 0.0;
	double e_s22 = 0.0;
	double e_s33 = 0.0;
	double e_s12 = 0.0;
	double e_s23 = 0.0;
	double e_s31 = 0.0;
	double en1_vm = 0.0;
	double en1_vmx = 0.0;
	double en1_vmy = 0.0;
	double en1_vmz = 0.0;
	double en2_vm = 0.0;
	double en2_vmx = 0.0;
	double en2_vmy = 0.0;
	double en2_vmz = 0.0;
	double en3_vm = 0.0;
	double en3_vmx = 0.0;
	double en3_vmy = 0.0;
	double en3_vmz = 0.0;
	double en4_vm = 0.0;
	double en4_vmx = 0.0;
	double en4_vmy = 0.0;
	double en4_vmz = 0.0;
	double en1_fx = 0.0;
	double en1_fy = 0.0;
	double en1_fz = 0.0;
	double en2_fx = 0.0;
	double en2_fy = 0.0;
	double en2_fz = 0.0;
	double en3_fx = 0.0;
	double en3_fy = 0.0;
	double en3_fz = 0.0;
	double en4_fx = 0.0;
	double en4_fy = 0.0;
	double en4_fz = 0.0;
	e_id = pcl_in_elem0[p_id0];
	assert(e_id < self.elem_num);
	size_t* const my_valid_elem_id = self.valid_elem_id + e_id;
	size_t* const my_node_has_elem = node_has_elem1 + e_id * 4;
	size_t* const my_node_elem_pair = node_elem_pair1 + e_id * 4;
	size_t my_valid_elem_num = 0;
	for (p_id = p_id0; p_id < p_id1; ++p_id)
	{
		// pcl index
		prev_p_id = prev_pcl_id0[p_id];
		assert(prev_p_id < self.prev_valid_pcl_num);
		ori_p_id = pcl_index1[prev_p_id];
		assert(ori_p_id < md.ori_pcl_num);
		pcl_index0[p_id] = ori_p_id;

		// map pcl mass
		p_m = pcl_m[ori_p_id];
		e_p_m += p_m;

		// map pcl volume
		p_vol = p_m / pcl_density1[prev_p_id];
		pcl_vol[p_id] = p_vol;
		e_p_vol += p_vol;

		// map pcl stress
		Stress& p_s1 = pcl_stress1[prev_p_id];
		Stress& p_s0 = pcl_stress0[p_id];
		p_s0.s11 = p_s1.s11;
		p_s0.s22 = p_s1.s22;
		p_s0.s33 = p_s1.s33;
		p_s0.s12 = p_s1.s12;
		p_s0.s23 = p_s1.s23;
		p_s0.s31 = p_s1.s31;
		e_s11 += p_s0.s11 * p_vol;
		e_s22 += p_s0.s22 * p_vol;
		e_s33 += p_s0.s33 * p_vol;
		e_s12 += p_s0.s12 * p_vol;
		e_s23 += p_s0.s23 * p_vol;
		e_s31 += p_s0.s31 * p_vol;

		// shape function
		ShapeFunc& p_N1 = pcl_N1[prev_p_id];
		ShapeFunc& p_N0 = pcl_N0[p_id];
		p_N0.N1 = p_N1.N1;
		p_N0.N2 = p_N1.N2;
		p_N0.N3 = p_N1.N3;
		p_N0.N4 = p_N1.N4;
		// map velocity
#define N_tol (1.0e-10)
		Velocity& p_v1 = pcl_v1[prev_p_id];
		Velocity& p_v0 = pcl_v0[p_id];
		p_v0.vx = p_v1.vx;
		p_v0.vy = p_v1.vy;
		p_v0.vz = p_v1.vz;
		p_N_m = (p_N0.N1 > N_tol ? p_N0.N1 : N_tol) * p_m;
		en1_vm += p_N_m;
		en1_vmx += p_N_m * p_v0.vx;
		en1_vmy += p_N_m * p_v0.vy;
		en1_vmz += p_N_m * p_v0.vz;
		p_N_m = (p_N0.N2 > N_tol ? p_N0.N2 : N_tol) * p_m;
		en2_vm += p_N_m;
		en2_vmx += p_N_m * p_v0.vx;
		en2_vmy += p_N_m * p_v0.vy;
		en2_vmz += p_N_m * p_v0.vz;
		p_N_m = (p_N0.N3 > N_tol ? p_N0.N3 : N_tol) * p_m;
		en3_vm += p_N_m;
		en3_vmx += p_N_m * p_v0.vx;
		en3_vmy += p_N_m * p_v0.vy;
		en3_vmz += p_N_m * p_v0.vz;
		p_N_m = (p_N0.N4 > N_tol ? p_N0.N4 : N_tol) * p_m;
		en4_vm += p_N_m;
		en4_vmx += p_N_m * p_v0.vx;
		en4_vmy += p_N_m * p_v0.vy;
		en4_vmz += p_N_m * p_v0.vz;

		// displacement (for contact)
		Displacement& p_d1 = pcl_disp1[prev_p_id];
		Displacement& p_d0 = pcl_disp0[p_id];
		p_d0.ux = p_d1.ux;
		p_d0.uy = p_d1.uy;
		p_d0.uz = p_d1.uz;

		// cal external load
		const Force& p_bf = pcl_bf[ori_p_id];
		one_fourth_bfx = one_fourth * p_bf.fx;
		one_fourth_bfy = one_fourth * p_bf.fy;
		one_fourth_bfz = one_fourth * p_bf.fz;
		const Force& p_t = pcl_t[ori_p_id];
		en1_fx += one_fourth_bfx + p_N0.N1 * p_t.fx;
		en1_fy += one_fourth_bfy + p_N0.N1 * p_t.fy;
		en1_fz += one_fourth_bfz + p_N0.N1 * p_t.fz;
		en2_fx += one_fourth_bfx + p_N0.N2 * p_t.fx;
		en2_fy += one_fourth_bfy + p_N0.N2 * p_t.fy;
		en2_fz += one_fourth_bfz + p_N0.N2 * p_t.fz;
		en3_fx += one_fourth_bfx + p_N0.N3 * p_t.fx;
		en3_fy += one_fourth_bfy + p_N0.N3 * p_t.fy;
		en3_fz += one_fourth_bfz + p_N0.N3 * p_t.fz;
		en4_fx += one_fourth_bfx + p_N0.N4 * p_t.fx;
		en4_fy += one_fourth_bfy + p_N0.N4 * p_t.fy;
		en4_fz += one_fourth_bfz + p_N0.N4 * p_t.fz;

		if (e_id != pcl_in_elem0[p_id + 1])
		{
			elem_pcl_m[e_id] = e_p_m;
			elem_density[e_id] = e_p_m / e_p_vol;

			ne_id = e_id * 4;
			ElemNodeVM& en1_v = elem_node_vm[ne_id];
			en1_v.vm = en1_vm;
			en1_v.vmx = en1_vmx;
			en1_v.vmy = en1_vmy;
			en1_v.vmz = en1_vmz;
			++ne_id;
			ElemNodeVM& en2_v = elem_node_vm[ne_id];
			en2_v.vm = en2_vm;
			en2_v.vmx = en2_vmx;
			en2_v.vmy = en2_vmy;
			en2_v.vmz = en2_vmz;
			++ne_id;
			ElemNodeVM& en3_v = elem_node_vm[ne_id];
			en3_v.vm = en3_vm;
			en3_v.vmx = en3_vmx;
			en3_v.vmy = en3_vmy;
			en3_v.vmz = en3_vmz;
			++ne_id;
			ElemNodeVM& en4_v = elem_node_vm[ne_id];
			en4_v.vm = en4_vm;
			en4_v.vmx = en4_vmx;
			en4_v.vmy = en4_vmy;
			en4_v.vmz = en4_vmz;

			e_s11 /= e_p_vol;
			e_s22 /= e_p_vol;
			e_s33 /= e_p_vol;
			e_s12 /= e_p_vol;
			e_s23 /= e_p_vol;
			e_s31 /= e_p_vol;
			if (e_p_vol > elem_vol[e_id])
				e_p_vol = elem_vol[e_id];
			DShapeFuncABC& e_dN = elem_dN_abc[e_id];
			// node 1
			ne_id = e_id * 4;
			Force& en1_f = elem_node_force[ne_id];
			en1_fx -= (e_dN.dN1_dx * e_s11 + e_dN.dN1_dy * e_s12 + e_dN.dN1_dz * e_s31) * e_p_vol;
			en1_f.fx = en1_fx;
			en1_fy -= (e_dN.dN1_dx * e_s12 + e_dN.dN1_dy * e_s22 + e_dN.dN1_dz * e_s23) * e_p_vol;
			en1_f.fy = en1_fy;
			en1_fz -= (e_dN.dN1_dx * e_s31 + e_dN.dN1_dy * e_s23 + e_dN.dN1_dz * e_s33) * e_p_vol;
			en1_f.fz = en1_fz;
			// node 2
			++ne_id;
			Force& en2_f = elem_node_force[ne_id];
			en2_fx -= (e_dN.dN2_dx * e_s11 + e_dN.dN2_dy * e_s12 + e_dN.dN2_dz * e_s31) * e_p_vol;
			en2_f.fx = en2_fx;
			en2_fy -= (e_dN.dN2_dx * e_s12 + e_dN.dN2_dy * e_s22 + e_dN.dN2_dz * e_s23) * e_p_vol;
			en2_f.fy = en2_fy;
			en2_fz -= (e_dN.dN2_dx * e_s31 + e_dN.dN2_dy * e_s23 + e_dN.dN2_dz * e_s33) * e_p_vol;
			en2_f.fz = en2_fz;
			// node 3
			++ne_id;
			Force& en3_f = elem_node_force[ne_id];
			en3_fx -= (e_dN.dN3_dx * e_s11 + e_dN.dN3_dy * e_s12 + e_dN.dN3_dz * e_s31) * e_p_vol;
			en3_f.fx = en3_fx;
			en3_fy -= (e_dN.dN3_dx * e_s12 + e_dN.dN3_dy * e_s22 + e_dN.dN3_dz * e_s23) * e_p_vol;
			en3_f.fy = en3_fy;
			en3_fz -= (e_dN.dN3_dx * e_s31 + e_dN.dN3_dy * e_s23 + e_dN.dN3_dz * e_s33) * e_p_vol;
			en3_f.fz = en3_fz;
			// node 4
			++ne_id;
			Force& en4_f = elem_node_force[ne_id];
			en4_fx -= (e_dN.dN4_dx * e_s11 + e_dN.dN4_dy * e_s12 + e_dN.dN4_dz * e_s31) * e_p_vol;
			en4_f.fx = en4_fx;
			en4_fy -= (e_dN.dN4_dx * e_s12 + e_dN.dN4_dy * e_s22 + e_dN.dN4_dz * e_s23) * e_p_vol;
			en4_f.fy = en4_fy;
			en4_fz -= (e_dN.dN4_dx * e_s31 + e_dN.dN4_dy * e_s23 + e_dN.dN4_dz * e_s33) * e_p_vol;
			en4_f.fz = en4_fz;

			ne_id = my_valid_elem_num * 4;
			my_valid_elem_id[my_valid_elem_num++] = e_id;

			ElemNodeIndex& eni = elem_node_id[e_id];
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

			e_p_m = 0.0;
			e_p_vol = 0.0;
			e_s11 = 0.0;
			e_s22 = 0.0;
			e_s33 = 0.0;
			e_s12 = 0.0;
			e_s23 = 0.0;
			e_s31 = 0.0;
			en1_vm = 0.0;
			en1_vmx = 0.0;
			en1_vmy = 0.0;
			en1_vmz = 0.0;
			en2_vm = 0.0;
			en2_vmx = 0.0;
			en2_vmy = 0.0;
			en2_vmz = 0.0;
			en3_vm = 0.0;
			en3_vmx = 0.0;
			en3_vmy = 0.0;
			en3_vmz = 0.0;
			en4_vm = 0.0;
			en4_vmx = 0.0;
			en4_vmy = 0.0;
			en4_vmz = 0.0;
			en1_fx = 0.0;
			en1_fy = 0.0;
			en1_fz = 0.0;
			en2_fx = 0.0;
			en2_fy = 0.0;
			en2_fz = 0.0;
			en3_fx = 0.0;
			en3_fy = 0.0;
			en3_fz = 0.0;
			en4_fx = 0.0;
			en4_fy = 0.0;
			en4_fz = 0.0;
		}
	}

#pragma omp critical
	self.valid_elem_num += my_valid_elem_num;

	// sort node-elem pair according to node id
	size_t ve_id;
	memset(my_cbin, 0, 0x100 * sizeof(size_t));
	for (ve_id = 0; ve_id < my_valid_elem_num; ++ve_id)
	{
		ne_id = ve_id * 4;
		++my_cbin[data_digit(my_node_has_elem[ne_id], 0)];
		++my_cbin[data_digit(my_node_has_elem[++ne_id], 0)];
		++my_cbin[data_digit(my_node_has_elem[++ne_id], 0)];
		++my_cbin[data_digit(my_node_has_elem[++ne_id], 0)];
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
		ne_id = ve_id * 4;
		pos_id = --my_sbin[data_digit(my_node_has_elem[ne_id], 0)];
		node_has_elem0[pos_id] = my_node_has_elem[ne_id];
		node_elem_pair0[pos_id] = my_node_elem_pair[ne_id];
		pos_id = --my_sbin[data_digit(my_node_has_elem[++ne_id], 0)];
		node_has_elem0[pos_id] = my_node_has_elem[ne_id];
		node_elem_pair0[pos_id] = my_node_elem_pair[ne_id];
		pos_id = --my_sbin[data_digit(my_node_has_elem[++ne_id], 0)];
		node_has_elem0[pos_id] = my_node_has_elem[ne_id];
		node_elem_pair0[pos_id] = my_node_elem_pair[ne_id];
		pos_id = --my_sbin[data_digit(my_node_has_elem[++ne_id], 0)];
		node_has_elem0[pos_id] = my_node_has_elem[ne_id];
		node_elem_pair0[pos_id] = my_node_elem_pair[ne_id];
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
	union
	{
		double n_fx;
		size_t in_fx;
	};
	union
	{
		double n_fy;
		size_t in_fy;
	};
	union
	{
		double n_fz;
		size_t in_fz;
	};
	n_fx = 0.0;
	n_fy = 0.0;
	n_fz = 0.0;
	double n_vm = 0.0;
	double n_vmx = 0.0;
	double n_vmy = 0.0;
	double n_vmz = 0.0;
	double f_ub = 0.0;
	double e_kin = 0.0;
	n_id = node_has_elem0[ve_id0];
	double vbc_len;
	assert(n_id < self.node_num);
	for (ve_id = ve_id0; ve_id < ve_id1; ++ve_id)
	{
		ne_id = node_elem_pair0[ve_id];
		assert(ne_id < self.elem_num * 4);
		e_id = ne_id / 4;
		n_am += elem_pcl_m[e_id];
		Force& nf = elem_node_force[ne_id];
		n_fx += nf.fx;
		n_fy += nf.fy;
		n_fz += nf.fz;
		ElemNodeVM& nvm = elem_node_vm[ne_id];
		n_vm += nvm.vm;
		n_vmx += nvm.vmx;
		n_vmy += nvm.vmy;
		n_vmz += nvm.vmz;

		if (n_id != node_has_elem0[ve_id + 1])
		{
			Acceleration& n_a = node_a[n_id];
			n_am *= one_fourth;
			node_am[n_id] = n_am;
			n_a.ax = n_fx / n_am;
			n_a.ay = n_fy / n_am;
			n_a.az = n_fz / n_am;
			Velocity& n_v = node_v[n_id];
			n_v.vx = n_vmx / n_vm + n_a.ax * dt;
			n_v.vy = n_vmy / n_vm + n_a.ay * dt;
			n_v.vz = n_vmz / n_vm + n_a.az * dt;
			NodeVBCVec& n_vbc_vec = node_vbc_vec[n_id];
			vbc_len = n_a.ax * n_vbc_vec.x + n_a.ay * n_vbc_vec.y + n_a.az * n_vbc_vec.z;
			n_a.ax -= vbc_len * n_vbc_vec.x;
			n_a.ay -= vbc_len * n_vbc_vec.y;
			n_a.az -= vbc_len * n_vbc_vec.z;
			vbc_len = n_v.vx * n_vbc_vec.x + n_v.vy * n_vbc_vec.y + n_v.vz * n_vbc_vec.z;
			n_v.vx -= vbc_len * n_vbc_vec.x;
			n_v.vy -= vbc_len * n_vbc_vec.y;
			n_v.vz -= vbc_len * n_vbc_vec.z;
			NodeHasVBC& n_has_vbc = node_has_vbc[n_id];
			bc_mask = size_t(n_has_vbc.has_vx_bc) + SIZE_MAX;
			n_a.iax &= bc_mask;
			n_v.ivx &= bc_mask;
			in_fx &= bc_mask;
			bc_mask = size_t(n_has_vbc.has_vy_bc) + SIZE_MAX;
			n_a.iay &= bc_mask;
			n_v.ivy &= bc_mask;
			in_fy &= bc_mask;
			bc_mask = size_t(n_has_vbc.has_vz_bc) + SIZE_MAX;
			n_a.iaz &= bc_mask;
			n_v.ivz &= bc_mask;
			in_fz &= bc_mask;

			f_ub += n_fx * n_fx + n_fy * n_fy + n_fz * n_fz;
			e_kin += n_am * (n_v.vx * n_v.vx + n_v.vy * n_v.vy + n_v.vz * n_v.vz);
			
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

#pragma omp critical
	{
		self.f_ub += f_ub;
		self.e_kin += e_kin;
	}
	
#pragma omp barrier

#pragma omp master
	{	
#ifdef _DEBUG
		self.prev_valid_pcl_num_tmp = self.prev_valid_pcl_num;
#endif
		self.prev_valid_pcl_num = self.valid_pcl_num;
		self.valid_pcl_num = 0;
		
		if (self.init_f_ub_is_init)
		{
			self.f_ub_ratio = sqrt(self.f_ub / self.init_f_ub);
		}
		else
		{
			self.init_f_ub_is_init = true;
			self.init_f_ub = self.f_ub;
			self.f_ub_ratio = 1.0;
		}

		if (self.e_kin_max_is_init)
			self.e_kin_ratio = sqrt(self.e_kin / self.e_kin_max);
		else
			self.e_kin_ratio = 1.0;

		//if (self.substep_index % 100 == 1)
		//	res_file_T3D_ME_mt_Geo << self.substep_index << ", "
		//		<< self.f_ub_ratio << ", "
		//		<< self.e_kin_ratio << "\n";

		if (self.e_kin < self.e_kin_prev)
		{
			if (!self.e_kin_max_is_init)
			{
				self.e_kin_max_is_init = true;
				self.e_kin_max = self.e_kin_prev;
			}
			self.e_kin_prev = 0.0;
			self.e_kin_ratio = 0.0;

			//if (self.f_ub_ratio < self.f_ub_ratio_limit &&
			//	self.e_kin_ratio < self.e_kin_ratio_limit)
			//{
			//	// exit
			//	self.cal_status = 2;
			//}
			//else
			//{
			//	// reset velocity
			//	self.cal_status = 1;
			//}
			self.cal_status = 1;
		}
		else // continue calculation
		{
			self.cal_status = 0;
			self.e_kin_prev = self.e_kin;
		}
	}

#pragma omp barrier
	if (self.cal_status == 1) // reset velocity
	{
		for (p_id = p_id0; p_id < p_id1; ++p_id)
		{
			Velocity& p_v0 = pcl_v0[p_id];
			p_v0.vx = 0.0;
			p_v0.vy = 0.0;
			p_v0.vz = 0.0;
		}

		n_id = node_has_elem0[ve_id0];
		assert(n_id < self.node_num);
		for (ve_id = ve_id0; ve_id < ve_id1; ++ve_id)
		{
			if (n_id != node_has_elem0[ve_id + 1])
			{
				Acceleration& n_a = node_a[n_id];
				n_a.ax = 0.0;
				n_a.ay = 0.0;
				n_a.az = 0.0;
				Velocity& n_v = node_v[n_id];
				n_v.vx = 0.0;
				n_v.vy = 0.0;
				n_v.vz = 0.0;
				n_id = node_has_elem0[ve_id + 1];
				assert(n_id < self.node_num || n_id == SIZE_MAX);
			}
		}
	}
	else if (self.cal_status == 2) // already converge
	{
#pragma omp master
		self.exit_calculation();
		return 1;
	}

	// cal element strain and "enhancement"
	double e_de_vol;
	for (ve_id = 0; ve_id < my_valid_elem_num; ++ve_id)
	{
		e_id = my_valid_elem_id[ve_id];
		assert(e_id < self.elem_num);

		ElemNodeIndex& eni = elem_node_id[e_id];
		Velocity& n_v1 = node_v[eni.n1];
		Velocity& n_v2 = node_v[eni.n2];
		Velocity& n_v3 = node_v[eni.n3];
		Velocity& n_v4 = node_v[eni.n4];
		DShapeFuncABC& e_dN = elem_dN_abc[e_id];
		StrainInc& e_de = elem_de[e_id];
		e_de.de11 = (e_dN.dN1_dx * n_v1.vx + e_dN.dN2_dx * n_v2.vx + e_dN.dN3_dx * n_v3.vx + e_dN.dN4_dx * n_v4.vx) * dt;
		e_de.de22 = (e_dN.dN1_dy * n_v1.vy + e_dN.dN2_dy * n_v2.vy + e_dN.dN3_dy * n_v3.vy + e_dN.dN4_dy * n_v4.vy) * dt;
		e_de.de33 = (e_dN.dN1_dz * n_v1.vz + e_dN.dN2_dz * n_v2.vz + e_dN.dN3_dz * n_v3.vz + e_dN.dN4_dz * n_v4.vz) * dt;
		e_de.de12 = (e_dN.dN1_dx * n_v1.vy + e_dN.dN2_dx * n_v2.vy + e_dN.dN3_dx * n_v3.vy + e_dN.dN4_dx * n_v4.vy
			+ e_dN.dN1_dy * n_v1.vx + e_dN.dN2_dy * n_v2.vx + e_dN.dN3_dy * n_v3.vx + e_dN.dN4_dy * n_v4.vx) * dt * 0.5;
		e_de.de23 = (e_dN.dN1_dy * n_v1.vz + e_dN.dN2_dy * n_v2.vz + e_dN.dN3_dy * n_v3.vz + e_dN.dN4_dy * n_v4.vz
			+ e_dN.dN1_dz * n_v1.vy + e_dN.dN2_dz * n_v2.vy + e_dN.dN3_dz * n_v3.vy + e_dN.dN4_dz * n_v4.vy) * dt * 0.5;
		e_de.de31 = (e_dN.dN1_dz * n_v1.vx + e_dN.dN2_dz * n_v2.vx + e_dN.dN3_dz * n_v3.vx + e_dN.dN4_dz * n_v4.vx
			+ e_dN.dN1_dx * n_v1.vz + e_dN.dN2_dx * n_v2.vz + e_dN.dN3_dx * n_v3.vz + e_dN.dN4_dx * n_v4.vz) * dt * 0.5;
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

	Acceleration* pn_a1, * pn_a2, * pn_a3, * pn_a4;
	StrainInc* pe_de, e_de_tmp;
	const double* estrain, * pstrain, * dstress;
	size_t p_e_id, pcl_in_mesh_num = 0;
	e_id = SIZE_MAX;
	thd.sorted_pcl_in_elem_id ^= 1;
	for (p_id = p_id0; p_id < p_id1; ++p_id)
	{
		if (e_id != pcl_in_elem0[p_id])
		{
			e_id = pcl_in_elem0[p_id];
			assert(e_id < self.elem_num);

			ElemNodeIndex& eni = elem_node_id[e_id];
			pn_a1 = node_a + eni.n1;
			pn_a2 = node_a + eni.n2;
			pn_a3 = node_a + eni.n3;
			pn_a4 = node_a + eni.n4;

			e_de_vol = one_fourth *
				(node_de_vol[eni.n1] + node_de_vol[eni.n2]
					+ node_de_vol[eni.n3] + node_de_vol[eni.n4]);

			pe_de = elem_de + e_id;
			e_de_vol *= one_third;
			pe_de->de11 += e_de_vol;
			pe_de->de22 += e_de_vol;
			pe_de->de33 += e_de_vol;
		}

		// update velocity
		ShapeFunc& p_N = pcl_N0[p_id];
		Velocity& p_v0 = pcl_v0[p_id];
		p_v0.vx += (p_N.N1 * pn_a1->ax + p_N.N2 * pn_a2->ax + p_N.N3 * pn_a3->ax + p_N.N4 * pn_a4->ax) * dt;
		p_v0.vy += (p_N.N1 * pn_a1->ay + p_N.N2 * pn_a2->ay + p_N.N3 * pn_a3->ay + p_N.N4 * pn_a4->ay) * dt;
		p_v0.vz += (p_N.N1 * pn_a1->az + p_N.N2 * pn_a2->az + p_N.N3 * pn_a3->az + p_N.N4 * pn_a4->az) * dt;

		// update location (in which element)
		ori_p_id = pcl_index0[p_id];
		assert(ori_p_id < md.ori_pcl_num);
		p_e_id = e_id;
		++pcl_in_mesh_num;
		pcl_in_elem1[p_id] = p_e_id;
		prev_pcl_id1[p_id] = p_id;
		assert(p_e_id < self.elem_num || p_e_id == SIZE_MAX);

		// update density
		pcl_density0[p_id] = elem_density[e_id];
		
		// update stress
		e_de_tmp.de11 = pe_de->de11;
		e_de_tmp.de22 = pe_de->de22;
		e_de_tmp.de33 = pe_de->de33;
		e_de_tmp.de12 = pe_de->de12;
		e_de_tmp.de23 = pe_de->de23;
		e_de_tmp.de31 = pe_de->de31;
		MatModel::MaterialModel& pcl_mm = *pcl_mat_model[ori_p_id];
		pcl_mm.integrate(e_de_tmp.de);
		dstress = pcl_mm.get_dstress();

		Stress& p_s = pcl_stress0[p_id];
		p_s.s11 += dstress[0];
		p_s.s22 += dstress[1];
		p_s.s33 += dstress[2];
		p_s.s12 += dstress[3];
		p_s.s23 += dstress[4];
		p_s.s31 += dstress[5];
				
		prev_p_id = prev_pcl_id0[p_id];
#ifdef _DEBUG
		assert(prev_p_id < self.prev_valid_pcl_num_tmp);
#endif
	}

#pragma omp critical
	self.valid_pcl_num += pcl_in_mesh_num;

#pragma omp master
	{
		pcl_in_elem0[self.prev_valid_pcl_num] = SIZE_MAX;
		pcl_in_elem1[self.prev_valid_pcl_num] = SIZE_MAX;
		self.valid_elem_num = 0;
		self.continue_calculation();
	}

#pragma omp barrier
	return 0;
}
