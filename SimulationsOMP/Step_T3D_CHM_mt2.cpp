#include "SimulationsOMP_pcp.h"

#include <fstream>
#include <iostream>
#include <omp.h>

#include "Step_T3D_CHM_mt.h"

#define one_fourth (0.25)
#define one_third (1.0/3.0)
#define N_min (1.0e-10)
#define Block_Low(th_id, th_num, data_num) ((th_id)*(data_num)/(th_num))

int substep_func_omp_T3D_CHM_mt2(
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
	typedef Step_T3D_CHM_mt::ThreadData ThreadData;

	Step_T3D_CHM_mt& self = *(Step_T3D_CHM_mt*)(_self);

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
	double* const pcl_vol = self.pcl_vol;
	MatModel::MaterialModel** const pcl_mat_model = self.pcl_mat_model;

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
	double* const elem_pcl_m_s = self.elem_pcl_m_s;
	double* const elem_pcl_m_f = self.elem_pcl_m_f;
	StrainInc* const elem_de = self.elem_de;
	double* const elem_p = self.elem_p;
	//double* const elem_n2_miu_div_k_vol = self.elem_n2_miu_div_k_vol;
	//Force* const elem_seep_force = self.elem_seep_force;
	double* const elem_m_de_vol_s = self.elem_m_de_vol_s;
	double* const elem_m_de_vol_f = self.elem_m_de_vol_f;

	ElemNodeVM* const elem_node_vm_s = self.elem_node_vm_s;
	ElemNodeVM* const elem_node_vm_f = self.elem_node_vm_f;
	Force* const elem_node_force_s = self.elem_node_force_s;
	Force* const elem_node_force_f = self.elem_node_force_f;

	Acceleration* const node_a_s = self.node_a_s;
	Acceleration* const node_a_f = self.node_a_f;
	Velocity* const node_v_s = self.node_v_s;
	Velocity* const node_v_f = self.node_v_f;
	NodeHasVBC* const node_has_vbc_s = self.node_has_vbc_s;
	NodeHasVBC* const node_has_vbc_f = self.node_has_vbc_f;
	double* const node_am_s = self.node_am_s;
	double* const node_am_f = self.node_am_f;
	double* const node_de_vol_s = self.node_de_vol_s;
	double* const node_de_vol_f = self.node_de_vol_f;

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
	double one_fourth_bfx_s, one_fourth_bfx_f;
	double one_fourth_bfy_s, one_fourth_bfy_f;
	double one_fourth_bfz_s, one_fourth_bfz_f;
	double en1_vm_s = 0.0;
	double en1_vmx_s = 0.0;
	double en1_vmy_s = 0.0;
	double en1_vmz_s = 0.0;
	double en1_vm_f = 0.0;
	double en1_vmx_f = 0.0;
	double en1_vmy_f = 0.0;
	double en1_vmz_f = 0.0;
	double en2_vm_s = 0.0;
	double en2_vmx_s = 0.0;
	double en2_vmy_s = 0.0;
	double en2_vmz_s = 0.0;
	double en2_vm_f = 0.0;
	double en2_vmx_f = 0.0;
	double en2_vmy_f = 0.0;
	double en2_vmz_f = 0.0;
	double en3_vm_s = 0.0;
	double en3_vmx_s = 0.0;
	double en3_vmy_s = 0.0;
	double en3_vmz_s = 0.0;
	double en3_vm_f = 0.0;
	double en3_vmx_f = 0.0;
	double en3_vmy_f = 0.0;
	double en3_vmz_f = 0.0;
	double en4_vm_s = 0.0;
	double en4_vmx_s = 0.0;
	double en4_vmy_s = 0.0;
	double en4_vmz_s = 0.0;
	double en4_vm_f = 0.0;
	double en4_vmx_f = 0.0;
	double en4_vmy_f = 0.0;
	double en4_vmz_f = 0.0;
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
	double e_fx_seep, e_fy_seep, e_fz_seep;
	double en1_fx_seep = 0.0;
	double en1_fy_seep = 0.0;
	double en1_fz_seep = 0.0;
	double en2_fx_seep = 0.0;
	double en2_fy_seep = 0.0;
	double en2_fz_seep = 0.0;
	double en3_fx_seep = 0.0;
	double en3_fy_seep = 0.0;
	double en3_fz_seep = 0.0;
	double en4_fx_seep = 0.0;
	double en4_fy_seep = 0.0;
	double en4_fz_seep = 0.0;
	double en1_fx_s = 0.0;
	double en1_fy_s = 0.0;
	double en1_fz_s = 0.0;
	double en2_fx_s = 0.0;
	double en2_fy_s = 0.0;
	double en2_fz_s = 0.0;
	double en3_fx_s = 0.0;
	double en3_fy_s = 0.0;
	double en3_fz_s = 0.0;
	double en4_fx_s = 0.0;
	double en4_fy_s = 0.0;
	double en4_fz_s = 0.0;
	double en1_fx_f = 0.0;
	double en1_fy_f = 0.0;
	double en1_fz_f = 0.0;
	double en2_fx_f = 0.0;
	double en2_fy_f = 0.0;
	double en2_fz_f = 0.0;
	double en3_fx_f = 0.0;
	double en3_fy_f = 0.0;
	double en3_fz_f = 0.0;
	double en4_fx_f = 0.0;
	double en4_fy_f = 0.0;
	double en4_fz_f = 0.0;
	e_id = pcl_in_elem0[p_id0];
	size_t* const my_valid_elem_id = self.valid_elem_id + e_id;
	size_t* const my_node_has_elem = node_has_elem1 + e_id * 4;
	size_t* const my_node_elem_pair = node_elem_pair1 + e_id * 4;
	size_t my_valid_elem_num = 0;
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
		p_N_m = (p_N0.N1 > N_tol ? p_N0.N1 : N_tol) * p_m_s;
		en1_vm_s += p_N_m;
		en1_vmx_s += p_N_m * p_v_s0.vx;
		en1_vmy_s += p_N_m * p_v_s0.vy;
		en1_vmz_s += p_N_m * p_v_s0.vz;
		p_N_m = (p_N0.N2 > N_tol ? p_N0.N2 : N_tol) * p_m_s;
		en2_vm_s += p_N_m;
		en2_vmx_s += p_N_m * p_v_s0.vx;
		en2_vmy_s += p_N_m * p_v_s0.vy;
		en2_vmz_s += p_N_m * p_v_s0.vz;
		p_N_m = (p_N0.N3 > N_tol ? p_N0.N3 : N_tol) * p_m_s;
		en3_vm_s += p_N_m;
		en3_vmx_s += p_N_m * p_v_s0.vx;
		en3_vmy_s += p_N_m * p_v_s0.vy;
		en3_vmz_s += p_N_m * p_v_s0.vz;
		p_N_m = (p_N0.N4 > N_tol ? p_N0.N4 : N_tol) * p_m_s;
		en4_vm_s += p_N_m;
		en4_vmx_s += p_N_m * p_v_s0.vx;
		en4_vmy_s += p_N_m * p_v_s0.vy;
		en4_vmz_s += p_N_m * p_v_s0.vz;
		// fluid phase
		Velocity& p_v_f1 = pcl_v_f1[prev_p_id];
		Velocity& p_v_f0 = pcl_v_f0[p_id];
		p_v_f0.vx = p_v_f1.vx;
		p_v_f0.vy = p_v_f1.vy;
		p_v_f0.vz = p_v_f1.vz;
		p_N_m = (p_N0.N1 > N_tol ? p_N0.N1 : N_tol) * p_m_f;
		en1_vm_f += p_N_m;
		en1_vmx_f += p_N_m * p_v_f0.vx;
		en1_vmy_f += p_N_m * p_v_f0.vy;
		en1_vmz_f += p_N_m * p_v_f0.vz;
		p_N_m = (p_N0.N2 > N_tol ? p_N0.N2 : N_tol) * p_m_f;
		en2_vm_f += p_N_m;
		en2_vmx_f += p_N_m * p_v_f0.vx;
		en2_vmy_f += p_N_m * p_v_f0.vy;
		en2_vmz_f += p_N_m * p_v_f0.vz;
		p_N_m = (p_N0.N3 > N_tol ? p_N0.N3 : N_tol) * p_m_f;
		en3_vm_f += p_N_m;
		en3_vmx_f += p_N_m * p_v_f0.vx;
		en3_vmy_f += p_N_m * p_v_f0.vy;
		en3_vmz_f += p_N_m * p_v_f0.vz;
		p_N_m = (p_N0.N4 > N_tol ? p_N0.N4 : N_tol) * p_m_f;
		en4_vm_f += p_N_m;
		en4_vmx_f += p_N_m * p_v_f0.vx;
		en4_vmy_f += p_N_m * p_v_f0.vy;
		en4_vmz_f += p_N_m * p_v_f0.vz;

		// displacement (for contact)
		Displacement& p_u_s1 = pcl_u_s1[prev_p_id];
		Displacement& p_u_s0 = pcl_u_s0[p_id];
		p_u_s0.ux = p_u_s1.ux;
		p_u_s0.uy = p_u_s1.uy;
		p_u_s0.uz = p_u_s1.uz;
		Displacement& p_u_f1 = pcl_u_f1[prev_p_id];
		Displacement& p_u_f0 = pcl_u_f0[p_id];
		p_u_f0.ux = p_u_f1.ux;
		p_u_f0.uy = p_u_f1.uy;
		p_u_f0.uz = p_u_f1.uz;

		// seepage force
		const double f_seep_tmp = md.miu / md.k * p_n * p_n * p_vol;
		e_fx_seep = (p_v_f0.vx - p_v_s0.vx) * f_seep_tmp;
		en1_fx_seep += p_N0.N1 * e_fx_seep;
		en2_fx_seep += p_N0.N2 * e_fx_seep;
		en3_fx_seep += p_N0.N3 * e_fx_seep;
		en4_fx_seep += p_N0.N4 * e_fx_seep;
		e_fy_seep = (p_v_f0.vy - p_v_s0.vy) * f_seep_tmp;
		en1_fy_seep += p_N0.N1 * e_fy_seep;
		en2_fy_seep += p_N0.N2 * e_fy_seep;
		en3_fy_seep += p_N0.N3 * e_fy_seep;
		en4_fy_seep += p_N0.N4 * e_fy_seep;
		e_fz_seep = (p_v_f0.vz - p_v_s0.vz) * f_seep_tmp;
		en1_fz_seep += p_N0.N1 * e_fz_seep;
		en2_fz_seep += p_N0.N2 * e_fz_seep;
		en3_fz_seep += p_N0.N3 * e_fz_seep;
		en4_fz_seep += p_N0.N4 * e_fz_seep;

		// solid external load
		const Force& p_bf_s = self.pcl_bf_s[ori_p_id];
		one_fourth_bfx_s = one_fourth * p_bf_s.fx;
		one_fourth_bfy_s = one_fourth * p_bf_s.fy;
		one_fourth_bfz_s = one_fourth * p_bf_s.fz;
		const Force& p_t = self.pcl_t[ori_p_id];
		en1_fx_s += one_fourth_bfx_s + p_N0.N1 * p_t.fx;
		en1_fy_s += one_fourth_bfy_s + p_N0.N1 * p_t.fy;
		en1_fz_s += one_fourth_bfz_s + p_N0.N1 * p_t.fz;
		en2_fx_s += one_fourth_bfx_s + p_N0.N2 * p_t.fx;
		en2_fy_s += one_fourth_bfy_s + p_N0.N2 * p_t.fy;
		en2_fz_s += one_fourth_bfz_s + p_N0.N2 * p_t.fz;
		en3_fx_s += one_fourth_bfx_s + p_N0.N3 * p_t.fx;
		en3_fy_s += one_fourth_bfy_s + p_N0.N3 * p_t.fy;
		en3_fz_s += one_fourth_bfz_s + p_N0.N3 * p_t.fz;
		en4_fx_s += one_fourth_bfx_s + p_N0.N4 * p_t.fx;
		en4_fy_s += one_fourth_bfy_s + p_N0.N4 * p_t.fy;
		en4_fz_s += one_fourth_bfz_s + p_N0.N4 * p_t.fz;

		// fluid external load
		const Force& p_bf_f = self.pcl_bf_f[ori_p_id];
		one_fourth_bfx_f = one_fourth * p_bf_f.fx;
		one_fourth_bfy_f = one_fourth * p_bf_f.fy;
		one_fourth_bfz_f = one_fourth * p_bf_f.fz;
		en1_fx_f += one_fourth_bfx_f;
		en1_fy_f += one_fourth_bfy_f;
		en1_fz_f += one_fourth_bfz_f;
		en2_fx_f += one_fourth_bfx_f;
		en2_fy_f += one_fourth_bfy_f;
		en2_fz_f += one_fourth_bfz_f;
		en3_fx_f += one_fourth_bfx_f;
		en3_fy_f += one_fourth_bfy_f;
		en3_fz_f += one_fourth_bfz_f;
		en4_fx_f += one_fourth_bfx_f;
		en4_fy_f += one_fourth_bfy_f;
		en4_fz_f += one_fourth_bfz_f;

		if (e_id != pcl_in_elem0[p_id + 1])
		{
			// v_s
			ElemNodeVM& en1_v_s = elem_node_vm_s[e_id * 4];
			en1_v_s.vm = en1_vm_s;
			en1_v_s.vmx = en1_vmx_s;
			en1_v_s.vmy = en1_vmy_s;
			en1_v_s.vmz = en1_vmz_s;
			ElemNodeVM& en2_v_s = elem_node_vm_s[e_id * 4 + 1];
			en2_v_s.vm = en2_vm_s;
			en2_v_s.vmx = en2_vmx_s;
			en2_v_s.vmy = en2_vmy_s;
			en2_v_s.vmz = en2_vmz_s;
			ElemNodeVM& en3_v_s = elem_node_vm_s[e_id * 4 + 2];
			en3_v_s.vm = en3_vm_s;
			en3_v_s.vmx = en3_vmx_s;
			en3_v_s.vmy = en3_vmy_s;
			en3_v_s.vmz = en3_vmz_s;
			ElemNodeVM& en4_v_s = elem_node_vm_s[e_id * 4 + 3];
			en4_v_s.vm = en4_vm_s;
			en4_v_s.vmx = en4_vmx_s;
			en4_v_s.vmy = en4_vmy_s;
			en4_v_s.vmz = en4_vmz_s;

			// v_f
			ElemNodeVM& en1_v_f = elem_node_vm_f[e_id * 4];
			en1_v_f.vm = en1_vm_f;
			en1_v_f.vmx = en1_vmx_f;
			en1_v_f.vmy = en1_vmy_f;
			en1_v_f.vmz = en1_vmz_f;
			ElemNodeVM& en2_v_f = elem_node_vm_f[e_id * 4 + 1];
			en2_v_f.vm = en2_vm_f;
			en2_v_f.vmx = en2_vmx_f;
			en2_v_f.vmy = en2_vmy_f;
			en2_v_f.vmz = en2_vmz_f;
			ElemNodeVM& en3_v_f = elem_node_vm_f[e_id * 4 + 2];
			en3_v_f.vm = en3_vm_f;
			en3_v_f.vmx = en3_vmx_f;
			en3_v_f.vmy = en3_vmy_f;
			en3_v_f.vmz = en3_vmz_f;
			ElemNodeVM& en4_v_f = elem_node_vm_f[e_id * 4 + 3];
			en4_v_f.vm = en4_vm_f;
			en4_v_f.vmx = en4_vmx_f;
			en4_v_f.vmy = en4_vmy_f;
			en4_v_f.vmz = en4_vmz_f;

			elem_pcl_m_s[e_id] = e_p_m_s;
			elem_pcl_m_f[e_id] = e_p_m_f;
			e_n = 1.0 - e_n / e_p_vol;
			elem_pcl_n[e_id] = e_n;
			elem_density_f[e_id] = e_p_m_f / e_p_vol_f;

			e_s11 /= e_p_vol;
			e_s22 /= e_p_vol;
			e_s33 /= e_p_vol;
			e_s12 /= e_p_vol;
			e_s23 /= e_p_vol;
			e_s31 /= e_p_vol;
			e_p /= e_p_vol;
			elem_p[e_id] = e_p;
			if (e_p_vol > elem_vol[e_id])
			{
				const double vol_ratio = elem_vol[e_id] / e_p_vol;
				en1_fx_seep *= vol_ratio;
				en2_fx_seep *= vol_ratio;
				en3_fx_seep *= vol_ratio;
				en4_fx_seep *= vol_ratio;
				en1_fy_seep *= vol_ratio;
				en2_fy_seep *= vol_ratio;
				en3_fy_seep *= vol_ratio;
				en4_fy_seep *= vol_ratio;
				en1_fz_seep *= vol_ratio;
				en2_fz_seep *= vol_ratio;
				en3_fz_seep *= vol_ratio;
				en4_fz_seep *= vol_ratio;
				e_p_vol = elem_vol[e_id];
			}

			const DShapeFuncABC& e_dN = elem_N_abc[e_id];
			// node 1
			Force& en1_f_s = elem_node_force_s[e_id * 4];
			en1_fx_s += en1_fx_seep;
			en1_fx_s -= (e_dN.dN1_dx * (e_s11 - (1.0 - e_n) * e_p) + e_dN.dN1_dy * e_s12 + e_dN.dN1_dz * e_s31) * e_p_vol;
			en1_f_s.fx = en1_fx_s;
			en1_fy_s += en1_fy_seep;
			en1_fy_s -= (e_dN.dN1_dx * e_s12 + e_dN.dN1_dy * (e_s22 - (1.0 - e_n) * e_p) + e_dN.dN1_dz * e_s23) * e_p_vol;
			en1_f_s.fy = en1_fy_s;
			en1_fz_s += en1_fz_seep;
			en1_fz_s -= (e_dN.dN1_dx * e_s31 + e_dN.dN1_dy * e_s23 + e_dN.dN1_dz * (e_s33 - (1.0 - e_n) * e_p)) * e_p_vol;
			en1_f_s.fz = en1_fz_s;
			// node 2
			Force& en2_f_s = elem_node_force_s[e_id * 4 + 1];
			en2_fx_s += en2_fx_seep;
			en2_fx_s -= (e_dN.dN2_dx * (e_s11 - (1.0 - e_n) * e_p) + e_dN.dN2_dy * e_s12 + e_dN.dN2_dz * e_s31) * e_p_vol;
			en2_f_s.fx = en2_fx_s;
			en2_fy_s += en2_fy_seep;
			en2_fy_s -= (e_dN.dN2_dx * e_s12 + e_dN.dN2_dy * (e_s22 - (1.0 - e_n) * e_p) + e_dN.dN2_dz * e_s23) * e_p_vol;
			en2_f_s.fy = en2_fy_s;
			en2_fz_s += en2_fz_seep;
			en2_fz_s -= (e_dN.dN2_dx * e_s31 + e_dN.dN2_dy * e_s23 + e_dN.dN2_dz * (e_s33 - (1.0 - e_n) * e_p)) * e_p_vol;
			en2_f_s.fz = en2_fz_s;
			// node 3
			Force& en3_f_s = elem_node_force_s[e_id * 4 + 2];
			en3_fx_s += en3_fx_seep;
			en3_fx_s -= (e_dN.dN3_dx * (e_s11 - (1.0 - e_n) * e_p) + e_dN.dN3_dy * e_s12 + e_dN.dN3_dz * e_s31) * e_p_vol;
			en3_f_s.fx = en3_fx_s;
			en3_fy_s += en3_fy_seep;
			en3_fy_s -= (e_dN.dN3_dx * e_s12 + e_dN.dN3_dy * (e_s22 - (1.0 - e_n) * e_p) + e_dN.dN3_dz * e_s23) * e_p_vol;
			en3_f_s.fy = en3_fy_s;
			en3_fz_s += en3_fz_seep;
			en3_fz_s -= (e_dN.dN3_dx * e_s31 + e_dN.dN3_dy * e_s23 + e_dN.dN3_dz * (e_s33 - (1.0 - e_n) * e_p)) * e_p_vol;
			en3_f_s.fz = en3_fz_s;
			// node 4
			Force& en4_f_s = elem_node_force_s[e_id * 4 + 3];
			en4_fx_s += en4_fx_seep;
			en4_fx_s -= (e_dN.dN4_dx * (e_s11 - (1.0 - e_n) * e_p) + e_dN.dN4_dy * e_s12 + e_dN.dN4_dz * e_s31) * e_p_vol;
			en4_f_s.fx = en4_fx_s;
			en4_fy_s += en4_fy_seep;
			en4_fy_s -= (e_dN.dN4_dx * e_s12 + e_dN.dN4_dy * (e_s22 - (1.0 - e_n) * e_p) + e_dN.dN4_dz * e_s23) * e_p_vol;
			en4_f_s.fy = en4_fy_s;
			en4_fz_s += en4_fz_seep;
			en4_fz_s -= (e_dN.dN4_dx * e_s31 + e_dN.dN4_dy * e_s23 + e_dN.dN4_dz * (e_s33 - (1.0 - e_n) * e_p)) * e_p_vol;
			en4_f_s.fz = en4_fz_s;
			// node 1
			Force& en1_f_f = elem_node_force_f[e_id * 4];
			en1_fx_f -= en1_fx_seep;
			en1_fx_f -= e_dN.dN1_dx * e_n * -e_p * e_p_vol;
			en1_f_f.fx = en1_fx_f;
			en1_fy_f -= en1_fy_seep;
			en1_fy_f -= e_dN.dN1_dy * e_n * -e_p * e_p_vol;
			en1_f_f.fy = en1_fy_f;
			en1_fz_f -= en1_fz_seep;
			en1_fz_f -= e_dN.dN1_dz * e_n * -e_p * e_p_vol;
			en1_f_f.fz = en1_fz_f;
			// node 2
			Force& en2_f_f = elem_node_force_f[e_id * 4 + 1];
			en2_fx_f -= en2_fx_seep;
			en2_fx_f -= e_dN.dN2_dx * e_n * -e_p * e_p_vol;
			en2_f_f.fx = en2_fx_f;
			en2_fy_f -= en2_fy_seep;
			en2_fy_f -= e_dN.dN2_dy * e_n * -e_p * e_p_vol;
			en2_f_f.fy = en2_fy_f;
			en2_fz_f -= en2_fz_seep;
			en2_fz_f -= e_dN.dN2_dz * e_n * -e_p * e_p_vol;
			en2_f_f.fz = en2_fz_f;
			// node 3
			Force& en3_f_f = elem_node_force_f[e_id * 4 + 2];
			en3_fx_f -= en3_fx_seep;
			en3_fx_f -= e_dN.dN3_dx * e_n * -e_p * e_p_vol;
			en3_f_f.fx = en3_fx_f;
			en3_fy_f -= en3_fy_seep;
			en3_fy_f -= e_dN.dN3_dy * e_n * -e_p * e_p_vol;
			en3_f_f.fy = en3_fy_f;
			en3_fz_f -= en3_fz_seep;
			en3_fz_f -= e_dN.dN3_dz * e_n * -e_p * e_p_vol;
			en3_f_f.fz = en3_fz_f;
			// node 4
			Force& en4_f_f = elem_node_force_f[e_id * 4 + 3];
			en4_fx_f -= en4_fx_seep;
			en4_fx_f -= e_dN.dN4_dx * e_n * -e_p * e_p_vol;
			en4_f_f.fx = en4_fx_f;
			en4_fy_f -= en4_fy_seep;
			en4_fy_f -= e_dN.dN4_dy * e_n * -e_p * e_p_vol;
			en4_f_f.fy = en4_fy_f;
			en4_fz_f -= en4_fz_seep;
			en4_fz_f -= e_dN.dN4_dz * e_n * -e_p * e_p_vol;
			en4_f_f.fz = en4_fz_f;

			ne_id = my_valid_elem_num * 4;
			my_valid_elem_id[my_valid_elem_num++] = e_id;

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
			en1_vm_f = 0.0;
			en1_vmx_f = 0.0;
			en1_vmy_f = 0.0;
			en1_vmz_f = 0.0;
			en2_vm_f = 0.0;
			en2_vmx_f = 0.0;
			en2_vmy_f = 0.0;
			en2_vmz_f = 0.0;
			en3_vm_f = 0.0;
			en3_vmx_f = 0.0;
			en3_vmy_f = 0.0;
			en3_vmz_f = 0.0;
			en4_vm_f = 0.0;
			en4_vmx_f = 0.0;
			en4_vmy_f = 0.0;
			en4_vmz_f = 0.0;
			en1_fx_seep = 0.0;
			en1_fy_seep = 0.0;
			en1_fz_seep = 0.0;
			en2_fx_seep = 0.0;
			en2_fy_seep = 0.0;
			en2_fz_seep = 0.0;
			en3_fx_seep = 0.0;
			en3_fy_seep = 0.0;
			en3_fz_seep = 0.0;
			en4_fx_seep = 0.0;
			en4_fy_seep = 0.0;
			en4_fz_seep = 0.0;
			en1_fx_s = 0.0;
			en1_fy_s = 0.0;
			en1_fz_s = 0.0;
			en2_fx_s = 0.0;
			en2_fy_s = 0.0;
			en2_fz_s = 0.0;
			en3_fx_s = 0.0;
			en3_fy_s = 0.0;
			en3_fz_s = 0.0;
			en4_fx_s = 0.0;
			en4_fy_s = 0.0;
			en4_fz_s = 0.0;
			en1_fx_f = 0.0;
			en1_fy_f = 0.0;
			en1_fz_f = 0.0;
			en2_fx_f = 0.0;
			en2_fy_f = 0.0;
			en2_fz_f = 0.0;
			en3_fx_f = 0.0;
			en3_fy_f = 0.0;
			en3_fz_f = 0.0;
			en4_fx_f = 0.0;
			en4_fy_f = 0.0;
			en4_fz_f = 0.0;
		}
	}

#pragma omp critical
	self.valid_elem_num += my_valid_elem_num;

	if (md.has_rigid_cylinder())
	{
		Force3D rc_force;
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
		Force3D rc_force;
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
	double n_vm_s = 0.0;
	double n_vmx_s = 0.0;
	double n_vmy_s = 0.0;
	double n_vmz_s = 0.0;
	double n_vm_f = 0.0;
	double n_vmx_f = 0.0;
	double n_vmy_f = 0.0;
	double n_vmz_f = 0.0;
	double n_am_s = 0.0;
	double n_am_f = 0.0;
	double n_fx_s = 0.0;
	double n_fy_s = 0.0;
	double n_fz_s = 0.0;
	double n_fx_f = 0.0;
	double n_fy_f = 0.0;
	double n_fz_f = 0.0;
	n_id = node_has_elem0[ve_id0];
	assert(n_id < self.node_num);
	for (ve_id = ve_id0; ve_id < ve_id1; ++ve_id)
	{
		ne_id = node_elem_pair0[ve_id];
		assert(ne_id < self.elem_num * 4);

		e_id = ne_id / 4;
		n_am_s += elem_pcl_m_s[e_id];
		n_am_f += elem_pcl_m_f[e_id];
		const Force& nf_s = elem_node_force_s[ne_id];
		n_fx_s += nf_s.fx;
		n_fy_s += nf_s.fy;
		n_fz_s += nf_s.fz;
		const Force& nf_f = elem_node_force_f[ne_id];
		n_fx_f += nf_f.fx;
		n_fy_f += nf_f.fy;
		n_fz_f += nf_f.fz;

		ElemNodeVM& nvm_s = elem_node_vm_s[ne_id];
		n_vm_s += nvm_s.vm;
		n_vmx_s += nvm_s.vmx;
		n_vmy_s += nvm_s.vmy;
		n_vmz_s += nvm_s.vmz;
		ElemNodeVM& nvm_f = elem_node_vm_f[ne_id];
		n_vm_f += nvm_f.vm;
		n_vmx_f += nvm_f.vmx;
		n_vmy_f += nvm_f.vmy;
		n_vmz_f += nvm_f.vmz;

		if (n_id != node_has_elem0[ve_id + 1])
		{
			// solid
			n_am_s *= one_fourth;
			node_am_s[n_id] = n_am_s;
			Acceleration& n_a_s = node_a_s[n_id];
			n_a_s.ax = n_fx_s / n_am_s;
			n_a_s.ay = n_fy_s / n_am_s;
			n_a_s.az = n_fz_s / n_am_s;
			Velocity& n_v_s = node_v_s[n_id];
			n_v_s.vx = n_vmx_s / n_vm_s + n_a_s.ax * dt;
			n_v_s.vy = n_vmy_s / n_vm_s + n_a_s.ay * dt;
			n_v_s.vz = n_vmz_s / n_vm_s + n_a_s.az * dt;
			NodeHasVBC& n_has_vbc_s = node_has_vbc_s[n_id];
			bc_mask = SIZE_MAX + size_t(n_has_vbc_s.has_vx_bc);
			n_a_s.iax &= bc_mask;
			n_v_s.ivx &= bc_mask;
			bc_mask = SIZE_MAX + size_t(n_has_vbc_s.has_vy_bc);
			n_a_s.iay &= bc_mask;
			n_v_s.ivy &= bc_mask;
			bc_mask = SIZE_MAX + size_t(n_has_vbc_s.has_vz_bc);
			n_a_s.iaz &= bc_mask;
			n_v_s.ivz &= bc_mask;
			// fluid
			n_am_f *= one_fourth;
			node_am_f[n_id] = n_am_f;
			Acceleration& n_a_f = node_a_f[n_id];
			n_a_f.ax = n_fx_f / n_am_f;
			n_a_f.ay = n_fy_f / n_am_f;
			n_a_f.az = n_fz_f / n_am_f;
			Velocity& n_v_f = node_v_f[n_id];
			n_v_f.vx = n_vmx_f / n_vm_f + n_a_f.ax * dt;
			n_v_f.vy = n_vmy_f / n_vm_f + n_a_f.ay * dt;
			n_v_f.vz = n_vmz_f / n_vm_f + n_a_f.az * dt;
			NodeHasVBC& n_has_vbc_f = node_has_vbc_f[n_id];
			bc_mask = SIZE_MAX + size_t(n_has_vbc_f.has_vx_bc);
			n_a_f.iax &= bc_mask;
			n_v_f.ivx &= bc_mask;
			bc_mask = SIZE_MAX + size_t(n_has_vbc_f.has_vy_bc);
			n_a_f.iay &= bc_mask;
			n_v_f.ivy &= bc_mask;
			bc_mask = SIZE_MAX + size_t(n_has_vbc_f.has_vz_bc);
			n_a_f.iaz &= bc_mask;
			n_v_f.ivz &= bc_mask;

			n_id = node_has_elem0[ve_id + 1];
			assert(n_id < self.node_num || n_id == SIZE_MAX);

			n_am_s = 0.0;
			n_am_f = 0.0;
			n_fx_s = 0.0;
			n_fy_s = 0.0;
			n_fz_s = 0.0;
			n_fx_f = 0.0;
			n_fy_f = 0.0;
			n_fz_f = 0.0;
			n_vm_s = 0.0;
			n_vmx_s = 0.0;
			n_vmy_s = 0.0;
			n_vmz_s = 0.0;
			n_vm_f = 0.0;
			n_vmx_f = 0.0;
			n_vmy_f = 0.0;
			n_vmz_f = 0.0;
		}
	}

#pragma omp master
	{
		if (md.has_rigid_cylinder())
		{
			RigidCylinder& rcy = md.get_rigid_cylinder();
			rcy.set_cont_force(self.cf_tmp);
			rcy.update_motion(dt);
			self.cf_tmp.reset();
		}

		if (md.has_t3d_rigid_mesh())
		{
			RigidObjectByT3DMesh& rb = md.get_t3d_rigid_mesh();
			rb.set_cont_force(self.cf_tmp);
			rb.update_motion(dt);
			self.cf_tmp.reset();
		}

#ifdef _DEBUG
		self.prev_valid_pcl_num_tmp = self.prev_valid_pcl_num;
#endif
		self.prev_valid_pcl_num = self.valid_pcl_num;
		self.valid_pcl_num = 0;
	}

#pragma omp barrier

	// cal element strain and "enhancement"
	double e_de_vol_s, e_de_vol_f;
	for (ve_id = 0; ve_id < my_valid_elem_num; ++ve_id)
	{
		e_id = my_valid_elem_id[ve_id];
		assert(e_id < self.elem_num);

		const ElemNodeIndex& eni = elem_node_id[e_id];
		const Velocity& n1_v_s = node_v_s[eni.n1];
		const Velocity& n2_v_s = node_v_s[eni.n2];
		const Velocity& n3_v_s = node_v_s[eni.n3];
		const Velocity& n4_v_s = node_v_s[eni.n4];
		const DShapeFuncABC& e_dN = elem_N_abc[e_id];
		StrainInc& e_de = elem_de[e_id];
		e_de.de11 = (e_dN.dN1_dx * n1_v_s.vx + e_dN.dN2_dx * n2_v_s.vx + e_dN.dN3_dx * n3_v_s.vx + e_dN.dN4_dx * n4_v_s.vx) * dt;
		e_de.de22 = (e_dN.dN1_dy * n1_v_s.vy + e_dN.dN2_dy * n2_v_s.vy + e_dN.dN3_dy * n3_v_s.vy + e_dN.dN4_dy * n4_v_s.vy) * dt;
		e_de.de33 = (e_dN.dN1_dz * n1_v_s.vz + e_dN.dN2_dz * n2_v_s.vz + e_dN.dN3_dz * n3_v_s.vz + e_dN.dN4_dz * n4_v_s.vz) * dt;
		e_de.de12 = (e_dN.dN1_dx * n1_v_s.vy + e_dN.dN2_dx * n2_v_s.vy + e_dN.dN3_dx * n3_v_s.vy + e_dN.dN4_dx * n4_v_s.vy
			+ e_dN.dN1_dy * n1_v_s.vx + e_dN.dN2_dy * n2_v_s.vx + e_dN.dN3_dy * n3_v_s.vx + e_dN.dN4_dy * n4_v_s.vx) * dt * 0.5;
		e_de.de23 = (e_dN.dN1_dy * n1_v_s.vz + e_dN.dN2_dy * n2_v_s.vz + e_dN.dN3_dy * n3_v_s.vz + e_dN.dN4_dy * n4_v_s.vz
			+ e_dN.dN1_dz * n1_v_s.vy + e_dN.dN2_dz * n2_v_s.vy + e_dN.dN3_dz * n3_v_s.vy + e_dN.dN4_dz * n4_v_s.vy) * dt * 0.5;
		e_de.de31 = (e_dN.dN1_dz * n1_v_s.vx + e_dN.dN2_dz * n2_v_s.vx + e_dN.dN3_dz * n3_v_s.vx + e_dN.dN4_dz * n4_v_s.vx
			+ e_dN.dN1_dx * n1_v_s.vz + e_dN.dN2_dx * n2_v_s.vz + e_dN.dN3_dx * n3_v_s.vz + e_dN.dN4_dx * n4_v_s.vz) * dt * 0.5;
		e_de_vol_s = e_de.de11 + e_de.de22 + e_de.de33;
		elem_m_de_vol_s[e_id] = elem_pcl_m_s[e_id] * e_de_vol_s;
		const Velocity& n1_v_f = node_v_f[eni.n1];
		const Velocity& n2_v_f = node_v_f[eni.n2];
		const Velocity& n3_v_f = node_v_f[eni.n3];
		const Velocity& n4_v_f = node_v_f[eni.n4];
		e_de_vol_f = (1.0 - elem_pcl_n[e_id]) / elem_pcl_n[e_id] * -e_de_vol_s
			- (e_dN.dN1_dx * n1_v_f.vx + e_dN.dN2_dx * n2_v_f.vx + e_dN.dN3_dx * n3_v_f.vx + e_dN.dN4_dx * n4_v_f.vx
				+ e_dN.dN1_dy * n1_v_f.vy + e_dN.dN2_dy * n2_v_f.vy + e_dN.dN3_dy * n3_v_f.vy + e_dN.dN4_dy * n4_v_f.vy
				+ e_dN.dN1_dz * n1_v_f.vz + e_dN.dN2_dz * n2_v_f.vz + e_dN.dN3_dz * n3_v_f.vz + e_dN.dN4_dz * n4_v_f.vz) * dt;
		elem_m_de_vol_f[e_id] = elem_pcl_m_f[e_id] * e_de_vol_f;
		e_de_vol_s *= one_third;
		e_de.de11 -= e_de_vol_s;
		e_de.de22 -= e_de_vol_s;
		e_de.de33 -= e_de_vol_s;
	}

#pragma omp barrier

	double n_am_de_vol_s = 0.0;
	double n_am_de_vol_f = 0.0;
	n_id = node_has_elem0[ve_id0];
	assert(n_id < self.node_num);
	for (ve_id = ve_id0; ve_id < ve_id1; ++ve_id)
	{
		e_id = node_elem_pair0[ve_id] / 4;
		assert(e_id < self.elem_num);
		n_am_de_vol_s += elem_m_de_vol_s[e_id];
		n_am_de_vol_f += elem_m_de_vol_f[e_id];
		if (n_id != node_has_elem0[ve_id + 1])
		{
			node_de_vol_s[n_id] = n_am_de_vol_s * one_fourth / node_am_s[n_id];
			node_de_vol_f[n_id] = n_am_de_vol_f * one_fourth / node_am_f[n_id];
			n_id = node_has_elem0[ve_id + 1];
			assert(n_id < self.node_num || n_id == SIZE_MAX);
			n_am_de_vol_s = 0.0;
			n_am_de_vol_f = 0.0;
		}
	}

#pragma omp barrier

	const Acceleration* pn1_a_s, * pn2_a_s, * pn3_a_s, * pn4_a_s;
	const Acceleration* pn1_a_f, * pn2_a_f, * pn3_a_f, * pn4_a_f;
	const Velocity* pn1_v_s, * pn2_v_s, * pn3_v_s, * pn4_v_s;
	const Velocity* pn1_v_f, * pn2_v_f, * pn3_v_f, * pn4_v_f;
	StrainInc* pe_de;
	const double* estrain, * pstrain, * dstress;
	double p_x, p_y, p_z, e_density_f;
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
			pn1_a_s = node_a_s + eni.n1;
			pn2_a_s = node_a_s + eni.n2;
			pn3_a_s = node_a_s + eni.n3;
			pn4_a_s = node_a_s + eni.n4;
			pn1_a_f = node_a_f + eni.n1;
			pn2_a_f = node_a_f + eni.n2;
			pn3_a_f = node_a_f + eni.n3;
			pn4_a_f = node_a_f + eni.n4;
			pn1_v_s = node_v_s + eni.n1;
			pn2_v_s = node_v_s + eni.n2;
			pn3_v_s = node_v_s + eni.n3;
			pn4_v_s = node_v_s + eni.n4;
			pn1_v_f = node_v_f + eni.n1;
			pn2_v_f = node_v_f + eni.n2;
			pn3_v_f = node_v_f + eni.n3;
			pn4_v_f = node_v_f + eni.n4;

			e_de_vol_s = (node_de_vol_s[eni.n1]
				+ node_de_vol_s[eni.n2]
				+ node_de_vol_s[eni.n3]
				+ node_de_vol_s[eni.n4]) * one_fourth;
			e_n = (e_de_vol_s + elem_pcl_n[e_id]) / (1.0 + e_de_vol_s);

			e_de_vol_f = (node_de_vol_f[eni.n1]
				+ node_de_vol_f[eni.n2]
				+ node_de_vol_f[eni.n3]
				+ node_de_vol_f[eni.n4]) * one_fourth;
			e_density_f = elem_density_f[e_id] / (1.0 - e_de_vol_f);
			e_p = elem_p[e_id] + self.Kf * e_de_vol_f;

			pe_de = elem_de + e_id;
			e_de_vol_s *= one_third;
			pe_de->de11 += e_de_vol_s;
			pe_de->de22 += e_de_vol_s;
			pe_de->de33 += e_de_vol_s;
		}

		// update velocity
		ShapeFunc& p_N = pcl_N0[p_id];
		Velocity& p_v_s0 = pcl_v_s0[p_id];
		p_v_s0.vx += (p_N.N1 * pn1_a_s->ax + p_N.N2 * pn2_a_s->ax + p_N.N3 * pn3_a_s->ax + p_N.N4 * pn4_a_s->ax) * dt;
		p_v_s0.vy += (p_N.N1 * pn1_a_s->ay + p_N.N2 * pn2_a_s->ay + p_N.N3 * pn3_a_s->ay + p_N.N4 * pn4_a_s->ay) * dt;
		p_v_s0.vz += (p_N.N1 * pn1_a_s->az + p_N.N2 * pn2_a_s->az + p_N.N3 * pn3_a_s->az + p_N.N4 * pn4_a_s->az) * dt;
		Velocity& p_v_f0 = pcl_v_f0[p_id];
		p_v_f0.vx += (p_N.N1 * pn1_a_f->ax + p_N.N2 * pn2_a_f->ax + p_N.N3 * pn3_a_f->ax + p_N.N4 * pn4_a_f->ax) * dt;
		p_v_f0.vy += (p_N.N1 * pn1_a_f->ay + p_N.N2 * pn2_a_f->ay + p_N.N3 * pn3_a_f->ay + p_N.N4 * pn4_a_f->ay) * dt;
		p_v_f0.vz += (p_N.N1 * pn1_a_f->az + p_N.N2 * pn2_a_f->az + p_N.N3 * pn3_a_f->az + p_N.N4 * pn4_a_f->az) * dt;

		// update displacement
		Displacement& p_u_s0 = pcl_u_s0[p_id];
		p_u_s0.ux += (p_N.N1 * pn1_v_s->vx + p_N.N2 * pn2_v_s->vx + p_N.N3 * pn3_v_s->vx + p_N.N4 * pn4_v_s->vx) * dt;
		p_u_s0.uy += (p_N.N1 * pn1_v_s->vy + p_N.N2 * pn2_v_s->vy + p_N.N3 * pn3_v_s->vy + p_N.N4 * pn4_v_s->vy) * dt;
		p_u_s0.uz += (p_N.N1 * pn1_v_s->vz + p_N.N2 * pn2_v_s->vz + p_N.N3 * pn3_v_s->vz + p_N.N4 * pn4_v_s->vz) * dt;
		Displacement& p_u_f0 = pcl_u_f0[p_id];
		p_u_f0.ux += (p_N.N1 * pn1_v_f->vx + p_N.N2 * pn2_v_f->vx + p_N.N3 * pn3_v_f->vx + p_N.N4 * pn4_v_f->vx) * dt;
		p_u_f0.uy += (p_N.N1 * pn1_v_f->vy + p_N.N2 * pn2_v_f->vy + p_N.N3 * pn3_v_f->vy + p_N.N4 * pn4_v_f->vy) * dt;
		p_u_f0.uz += (p_N.N1 * pn1_v_f->vz + p_N.N2 * pn2_v_f->vz + p_N.N3 * pn3_v_f->vz + p_N.N4 * pn4_v_f->vz) * dt;

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

		// update stress
		MatModel::MaterialModel& pcl_mm = *pcl_mat_model[ori_p_id];
		pcl_mm.integrate(pe_de->de);
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
		Strain& p_e1 = pcl_strain1[prev_p_id];
		Strain& p_e0 = pcl_strain0[p_id];
		p_e0.e11 = p_e1.e11 + pe_de->de11;
		p_e0.e22 = p_e1.e22 + pe_de->de22;
		p_e0.e33 = p_e1.e33 + pe_de->de33;
		p_e0.e12 = p_e1.e12 + pe_de->de12;
		p_e0.e23 = p_e1.e23 + pe_de->de23;
		p_e0.e31 = p_e1.e31 + pe_de->de31;

		estrain = pcl_mm.get_dstrain_e();
		Strain& p_ee1 = pcl_estrain1[prev_p_id];
		Strain& p_ee0 = pcl_estrain0[p_id];
		p_ee0.e11 = p_ee1.e11 + estrain[0];
		p_ee0.e22 = p_ee1.e22 + estrain[1];
		p_ee0.e33 = p_ee1.e33 + estrain[2];
		p_ee0.e12 = p_ee1.e12 + estrain[3];
		p_ee0.e23 = p_ee1.e23 + estrain[4];
		p_ee0.e31 = p_ee1.e31 + estrain[5];

		pstrain = pcl_mm.get_dstrain_p();
		Strain& p_pe1 = pcl_pstrain1[prev_p_id];
		Strain& p_pe0 = pcl_pstrain0[p_id];
		p_pe0.e11 = p_pe1.e11 + pstrain[0];
		p_pe0.e22 = p_pe1.e22 + pstrain[1];
		p_pe0.e33 = p_pe1.e33 + pstrain[2];
		p_pe0.e12 = p_pe1.e12 + pstrain[3];
		p_pe0.e23 = p_pe1.e23 + pstrain[4];
		p_pe0.e31 = p_pe1.e31 + pstrain[5];
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
