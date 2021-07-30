#include "SimulationsOMP_pcp.h"

#include <fstream>
#include <iostream>
#include <omp.h>

#include "Step_T3D_CHM_mt_Geo.h"

#define one_fourth (0.25)
#define one_third (1.0/3.0)
#define N_min (1.0e-10)
#define Block_Low(th_id, th_num, data_num) ((th_id)*(data_num)/(th_num))

#ifdef _DEBUG
static std::fstream res_file_t2d_me_mt;
#endif

Step_T3D_CHM_mt_Geo::Step_T3D_CHM_mt_Geo(const char* _name) : 
	Step_OMP(_name, "Step_T3D_CHM_mt_Geo", &substep_func_omp_T3D_CHM_mt_Geo) {}

Step_T3D_CHM_mt_Geo::~Step_T3D_CHM_mt_Geo() {}

int Step_T3D_CHM_mt_Geo::init_calculation()
{
#ifdef _DEBUG
	res_file_t2d_me_mt.open("t3d_chm_mt_res.txt", std::ios::out | std::ios::binary);
#endif

	Model_T3D_CHM_mt &md = *(Model_T3D_CHM_mt *)model;

	omp_set_num_threads(thread_num);

	pcl_m_s = md.pcl_m_s;
	pcl_density_s = md.pcl_density_s;
	pcl_vol_s = md.pcl_vol_s;
	pcl_bf_s = md.pcl_bf_s;
	pcl_bf_f = md.pcl_bf_f;
	pcl_t = md.pcl_t;
	pcl_pos = md.pcl_pos;
	pcl_vol = md.pcl_vol;
	pcl_mat_model = md.pcl_mat_model;

	Model_T3D_CHM_mt::SortedPclVarArrays& md_spva0
		= md.sorted_pcl_var_arrays[0];
	SortedPclVarArrays& spva0 = sorted_pcl_var_arrays[0];
	spva0.pcl_index = md_spva0.pcl_index;
	spva0.pcl_n = md_spva0.pcl_n;
	spva0.pcl_v_s = md_spva0.pcl_v_s;
	spva0.pcl_u_s = md_spva0.pcl_u_s;
	spva0.pcl_stress = md_spva0.pcl_stress;
	spva0.pcl_strain = md_spva0.pcl_strain;
	spva0.pcl_estrain = md_spva0.pcl_estrain;
	spva0.pcl_pstrain = md_spva0.pcl_pstrain;
	spva0.pcl_N = md_spva0.pcl_N;

	elem_num = md.elem_num;
	node_num = md.node_num;

	elem_node_id = md.elem_node_id;
	elem_N_abc = md.elem_N_abc;
	elem_N_d = md.elem_N_d;
	elem_vol = md.elem_vol;

	elem_pcl_n = md.elem_pcl_n;
	elem_pcl_m_s = md.elem_pcl_m_s;
	elem_de = md.elem_de;

	elem_m_de_vol_s = md.elem_m_de_vol_s;

	// element-node data
	elem_node_vm_s = md.elem_node_vm_s;
	elem_node_force_s = md.elem_node_force_s;

	node_a_s = md.node_a_s;
	node_v_s = md.node_v_s;
	node_has_vbc_s = md.node_has_vbc_s;
	node_am_s = md.node_am_s;
	node_de_vol_s = md.node_de_vol_s;

	Kf = md.Kf; k = md.k; miu = md.miu;

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
	pcl_in_elems[0][md.pcl_num] = SIZE_MAX;
	pcl_in_elems[1][md.pcl_num] = SIZE_MAX;
	node_has_elems[0][-1] = SIZE_MAX;
	node_has_elems[1][-1] = SIZE_MAX;

	valid_pcl_num = 0;
	valid_elem_num = 0;
#pragma omp parallel
	{
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

		prev_pcl_id0 = prev_pcl_ids[0];
		prev_pcl_id1 = prev_pcl_ids[1];
		pcl_in_elem0 = pcl_in_elems[0];
		pcl_in_elem1 = pcl_in_elems[1];
		node_has_elem0 = node_has_elems[0];
		node_has_elem1 = node_has_elems[1];
		node_elem_pair0 = node_elem_pairs[0];
		node_elem_pair1 = node_elem_pairs[1];
		
		size_t my_th_id = size_t(omp_get_thread_num());

		ThreadData& thd = thread_datas[my_th_id];
		new (&thd) ThreadData;

		size_t p_id, ori_p_id, e_id;
		const size_t p_id0 = Block_Low(my_th_id, thread_num, md.pcl_num);
		const size_t p_id1 = Block_Low(my_th_id + 1, thread_num, md.pcl_num);
		size_t pcl_in_mesh_num = 0;
		for (p_id = p_id0; p_id < p_id1; ++p_id)
		{
			ori_p_id = spva0.pcl_index[p_id];
			const double p_vol = pcl_m_s[ori_p_id]
				/ (pcl_density_s[ori_p_id] * (1.0 - spva0.pcl_n[p_id]));
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
		
		size_t  bin_id, th_id, pos_id;
		size_t digit_disp, elem_num_tmp, * other_cbin;
		size_t* const my_cbin = elem_count_bin + my_th_id * 0x100;
		size_t* const my_sbin = elem_sum_bin + my_th_id * 0x100;
#define data_digit(num, disp) (((num) >> (disp)) & 0xFF)
#define swap(a, b) \
		(a) = (a) ^ (b); \
		(b) = (a) ^ (b); \
		(a) = (a) ^ (b)
		for (digit_disp = 0, elem_num_tmp = elem_num; 
			elem_num_tmp; digit_disp += 8, elem_num_tmp >>= 8)
		{
			memset(my_cbin, 0, 0x100 * sizeof(size_t));

			for (p_id = p_id0; p_id < p_id1; ++p_id)
			{
				++my_cbin[data_digit(pcl_in_elem0[p_id], digit_disp)];
				assert(pcl_in_elem0[p_id] < elem_num || pcl_in_elem0[p_id] == SIZE_MAX);
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
				assert((pcl_in_elem0[p_id] < elem_num ||
						pcl_in_elem0[p_id] == SIZE_MAX) &&
					   (prev_pcl_id0[p_id] < md.pcl_num));
			}

			swap(pcl_in_elem_ui0, pcl_in_elem_ui1);
			swap(prev_pcl_id_ui0, prev_pcl_id_ui1);
#pragma omp barrier
		}
		
#pragma omp single
		{
			pcl_in_elems[0] = pcl_in_elem0;
			pcl_in_elems[1] = pcl_in_elem1;
			prev_pcl_ids[0] = prev_pcl_id0;
			prev_pcl_ids[1] = prev_pcl_id1;
			md.pcl_num = valid_pcl_num;
		}

		thd.p_id0 = Block_Low(my_th_id, thread_num, valid_pcl_num);
		e_id = pcl_in_elem0[thd.p_id0];
		while (thd.p_id0 != SIZE_MAX && e_id == pcl_in_elem0[--thd.p_id0]);
		++thd.p_id0;
		assert(thd.p_id0 <= valid_pcl_num);
		thd.p_id1 = Block_Low(my_th_id + 1, thread_num, valid_pcl_num);
		e_id = pcl_in_elem0[thd.p_id1];
		while (thd.p_id1 != SIZE_MAX && e_id == pcl_in_elem0[--thd.p_id1]);
		++thd.p_id1;
		assert(thd.p_id1 <= valid_pcl_num);

		size_t ne_id;
		e_id = pcl_in_elem0[p_id0];
		size_t* const my_valid_elem_id = valid_elem_id + e_id;
		size_t* const my_node_has_elem = node_has_elem1 + e_id * 4;
		size_t* const my_node_elem_pair = node_elem_pair1 + e_id * 4;
		size_t my_valid_elem_num = 0;
		for (p_id = p_id0; p_id < p_id1; ++p_id)
		{
			if (e_id != pcl_in_elem0[p_id + 1])
			{
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
				assert(e_id < elem_num || e_id == SIZE_MAX);
			}
		}

		thd.valid_elem_num = my_valid_elem_num;
		thd.valid_elem_id = my_valid_elem_id;

#pragma omp critical
		valid_elem_num += my_valid_elem_num;

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
		size_t ve_id0 = Block_Low(my_th_id, thread_num, valid_elem_num * 4);
		size_t ve_id1 = Block_Low(my_th_id + 1, thread_num, valid_elem_num * 4);
		size_t node_num_tmp = node_num >> 8;
		for (digit_disp = 8; node_num_tmp; digit_disp += 8, node_num_tmp >>= 8)
		{
			memset(my_cbin, 0, sizeof(size_t) * 0x100);

			for (ve_id = ve_id0; ve_id < ve_id1; ++ve_id)
			{
				++my_cbin[data_digit(node_has_elem0[ve_id], digit_disp)];
				assert(node_has_elem0[ve_id] < node_num);
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
				assert(node_has_elem0[ve_id] < node_num);
				assert(node_elem_pair0[ve_id] < elem_num * 4);
			}

			swap(node_has_elem_ui0, node_has_elem_ui1);
			swap(node_elem_pair_ui0, node_elem_pair_ui1);
#pragma omp barrier
		}

#pragma omp single
		{
			node_has_elems[0] = node_has_elem0;
			node_has_elems[1] = node_has_elem1;
			node_elem_pairs[0] = node_elem_pair0;
			node_elem_pairs[1] = node_elem_pair1;
			node_has_elems[0][valid_elem_num * 4] = SIZE_MAX;
			node_has_elems[1][valid_elem_num * 4] = SIZE_MAX;
		}

		// modify ne_id0, ne_id1
		size_t n_id;
		n_id = node_has_elem0[ve_id0];
		while (ve_id0 != SIZE_MAX && n_id == node_has_elem0[--ve_id0]);
		++ve_id0;
		assert(ve_id0 <= valid_elem_num * 4);
		n_id = node_has_elem0[ve_id1];
		while (ve_id1 != SIZE_MAX && n_id == node_has_elem0[--ve_id1]);
		++ve_id1;
		assert(ve_id1 <= valid_elem_num * 4);
		thd.ve_id0 = ve_id0;
		thd.ve_id1 = ve_id1;
	}

	return 0;
}

int Step_T3D_CHM_mt_Geo::finalize_calculation()
{
	//Model_T3D_CHM_mt &md = *(Model_T3D_CHM_mt *)model;
	for (size_t t_id = 0; t_id < thread_num; ++t_id)
		thread_datas[t_id].~ThreadData();
	return 0;
}

int substep_func_omp_T3D_CHM_mt_Geo(
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
	typedef Step_T3D_CHM_mt_Geo::ThreadData ThreadData;

	Step_T3D_CHM_mt_Geo& self = *(Step_T3D_CHM_mt_Geo*)(_self);
	
	if (self.valid_pcl_num == 0)
	{
#pragma omp master
		self.abort_calculation();

#pragma omp barrier
		return 0;
	}
	
	Model_T3D_CHM_mt& md = *(Model_T3D_CHM_mt *)(self.model);

	const double* const pcl_m_s = self.pcl_m_s;
	const double* const pcl_density_s = self.pcl_density_s;
	const double* const pcl_vol_s = self.pcl_vol_s;
	const Force* const pcl_bf_s = self.pcl_bf_s;
	const Force* const pcl_t = self.pcl_t;
	const Position* const pcl_pos = self.pcl_pos;
	double* const pcl_vol = self.pcl_vol;
	MatModel::MaterialModel** const pcl_mat_model = self.pcl_mat_model;
	
	const size_t thread_num = self.thread_num;
	ThreadData &thd = self.thread_datas[my_th_id];
	SortedPclVarArrays &spva0 = self.sorted_pcl_var_arrays[0];
	size_t* const pcl_index0 = spva0.pcl_index;
	double* const pcl_n0 = spva0.pcl_n;
	Displacement* const pcl_u_s0 = spva0.pcl_u_s;
	Velocity* const pcl_v_s0 = spva0.pcl_v_s;
	Stress* const pcl_stress0 = spva0.pcl_stress;
	Strain* const pcl_strain0 = spva0.pcl_strain;
	Strain* const pcl_estrain0 = spva0.pcl_estrain;
	Strain* const pcl_pstrain0 = spva0.pcl_pstrain;
	ShapeFunc* const pcl_N0 = spva0.pcl_N;

	const ElemNodeIndex* const elem_node_id = self.elem_node_id;
	const DShapeFuncABC* const elem_N_abc = self.elem_N_abc;
	const DShapeFuncD* const elem_N_d = self.elem_N_d;
	const double* const elem_vol = self.elem_vol;

	double* const elem_pcl_n = self.elem_pcl_n;
	double* const elem_pcl_m_s = self.elem_pcl_m_s;
	StrainInc* const elem_de = self.elem_de;
	double* const elem_m_de_vol_s = self.elem_m_de_vol_s;

	ElemNodeVM* const elem_node_vm_s = self.elem_node_vm_s;
	Force* const elem_node_force_s = self.elem_node_force_s;
	
	Acceleration* const node_a_s = self.node_a_s;
	Velocity* const node_v_s = self.node_v_s;
	NodeHasVBC* const node_has_vbc_s = self.node_has_vbc_s;
	double* const node_am_s = self.node_am_s;
	double* const node_de_vol_s = self.node_de_vol_s;
	
	const size_t* prev_pcl_id0 = self.prev_pcl_ids[0];
	const size_t* pcl_in_elem0 = self.pcl_in_elems[0];
	const size_t* node_has_elem0 = self.node_has_elems[0];
	const size_t* node_elem_pair0 = self.node_elem_pairs[0];
	
#pragma omp master
	{
		self.f_ub = 0.0;
		self.e_kin = 0.0;
	}
	
	// update p_id0, p_id1
	const size_t p_id0 = thd.p_id0;
	const size_t p_id1 = thd.p_id1;
	size_t p_id, ori_p_id, prev_p_id, e_id, ne_id;
	double p_n, p_m_s, p_vol, p_N_m;
	double one_fourth_bfx_s, one_fourth_bfy_s, one_fourth_bfz_s;
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
	double e_n = 0.0;
	double e_p_vol = 0.0;
	double e_s11 = 0.0;
	double e_s22 = 0.0;
	double e_s33 = 0.0;
	double e_s12 = 0.0;
	double e_s23 = 0.0;
	double e_s31 = 0.0;
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
	e_id = pcl_in_elem0[p_id0];
	for (p_id = p_id0; p_id < p_id1; ++p_id)
	{
		prev_p_id = prev_pcl_id0[p_id];
		assert(prev_p_id < md.pcl_num);
		
		// ori_p_id
		ori_p_id = pcl_index0[prev_p_id];
		assert(ori_p_id < md.ori_pcl_num);

		// map pcl mass and volume
		// m_s
		p_m_s = pcl_m_s[ori_p_id];
		e_p_m_s += p_m_s;
		// e_n
		e_n += pcl_vol_s[ori_p_id];
		// vol
		p_n = pcl_n0[prev_p_id];
		p_vol = pcl_vol_s[ori_p_id] / (1.0 - p_n);
		pcl_vol[p_id] = p_vol;
		e_p_vol += p_vol;

		// map stress
		Stress& p_s0 = pcl_stress0[p_id];
		e_s11 += p_s0.s11 * p_vol;
		e_s22 += p_s0.s22 * p_vol;
		e_s33 += p_s0.s33 * p_vol;
		e_s12 += p_s0.s12 * p_vol;
		e_s23 += p_s0.s23 * p_vol;
		e_s31 += p_s0.s31 * p_vol;

		// map velocity
		ShapeFunc& p_N0 = pcl_N0[p_id];
		// solid velocity
		Velocity &p_v_s0 = pcl_v_s0[p_id];
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

		// solid external load
		const Force &p_bf_s = self.pcl_bf_s[ori_p_id];
		one_fourth_bfx_s = one_fourth * p_bf_s.fx;
		one_fourth_bfy_s = one_fourth * p_bf_s.fy;
		one_fourth_bfz_s = one_fourth * p_bf_s.fz;
		const Force &p_t = self.pcl_t[ori_p_id];
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
			
			elem_pcl_m_s[e_id] = e_p_m_s;
			e_n = 1.0 - e_n / e_p_vol;
			elem_pcl_n[e_id] = e_n;
			
			e_s11 /= e_p_vol;
			e_s22 /= e_p_vol;
			e_s33 /= e_p_vol;
			e_s12 /= e_p_vol;
			e_s23 /= e_p_vol;
			e_s31 /= e_p_vol;
			if (e_p_vol > elem_vol[e_id])
				e_p_vol = elem_vol[e_id];

			const DShapeFuncABC& e_dN = elem_N_abc[e_id];
			// node 1
			Force& en1_f_s = elem_node_force_s[e_id * 4];
			en1_fx_s -= (e_dN.dN1_dx * e_s11 + e_dN.dN1_dy * e_s12 + e_dN.dN1_dz * e_s31) * e_p_vol;
			en1_f_s.fx = en1_fx_s;
			en1_fy_s -= (e_dN.dN1_dx * e_s12 + e_dN.dN1_dy * e_s22 + e_dN.dN1_dz * e_s23) * e_p_vol;
			en1_f_s.fy = en1_fy_s;
			en1_fz_s -= (e_dN.dN1_dx * e_s31 + e_dN.dN1_dy * e_s23 + e_dN.dN1_dz * e_s33) * e_p_vol;
			en1_f_s.fz = en1_fz_s;
			// node 2
			Force& en2_f_s = elem_node_force_s[e_id * 4 + 1];
			en2_fx_s -= (e_dN.dN2_dx * e_s11 + e_dN.dN2_dy * e_s12 + e_dN.dN2_dz * e_s31) * e_p_vol;
			en2_f_s.fx = en2_fx_s;
			en2_fy_s -= (e_dN.dN2_dx * e_s12 + e_dN.dN2_dy * e_s22 + e_dN.dN2_dz * e_s23) * e_p_vol;
			en2_f_s.fy = en2_fy_s;
			en2_fz_s -= (e_dN.dN2_dx * e_s31 + e_dN.dN2_dy * e_s23 + e_dN.dN2_dz * e_s33) * e_p_vol;
			en2_f_s.fz = en2_fz_s;
			// node 3
			Force& en3_f_s = elem_node_force_s[e_id * 4 + 2];
			en3_fx_s -= (e_dN.dN3_dx * e_s11 + e_dN.dN3_dy * e_s12 + e_dN.dN3_dz * e_s31) * e_p_vol;
			en3_f_s.fx = en3_fx_s;
			en3_fy_s -= (e_dN.dN3_dx * e_s12 + e_dN.dN3_dy * e_s22 + e_dN.dN3_dz * e_s23) * e_p_vol;
			en3_f_s.fy = en3_fy_s;
			en3_fz_s -= (e_dN.dN3_dx * e_s31 + e_dN.dN3_dy * e_s23 + e_dN.dN3_dz * e_s33) * e_p_vol;
			en3_f_s.fz = en3_fz_s;
			// node 4
			Force& en4_f_s = elem_node_force_s[e_id * 4 + 3];
			en4_fx_s -= (e_dN.dN4_dx * e_s11 + e_dN.dN4_dy * e_s12 + e_dN.dN4_dz * e_s31) * e_p_vol;
			en4_f_s.fx = en4_fx_s;
			en4_fy_s -= (e_dN.dN4_dx * e_s12 + e_dN.dN4_dy * e_s22 + e_dN.dN4_dz * e_s23) * e_p_vol;
			en4_f_s.fy = en4_fy_s;
			en4_fz_s -= (e_dN.dN4_dx * e_s31 + e_dN.dN4_dy * e_s23 + e_dN.dN4_dz * e_s33) * e_p_vol;
			en4_f_s.fz = en4_fz_s;

			e_id = pcl_in_elem0[p_id + 1];
			assert(e_id < self.elem_num || e_id == SIZE_MAX);

			e_p_m_s = 0.0;
			e_n = 0.0;
			e_p_vol = 0.0;
			e_s11 = 0.0;
			e_s22 = 0.0;
			e_s33 = 0.0;
			e_s12 = 0.0;
			e_s23 = 0.0;
			e_s31 = 0.0;
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
		}
	}
	
#pragma omp barrier
	// update node variables
	size_t ve_id, n_id, bc_mask;
	const size_t ve_id0 = thd.ve_id0;
	const size_t ve_id1 = thd.ve_id1;
	double n_vm_s = 0.0;
	double n_vmx_s = 0.0;
	double n_vmy_s = 0.0;
	double n_vmz_s = 0.0;
	double n_am_s = 0.0;
	union
	{
		struct { double n_fx_s, n_fy_s, n_fz_s; };
		struct { size_t in_fx_s, in_fy_s, in_fz_s; };
	};
	n_fx_s = 0.0;
	n_fy_s = 0.0;
	n_fz_s = 0.0;
	double f_ub = 0.0;
	double e_kin = 0.0;
	n_id = node_has_elem0[ve_id0];
	assert(n_id < self.node_num);
	for (ve_id = ve_id0; ve_id < ve_id1; ++ve_id)
	{
		ne_id = node_elem_pair0[ve_id];
		assert(ne_id < self.elem_num * 4);
		e_id = ne_id / 4;

		ElemNodeVM& nvm_s = elem_node_vm_s[ne_id];
		n_vm_s += nvm_s.vm;
		n_vmx_s += nvm_s.vmx;
		n_vmy_s += nvm_s.vmy;
		n_vmz_s += nvm_s.vmz;

		n_am_s += elem_pcl_m_s[e_id];
		const Force &nf_s = elem_node_force_s[ne_id];
		n_fx_s += nf_s.fx;
		n_fy_s += nf_s.fy;
		n_fz_s += nf_s.fz;

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
			bc_mask = size_t(n_has_vbc_s.has_vx_bc) + SIZE_MAX;
			n_a_s.iax &= bc_mask;
			n_v_s.ivx &= bc_mask;
			in_fx_s &= bc_mask;
			bc_mask = size_t(n_has_vbc_s.has_vy_bc) + SIZE_MAX;
			n_a_s.iay &= bc_mask;
			n_v_s.ivy &= bc_mask;
			in_fy_s &= bc_mask;
			bc_mask = size_t(n_has_vbc_s.has_vz_bc) + SIZE_MAX;
			n_a_s.iaz &= bc_mask;
			n_v_s.ivz &= bc_mask;
			in_fz_s &= bc_mask;

			f_ub += n_fx_s * n_fx_s + n_fy_s * n_fy_s + n_fz_s * n_fz_s;
			e_kin += n_am_s * (n_v_s.vx * n_v_s.vx + n_v_s.vy * n_v_s.vy + n_v_s.vz * n_v_s.vz);
						
			n_id = node_has_elem0[ve_id + 1];
			assert(n_id < self.node_num || n_id == SIZE_MAX);

			n_am_s = 0.0;
			n_fx_s = 0.0;
			n_fy_s = 0.0;
			n_fz_s = 0.0;
			n_vm_s = 0.0;
			n_vmx_s = 0.0;
			n_vmy_s = 0.0;
			n_vmz_s = 0.0;
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

		if (self.e_kin < self.e_kin_prev)
		{
			if (!self.e_kin_max_is_init)
			{
				self.e_kin_max_is_init = true;
				self.e_kin_max = self.e_kin_prev;
			}
			self.e_kin_prev = 0.0;
			self.e_kin_ratio = 0.0;

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
			Velocity& p_v0 = pcl_v_s0[p_id];
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
				Acceleration& n_a = node_a_s[n_id];
				n_a.ax = 0.0;
				n_a.ay = 0.0;
				n_a.az = 0.0;
				Velocity& n_v = node_v_s[n_id];
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
	const size_t* my_valid_elem_id = thd.valid_elem_id;
	const size_t my_valid_elem_num = thd.valid_elem_num;
	double e_de_vol_s;
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
		e_de_vol_s *= one_third;
		e_de.de11 -= e_de_vol_s;
		e_de.de22 -= e_de_vol_s;
		e_de.de33 -= e_de_vol_s;
	}

#pragma omp barrier
	double n_am_de_vol_s = 0.0;
	n_id = node_has_elem0[ve_id0];
	assert(n_id < self.node_num);
	for (ve_id = ve_id0; ve_id < ve_id1; ++ve_id)
	{
		e_id = node_elem_pair0[ve_id] / 4;
		assert(e_id < self.elem_num);
		n_am_de_vol_s += elem_m_de_vol_s[e_id];
		if (n_id != node_has_elem0[ve_id + 1])
		{
			node_de_vol_s[n_id] = n_am_de_vol_s * one_fourth / node_am_s[n_id];
			n_id = node_has_elem0[ve_id + 1];
			assert(n_id < self.node_num || n_id == SIZE_MAX);
			n_am_de_vol_s = 0.0;
		}
	}

#pragma omp barrier

	const Acceleration* pn1_a_s, *pn2_a_s, *pn3_a_s, *pn4_a_s;
	const Velocity* pn1_v_s, *pn2_v_s, *pn3_v_s, *pn4_v_s;
	StrainInc* pe_de;
	const double *estrain, *pstrain, *dstress;
	double p_x, p_y, p_z;
	size_t p_e_id, pcl_in_mesh_num = 0;
	e_id = SIZE_MAX;
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
			pn1_v_s = node_v_s + eni.n1;
			pn2_v_s = node_v_s + eni.n2;
			pn3_v_s = node_v_s + eni.n3;
			pn4_v_s = node_v_s + eni.n4;

			e_de_vol_s = (node_de_vol_s[eni.n1]
						+ node_de_vol_s[eni.n2]
						+ node_de_vol_s[eni.n3]
						+ node_de_vol_s[eni.n4]) * one_fourth;

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
		
		// update stress
		assert(ori_p_id < md.ori_pcl_num);
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
	}

#pragma omp critical
	self.valid_pcl_num += pcl_in_mesh_num;

#pragma omp master
	{
		self.continue_calculation();
	}

#pragma omp barrier
	return 0;
}
