#include "SimulationsOMP_pcp.h"

#include <fstream>
#include <omp.h>

#include "Step_T2D_ME_mt_Geo.h"

#define one_third (1.0/3.0)
#define N_min (1.0e-10)
#define Block_Low(th_id, th_num, data_num) ((th_id)*(data_num)/(th_num))

#ifdef _DEBUG
static std::fstream res_file_t2d_me_mt;
#endif

Step_T2D_ME_mt_Geo::Step_T2D_ME_mt_Geo(const char* _name) : 
	Step_OMP(_name, "Step_T2D_ME_mt_Geo", &substep_func_omp_T2D_ME_mt_Geo) {}

Step_T2D_ME_mt_Geo::~Step_T2D_ME_mt_Geo() {}

int Step_T2D_ME_mt_Geo::init_calculation()
{
#ifdef _DEBUG
	res_file_t2d_me_mt.open("t3d_chm_mt_res.txt", std::ios::out | std::ios::binary);
#endif

	Model_T2D_ME_mt &md = *(Model_T2D_ME_mt *)model;

	omp_set_num_threads(thread_num);

	pcl_m = md.pcl_m;
	pcl_bf = md.pcl_bf;
	pcl_t = md.pcl_t;
	pcl_pos = md.pcl_pos;
	pcl_vol = md.pcl_vol;
	pcl_mat_model = md.pcl_mat_model;

	Model_T2D_ME_mt::SortedPclVarArrays& md_spva0
		= md.sorted_pcl_var_arrays[0];
	SortedPclVarArrays& spva0 = sorted_pcl_var_arrays[0];
	spva0.pcl_index = md_spva0.pcl_index;
	spva0.pcl_v = md_spva0.pcl_v;
	spva0.pcl_disp = md_spva0.pcl_disp;
	spva0.pcl_density = md_spva0.pcl_density;
	spva0.pcl_stress = md_spva0.pcl_stress;
	spva0.pcl_N = md_spva0.pcl_N;

	elem_num = md.elem_num;
	node_num = md.node_num;

	elem_node_id = md.elem_node_id;
	elem_N_ab = md.elem_dN_ab;
	elem_N_c = md.elem_dN_c;
	elem_area = md.elem_area;

	elem_pcl_m = md.elem_pcl_m;
	elem_de = md.elem_de;

	elem_m_de_vol = md.elem_m_de_vol;

	// element-node data
	elem_node_vm = md.elem_node_vm;
	elem_node_force = md.elem_node_force;

	node_a = md.node_a;
	node_v = md.node_v;
	node_has_vbc = md.node_has_vbc;
	node_am = md.node_am;
	node_de_vol = md.node_de_vol;

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
			const double p_vol = pcl_m[ori_p_id] / spva0.pcl_density[ori_p_id];
			Position& p_p = pcl_pos[ori_p_id];
			Displacement& p_disp = spva0.pcl_disp[p_id];
			p_p.x += p_disp.ux;
			p_p.y += p_disp.uy;
			p_disp.ux = 0.0;
			p_disp.uy = 0.0;
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
		size_t* const my_node_has_elem = node_has_elem1 + e_id * 3;
		size_t* const my_node_elem_pair = node_elem_pair1 + e_id * 3;
		size_t my_valid_elem_num = 0;
		for (p_id = p_id0; p_id < p_id1; ++p_id)
		{
			if (e_id != pcl_in_elem0[p_id + 1])
			{
				ne_id = my_valid_elem_num * 3;
				my_valid_elem_id[my_valid_elem_num++] = e_id;

				const ElemNodeIndex& eni = elem_node_id[e_id];
				my_node_has_elem[ne_id] = eni.n1;
				my_node_elem_pair[ne_id] = e_id * 3;
				my_node_has_elem[++ne_id] = eni.n2;
				my_node_elem_pair[ne_id] = e_id * 3 + 1;
				my_node_has_elem[++ne_id] = eni.n3;
				my_node_elem_pair[ne_id] = e_id * 3 + 2;

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
			++my_cbin[data_digit(my_node_has_elem[ve_id * 3], 0)];
			++my_cbin[data_digit(my_node_has_elem[ve_id * 3 + 1], 0)];
			++my_cbin[data_digit(my_node_has_elem[ve_id * 3 + 2], 0)];
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
			pos_id = --my_sbin[data_digit(my_node_has_elem[ve_id * 3], 0)];
			node_has_elem0[pos_id] = my_node_has_elem[ve_id * 3];
			node_elem_pair0[pos_id] = my_node_elem_pair[ve_id * 3];
			pos_id = --my_sbin[data_digit(my_node_has_elem[ve_id * 3 + 1], 0)];
			node_has_elem0[pos_id] = my_node_has_elem[ve_id * 3 + 1];
			node_elem_pair0[pos_id] = my_node_elem_pair[ve_id * 3 + 1];
			pos_id = --my_sbin[data_digit(my_node_has_elem[ve_id * 3 + 2], 0)];
			node_has_elem0[pos_id] = my_node_has_elem[ve_id * 3 + 2];
			node_elem_pair0[pos_id] = my_node_elem_pair[ve_id * 3 + 2];
		}

#pragma omp barrier
		size_t ve_id0 = Block_Low(my_th_id, thread_num, valid_elem_num * 3);
		size_t ve_id1 = Block_Low(my_th_id + 1, thread_num, valid_elem_num * 3);
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
				assert(node_elem_pair0[ve_id] < elem_num * 3);
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
			node_has_elems[0][valid_elem_num * 3] = SIZE_MAX;
			node_has_elems[1][valid_elem_num * 3] = SIZE_MAX;
		}

		// modify ne_id0, ne_id1
		size_t n_id;
		n_id = node_has_elem0[ve_id0];
		while (ve_id0 != SIZE_MAX && n_id == node_has_elem0[--ve_id0]);
		++ve_id0;
		assert(ve_id0 <= valid_elem_num * 3);
		n_id = node_has_elem0[ve_id1];
		while (ve_id1 != SIZE_MAX && n_id == node_has_elem0[--ve_id1]);
		++ve_id1;
		assert(ve_id1 <= valid_elem_num * 3);
		thd.ve_id0 = ve_id0;
		thd.ve_id1 = ve_id1;
	}

	reorder(md, valid_pcl_num, prev_pcl_ids[0]);

	init_f_ub_is_init = false;
	e_kin_max_is_init = false;
	e_kin_prev = 0.0;
	return 0;
}

int Step_T2D_ME_mt_Geo::finalize_calculation()
{
	for (size_t t_id = 0; t_id < thread_num; ++t_id)
		thread_datas[t_id].~ThreadData();
	return 0;
}

int substep_func_omp_T2D_ME_mt_Geo(
	void* _self,
	size_t my_th_id,
	double dt,
	double cur_time,
	size_t substp_id)
{
	typedef Model_T2D_ME_mt::ShapeFunc ShapeFunc;
	typedef Model_T2D_ME_mt::ShapeFuncAB DShapeFuncAB;
	typedef Model_T2D_ME_mt::ShapeFuncC DShapeFuncC;
	typedef Model_T2D_ME_mt::Force Force;
	typedef Model_T2D_ME_mt::Position Position;
	typedef Model_T2D_ME_mt::Displacement Displacement;
	typedef Model_T2D_ME_mt::Velocity Velocity;
	typedef Model_T2D_ME_mt::Acceleration Acceleration;
	typedef Model_T2D_ME_mt::Stress Stress;
	typedef Model_T2D_ME_mt::Strain Strain;
	typedef Model_T2D_ME_mt::StrainInc StrainInc;
	typedef Model_T2D_ME_mt::SortedPclVarArrays SortedPclVarArrays;
	typedef Model_T2D_ME_mt::ElemNodeIndex ElemNodeIndex;
	typedef Model_T2D_ME_mt::ElemNodeVM ElemNodeVM;
	typedef Model_T2D_ME_mt::NodeHasVBC NodeHasVBC;
	typedef Step_T2D_ME_mt_Geo::ThreadData ThreadData;

	Step_T2D_ME_mt_Geo& self = *(Step_T2D_ME_mt_Geo*)(_self);
	
	if (self.valid_pcl_num == 0)
	{
#pragma omp master
		self.abort_calculation();

#pragma omp barrier
		return 0;
	}
	
	Model_T2D_ME_mt& md = *(Model_T2D_ME_mt *)(self.model);

	const double* const pcl_m = self.pcl_m;
	const Force* const pcl_bf = self.pcl_bf;
	const Force* const pcl_t = self.pcl_t;
	const Position* const pcl_pos = self.pcl_pos;
	double* const pcl_vol = self.pcl_vol;
	MatModel::MaterialModel** const pcl_mat_model = self.pcl_mat_model;
	
	const size_t thread_num = self.thread_num;
	ThreadData &thd = self.thread_datas[my_th_id];

	SortedPclVarArrays &spva0 = self.sorted_pcl_var_arrays[0];
	size_t* const pcl_index0 = spva0.pcl_index;
	const double* const pcl_density0 = spva0.pcl_density;
	Velocity* const pcl_v0 = spva0.pcl_v;
	Displacement* const pcl_disp = spva0.pcl_disp;
	ShapeFunc* const pcl_N0 = spva0.pcl_N;
	Stress* const pcl_stress0 = spva0.pcl_stress;

	const ElemNodeIndex* const elem_node_id = self.elem_node_id;
	const DShapeFuncAB *const elem_N_ab = self.elem_N_ab;
	const DShapeFuncC *const elem_N_c = self.elem_N_c;
	const double* const elem_area = self.elem_area;

	double* const elem_pcl_m = self.elem_pcl_m;
	StrainInc* const elem_de = self.elem_de;
	double* const elem_m_de_vol = self.elem_m_de_vol;

	ElemNodeVM* const elem_node_vm = self.elem_node_vm;
	Force* const elem_node_force = self.elem_node_force;
	
	Acceleration* const node_a = self.node_a;
	Velocity* const node_v = self.node_v;
	NodeHasVBC* const node_has_vbc = self.node_has_vbc;
	double* const node_am = self.node_am;
	double* const node_de_vol = self.node_de_vol;
	
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
	size_t p_id, ori_p_id, e_id, ne_id;
	double p_m, p_vol, p_N_m, one_third_bfx, one_third_bfy;
	double en1_vm = 0.0;
	double en1_vmx = 0.0;
	double en1_vmy = 0.0;
	double en2_vm = 0.0;
	double en2_vmx = 0.0;
	double en2_vmy = 0.0;
	double en3_vm = 0.0;
	double en3_vmx = 0.0;
	double en3_vmy = 0.0;
	double e_p_m = 0.0;
	double e_p_vol = 0.0;
	double e_s11 = 0.0;
	double e_s22 = 0.0;
	double e_s12 = 0.0;
	double en1_fx = 0.0;
	double en1_fy = 0.0;
	double en2_fx = 0.0;
	double en2_fy = 0.0;
	double en3_fx = 0.0;
	double en3_fy = 0.0;
	e_id = pcl_in_elem0[p_id0];
	for (p_id = p_id0; p_id < p_id1; ++p_id)
	{
		// ori_p_id
		ori_p_id = pcl_index0[p_id];
		assert(ori_p_id < md.ori_pcl_num);

		// map pcl mass and volume
		// m_s
		p_m = pcl_m[ori_p_id];
		e_p_m += p_m;

		// map pcl volume
		p_vol = p_m / pcl_density0[p_id];
		pcl_vol[p_id] = p_vol;
		e_p_vol += p_vol;
		
		// map stress
		Stress& p_s0 = pcl_stress0[p_id];
		e_s11 += p_s0.s11 * p_vol;
		e_s22 += p_s0.s22 * p_vol;
		e_s12 += p_s0.s12 * p_vol;

		// map velocity
		ShapeFunc& p_N0 = pcl_N0[p_id];
		// solid velocity
		Velocity &p_v0 = pcl_v0[p_id];
		p_N_m = (p_N0.N1 > N_tol ? p_N0.N1 : N_tol) * p_m;
		en1_vm += p_N_m;
		en1_vmx += p_N_m * p_v0.vx;
		en1_vmy += p_N_m * p_v0.vy;
		p_N_m = (p_N0.N2 > N_tol ? p_N0.N2 : N_tol) * p_m;
		en2_vm += p_N_m;
		en2_vmx += p_N_m * p_v0.vx;
		en2_vmy += p_N_m * p_v0.vy;
		p_N_m = (p_N0.N3 > N_tol ? p_N0.N3 : N_tol) * p_m;
		en3_vm += p_N_m;
		en3_vmx += p_N_m * p_v0.vx;
		en3_vmy += p_N_m * p_v0.vy;

		// solid external load
		const Force &p_bf_s = self.pcl_bf[ori_p_id];
		one_third_bfx = one_third * p_bf_s.fx;
		one_third_bfy = one_third * p_bf_s.fy;
		const Force &p_t = self.pcl_t[ori_p_id];
		en1_fx += one_third_bfx + p_N0.N1 * p_t.fx;
		en1_fy += one_third_bfy + p_N0.N1 * p_t.fy;
		en2_fx += one_third_bfx + p_N0.N2 * p_t.fx;
		en2_fy += one_third_bfy + p_N0.N2 * p_t.fy;
		en3_fx += one_third_bfx + p_N0.N3 * p_t.fx;
		en3_fy += one_third_bfy + p_N0.N3 * p_t.fy;

		if (e_id != pcl_in_elem0[p_id + 1])
		{
			elem_pcl_m[e_id] = e_p_m;

			// v_s
			ElemNodeVM& en1_v_s = elem_node_vm[e_id * 3];
			en1_v_s.vm = en1_vm;
			en1_v_s.vmx = en1_vmx;
			en1_v_s.vmy = en1_vmy;
			ElemNodeVM& en2_v_s = elem_node_vm[e_id * 3 + 1];
			en2_v_s.vm = en2_vm;
			en2_v_s.vmx = en2_vmx;
			en2_v_s.vmy = en2_vmy;
			ElemNodeVM& en3_v_s = elem_node_vm[e_id * 3 + 2];
			en3_v_s.vm = en3_vm;
			en3_v_s.vmx = en3_vmx;
			en3_v_s.vmy = en3_vmy;
			
			e_s11 /= e_p_vol;
			e_s22 /= e_p_vol;
			e_s12 /= e_p_vol;
			if (e_p_vol > elem_area[e_id])
				e_p_vol = elem_area[e_id];

			const DShapeFuncAB& e_dN = elem_N_ab[e_id];
			// node 1
			Force& en1_f = elem_node_force[e_id * 3];
			en1_fx -= (e_dN.dN1_dx * e_s11 + e_dN.dN1_dy * e_s12) * e_p_vol;
			en1_f.fx = en1_fx;
			en1_fy -= (e_dN.dN1_dx * e_s12 + e_dN.dN1_dy * e_s22) * e_p_vol;
			en1_f.fy = en1_fy;
			// node 2
			Force& en2_f = elem_node_force[e_id * 3 + 1];
			en2_fx -= (e_dN.dN2_dx * e_s11 + e_dN.dN2_dy * e_s12) * e_p_vol;
			en2_f.fx = en2_fx;
			en2_fy -= (e_dN.dN2_dx * e_s12 + e_dN.dN2_dy * e_s22) * e_p_vol;
			en2_f.fy = en2_fy;
			// node 3
			Force& en3_f = elem_node_force[e_id * 3 + 2];
			en3_fx -= (e_dN.dN3_dx * e_s11 + e_dN.dN3_dy * e_s12) * e_p_vol;
			en3_f.fx = en3_fx;
			en3_fy -= (e_dN.dN3_dx * e_s12 + e_dN.dN3_dy * e_s22) * e_p_vol;
			en3_f.fy = en3_fy;

			e_id = pcl_in_elem0[p_id + 1];
			assert(e_id < self.elem_num || e_id == SIZE_MAX);

			e_p_m = 0.0;
			e_p_vol = 0.0;
			e_s11 = 0.0;
			e_s22 = 0.0;
			e_s12 = 0.0;
			en1_vm = 0.0;
			en1_vmx = 0.0;
			en1_vmy = 0.0;
			en2_vm = 0.0;
			en2_vmx = 0.0;
			en2_vmy = 0.0;
			en3_vm = 0.0;
			en3_vmx = 0.0;
			en3_vmy = 0.0;
			en1_fx = 0.0;
			en1_fy = 0.0;
			en2_fx = 0.0;
			en2_fy = 0.0;
			en3_fx = 0.0;
			en3_fy = 0.0;
		}
	}
	
#pragma omp barrier
	// update node variables
	size_t ve_id, n_id, bc_mask;
	const size_t ve_id0 = thd.ve_id0;
	const size_t ve_id1 = thd.ve_id1;
	double n_vm = 0.0;
	double n_vmx = 0.0;
	double n_vmy = 0.0;
	double n_am = 0.0;
	union
	{
		struct { double n_fx, n_fy; };
		struct { size_t in_fx, in_fy; };
	};
	n_fx = 0.0;
	n_fy = 0.0;
	double f_ub = 0.0;
	double e_kin = 0.0;
	n_id = node_has_elem0[ve_id0];
	assert(n_id < self.node_num);
	for (ve_id = ve_id0; ve_id < ve_id1; ++ve_id)
	{
		ne_id = node_elem_pair0[ve_id];
		assert(ne_id < self.elem_num * 3);
		e_id = ne_id / 3;

		ElemNodeVM& nvm = elem_node_vm[ne_id];
		n_vm += nvm.vm;
		n_vmx += nvm.vmx;
		n_vmy += nvm.vmy;

		n_am += elem_pcl_m[e_id];
		const Force &nf_s = elem_node_force[ne_id];
		n_fx += nf_s.fx;
		n_fy += nf_s.fy;

		if (n_id != node_has_elem0[ve_id + 1])
		{
			// solid
			n_am *= one_third;
			node_am[n_id] = n_am;
			Acceleration& n_a = node_a[n_id];
			n_a.ax = n_fx / n_am;
			n_a.ay = n_fy / n_am;
			Velocity& n_v = node_v[n_id];
			n_v.vx = n_vmx / n_vm + n_a.ax * dt;
			n_v.vy = n_vmy / n_vm + n_a.ay * dt;
			NodeHasVBC& n_has_vbc = node_has_vbc[n_id];
			bc_mask = size_t(n_has_vbc.has_vx_bc) + SIZE_MAX;
			n_a.iax &= bc_mask;
			n_v.ivx &= bc_mask;
			in_fx &= bc_mask;
			bc_mask = size_t(n_has_vbc.has_vy_bc) + SIZE_MAX;
			n_a.iay &= bc_mask;
			n_v.ivy &= bc_mask;
			in_fy &= bc_mask;

			f_ub += n_fx * n_fx + n_fy * n_fy;
			e_kin += n_am * (n_v.vx * n_v.vx + n_v.vy * n_v.vy);
						
			n_id = node_has_elem0[ve_id + 1];
			assert(n_id < self.node_num || n_id == SIZE_MAX);

			n_am = 0.0;
			n_fx = 0.0;
			n_fy = 0.0;
			n_vm = 0.0;
			n_vmx = 0.0;
			n_vmy = 0.0;
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
			Velocity& p_v0 = pcl_v0[p_id];
			p_v0.vx = 0.0;
			p_v0.vy = 0.0;
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
				Velocity& n_v = node_v[n_id];
				n_v.vx = 0.0;
				n_v.vy = 0.0;
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
	double e_de_vol;
	for (ve_id = 0; ve_id < my_valid_elem_num; ++ve_id)
	{
		e_id = my_valid_elem_id[ve_id];
		assert(e_id < self.elem_num);

		const ElemNodeIndex& eni = elem_node_id[e_id];
		const Velocity& n1_v = node_v[eni.n1];
		const Velocity& n2_v = node_v[eni.n2];
		const Velocity& n3_v = node_v[eni.n3];
		const DShapeFuncAB& e_dN = elem_N_ab[e_id];
		StrainInc& e_de = elem_de[e_id];
		e_de.de11 = (e_dN.dN1_dx * n1_v.vx + e_dN.dN2_dx * n2_v.vx + e_dN.dN3_dx * n3_v.vx) * dt;
		e_de.de22 = (e_dN.dN1_dy * n1_v.vy + e_dN.dN2_dy * n2_v.vy + e_dN.dN3_dy * n3_v.vy) * dt;
		e_de.de12 = (e_dN.dN1_dx * n1_v.vy + e_dN.dN2_dx * n2_v.vy + e_dN.dN3_dx * n3_v.vy
				   + e_dN.dN1_dy * n1_v.vx + e_dN.dN2_dy * n2_v.vx + e_dN.dN3_dy * n3_v.vx) * dt * 0.5;
		e_de_vol = e_de.de11 + e_de.de22;
		elem_m_de_vol[e_id] = elem_pcl_m[e_id] * e_de_vol;
		e_de_vol *= one_third;
		e_de.de11 -= e_de_vol;
		e_de.de22 -= e_de_vol;
	}

#pragma omp barrier
	double n_am_de_vol = 0.0;
	n_id = node_has_elem0[ve_id0];
	assert(n_id < self.node_num);
	for (ve_id = ve_id0; ve_id < ve_id1; ++ve_id)
	{
		e_id = node_elem_pair0[ve_id] / 3;
		assert(e_id < self.elem_num);
		n_am_de_vol += elem_m_de_vol[e_id];
		if (n_id != node_has_elem0[ve_id + 1])
		{
			node_de_vol[n_id] = n_am_de_vol * one_third / node_am[n_id];
			n_id = node_has_elem0[ve_id + 1];
			assert(n_id < self.node_num || n_id == SIZE_MAX);
			n_am_de_vol = 0.0;
		}
	}

#pragma omp barrier

	const Acceleration* pn1_a_s, *pn2_a_s, *pn3_a_s;
	const Velocity* pn1_v_s, *pn2_v_s, *pn3_v_s;
	StrainInc* pe_de;
	const double *estrain, *pstrain, *dstress;
	double p_x, p_y;
	size_t p_e_id;
	e_id = SIZE_MAX;
	for (p_id = p_id0; p_id < p_id1; ++p_id)
	{
		if (e_id != pcl_in_elem0[p_id])
		{
			e_id = pcl_in_elem0[p_id];
			assert(e_id < self.elem_num);

			const ElemNodeIndex& eni = elem_node_id[e_id];
			pn1_a_s = node_a + eni.n1;
			pn2_a_s = node_a + eni.n2;
			pn3_a_s = node_a + eni.n3;
			pn1_v_s = node_v + eni.n1;
			pn2_v_s = node_v + eni.n2;
			pn3_v_s = node_v + eni.n3;

			e_de_vol = (node_de_vol[eni.n1]
					  + node_de_vol[eni.n2]
					  + node_de_vol[eni.n3]) * one_third;

			pe_de = elem_de + e_id;
			e_de_vol *= one_third;
			pe_de->de11 += e_de_vol;
			pe_de->de22 += e_de_vol;
		}

		// update velocity
		ShapeFunc& p_N = pcl_N0[p_id];
		Velocity& p_v0 = pcl_v0[p_id];
		p_v0.vx += (p_N.N1 * pn1_a_s->ax + p_N.N2 * pn2_a_s->ax + p_N.N3 * pn3_a_s->ax) * dt;
		p_v0.vy += (p_N.N1 * pn1_a_s->ay + p_N.N2 * pn2_a_s->ay + p_N.N3 * pn3_a_s->ay) * dt;
		
		// update stress
		ori_p_id = spva0.pcl_index[p_id];
		assert(ori_p_id < md.ori_pcl_num);
		MatModel::MaterialModel& pcl_mm = *pcl_mat_model[ori_p_id];
		double dstrain[6] = { pe_de->de11, pe_de->de22, 0.0, pe_de->de12, 0.0, 0.0 };
		pcl_mm.integrate(dstrain);
		dstress = pcl_mm.get_dstress();
		Stress& p_s = pcl_stress0[p_id];
		p_s.s11 += dstress[0];
		p_s.s22 += dstress[1];
		p_s.s12 += dstress[3];
	}

#pragma omp master
	{
		self.continue_calculation();
	}

#pragma omp barrier
	return 0;
}

void Step_T2D_ME_mt_Geo::reorder(
	Model_T2D_ME_mt& md,
	size_t valid_pcl_num,
	const size_t* prev_pcl_ids)
{
	MemoryUtils::ItemArrayFast<char> tmp_mem;
	size_t tmp_mem_size;

	SortedPclVarArrays& spva = md.sorted_pcl_var_arrays[0];

	size_t* pcl_index = spva.pcl_index;
	tmp_mem_size = sizeof(size_t) * valid_pcl_num;
	size_t* pcl_index_tmp = (size_t *)tmp_mem.resize(tmp_mem_size);
	memcpy(pcl_index_tmp, pcl_index, tmp_mem_size);
	for (size_t p_id = 0; p_id < valid_pcl_num; ++p_id)
		pcl_index[p_id] = pcl_index_tmp[prev_pcl_ids[p_id]];

	double* pcl_density = spva.pcl_density;
	tmp_mem_size = sizeof(double) * valid_pcl_num;
	double* pcl_density_tmp = (double *)tmp_mem.resize(tmp_mem_size);
	memcpy(pcl_density_tmp, pcl_density, tmp_mem_size);
	for (size_t p_id = 0; p_id < valid_pcl_num; ++p_id)
		pcl_density[p_id] = pcl_density_tmp[prev_pcl_ids[p_id]];

	Velocity* pcl_v = spva.pcl_v;
	tmp_mem_size = sizeof(Velocity) * valid_pcl_num;
	Velocity* pcl_v_tmp = (Velocity*)tmp_mem.resize(tmp_mem_size);
	memcpy(pcl_v_tmp, pcl_v, tmp_mem_size);
	for (size_t p_id = 0; p_id < valid_pcl_num; ++p_id)
	{
		Velocity& p_v = pcl_v[p_id];
		Velocity& p_v_tmp = pcl_v_tmp[prev_pcl_ids[p_id]];
		p_v.vx = p_v_tmp.vx;
		p_v.vy = p_v_tmp.vy;
	}

	Displacement* pcl_disp = spva.pcl_disp;
	tmp_mem_size = sizeof(Displacement) * valid_pcl_num;
	Displacement* pcl_disp_tmp = (Displacement*)tmp_mem.resize(tmp_mem_size);
	memcpy(pcl_disp_tmp, pcl_disp, tmp_mem_size);
	for (size_t p_id = 0; p_id < valid_pcl_num; ++p_id)
	{
		Displacement& p_u = pcl_disp[p_id];
		Displacement& p_u_tmp = pcl_disp_tmp[prev_pcl_ids[p_id]];
		p_u.ux = p_u_tmp.ux;
		p_u.uy = p_u_tmp.uy;
	}

	ShapeFunc *pcl_N = spva.pcl_N;
	tmp_mem_size = sizeof(ShapeFunc) * valid_pcl_num;
	ShapeFunc* pcl_N_tmp = (ShapeFunc *)tmp_mem.resize(tmp_mem_size);
	memcpy(pcl_N_tmp, pcl_N, tmp_mem_size);
	for (size_t p_id = 0; p_id < valid_pcl_num; ++p_id)
	{
		ShapeFunc& p_N = pcl_N[p_id];
		ShapeFunc& p_N_tmp = pcl_N_tmp[prev_pcl_ids[p_id]];
		p_N.N1 = p_N_tmp.N1;
		p_N.N2 = p_N_tmp.N2;
		p_N.N3 = p_N_tmp.N3;
	}

	Stress* pcl_stress = spva.pcl_stress;
	tmp_mem_size = sizeof(Stress) * valid_pcl_num;
	Stress* pcl_stress_tmp = (Stress*)tmp_mem.resize(tmp_mem_size);
	memcpy(pcl_stress_tmp, pcl_stress, tmp_mem_size);
	for (size_t p_id = 0; p_id < valid_pcl_num; ++p_id)
	{
		Stress& p_s = pcl_stress[p_id];
		Stress& p_s_tmp = pcl_stress_tmp[prev_pcl_ids[p_id]];
		p_s.s11 = p_s_tmp.s11;
		p_s.s22 = p_s_tmp.s22;
		p_s.s12 = p_s_tmp.s12;
	}

	Strain* pcl_strain = spva.pcl_strain;
	tmp_mem_size = sizeof(Strain) * valid_pcl_num;
	Strain* pcl_strain_tmp = (Strain*)tmp_mem.resize(tmp_mem_size);
	memcpy(pcl_strain_tmp, pcl_strain, tmp_mem_size);
	for (size_t p_id = 0; p_id < valid_pcl_num; ++p_id)
	{
		Strain& p_e = pcl_strain[p_id];
		Strain& p_e_tmp = pcl_strain_tmp[prev_pcl_ids[p_id]];
		p_e.e11 = p_e_tmp.e11;
		p_e.e22 = p_e_tmp.e22;
		p_e.e12 = p_e_tmp.e12;
	}

	Strain* pcl_estrain = spva.pcl_estrain;
	tmp_mem_size = sizeof(Strain) * valid_pcl_num;
	Strain* pcl_estrain_tmp = (Strain*)tmp_mem.resize(tmp_mem_size);
	memcpy(pcl_estrain_tmp, pcl_estrain, tmp_mem_size);
	for (size_t p_id = 0; p_id < valid_pcl_num; ++p_id)
	{
		Strain& p_e = pcl_estrain[p_id];
		Strain& p_e_tmp = pcl_estrain_tmp[prev_pcl_ids[p_id]];
		p_e.e11 = p_e_tmp.e11;
		p_e.e22 = p_e_tmp.e22;
		p_e.e12 = p_e_tmp.e12;
	}

	Strain* pcl_pstrain = spva.pcl_pstrain;
	tmp_mem_size = sizeof(Strain) * valid_pcl_num;
	Strain* pcl_pstrain_tmp = (Strain*)tmp_mem.resize(tmp_mem_size);
	memcpy(pcl_pstrain_tmp, pcl_pstrain, tmp_mem_size);
	for (size_t p_id = 0; p_id < valid_pcl_num; ++p_id)
	{
		Strain& p_e = pcl_pstrain[p_id];
		Strain& p_e_tmp = pcl_pstrain_tmp[prev_pcl_ids[p_id]];
		p_e.e11 = p_e_tmp.e11;
		p_e.e22 = p_e_tmp.e22;
		p_e.e12 = p_e_tmp.e12;
	}
}
