#include "SimulationsOMP_pcp.h"

#include <iostream>
#include <omp.h>

#include "Step_T3D_ME_mt.h"

#define one_third (1.0/3.0)
#define one_fourth (0.25)
#define Block_Low(th_id, th_num, data_num) ((th_id)*(data_num)/(th_num))

static std::fstream res_file_t3d_me_mt;

namespace
{
	inline void divide_task_to_thread(
		size_t thread_num,
		size_t data_num,
		size_t data_range[] // len = thread_num+1
	)
	{
		data_range[0] = 0;
		data_range[thread_num] = data_num;
		for (register size_t i = 1; i < thread_num; ++i)
			data_range[i] = Block_Low(i, thread_num, data_num);
	}

	size_t find_cloest_num_in_sorted_list(
		size_t num,
		size_t* num_list,
		size_t list_len
		)
	{
		size_t lid = 0;
		size_t uid = list_len - 1;
		size_t mid;
		while ((uid - lid) > 1)
		{
			mid = (lid + uid) / 2;
			if (num > num_list[mid])
				lid = mid;
			else
				uid = mid;
		}
		return (num - num_list[lid]) > (num_list[uid] - num) ? uid : lid;
	}
}

Step_T3D_ME_mt::Step_T3D_ME_mt(const char* _name) : 
	Step_OMP(_name, "Step_T3D_ME_mt", &substep_func_omp_T3D_ME_mt) {}

Step_T3D_ME_mt::~Step_T3D_ME_mt() {}

int Step_T3D_ME_mt::init_calculation()
{
	res_file_t3d_me_mt.open("t3d_stp_mt.txt", std::ios::out | std::ios::binary);
	
	Model_T3D_ME_mt &md = *(Model_T3D_ME_mt *)model;

	omp_set_num_threads(thread_num);

	pcl_m = md.pcl_m;
	pcl_bf = md.pcl_bf;
	pcl_t = md.pcl_t;
	pcl_pos = md.pcl_pos;
	pcl_vol = md.pcl_vol;
	pcl_N = md.pcl_N;
	pcl_mat_model = md.pcl_mat_model;

	size_t* radix_mem = (size_t *)radix_sort_var_mem.alloc(sizeof(size_t) * ((md.pcl_num + 2) * 2 + md.pcl_num) + Cache_Alignment);

	Model_T3D_ME_mt::SortedPclVarArrays &md_spva0 = md.sorted_pcl_var_arrays[0];
	SortedPclVarArrays &spva0 = sorted_pcl_var_arrays[0];
	spva0.pcl_index = md_spva0.pcl_index;
	spva0.pcl_density = md_spva0.pcl_density;
	spva0.pcl_disp = md_spva0.pcl_disp;
	spva0.pcl_v = md_spva0.pcl_v;
	spva0.pcl_stress = md_spva0.pcl_stress;
	spva0.pcl_strain = md_spva0.pcl_strain;
	spva0.pcl_estrain = md_spva0.pcl_estrain;
	spva0.pcl_pstrain = md_spva0.pcl_pstrain;
	spva0.pcl_in_elem = radix_mem + 1;
	spva0.pcl_in_elem[-1] = md.elem_num;
	spva0.pcl_in_elem[md.pcl_num] = md.elem_num;

	Model_T3D_ME_mt::SortedPclVarArrays &md_spva1 = md.sorted_pcl_var_arrays[1];
	SortedPclVarArrays& spva1 = sorted_pcl_var_arrays[1];
	spva1.pcl_index = md_spva1.pcl_index;
	spva1.pcl_density = md_spva1.pcl_density;
	spva1.pcl_disp = md_spva1.pcl_disp;
	spva1.pcl_v = md_spva1.pcl_v;
	spva1.pcl_stress = md_spva1.pcl_stress;
	spva1.pcl_strain = md_spva1.pcl_strain;
	spva1.pcl_estrain = md_spva1.pcl_estrain;
	spva1.pcl_pstrain = md_spva1.pcl_pstrain;
	spva1.pcl_in_elem = cache_aligned(radix_mem + md.pcl_num + 2);
	spva1.pcl_in_elem[-1] = md.elem_num;
	spva1.pcl_in_elem[md.pcl_num] = md.elem_num;

	pcl_in_elem_tmp = spva1.pcl_in_elem + md.pcl_num + 1;
	
	elem_num = md.elem_num;
	
	elem_node_id = md.elem_node_id;
	elem_id_array = md.elem_id_array;
	node_elem_id_array = md.node_elem_id_array;
	node_elem_list = md.node_elem_list;
	elem_dN_abc = md.elem_dN_abc;
	elem_dN_d = md.elem_dN_d;
	elem_vol = md.elem_vol;

	elem_substep_id = md.elem_substep_id;
	elem_density = md.elem_density;
	elem_pcl_m = md.elem_pcl_m;
	elem_pcl_vol = md.elem_pcl_vol;
	elem_de = md.elem_de;
	elem_m_de_vol = md.elem_m_de_vol;

	elem_node_vm = md.elem_node_vm;
	elem_node_force = md.elem_node_force;

	node_substep_id = md.node_substep_id;
	node_a = md.node_a;
	node_v = md.node_v;
	node_has_vbc = md.node_has_vbc;
	node_am = md.node_am;
	node_de_vol = md.node_de_vol;
	
	memset(elem_substep_id, 0xFF, sizeof(size_t) * md.elem_num);
	memset(node_substep_id, 0xFF, sizeof(size_t) * md.node_num);
	//std::cout << (elem_substep_id[0] == SIZE_MAX) << ", "
	//	<< (elem_substep_id[1] == SIZE_MAX) << ", "
	//	<< (elem_substep_id[md.elem_num - 1] == SIZE_MAX) << ",\n"
	//	<< (node_substep_id[0] == SIZE_MAX) << ", "
	//	<< (node_substep_id[1] == SIZE_MAX) << ", "
	//	<< (node_substep_id[md.node_num - 1] == SIZE_MAX) << ",\n";

	size_t th_num_1 = thread_num + 1;
	char* mem_range = (char*)task_range_mem.alloc((sizeof(ThreadData) + sizeof(size_t) * 2) * th_num_1);
	thread_datas = (ThreadData *)mem_range;
	mem_range += sizeof(ThreadData) * th_num_1;
	node_range = (size_t*)mem_range;
	mem_range += sizeof(size_t) * th_num_1;
	node_elem_range = (size_t*)mem_range;
	divide_task_to_thread(thread_num, md.elem_num * 4, node_elem_range);
	//for (size_t th_id = 0; th_id < thread_num; ++th_id)
	//	std::cout << th_id << ":\n"
	//		<< node_elem_range[th_id] << ", " << node_elem_range[th_id+1] << ",\n";

	node_range[0] = 0;
	node_range[thread_num] = md.node_num;
	size_t* node_elem_list_1 = md.node_elem_list - 1;
	for (size_t th_id = 1; th_id < thread_num; ++th_id)
	{
		node_range[th_id] = find_cloest_num_in_sorted_list(
			node_elem_range[th_id],
			node_elem_list_1,
			md.node_num + 1
			);
		node_elem_range[th_id] = node_elem_list_1[node_range[th_id]];
	}

	//std::cout << "\n" << node_elem_list_1[33] << ", "
	//	<< node_elem_list_1[34] << ", "
	//	<< node_elem_list_1[35] << ",\n";
	//std::cout << "\n" << node_elem_list_1[64] << ", "
	//	<< node_elem_list_1[65] << ", "
	//	<< node_elem_list_1[66] << ",\n\n";
	//for (size_t th_id = 0; th_id < thread_num; ++th_id)
	//{
	//	std::cout << th_id << ":\n"
	//		<< node_range[th_id] << ", "
	//		<< node_range[th_id + 1] << ",\n"
	//		<< node_elem_range[th_id] << ", "
	//		<< node_elem_range[th_id + 1] << ",\n";
	//}

	if (md.has_rigid_cylinder())
	{
		prc = &md.get_rigid_cylinder();
		prc->reset_f_cont();
		Kn_cont = md.Kn_cont;
		Kt_cont = md.Kt_cont;
		rc_fx_cont = 0.0;
		rc_fy_cont = 0.0;
		rc_fz_cont = 0.0;
		rc_mx_cont = 0.0;
		rc_my_cont = 0.0;
		rc_mz_cont = 0.0;
		contact_substep_id = md.contact_substep_id;
		prev_contact_pos = md.prev_contact_pos;
		prev_contact_force = md.prev_contact_force;
		memset(contact_substep_id, 0xFF, sizeof(size_t) * md.ori_pcl_num);
	}
	
	elem_count_bin = (size_t *)elem_bin_mem.alloc(sizeof(size_t) * thread_num * 0x100 * 2);
	elem_sum_bin = elem_count_bin + thread_num * 0x100;

	//for (size_t p_id = 0; p_id < md.pcl_num; ++p_id)
	//{
	//	Position& p_p = pcl_pos[p_id];
	//	size_t e_id = md.find_pcl_in_which_elem(p_p.x, p_p.y, p_p.z);
	//	res_file_t3d_me_mt << p_id << ", " << e_id << ",\n";
	//}

	pcl_num = 0;
#pragma omp parallel
	{
		size_t my_th_id = size_t(omp_get_thread_num());
		
		ThreadData &thd = thread_datas[my_th_id];
		thd.sorted_pcl_var_id = 0;

		size_t p_id0 = Block_Low(my_th_id, thread_num, md.pcl_num);
		size_t p_id1 = Block_Low(my_th_id + 1, thread_num, md.pcl_num);
		size_t pcl_in_mesh_num = 0;
		size_t p_id, ori_p_id, e_id;
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

			e_id = md.find_pcl_in_which_elem(p_p.x, p_p.y, p_p.z);
			spva0.pcl_in_elem[p_id] = e_id;
			if (e_id != elem_num)
				++pcl_in_mesh_num;
		}

#pragma omp critical
		{
			pcl_num += pcl_in_mesh_num;
		}

#pragma omp barrier

		//for (size_t p_id = 0; p_id < md.pcl_num; ++p_id)
		//{
		//	Position& p_p = pcl_pos[p_id];
		//	res_file_t3d_me_mt << spva0.pcl_index[p_id] << ", "
		//		<< spva0.pcl_in_elem[p_id] << ", "
		//		<< md.pcl_pos[spva0.pcl_index[p_id]].x << ",\n";
		//}

		// sort particle variables
		size_t* pcl_index0 = spva0.pcl_index;
		double* pcl_density0 = spva0.pcl_density;
		Displacement* pcl_disp0 = spva0.pcl_disp;
		Velocity* pcl_v0 = spva0.pcl_v;
		Stress* pcl_stress0 = spva0.pcl_stress;
		Strain* pcl_strain0 = spva0.pcl_strain;
		Strain* pcl_estrain0 = spva0.pcl_estrain;
		Strain* pcl_pstrain0 = spva0.pcl_pstrain;
		size_t* pcl_in_elem0 = spva0.pcl_in_elem;
		size_t* pcl_index1;
		double* pcl_density1;
		Displacement* pcl_disp1;
		Velocity* pcl_v1;
		Stress* pcl_stress1;
		Strain* pcl_strain1;
		Strain* pcl_estrain1;
		Strain* pcl_pstrain1;
		size_t* pcl_in_elem1;
		size_t data_digit, bin_id, th_id, pos_id;
		size_t* other_cbin;
		size_t* my_cbin = elem_count_bin + my_th_id * 0x100;
		size_t* my_sbin = elem_sum_bin + my_th_id * 0x100;
		for (size_t digit_disp = 0, elem_num_tmp = elem_num;
			elem_num_tmp; digit_disp += 8, elem_num_tmp >>= 8)
		{
			thd.sorted_pcl_var_id ^= 1;
			SortedPclVarArrays& spva1 = sorted_pcl_var_arrays[thd.sorted_pcl_var_id];
			pcl_index1 = spva1.pcl_index;
			pcl_density1 = spva1.pcl_density;
			pcl_disp1 = spva1.pcl_disp;
			pcl_v1 = spva1.pcl_v;
			pcl_stress1 = spva1.pcl_stress;
			pcl_strain1 = spva1.pcl_strain;
			pcl_estrain1 = spva1.pcl_estrain;
			pcl_pstrain1 = spva1.pcl_pstrain;
			pcl_in_elem1 = spva1.pcl_in_elem;
			memset(my_cbin, 0, 0x100 * sizeof(size_t));

			for (p_id = p_id0; p_id < p_id1; ++p_id)
			{
				data_digit = (pcl_in_elem0[p_id] >> digit_disp) & 0xFF;
				++my_cbin[data_digit];
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
				data_digit = (pcl_in_elem0[p_id] >> digit_disp) & 0xFF;
				pos_id = --my_sbin[data_digit];
				pcl_index1[pos_id] = pcl_index0[p_id];
				pcl_density1[pos_id] = pcl_density0[p_id];
				Displacement &p_d1 = pcl_disp1[pos_id];
				Displacement &p_d0 = pcl_disp0[p_id];
				p_d1.ux = p_d0.ux;
				p_d1.uy = p_d0.uy;
				p_d1.uz = p_d0.uz;
				Velocity& p_v1 = pcl_v1[pos_id];
				Velocity& p_v0 = pcl_v0[p_id];
				p_v1.vx = p_v0.vx;
				p_v1.vy = p_v0.vy;
				p_v1.vz = p_v0.vz;
				Stress& p_s1 = pcl_stress1[pos_id];
				Stress& p_s0 = pcl_stress0[p_id];
				p_s1.s11 = p_s0.s11;
				p_s1.s22 = p_s0.s22;
				p_s1.s33 = p_s0.s33;
				p_s1.s12 = p_s0.s12;
				p_s1.s23 = p_s0.s23;
				p_s1.s31 = p_s0.s31;
				Strain& p_e1 = pcl_strain1[pos_id];
				Strain &p_e0 = pcl_strain0[p_id];
				p_e1.e11 = p_e0.e11;
				p_e1.e22 = p_e0.e22;
				p_e1.e33 = p_e0.e33;
				p_e1.e12 = p_e0.e12;
				p_e1.e23 = p_e0.e23;
				p_e1.e31 = p_e0.e31;
				Strain& p_ee1 = pcl_estrain1[pos_id];
				Strain& p_ee0 = pcl_estrain0[p_id];
				p_ee1.e11 = p_ee0.e11;
				p_ee1.e22 = p_ee0.e22;
				p_ee1.e33 = p_ee0.e33;
				p_ee1.e12 = p_ee0.e12;
				p_ee1.e23 = p_ee0.e23;
				p_ee1.e31 = p_ee0.e31;
				Strain& p_pe1 = pcl_pstrain1[pos_id];
				Strain& p_pe0 = pcl_pstrain0[p_id];
				p_pe1.e11 = p_pe0.e11;
				p_pe1.e22 = p_pe0.e22;
				p_pe1.e33 = p_pe0.e33;
				p_pe1.e12 = p_pe0.e12;
				p_pe1.e23 = p_pe0.e23;
				p_pe1.e31 = p_pe0.e31;
				pcl_in_elem1[pos_id] = pcl_in_elem0[p_id];
			}

			pcl_index0 = pcl_index1;
			pcl_density0 = pcl_density1;
			pcl_disp0 = pcl_disp1;
			pcl_v0 = pcl_v1;
			pcl_stress0 = pcl_stress1;
			pcl_strain0 = pcl_strain1;
			pcl_estrain0 = pcl_estrain1;
			pcl_pstrain0 = pcl_pstrain1;
			pcl_in_elem0 = pcl_in_elem1;
#pragma omp barrier
		}
	}
	
	//ThreadData& thd = thread_datas[0];
	//SortedPclVarArrays& spva2 = sorted_pcl_var_arrays[thd.sorted_pcl_var_id];
	//for (size_t p_id = 0; p_id < md.pcl_num; ++p_id)
	//{
	//	Position& p_p = pcl_pos[p_id];
	//	res_file_t3d_me_mt << spva2.pcl_index[p_id] << ", "
	//		<< spva2.pcl_in_elem[p_id] << ", "
	//		<< md.pcl_pos[spva2.pcl_index[p_id]].x << ",\n";
	//}

	return 0;
}

int Step_T3D_ME_mt::finalize_calculation()
{
	Model_T3D_ME_mt &md = *(Model_T3D_ME_mt *)model;
	md.pcl_num = pcl_num;
	return 0;
}

int substep_func_omp_T3D_ME_mt(
	void* _self,
	size_t my_th_id,
	double dt,
	double cur_time,
	size_t substp_id
	)
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
	typedef Model_T3D_ME_mt::ElemNodeForce ElemNodeForce;
	typedef Model_T3D_ME_mt::NodeHasVBC NodeHasVBC;
	typedef Step_T3D_ME_mt::ThreadData ThreadData;
	typedef Step_T3D_ME_mt::SortedPclVarArrays SortedPclVarArrays;

	Step_T3D_ME_mt& self = *(Step_T3D_ME_mt*)(_self);
	Model_T3D_ME_mt& md = *(Model_T3D_ME_mt*)(self.model);
	
	double* pcl_m = self.pcl_m;
	Force *pcl_bf = self.pcl_bf;
	Force *pcl_t = self.pcl_t;
	Position* pcl_pos = self.pcl_pos;
	double* pcl_vol = self.pcl_vol;
	ShapeFunc* pcl_N = self.pcl_N;
	MatModel::MaterialModel** pcl_mat_model = self.pcl_mat_model;

	ThreadData& thd = self.thread_datas[my_th_id];
	SortedPclVarArrays *sorted_pcl_var_arrays = self.sorted_pcl_var_arrays;
	SortedPclVarArrays &spva0 = sorted_pcl_var_arrays[thd.sorted_pcl_var_id];
	size_t* pcl_index0 = spva0.pcl_index;
	double* pcl_density0 = spva0.pcl_density;
	Displacement* pcl_disp0 = spva0.pcl_disp;
	Velocity* pcl_v0 = spva0.pcl_v;
	Stress* pcl_stress0 = spva0.pcl_stress;
	Strain* pcl_strain0 = spva0.pcl_strain;
	Strain* pcl_estrain0 = spva0.pcl_estrain;
	Strain* pcl_pstrain0 = spva0.pcl_pstrain;
	size_t* pcl_in_elem0 = spva0.pcl_in_elem;
	
	size_t* pcl_in_elem_tmp = self.pcl_in_elem_tmp;

	ElemNodeIndex* elem_node_id = self.elem_node_id;
	size_t* elem_id_array = self.elem_id_array;
	size_t* node_elem_id_array = self.node_elem_id_array;
	size_t* node_elem_list = self.node_elem_list;
	DShapeFuncABC* elem_dN_abc = self.elem_dN_abc;
	DShapeFuncD* elem_dN_d = self.elem_dN_d;
	double* elem_vol = self.elem_vol;
	
	size_t* elem_substep_id = self.elem_substep_id;
	double* elem_density = self.elem_density;
	double* elem_pcl_m = self.elem_pcl_m;
	double* elem_pcl_vol = self.elem_pcl_vol;
	StrainInc* elem_de = self.elem_de;
	double* elem_m_de_vol = self.elem_m_de_vol;

	ElemNodeVM* elem_node_vm = self.elem_node_vm;
	ElemNodeForce* elem_node_force = self.elem_node_force;
	
	size_t* node_substep_id = self.node_substep_id;
	Acceleration* node_a = self.node_a;
	Velocity* node_v = self.node_v;
	NodeHasVBC* node_has_vbc = self.node_has_vbc;
	double* node_am = self.node_am;
	double* node_de_vol = self.node_de_vol;

	// cal p_id0, p_id1
	size_t pcl_num = self.pcl_num;
	size_t thread_num = self.thread_num;
	size_t p_id0, p_id1, e_id;
	p_id0 = Block_Low(my_th_id, thread_num, pcl_num);
	if (p_id0 != 0)
	{
		e_id = pcl_in_elem0[p_id0];
		while (e_id == pcl_in_elem0[p_id0 - 1])
			--p_id0;
	}
	p_id1 = Block_Low(my_th_id + 1, thread_num, pcl_num);
	if (p_id1 < pcl_num)
	{
		e_id = pcl_in_elem0[p_id1];
		while (e_id == pcl_in_elem0[p_id1 - 1])
			--p_id1;
	}

	// map pcl to mesh
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
	size_t p_id, ori_p_id;
	double p_m, p_vol, p_N_m;
	double one_fourth_bfx, one_fourth_bfy, one_fourth_bfz;
	e_id = pcl_in_elem0[p_id0];
	for (p_id = p_id0; p_id < p_id1; ++p_id)
	{
		// map pcl mass
		ori_p_id = pcl_index0[p_id];
		p_m = pcl_m[ori_p_id];
		e_p_m += p_m;

		// map pcl volume
		p_vol = p_m / pcl_density0[p_id];
		pcl_vol[p_id] = p_vol;
		e_p_vol += p_vol;

		// map stress
		Stress &p_s = pcl_stress0[p_id];
		e_s11 += p_s.s11 * p_vol;
		e_s22 += p_s.s22 * p_vol;
		e_s33 += p_s.s33 * p_vol;
		e_s12 += p_s.s12 * p_vol;
		e_s23 += p_s.s23 * p_vol;
		e_s31 += p_s.s31 * p_vol;

		// cal shape function
		Position& p_p = pcl_pos[ori_p_id];
		Displacement& p_d = pcl_disp0[p_id];
		ShapeFunc &p_N = pcl_N[p_id];
		md.cal_N(p_p.x + p_d.ux, p_p.y + p_d.uy, p_p.z + p_d.uz, e_id, p_N);
		// map velocity
		Velocity &p_v = pcl_v0[p_id];
		p_N_m = p_N.N1 * p_m;
		en1_vm += p_N_m;
		en1_vmx += p_N_m * p_v.vx;
		en1_vmy += p_N_m * p_v.vy;
		en1_vmz += p_N_m * p_v.vz;
		p_N_m = p_N.N2 * p_m;
		en2_vm += p_N_m;
		en2_vmx += p_N_m * p_v.vx;
		en2_vmy += p_N_m * p_v.vy;
		en2_vmz += p_N_m * p_v.vz;
		p_N_m = p_N.N3 * p_m;
		en3_vm += p_N_m;
		en3_vmx += p_N_m * p_v.vx;
		en3_vmy += p_N_m * p_v.vy;
		en3_vmz += p_N_m * p_v.vz;
		p_N_m = p_N.N4 * p_m;
		en4_vm += p_N_m;
		en4_vmx += p_N_m * p_v.vx;
		en4_vmy += p_N_m * p_v.vy;
		en4_vmz += p_N_m * p_v.vz;

		// cal external load
		Force& p_bf = pcl_bf[ori_p_id];
		one_fourth_bfx = one_fourth * p_bf.fx;
		one_fourth_bfy = one_fourth * p_bf.fy;
		one_fourth_bfz = one_fourth * p_bf.fz;
		Force& p_t = pcl_t[ori_p_id];
		en1_fx += one_fourth_bfx + p_N.N1 * p_t.fx;
		en1_fy += one_fourth_bfy + p_N.N1 * p_t.fy;
		en1_fz += one_fourth_bfz + p_N.N1 * p_t.fz;
		en2_fx += one_fourth_bfx + p_N.N2 * p_t.fx;
		en2_fy += one_fourth_bfy + p_N.N2 * p_t.fy;
		en2_fz += one_fourth_bfz + p_N.N2 * p_t.fz;
		en3_fx += one_fourth_bfx + p_N.N3 * p_t.fx;
		en3_fy += one_fourth_bfy + p_N.N3 * p_t.fy;
		en3_fz += one_fourth_bfz + p_N.N3 * p_t.fz;
		en4_fx += one_fourth_bfx + p_N.N4 * p_t.fx;
		en4_fy += one_fourth_bfy + p_N.N4 * p_t.fy;
		en4_fz += one_fourth_bfz + p_N.N4 * p_t.fz;
	
		if (e_id != pcl_in_elem0[p_id + 1])
		{
			elem_substep_id[e_id] = substp_id;

			elem_pcl_m[e_id] = e_p_m;
			elem_density[e_id] = e_p_m / e_p_vol;
			
			ElemNodeVM& en1_v = elem_node_vm[e_id * 4];
			en1_v.vm = en1_vm;
			en1_v.vmx = en1_vmx;
			en1_v.vmy = en1_vmy;
			en1_v.vmz = en1_vmz;
			ElemNodeVM& en2_v = *(&en1_v + 1);
			en2_v.vm = en2_vm;
			en2_v.vmx = en2_vmx;
			en2_v.vmy = en2_vmy;
			en2_v.vmz = en2_vmz;
			ElemNodeVM& en3_v = *(&en2_v + 1);
			en3_v.vm = en3_vm;
			en3_v.vmx = en3_vmx;
			en3_v.vmy = en3_vmy;
			en3_v.vmz = en3_vmz;
			ElemNodeVM& en4_v = *(&en3_v + 1);
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
			DShapeFuncABC &e_dN = elem_dN_abc[e_id];
			// node 1
			ElemNodeForce& en1_f = elem_node_force[e_id * 4];
			en1_fx -= (e_dN.dN1_dx * e_s11 + e_dN.dN1_dy * e_s12 + e_dN.dN1_dz * e_s31) * e_p_vol;
			en1_f.fx = en1_fx;
			en1_fy -= (e_dN.dN1_dx * e_s12 + e_dN.dN1_dy * e_s22 + e_dN.dN1_dz * e_s23) * e_p_vol;
			en1_f.fy = en1_fy;
			en1_fz -= (e_dN.dN1_dx * e_s31 + e_dN.dN1_dy * e_s23 + e_dN.dN1_dz * e_s33) * e_p_vol;
			en1_f.fz = en1_fz;
			// node 2
			ElemNodeForce& en2_f = *(&en1_f + 1);
			en2_fx -= (e_dN.dN2_dx * e_s11 + e_dN.dN2_dy * e_s12 + e_dN.dN2_dz * e_s31) * e_p_vol;
			en2_f.fx = en2_fx;
			en2_fy -= (e_dN.dN2_dx * e_s12 + e_dN.dN2_dy * e_s22 + e_dN.dN2_dz * e_s23) * e_p_vol;
			en2_f.fy = en2_fy;
			en2_fz -= (e_dN.dN2_dx * e_s31 + e_dN.dN2_dy * e_s23 + e_dN.dN2_dz * e_s33) * e_p_vol;
			en2_f.fz = en2_fz;
			// node 3
			ElemNodeForce& en3_f = *(&en2_f + 1);
			en3_fx -= (e_dN.dN3_dx * e_s11 + e_dN.dN3_dy * e_s12 + e_dN.dN3_dz * e_s31) * e_p_vol;
			en3_f.fx = en3_fx;
			en3_fy -= (e_dN.dN3_dx * e_s12 + e_dN.dN3_dy * e_s22 + e_dN.dN3_dz * e_s23) * e_p_vol;
			en3_f.fy = en3_fy;
			en3_fz -= (e_dN.dN3_dx * e_s31 + e_dN.dN3_dy * e_s23 + e_dN.dN3_dz * e_s33) * e_p_vol;
			en3_f.fz = en3_fz;
			// node 4
			ElemNodeForce& en4_f = *(&en3_f + 1);
			en4_fx -= (e_dN.dN4_dx * e_s11 + e_dN.dN4_dy * e_s12 + e_dN.dN4_dz * e_s31) * e_p_vol;
			en4_f.fx = en4_fx;
			en4_fy -= (e_dN.dN4_dx * e_s12 + e_dN.dN4_dy * e_s22 + e_dN.dN4_dz * e_s23) * e_p_vol;
			en4_f.fy = en4_fy;
			en4_fz -= (e_dN.dN4_dx * e_s31 + e_dN.dN4_dy * e_s23 + e_dN.dN4_dz * e_s33) * e_p_vol;
			en4_f.fz = en4_fz;

			// reset
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
			
			e_id = pcl_in_elem0[p_id + 1];
		}
	}

	if (md.has_rigid_cylinder())
	{
		ContactForce3D cf;
		cf.reset();
		self.apply_rigid_cylinder(
			p_id0, p_id1,
			pcl_in_elem0,
			spva0, cf,
			substp_id
			);
#pragma omp critical
		{
			self.prc->combine_f_cont(cf);
		}
	}

#pragma omp barrier

#pragma omp master
	{
		if (md.has_rigid_cylinder())
		{
			RigidCylinder& rc = *(self.prc);
			rc.update_motion(dt);
			const ContactForce3D& cf = rc.get_cont_force();
			self.rc_fx_cont = cf.fx;
			self.rc_fy_cont = cf.fy;
			self.rc_fz_cont = cf.fz;
			self.rc_mx_cont = cf.mx;
			self.rc_my_cont = cf.my;
			self.rc_mz_cont = cf.mz;
		}
	}

	// update node variables
	size_t* node_range = self.node_range;
	size_t n_id0 = node_range[my_th_id];
	size_t n_id1 = node_range[my_th_id + 1];
	size_t* node_elem_range = self.node_elem_range;
	size_t ne_id = node_elem_range[my_th_id];
	size_t n_id, ne_id1, n_var_id, bc_mask;
	double n_am, n_fx, n_fy, n_fz;
	double n_vm, n_vmx, n_vmy, n_vmz;
	size_t* pn_substp_id;
	for (n_id = n_id0; n_id < n_id1; ++n_id)
	{
		n_am = 0.0;
		n_fx = 0.0;
		n_fy = 0.0;
		n_fz = 0.0;
		n_vm = 0.0;
		n_vmx = 0.0;
		n_vmy = 0.0;
		n_vmz = 0.0;
		pn_substp_id = &node_substep_id[n_id];
		ne_id1 = node_elem_list[n_id];
		for (; ne_id < ne_id1; ++ne_id)
		{
			e_id = elem_id_array[ne_id];
			if (elem_substep_id[e_id] == substp_id)
			{
				*pn_substp_id = substp_id;
				n_am += elem_pcl_m[e_id];
				n_var_id = node_elem_id_array[ne_id];
				ElemNodeForce& nf = elem_node_force[n_var_id];
				n_fx += nf.fx;
				n_fy += nf.fy;
				n_fz += nf.fz;
				ElemNodeVM &nvm = elem_node_vm[n_var_id];
				n_vm += nvm.vm;
				n_vmx += nvm.vmx;
				n_vmy += nvm.vmy;
				n_vmz += nvm.vmz;
			}
		}
		if (*pn_substp_id == substp_id)
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
			NodeHasVBC& n_has_vbc = node_has_vbc[n_id];
			bc_mask = size_t(n_has_vbc.has_vx_bc) + SIZE_MAX;
			n_a.iax &= bc_mask;
			n_v.ivx &= bc_mask;
			bc_mask = size_t(n_has_vbc.has_vy_bc) + SIZE_MAX;
			n_a.iay &= bc_mask;
			n_v.ivy &= bc_mask;
			bc_mask = size_t(n_has_vbc.has_vz_bc) + SIZE_MAX;
			n_a.iaz &= bc_mask;
			n_v.ivz &= bc_mask;
		}
	}

#pragma omp master
	{
//		if (md.has_rigid_rect())
//		{
//			RigidRect& rr = md.get_rigid_rect();
//			rr.update_motion(dt);
//			self.rr_fx_cont = rr.get_fx_contact();
//			self.rr_fy_cont = rr.get_fy_contact();
//			self.rr_m_cont = rr.get_m_contact();
//			rr.reset_f_contact();
//		}

		self.pcl_num = 0;
	}

#pragma omp barrier
	
	// update particle variables
	Acceleration *pn_a1, *pn_a2, *pn_a3, *pn_a4;
	Velocity *pn_v1, *pn_v2, *pn_v3, *pn_v4;
	double p_x, p_y, p_z, e_de_vol;
	size_t p_e_id;
	size_t elem_num = self.elem_num;
	e_id = elem_num;
	size_t pcl_in_mesh_num = 0;
	for (p_id = p_id0; p_id < p_id1; ++p_id)
	{
		if (e_id != pcl_in_elem0[p_id])
		{
			e_id = pcl_in_elem0[p_id];

			ElemNodeIndex& e_n_id = elem_node_id[e_id];
			pn_a1 = node_a + e_n_id.n1;
			pn_a2 = node_a + e_n_id.n2;
			pn_a3 = node_a + e_n_id.n3;
			pn_a4 = node_a + e_n_id.n4;

			pn_v1 = node_v + e_n_id.n1;
			pn_v2 = node_v + e_n_id.n2;
			pn_v3 = node_v + e_n_id.n3;
			pn_v4 = node_v + e_n_id.n4;

			DShapeFuncABC& e_dN = elem_dN_abc[e_id];
			StrainInc& e_de = elem_de[e_id];
			e_de.de11 = (e_dN.dN1_dx * pn_v1->vx + e_dN.dN2_dx * pn_v2->vx + e_dN.dN3_dx * pn_v3->vx + e_dN.dN4_dx * pn_v4->vx) * dt;
			e_de.de22 = (e_dN.dN1_dy * pn_v1->vy + e_dN.dN2_dy * pn_v2->vy + e_dN.dN3_dy * pn_v3->vy + e_dN.dN4_dy * pn_v4->vy) * dt;
			e_de.de33 = (e_dN.dN1_dz * pn_v1->vz + e_dN.dN2_dz * pn_v2->vz + e_dN.dN3_dz * pn_v3->vz + e_dN.dN4_dz * pn_v4->vz) * dt;
			e_de.de12 = (e_dN.dN1_dx * pn_v1->vy + e_dN.dN2_dx * pn_v2->vy + e_dN.dN3_dx * pn_v3->vy + e_dN.dN4_dx * pn_v4->vy
					   + e_dN.dN1_dy * pn_v1->vx + e_dN.dN2_dy * pn_v2->vx + e_dN.dN3_dy * pn_v3->vx + e_dN.dN4_dy * pn_v4->vx) * dt * 0.5;
			e_de.de23 = (e_dN.dN1_dy * pn_v1->vz + e_dN.dN2_dy * pn_v2->vz + e_dN.dN3_dy * pn_v3->vz + e_dN.dN4_dy * pn_v4->vz
					   + e_dN.dN1_dz * pn_v1->vy + e_dN.dN2_dz * pn_v2->vy + e_dN.dN3_dz * pn_v3->vy + e_dN.dN4_dz * pn_v4->vy) * dt * 0.5;
			e_de.de31 = (e_dN.dN1_dz * pn_v1->vx + e_dN.dN2_dz * pn_v2->vx + e_dN.dN3_dz * pn_v3->vx + e_dN.dN4_dz * pn_v4->vx
					   + e_dN.dN1_dx * pn_v1->vz + e_dN.dN2_dx * pn_v2->vz + e_dN.dN3_dx * pn_v3->vz + e_dN.dN4_dx * pn_v4->vz) * dt * 0.5;
			e_de_vol = e_de.de11 + e_de.de22 + e_de.de33;
			elem_m_de_vol[e_id] = elem_pcl_m[e_id] * e_de_vol;
			e_de_vol *= one_third;
			e_de.de11 -= e_de_vol;
			e_de.de22 -= e_de_vol;
			e_de.de33 -= e_de_vol;
		}

		// update velocity
		ShapeFunc& p_N = pcl_N[p_id];
		Velocity& p_v = pcl_v0[p_id];
		p_v.vx += (p_N.N1 * pn_a1->ax + p_N.N2 * pn_a2->ax + p_N.N3 * pn_a3->ax + p_N.N4 * pn_a4->ax) * dt;
		p_v.vy += (p_N.N1 * pn_a1->ay + p_N.N2 * pn_a2->ay + p_N.N3 * pn_a3->ay + p_N.N4 * pn_a4->ay) * dt;
		p_v.vz += (p_N.N1 * pn_a1->az + p_N.N2 * pn_a2->az + p_N.N3 * pn_a3->az + p_N.N4 * pn_a4->az) * dt;

		// update displacement
		Displacement& p_d = pcl_disp0[p_id];
		p_d.ux += (p_N.N1 * pn_v1->vx + p_N.N2 * pn_v2->vx + p_N.N3 * pn_v3->vx + p_N.N4 * pn_v4->vx) * dt;
		p_d.uy += (p_N.N1 * pn_v1->vy + p_N.N2 * pn_v2->vy + p_N.N3 * pn_v3->vy + p_N.N4 * pn_v4->vy) * dt;
		p_d.uz += (p_N.N1 * pn_v1->vz + p_N.N2 * pn_v2->vz + p_N.N3 * pn_v3->vz + p_N.N4 * pn_v4->vz) * dt;

		// update location (in which element)
		Position& p_p = pcl_pos[ori_p_id];
		p_x = p_p.x + p_d.ux;
		p_y = p_p.y + p_d.uy;
		p_z = p_p.z + p_d.uz;
		p_e_id = e_id;
		if (!md.is_in_element(p_x, p_y, p_z, e_id))
			p_e_id = md.find_pcl_in_which_elem(p_x, p_y, p_z);
		if (p_e_id != elem_num)
			++pcl_in_mesh_num;
		pcl_in_elem_tmp[p_id] = p_e_id;
	}

#pragma omp critical
	{
		self.pcl_num += pcl_in_mesh_num;
	}

#pragma omp barrier

	double n_am_de_vol;
	ne_id = node_elem_range[my_th_id];
	for (n_id = n_id0; n_id < n_id1; ++n_id)
	{
		if (node_substep_id[n_id] == substp_id)
		{
			n_am_de_vol = 0.0;
			ne_id1 = node_elem_list[n_id];
			for (; ne_id < ne_id1; ++ne_id)
			{
				e_id = elem_id_array[ne_id];
				if (elem_substep_id[e_id] == substp_id)
					n_am_de_vol += elem_m_de_vol[e_id];
			}
			node_de_vol[n_id] = n_am_de_vol * one_fourth / node_am[n_id];
		}
	}

#pragma omp barrier

	StrainInc* pe_de;
	int mm_res;
	const double* dstress;
	e_id = elem_num;
	for (p_id = p_id0; p_id < p_id1; ++p_id)
	{
		if (e_id != pcl_in_elem0[p_id])
		{
			e_id = pcl_in_elem0[p_id];

			ElemNodeIndex &e_n_id = elem_node_id[e_id];
			e_de_vol = one_fourth *
				 (node_de_vol[e_n_id.n1] + node_de_vol[e_n_id.n2]
				+ node_de_vol[e_n_id.n3] + node_de_vol[e_n_id.n4]);
			elem_density[e_id] /= 1.0 + e_de_vol;
			
			pe_de = elem_de + e_id;
			e_de_vol *= one_third;
			pe_de->de11 += e_de_vol;
			pe_de->de22 += e_de_vol;
			pe_de->de33 += e_de_vol;
		}

		// update density
		pcl_density0[p_id] = elem_density[e_id];

		// update stress
		MatModel::MaterialModel& pcl_mm = *pcl_mat_model[pcl_index0[p_id]];
		mm_res = pcl_mm.integrate(pe_de->de);
		dstress = pcl_mm.get_dstress();
		Stress& p_s = pcl_stress0[p_id];
		p_s.s11 += dstress[0];
		p_s.s22 += dstress[1];
		p_s.s33 += dstress[2];
		p_s.s12 += dstress[3];
		p_s.s23 += dstress[4];
		p_s.s31 += dstress[5];
	}

#pragma omp barrier

#pragma omp master
	{
		if (self.pcl_num)
			self.continue_calculation();
		else
			self.exit_calculation();
	}

	pcl_in_elem0 = pcl_in_elem_tmp;
	
	// sort particle variables
	size_t *pcl_index1;
	double *pcl_density1;
	Displacement* pcl_disp1;
	Velocity* pcl_v1;
	Stress* pcl_stress1;
	Strain* pcl_strain1;
	Strain* pcl_estrain1;
	Strain* pcl_pstrain1;
	size_t* pcl_in_elem1;
	size_t data_digit, bin_id, th_id, pos_id;
	size_t *other_cbin;
	size_t *elem_count_bin = self.elem_count_bin;
	size_t *my_cbin = elem_count_bin + my_th_id * 0x100;
	size_t *elem_sum_bin = self.elem_sum_bin;
	size_t *my_sbin = elem_sum_bin + my_th_id * 0x100;
	for (size_t digit_disp = 0, elem_num_tmp = elem_num;
		 elem_num_tmp; digit_disp += 8, elem_num_tmp >>= 8)
	{
		thd.sorted_pcl_var_id ^= 1;
		SortedPclVarArrays& spva1 = sorted_pcl_var_arrays[thd.sorted_pcl_var_id];
		pcl_index1 = spva1.pcl_index;
		pcl_density1 = spva1.pcl_density;
		pcl_disp1 = spva1.pcl_disp;
		pcl_v1 = spva1.pcl_v;
		pcl_stress1 = spva1.pcl_stress;
		pcl_strain1 = spva1.pcl_strain;
		pcl_estrain1 = spva1.pcl_estrain;
		pcl_pstrain1 = spva1.pcl_pstrain;
		pcl_in_elem1 = spva1.pcl_in_elem;
		memset(my_cbin, 0, 0x100 * sizeof(size_t));

		for (p_id = p_id0; p_id < p_id1; ++p_id)
		{
			data_digit = (pcl_in_elem0[p_id] >> digit_disp) & 0xFF;
			++my_cbin[data_digit];
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
			data_digit = (pcl_in_elem0[p_id] >> digit_disp) & 0xFF;
			pos_id = --my_sbin[data_digit];
			pcl_index1[pos_id] = pcl_index0[p_id];
			pcl_density1[pos_id] = pcl_density0[p_id];
			Displacement& p_d1 = pcl_disp1[pos_id];
			Displacement& p_d0 = pcl_disp0[p_id];
			p_d1.ux = p_d0.ux;
			p_d1.uy = p_d0.uy;
			p_d1.uz = p_d0.uz;
			Velocity& p_v1 = pcl_v1[pos_id];
			Velocity& p_v0 = pcl_v0[p_id];
			p_v1.vx = p_v0.vx;
			p_v1.vy = p_v0.vy;
			p_v1.vz = p_v0.vz;
			Stress& p_s1 = pcl_stress1[pos_id];
			Stress& p_s0 = pcl_stress0[p_id];
			p_s1.s11 = p_s0.s11;
			p_s1.s22 = p_s0.s22;
			p_s1.s33 = p_s0.s33;
			p_s1.s12 = p_s0.s12;
			p_s1.s23 = p_s0.s23;
			p_s1.s31 = p_s0.s31;
			Strain& p_e1 = pcl_strain1[pos_id];
			Strain& p_e0 = pcl_strain0[p_id];
			p_e1.e11 = p_e0.e11;
			p_e1.e22 = p_e0.e22;
			p_e1.e33 = p_e0.e33;
			p_e1.e12 = p_e0.e12;
			p_e1.e23 = p_e0.e23;
			p_e1.e31 = p_e0.e31;
			Strain& p_ee1 = pcl_estrain1[pos_id];
			Strain& p_ee0 = pcl_estrain0[p_id];
			p_ee1.e11 = p_ee0.e11;
			p_ee1.e22 = p_ee0.e22;
			p_ee1.e33 = p_ee0.e33;
			p_ee1.e12 = p_ee0.e12;
			p_ee1.e23 = p_ee0.e23;
			p_ee1.e31 = p_ee0.e31;
			Strain& p_pe1 = pcl_pstrain1[pos_id];
			Strain& p_pe0 = pcl_pstrain0[p_id];
			p_pe1.e11 = p_pe0.e11;
			p_pe1.e22 = p_pe0.e22;
			p_pe1.e33 = p_pe0.e33;
			p_pe1.e12 = p_pe0.e12;
			p_pe1.e23 = p_pe0.e23;
			p_pe1.e31 = p_pe0.e31;
			pcl_in_elem1[pos_id] = pcl_in_elem0[p_id];
		}

		pcl_index0 = pcl_index1;
		pcl_density0 = pcl_density1;
		pcl_disp0 = pcl_disp1;
		pcl_v0 = pcl_v1;
		pcl_stress0 = pcl_stress1;
		pcl_strain0 = pcl_strain1;
		pcl_estrain0 = pcl_estrain1;
		pcl_pstrain0 = pcl_pstrain1;
		pcl_in_elem0 = pcl_in_elem1;
#pragma omp barrier
	}

	return 0;
}

int Step_T3D_ME_mt::apply_rigid_cylinder(
	size_t p_id0,
	size_t p_id1,
	size_t* pcl_in_elem,
	SortedPclVarArrays &cur_spva,
	ContactForce3D& rc_cf,
	size_t substp_id
	) noexcept
{
	double p_x, p_y, p_z;
	double f_cont, fx_cont, fy_cont, fz_cont;
	size_t e_id, ori_p_id;
	double dist;
	Vector3D lnorm, gnorm, ft_cont;
	Point3D lcontpos;
	for (size_t p_id = p_id0; p_id < p_id1; ++p_id)
	{
		ori_p_id = cur_spva.pcl_index[p_id];
		Position &p_p = pcl_pos[ori_p_id];
		Displacement &p_d = cur_spva.pcl_disp[p_id];
		p_x = p_p.x + p_d.ux;
		p_y = p_p.y + p_d.uy;
		p_z = p_p.z + p_d.uz;
		if (prc->detect_collision_with_point(
			p_x, p_y, p_z, pcl_vol[p_id],
			dist, lnorm, lcontpos))
		{
			// normal force
			prc->get_global_vector(lnorm, gnorm);
			f_cont = Kn_cont * dist;
			fx_cont = f_cont * gnorm.x;
			fy_cont = f_cont * gnorm.y;
			fz_cont = f_cont * gnorm.z;
			// tangential force
			if (contact_substep_id[ori_p_id] == substp_id)
			{
				Position& cp = prev_contact_pos[ori_p_id];
				cp.x = lcontpos.x;
				cp.y = lcontpos.y;
				cp.z = lcontpos.z;
				Force& cf = prev_contact_force[ori_p_id];
				cf.fx = 0.0;
				cf.fy = 0.0;
				cf.fz = 0.0;
			}
			else
			{
				Position& cp = prev_contact_pos[ori_p_id];
				double rx = lcontpos.x - cp.x;
				double ry = lcontpos.y - cp.y;
				double rz = lcontpos.z - cp.z;
				cp.x = lcontpos.x;
				cp.y = lcontpos.y;
				cp.z = lcontpos.z;
				double norm_len = rx * gnorm.x + ry * gnorm.y + rz * gnorm.z;
				rx -= norm_len * gnorm.x;
				ry -= norm_len * gnorm.y;
				rz -= norm_len * gnorm.z;
				Force& cf = prev_contact_force[ori_p_id];
				cf.fx += rx * Kt_cont;
				cf.fy += ry * Kt_cont;
				cf.fz += rz * Kt_cont;
				// apply contact consitutive model here
				fx_cont += cf.fx;
				fy_cont += cf.fy;
				fz_cont += cf.fz;
			}
			contact_substep_id[ori_p_id] = substp_id + 1;
			// apply contact force to rigid body
			const Point3D& rc_cen = prc->get_centre();
			rc_cf.add_force(
				p_x, p_y, p_z,
				-fx_cont, -fy_cont, -fz_cont,
				rc_cen.x, rc_cen.y, rc_cen.z 
				);
			// apply contact force to mesh
			ShapeFunc& p_N = pcl_N[p_id];
			e_id = pcl_in_elem[p_id];
			ElemNodeForce& en_f1 = elem_node_force[e_id * 4];
			en_f1.fx += p_N.N1 * fx_cont;
			en_f1.fy += p_N.N1 * fy_cont;
			en_f1.fz += p_N.N1 * fz_cont;
			ElemNodeForce& en_f2 = elem_node_force[e_id * 4 + 1];
			en_f2.fx += p_N.N2 * fx_cont;
			en_f2.fy += p_N.N2 * fy_cont;
			en_f2.fz += p_N.N2 * fz_cont;
			ElemNodeForce& en_f3 = elem_node_force[e_id * 4 + 2];
			en_f3.fx += p_N.N3 * fx_cont;
			en_f3.fy += p_N.N3 * fy_cont;
			en_f3.fz += p_N.N3 * fz_cont;
			ElemNodeForce& en_f4 = elem_node_force[e_id * 4 + 3];
			en_f4.fx += p_N.N4 * fx_cont;
			en_f4.fy += p_N.N4 * fy_cont;
			en_f4.fz += p_N.N4 * fz_cont;
		}
	}

	return 0;
}
