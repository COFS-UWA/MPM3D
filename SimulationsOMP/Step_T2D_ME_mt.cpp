#include "SimulationsOMP_pcp.h"

#include <iostream>
#include <omp.h>

#include "Step_T2D_ME_mt.h"

#define one_third (1.0f/3.0f)
#define N_min (1.0e-10f)
#define Block_Low(th_id, th_num, data_num) ((th_id)*(data_num)/(th_num))

Step_T2D_ME_mt::Step_T2D_ME_mt(const char* _name) : 
	Step_OMP(_name, "Step_T2D_ME_mt", &substep_func_omp_T2D_ME_mt) {}

Step_T2D_ME_mt::~Step_T2D_ME_mt() {}

namespace
{
	inline void divide_task_to_thread(
		uint32_t thread_num,
		uint32_t data_num,
		uint32_t data_range[] // len = thread_num+1
		)
	{
		data_range[0] = 0;
		data_range[thread_num] = data_num;
		register uint32_t i;
		for (i = 1; i < thread_num; ++i)
			data_range[i] = Block_Low(i, thread_num, data_num);
	}

	inline void swap(size_t &a, size_t &b)
	{
		a = a ^ b;
		b = a ^ b;
		a = a ^ b;
	}
}

int Step_T2D_ME_mt::init_calculation()
{
	Model_T2D_ME_mt& md = *(Model_T2D_ME_mt*)model;

	omp_set_num_threads(thread_num);

	uint32_t th_num_1 = thread_num + 1;
	char *mem_range = (char *)task_range_mem.alloc((sizeof(PclRange) + sizeof(ElemRange) + sizeof(uint32_t) * 2) * size_t(th_num_1));
	pcl_range = (PclRange *)mem_range;
	mem_range += sizeof(PclRange) * th_num_1;
	elem_range = (ElemRange *)mem_range;
	mem_range += sizeof(ElemRange) * th_num_1;
	node_elem_range = (uint32_t *)mem_range;
	mem_range += sizeof(uint32_t) * th_num_1;
	node_range = (uint32_t*)mem_range;

	pcl_range[0].id = 0;
	elem_range[0].id = 0;
	divide_task_to_thread(thread_num, md.node_num, node_range);
	node_elem_range[0] = 0;
	for (uint32_t th_id = 0; th_id < thread_num; ++th_id)
		node_elem_range[th_id + 1] = md.node_elem_list[node_range[th_id]];

	new_to_ori_pcl_map0 = (uint32_t *)sort_var_mem.alloc(sizeof(uint32_t) * (md.pcl_num * 4 + 2) + Cache_Alignment * 3);
	new_to_ori_pcl_map1 = cache_aligned(new_to_ori_pcl_map0 + md.pcl_num);
	pcl_in_elem_array0 = cache_aligned(new_to_ori_pcl_map1 + md.pcl_num + 1);
	pcl_in_elem_array1 = cache_aligned(pcl_in_elem_array0 + md.pcl_num + 1);
	pcl_in_elem_array0[-1] = UINT32_MAX;
	pcl_in_elem_array1[-1] = UINT32_MAX;

	elem_count_bin = (uint32_t *)elem_bin_mem.alloc(sizeof(uint32_t) * thread_num * 0x100 * 2);
	elem_sum_bin = elem_count_bin + size_t(thread_num) * 0x100;

	// accelerate substep function
	pcl_num = md.pcl_num;
	elem_num = md.elem_num;
	node_num = md.node_num;

	pcl_m = md.pcl_m;
	pcl_bf = md.pcl_bf;
	pcl_t = md.pcl_t;
	pcl_pos = md.pcl_pos;
	pcl_mat_model = md.pcl_mat_model;

	pcl_sorted_var_id = 0;

	PclSortedVarArray& md_pscv0 = md.pcl_sorted_var_array[0];
	PclSortedVarArray& pscv0 = pcl_sorted_var_array[0];
	pscv0.pcl_index = md_pscv0.pcl_index;
	pscv0.pcl_density = md_pscv0.pcl_density;
	pscv0.pcl_disp = md_pscv0.pcl_disp;
	pscv0.pcl_v = md_pscv0.pcl_v;
	pscv0.pcl_N = md_pscv0.pcl_N;
	pscv0.pcl_stress = md_pscv0.pcl_stress;
	pscv0.elem_has_pcl_num = md_pscv0.elem_has_pcl_num;

	PclSortedVarArray& md_pscv1 = md.pcl_sorted_var_array[1];
	PclSortedVarArray& pscv1 = pcl_sorted_var_array[1];
	pscv1.pcl_index = md_pscv1.pcl_index;
	pscv1.pcl_density = md_pscv1.pcl_density;
	pscv1.pcl_disp = md_pscv1.pcl_disp;
	pscv1.pcl_v = md_pscv1.pcl_v;
	pscv1.pcl_N = md_pscv1.pcl_N;
	pscv1.pcl_stress = md_pscv1.pcl_stress;
	pscv1.elem_has_pcl_num = md_pscv1.elem_has_pcl_num;

	elem_node_id = md.elem_node_id;
	elem_area = md.elem_area;
	elem_sf_ab = md.elem_sf_ab;
	elem_sf_c = md.elem_sf_c;

	elem_density = md.elem_density;
	elem_de = md.elem_de;
	elem_stress = md.elem_stress;
	elem_am = md.elem_am;
	elem_am_de_vol = md.elem_am_de_vol;

	elem_node_vm = md.elem_node_vm;
	elem_node_force = md.elem_node_force;

	elem_id_array = md.elem_id_array;
	node_elem_id_array = md.node_elem_id_array;
	node_elem_list = md.node_elem_list;
	node_a = md.node_a;
	node_v = md.node_v;
	node_has_vbc = md.node_has_vbc;
	node_am = md.node_am;
	node_de_vol = md.node_de_vol;

	union
	{
		struct
		{
			uint32_t* new_to_ori_pcl_map;
			uint32_t* new_to_ori_pcl_map_tmp;
			uint32_t* pcl_in_elem_array;
			uint32_t* pcl_in_elem_array_tmp;
		};
		struct
		{
			size_t new_to_ori_pcl_map_ui;
			size_t new_to_ori_pcl_map_tmp_ui;
			size_t pcl_in_elem_array_ui;
			size_t pcl_in_elem_array_tmp_ui;
		};
	};
	new_to_ori_pcl_map = new_to_ori_pcl_map0;
	new_to_ori_pcl_map_tmp = new_to_ori_pcl_map1;
	pcl_in_elem_array = pcl_in_elem_array0;
	pcl_in_elem_array_tmp = pcl_in_elem_array1;
	
	uint32_t th_id;
	for (th_id = 0; th_id < thread_num; ++th_id)
		pcl_range[th_id+1].id = Block_Low(th_id, thread_num, pcl_num);
	PclSortedVarArray &psva = md.pcl_sorted_var_array[0];
	uint32_t pcl_in_elem_id;
	uint32_t pcl_num1 = 0;
#pragma omp parallel reduction(+:pcl_num1)
	{
		uint32_t my_th_id = omp_get_thread_num();
		uint32_t p_id0 = pcl_range[my_th_id].id;
		uint32_t p_id1 = pcl_range[my_th_id + 1].id;
		uint32_t pcl_in_mesh_num = 0;
		for (uint32_t p_id = p_id0; p_id < p_id1; ++p_id)
		{
			PclDisp& p_d = psva.pcl_disp[p_id];
			p_d.ux = 0.0f;
			p_d.uy = 0.0f;

			PclPos& p_p = md.pcl_pos[p_id];
			PclShapeFunc& p_N = psva.pcl_N[p_id];
			pcl_in_elem_id = md.find_pcl_in_which_elem(p_p.x, p_p.y, p_N);
			if (pcl_in_elem_id != elem_num)
			{
				if (p_N.N1 < N_min)
					p_N.N1 = N_min;
				if (p_N.N2 < N_min)
					p_N.N2 = N_min;
				if (p_N.N3 < N_min)
					p_N.N3 = N_min;
				++pcl_in_mesh_num;
			}
			new_to_ori_pcl_map[p_id] = p_id;
			pcl_in_elem_array[p_id] = pcl_in_elem_id;
		}
		pcl_num1 += pcl_in_mesh_num;

		uint32_t* my_cbin = elem_count_bin + size_t(my_th_id) * 0x100;
		uint32_t* my_sbin = elem_sum_bin + size_t(my_th_id) * 0x100;
		for (uint32_t digit_disp = 0; digit_disp < 32; digit_disp += 8)
		{
			// cal historgram
			memset(my_cbin, 0, sizeof(uint32_t) * 0x100);

			uint32_t data_digit;
			for (uint32_t p_id = p_id0; p_id < p_id1; ++p_id)
			{
				data_digit = (pcl_in_elem_array[p_id] >> digit_disp) & 0xFF;
				++my_cbin[data_digit];
			}

			uint32_t bin_id;
			my_sbin[0] = my_cbin[0];
			for (bin_id = 1; bin_id < 0x100; ++bin_id)
			{
				my_cbin[bin_id] += my_cbin[bin_id - 1];
				my_sbin[bin_id] = my_cbin[bin_id];
			}
#pragma omp barrier

			uint32_t* other_cbin;
			uint32_t th_id;
			for (th_id = 0; th_id < my_th_id; ++th_id)
			{
				other_cbin = elem_count_bin + (size_t(th_id) << 8);
				for (bin_id = 0; bin_id < 0x100; ++bin_id)
					my_sbin[bin_id] += other_cbin[bin_id];
			}

			for (th_id = my_th_id + 1; th_id < thread_num; ++th_id)
			{
				other_cbin = elem_count_bin + (size_t(th_id) << 8);
				for (bin_id = 1; bin_id < 0x100; ++bin_id)
					my_sbin[bin_id] += other_cbin[bin_id - 1];
			}

			// reorder memory
			for (uint32_t p_id = p_id1; p_id-- > p_id0;)
			{
				data_digit = (pcl_in_elem_array[p_id] >> digit_disp) & 0xFF;
				pcl_in_elem_array_tmp[--my_cbin[data_digit]] = pcl_in_elem_array[p_id];
				new_to_ori_pcl_map_tmp[my_cbin[data_digit]] = new_to_ori_pcl_map[p_id];
			}

			swap(new_to_ori_pcl_map_ui, new_to_ori_pcl_map_tmp_ui);
			swap(pcl_in_elem_array_ui, pcl_in_elem_array_tmp_ui);
#pragma omp barrier
		}

		// divide pcl and element task
		p_id1 = Block_Low(my_th_id + 1, thread_num, pcl_num1);
		pcl_in_elem_id = pcl_in_elem_array[int64_t(p_id1) - 1];
		while (p_id1 < pcl_num1 && pcl_in_elem_id == pcl_in_elem_array[p_id1])
			++p_id1;
		pcl_range[my_th_id + 1].id = p_id1;
		elem_range[my_th_id + 1].id = pcl_in_elem_id + 1;
	}

	pcl_num = pcl_num1;
	//for (uint32_t th_id = 0; th_id < th_num_1; ++th_id)
	//	std::cout << pcl_range[th_id] << ", ";
	//std::cout << "\n";
	//for (uint32_t th_id = 0; th_id < th_num_1; ++th_id)
	//	std::cout << elem_range[th_id] << ", ";
	//std::cout << "\n";

	return 0;
}

int Step_T2D_ME_mt::finalize_calculation()
{
	Model_T2D_ME_mt &md = *(Model_T2D_ME_mt *)model;
	md.pcl_num = pcl_num;
	return 0;
}

int substep_func_omp_T2D_ME_mt(
	void* _self,
	uint32_t my_th_id,
	float dt,
	float cur_time,
	uint32_t substp_id
	)
{
	typedef Model_T2D_ME_mt::PclBodyForce PclBodyForce;
	typedef Model_T2D_ME_mt::PclTraction PclTraction;
	typedef Model_T2D_ME_mt::PclPos PclPos;
	typedef Model_T2D_ME_mt::PclSortedVarArray PclSortedVarArray;
	typedef Model_T2D_ME_mt::PclDisp PclDisp;
	typedef Model_T2D_ME_mt::PclV PclV;
	typedef Model_T2D_ME_mt::PclShapeFunc PclShapeFunc;
	typedef Model_T2D_ME_mt::PclStress PclStress;
	typedef Model_T2D_ME_mt::ElemNodeIndex ElemNodeIndex;
	typedef Model_T2D_ME_mt::ElemShapeFuncAB ElemShapeFuncAB;
	typedef Model_T2D_ME_mt::ElemShapeFuncC ElemShapeFuncC;
	typedef Model_T2D_ME_mt::ElemStrainInc ElemStrainInc;
	typedef Model_T2D_ME_mt::ElemStress ElemStress;
	typedef Model_T2D_ME_mt::ElemNodeVM ElemNodeVM;
	typedef Model_T2D_ME_mt::ElemNodeForce ElemNodeForce;
	typedef Model_T2D_ME_mt::NodeA NodeA;
	typedef Model_T2D_ME_mt::NodeV NodeV;
	typedef Model_T2D_ME_mt::NodeHasVBC NodeHasVBC;
	typedef Model_T2D_ME_mt::NodePos NodePos;

	Step_T2D_ME_mt &self = *(Step_T2D_ME_mt *)(_self);
	Model_T2D_ME_mt &md = *(Model_T2D_ME_mt*)(self.model);

	PclSortedVarArray& pscv0 = self.pcl_sorted_var_array[self.pcl_sorted_var_id];
	uint32_t* pcl_index0 = pscv0.pcl_index;
	float* pcl_density0 = pscv0.pcl_density;
	PclDisp* pcl_disp0 = pscv0.pcl_disp; 
	PclV* pcl_v0 = pscv0.pcl_v;
	PclShapeFunc* pcl_N0 = pscv0.pcl_N;
	PclStress* pcl_stress0 = pscv0.pcl_stress;
	uint32_t* elem_has_pcl_num0 = pscv0.elem_has_pcl_num;
	self.pcl_sorted_var_id ^= 1;
	PclSortedVarArray& pscv1 = self.pcl_sorted_var_array[self.pcl_sorted_var_id];
	uint32_t* pcl_index1 = pscv1.pcl_index;
	float* pcl_density1 = pscv1.pcl_density;
	PclDisp* pcl_disp1 = pscv1.pcl_disp;
	PclV* pcl_v1 = pscv1.pcl_v;
	PclShapeFunc* pcl_N1 = pscv1.pcl_N;
	PclStress* pcl_stress1 = pscv1.pcl_stress;
	uint32_t* elem_has_pcl_num1 = pscv1.elem_has_pcl_num;

	union
	{
		struct
		{
			uint32_t* new_to_ori_pcl_map;
			uint32_t* new_to_ori_pcl_map_tmp;
			uint32_t* pcl_in_elem_array;
			uint32_t* pcl_in_elem_array_tmp;
		};
		struct
		{
			size_t new_to_ori_pcl_map_ui;
			size_t new_to_ori_pcl_map_tmp_ui;
			size_t pcl_in_elem_array_ui;
			size_t pcl_in_elem_array_tmp_ui;
		};
	};
	new_to_ori_pcl_map = self.new_to_ori_pcl_map0;
	new_to_ori_pcl_map_tmp = self.new_to_ori_pcl_map1;
	pcl_in_elem_array = self.pcl_in_elem_array0;
	pcl_in_elem_array_tmp = self.pcl_in_elem_array1;

	// cal elem_has_pcl_num
	uint32_t e_id, e_id0, e_id1;
	e_id0 = self.elem_range[my_th_id].id;
	e_id1 = self.elem_range[my_th_id + 1].id;
	memset(elem_has_pcl_num1 + e_id0, 0, (e_id1 - e_id0) * sizeof(uint32_t));

	uint32_t p_id, p_id0, p_id1;
	p_id0 = self.pcl_range[my_th_id].id;
	p_id1 = self.pcl_range[my_th_id + 1].id;
	for (p_id = p_id0; p_id < p_id1; ++p_id)
		++elem_has_pcl_num1[pcl_in_elem_array[p_id]];

	// velocity mapping
	// nodal force integration
	uint32_t g_p_id = p_id0;
	for (e_id = e_id0; e_id < e_id1; ++e_id)
	{
		uint32_t ne_var_start_id = e_id * 3;
		ElemNodeVM& en_vm1 = self.elem_node_vm[ne_var_start_id];
		ElemNodeVM& en_vm2 = self.elem_node_vm[ne_var_start_id + 1];
		ElemNodeVM& en_vm3 = self.elem_node_vm[ne_var_start_id + 2];
		en_vm1.vm = 0.0f;
		en_vm1.vmx = 0.0f;
		en_vm1.vmy = 0.0f;
		en_vm2.vm = 0.0f;
		en_vm2.vmx = 0.0f;
		en_vm2.vmy = 0.0f;
		en_vm3.vm = 0.0f;
		en_vm3.vmx = 0.0f;
		en_vm3.vmy = 0.0f;
		ElemNodeForce& en_f1 = self.elem_node_force[ne_var_start_id];
		ElemNodeForce& en_f2 = self.elem_node_force[ne_var_start_id + 1];
		ElemNodeForce& en_f3 = self.elem_node_force[ne_var_start_id + 2];
		en_f1.fx = 0.0f;
		en_f1.fy = 0.0f;
		en_f2.fx = 0.0f;
		en_f2.fy = 0.0f;
		en_f3.fx = 0.0f;
		en_f3.fy = 0.0f;
		self.elem_am[e_id] = 0.0f;
		uint32_t p_num = elem_has_pcl_num0[e_id];
		if (p_num == 0)
			continue; // no particle
		float e_pcl_m = 0.0f;
		float e_pcl_bfx = 0.0f;
		float e_pcl_bfy = 0.0f;
		float e_pcl_vol = 0.0f;
		float e_s11 = 0.0f;
		float e_s22 = 0.0f;
		float e_s12 = 0.0f;
		float p_m, p_vol, p_N_m;
		for (p_id = 0; p_id < p_num; ++p_id, ++g_p_id)
		{
			uint32_t old_pcl_id = new_to_ori_pcl_map[g_p_id];
			
			uint32_t ori_pcl_id = pcl_index1[old_pcl_id];
			pcl_index0[g_p_id] = ori_pcl_id;

			p_m = self.pcl_m[ori_pcl_id];
			e_pcl_m += p_m;

			PclBodyForce& p_bf = self.pcl_bf[ori_pcl_id];
			e_pcl_bfx += p_bf.bfx;
			e_pcl_bfy += p_bf.bfy;

			p_vol = p_m / pcl_density1[old_pcl_id];
			e_pcl_vol += p_vol;
			PclStress& p_s = pcl_stress1[old_pcl_id];
			e_s11 += p_s.s11 * p_vol;
			e_s22 += p_s.s22 * p_vol;
			e_s12 += p_s.s12 * p_vol;
			PclShapeFunc& p_N = pcl_N1[old_pcl_id];
			PclShapeFunc& p_N0 = pcl_N0[g_p_id];
			p_N0.N1 = p_N.N1;
			p_N0.N2 = p_N.N2;
			p_N0.N3 = p_N.N3;
			PclV& p_v = pcl_v1[old_pcl_id];
			PclV& p_v0 = pcl_v0[g_p_id];
			p_v0.vx = p_v.vx;
			p_v0.vy = p_v.vy;
			p_N_m = p_N.N1 * p_m;
			en_vm1.vm += p_N_m;
			en_vm1.vmx += p_N_m * p_v.vx;
			en_vm1.vmy += p_N_m * p_v.vy;
			p_N_m = p_N.N2 * p_m;
			en_vm2.vm += p_N_m;
			en_vm2.vmx += p_N_m * p_v.vx;
			en_vm2.vmy += p_N_m * p_v.vy;
			p_N_m = p_N.N3 * p_m;
			en_vm3.vm += p_N_m;
			en_vm3.vmx += p_N_m * p_v.vx;
			en_vm3.vmy += p_N_m * p_v.vy;
			PclTraction& p_t = self.pcl_t[ori_pcl_id];
			en_f1.fx += p_N.N1 * p_t.tx;
			en_f1.fy += p_N.N1 * p_t.ty;
			en_f2.fx += p_N.N2 * p_t.tx;
			en_f2.fy += p_N.N2 * p_t.ty;
			en_f3.fx += p_N.N3 * p_t.tx;
			en_f3.fy += p_N.N3 * p_t.ty;
		}
		self.elem_density[e_id] = e_pcl_m / e_pcl_vol;
		e_s11 /= e_pcl_vol;
		e_s22 /= e_pcl_vol;
		e_s12 /= e_pcl_vol;
		if (e_pcl_vol > self.elem_area[e_id])
			e_pcl_vol = self.elem_area[e_id];
		e_pcl_m *= one_third;
		self.elem_am[e_id] = e_pcl_m;
		e_pcl_bfx *= one_third;
		e_pcl_bfy *= one_third;
		ElemShapeFuncAB& e_sf = self.elem_sf_ab[e_id];
		en_f1.fx += e_pcl_bfx;
		en_f1.fx -= (e_sf.dN1_dx * e_s11 + e_sf.dN1_dy * e_s12) * e_pcl_vol;
		en_f1.fy += e_pcl_bfy;
		en_f1.fy -= (e_sf.dN1_dx * e_s12 + e_sf.dN1_dy * e_s22) * e_pcl_vol;
		en_f2.fx += e_pcl_bfx;
		en_f2.fx -= (e_sf.dN2_dx * e_s11 + e_sf.dN2_dy * e_s12) * e_pcl_vol;
		en_f2.fy += e_pcl_bfy;
		en_f2.fy -= (e_sf.dN2_dx * e_s12 + e_sf.dN2_dy * e_s22) * e_pcl_vol;
		en_f3.fx += e_pcl_bfx;
		en_f3.fx -= (e_sf.dN3_dx * e_s11 + e_sf.dN3_dy * e_s12) * e_pcl_vol;
		en_f3.fy += e_pcl_bfy;
		en_f3.fy -= (e_sf.dN3_dx * e_s12 + e_sf.dN3_dy * e_s22) * e_pcl_vol;
	}
#pragma omp barrier

	// update nodal variables
	uint32_t n_id, n_id0, n_id1;
	uint32_t ne_id, ne_id1;
	n_id0 = self.node_range[my_th_id];
	n_id1 = self.node_range[my_th_id + 1];
	ne_id = self.node_elem_range[my_th_id];
	for (n_id = n_id0; n_id < n_id1; ++n_id)
	{
		float n_am = 0.0f;
		float n_fx = 0.0f;
		float n_fy = 0.0f;
		float n_vm = 0.0f;
		float n_vmx = 0.0f;
		float n_vmy = 0.0f;
		ne_id1 = self.node_elem_list[n_id];
		for (; ne_id < ne_id1; ++ne_id)
		{
			n_am += self.elem_am[self.elem_id_array[ne_id]];
			uint32_t node_var_id = self.node_elem_id_array[ne_id];
			ElemNodeForce& nf = self.elem_node_force[node_var_id];
			n_fx += nf.fx;
			n_fy += nf.fy;
			ElemNodeVM& nvm = self.elem_node_vm[node_var_id];
			n_vm += nvm.vm;
			n_vmx += nvm.vmx;
			n_vmy += nvm.vmy;
		}
		NodeA& node_a = self.node_a[n_id];
		if (n_am != 0.0f)
		{
			self.node_am[n_id] = n_am;
			node_a.ax = n_fx / n_am;
			node_a.ay = n_fy / n_am;
		}
		if (n_vm != 0.0f)
		{
			NodeV& node_v = self.node_v[n_id];
			node_v.vx = n_vmx / n_vm + node_a.ax * dt;
			node_v.vy = n_vmy / n_vm + node_a.ay * dt;
			NodeHasVBC& node_has_vbc = self.node_has_vbc[n_id];
			node_v.vx_ui &= (uint32_t(node_has_vbc.has_vx_bc) + UINT32_MAX - 1);
			node_v.vy_ui &= (uint32_t(node_has_vbc.has_vy_bc) + UINT32_MAX - 1);
		}
	}
#pragma omp barrier

	// cal element strain and strain enhancement
	for (e_id = e_id0; e_id < e_id1; ++e_id)
	{
		ElemNodeIndex& e_n_id = self.elem_node_id[e_id];
		ElemShapeFuncAB& e_sf = self.elem_sf_ab[e_id];
		ElemStrainInc& e_de = self.elem_de[e_id];
		NodeV& n_v1 = self.node_v[e_n_id.n1];
		NodeV& n_v2 = self.node_v[e_n_id.n2];
		NodeV& n_v3 = self.node_v[e_n_id.n3];
		e_de.de11 = (e_sf.dN1_dx * n_v1.vx + e_sf.dN2_dx * n_v2.vx + e_sf.dN3_dx * n_v3.vx) * dt;
		e_de.de22 = (e_sf.dN1_dy * n_v1.vy + e_sf.dN2_dy * n_v2.vy + e_sf.dN3_dy * n_v3.vy) * dt;
		e_de.de12 = (e_sf.dN1_dx * n_v1.vy + e_sf.dN2_dx * n_v2.vy + e_sf.dN3_dx * n_v3.vy
				   + e_sf.dN1_dy * n_v1.vx + e_sf.dN2_dy * n_v2.vx + e_sf.dN3_dy * n_v3.vx) * dt * 0.5f;
		float e_de_vol = e_de.de11 + e_de.de22;
		self.elem_am_de_vol[e_id] = self.elem_am[e_id] * e_de_vol;
		e_de_vol *= one_third;
		e_de.de11 -= e_de_vol;
		e_de.de22 -= e_de_vol;
	}
#pragma omp barrier

	ne_id = self.node_elem_range[my_th_id];
	for (n_id = n_id0; n_id < n_id1; ++n_id)
	{
		if (self.node_am[n_id] != 0.0f)
		{
			ne_id1 = self.node_elem_list[n_id];
			float n_am_de_vol = 0.0f;
			for (; ne_id < ne_id1; ++ne_id)
				n_am_de_vol += self.elem_am_de_vol[self.elem_id_array[ne_id]];
			self.node_de_vol[n_id] = n_am_de_vol / self.node_am[n_id];
		}
		else
		{
			ne_id = self.node_elem_list[n_id];
			self.node_de_vol[n_id] = 0.0f;
		}
	}

#pragma omp master
	self.new_pcl_num = 0;

#pragma omp barrier

	// update particle variables
	double dstrain[6] = { 0.0 };
	uint32_t pcl_in_elem_id;
	uint32_t pcl_in_mesh_num = 0;
	g_p_id = p_id0;
	for (e_id = e_id0; e_id < e_id1; ++e_id)
	{
		ElemNodeIndex& e_n_id = self.elem_node_id[e_id];
		float e_de_vol = (self.node_de_vol[e_n_id.n1]
						+ self.node_de_vol[e_n_id.n2]
						+ self.node_de_vol[e_n_id.n3]) * one_third;
		self.elem_density[e_id] /= (1.0f + e_de_vol);
		e_de_vol *= one_third;
		ElemStrainInc& e_de = self.elem_de[e_id];
		e_de.de11 += e_de_vol;
		e_de.de22 += e_de_vol;
		dstrain[0] = double(e_de.de11);
		dstrain[1] = double(e_de.de22);
		dstrain[3] = double(e_de.de12);
		uint32_t p_num = elem_has_pcl_num0[e_id];
		for (uint32_t p_id = 0; p_id < p_num; ++p_id, ++g_p_id)
		{
			uint32_t ori_pcl_id = pcl_index0[g_p_id];

			pcl_density0[g_p_id] = self.elem_density[e_id];

			MatModel::MaterialModel &pcl_mm = *self.pcl_mat_model[ori_pcl_id];
			int32_t mm_res = pcl_mm.integrate(dstrain);
			const double *dstress = pcl_mm.get_dstress();
			PclStress& p_s = pcl_stress0[g_p_id];
			p_s.s11 += float(dstress[0]);
			p_s.s22 += float(dstress[1]);
			p_s.s12 += float(dstress[3]);

			PclShapeFunc& p_N = pcl_N0[g_p_id];

			NodeA& n_a1 = self.node_a[e_n_id.n1];
			NodeA& n_a2 = self.node_a[e_n_id.n2];
			NodeA& n_a3 = self.node_a[e_n_id.n3];
			PclV& p_v = pcl_v0[g_p_id];
			p_v.vx += (p_N.N1 * n_a1.ax + p_N.N2 * n_a2.ax + p_N.N3 * n_a3.ax) * dt;
			p_v.vy += (p_N.N1 * n_a1.ay + p_N.N2 * n_a2.ay + p_N.N3 * n_a3.ay) * dt;

			NodeV& n_v1 = self.node_v[e_n_id.n1];
			NodeV& n_v2 = self.node_v[e_n_id.n2];
			NodeV& n_v3 = self.node_v[e_n_id.n3];
			PclDisp& p_d = pcl_disp0[g_p_id];
			p_d.ux += (p_N.N1 * n_v1.vx + p_N.N2 * n_v2.vx + p_N.N3 * n_v3.vx) * dt;
			p_d.uy += (p_N.N1 * n_v1.vy + p_N.N2 * n_v2.vy + p_N.N3 * n_v3.vy) * dt;

			PclPos& p_p = self.pcl_pos[ori_pcl_id];
			float pcl_x = p_p.x + p_d.ux;
			float pcl_y = p_p.y + p_d.uy;

			pcl_in_elem_id = e_id;
			if (!md.is_in_element(pcl_x, pcl_y, e_id, p_N))
				pcl_in_elem_id = md.find_pcl_in_which_elem(pcl_x, pcl_y, p_N);
			if (pcl_in_elem_id != self.elem_num)
			{
				if (p_N.N1 < N_min)
					p_N.N1 = N_min;
				if (p_N.N2 < N_min)
					p_N.N2 = N_min;
				if (p_N.N3 < N_min)
					p_N.N3 = N_min;
				++pcl_in_mesh_num;
			}
			new_to_ori_pcl_map[g_p_id] = g_p_id;
			pcl_in_elem_array[g_p_id] = pcl_in_elem_id;
		}
	}

#pragma omp critical
	self.new_pcl_num += pcl_in_mesh_num;

#pragma omp barrier

#pragma omp master
	self.pcl_num = self.new_pcl_num;

	// sort particle variables
	uint32_t* my_cbin = self.elem_count_bin + size_t(my_th_id) * 0x100;
	uint32_t* my_sbin = self.elem_sum_bin + size_t(my_th_id) * 0x100;
	for (uint32_t digit_disp = 0, elem_num_tmp = self.elem_num;
		 elem_num_tmp; digit_disp += 8, elem_num_tmp = elem_num_tmp >> 8)
	{
		memset(my_cbin, 0, 0x100 * sizeof(uint32_t));

		uint32_t data_digit;
		for (p_id = p_id0; p_id < p_id1; ++p_id)
		{
			data_digit = (pcl_in_elem_array[p_id] >> digit_disp) & 0xFF;
			++my_cbin[data_digit];
		}

		uint32_t bin_id;
		my_sbin[0] = my_cbin[0];
		for (bin_id = 1; bin_id < 0x100; ++bin_id)
		{
			my_cbin[bin_id] += my_cbin[bin_id - 1];
			my_sbin[bin_id] = my_cbin[bin_id];
		}
#pragma omp barrier

		uint32_t* other_cbin;
		uint32_t th_id;
		for (th_id = 0; th_id < my_th_id; ++th_id)
		{
			other_cbin = self.elem_count_bin + size_t(th_id) * 0x100;
			for (bin_id = 0; bin_id < 0x100; ++bin_id)
				my_sbin[bin_id] += other_cbin[bin_id];
		}

		uint32_t thread_num = self.thread_num;
		for (th_id = my_th_id + 1; th_id < thread_num; ++th_id)
		{
			other_cbin = self.elem_count_bin + size_t(th_id) * 0x100;
			for (bin_id = 1; bin_id < 0x100; ++bin_id)
				my_sbin[bin_id] += other_cbin[bin_id - 1];
		}

		// reorder memory
		for (p_id = p_id1; p_id-- > p_id0;)
		{
			data_digit = (pcl_in_elem_array[p_id] >> digit_disp) & 0xFF;
			pcl_in_elem_array_tmp[--my_cbin[data_digit]] = pcl_in_elem_array[p_id];
			new_to_ori_pcl_map_tmp[my_cbin[data_digit]] = new_to_ori_pcl_map[p_id];
		}

		swap(new_to_ori_pcl_map_ui, new_to_ori_pcl_map_tmp_ui);
		swap(pcl_in_elem_array_ui, pcl_in_elem_array_tmp_ui);
#pragma omp barrier
	}

	// divide pcl and element task
	p_id1 = Block_Low(my_th_id+1, self.thread_num, self.pcl_num);
	pcl_in_elem_id = pcl_in_elem_array[int64_t(p_id1) - 1];
	while (p_id1 < self.pcl_num && pcl_in_elem_id == pcl_in_elem_array[p_id1])
		++p_id1;
	self.pcl_range[my_th_id + 1].id = p_id1;
	self.elem_range[my_th_id + 1].id = pcl_in_elem_id + 1;
	//for (uint32_t th_id = 0; th_id < thread_num + 1; ++th_id)
	//	std::cout << pcl_range[th_id] << ", ";
	//std::cout << "\n";
	//for (uint32_t th_id = 0; th_id < thread_num + 1; ++th_id)
	//	std::cout << elem_range[th_id] << ", ";
	//std::cout << "\n";

	return 0;
}
