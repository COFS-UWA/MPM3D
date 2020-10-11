#include "SimulationsOMP_pcp.h"

#include <omp.h>

#include "Step_T2D_ME_mt.h"

Step_T2D_ME_mt::Step_T2D_ME_mt(const char* _name) : 
	Step(_name, "Step_T2D_ME_mt", &solve_substep_T2D_ME_mt),
	thread_num(1), pcl_range(nullptr),
	new_to_ori_pcl_map0(nullptr), elem_bin(nullptr) {}

Step_T2D_ME_mt::~Step_T2D_ME_mt() { clear_mem(); }

void Step_T2D_ME_mt::clear_mem()
{
	if (pcl_range)
	{
		delete[] pcl_range;
		pcl_range = nullptr;
	}
	if (new_to_ori_pcl_map0)
	{
		delete[] new_to_ori_pcl_map0;
		new_to_ori_pcl_map0 = nullptr;
	}
	if (elem_bin)
	{
		delete[] elem_bin;
		elem_bin = nullptr;
	}
}

#define one_third (1.0f/3.0f)
#define N_min (1.0e-10f)

namespace Step_T2D_ME_mt_Internal
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
	
	inline void divide_task_to_thread(
		uint32_t thread_num,
		uint32_t data_num,
		uint32_t data_range[] // len = thread_num+1
		)
	{
		register uint32_t i;
		uint32_t residue = data_num % thread_num;
		data_range[0] = 0;
		data_range[thread_num] = data_num;
		uint32_t th_min_resi = thread_num - residue;
		uint32_t data_inv = data_num / thread_num;
		for (i = 0; i < th_min_resi; ++i)
			data_range[i + 1] = data_range[i] + data_inv;
		--thread_num;
		++data_inv;
		for (i = th_min_resi; i < thread_num; ++i)
			data_range[i + 1] = data_range[i] + data_inv;
	}

	inline void swap(size_t &a, size_t &b)
	{
		a = a ^ b;
		b = a ^ b;
		a = a ^ b;
	}

	uint32_t pcl_num;
	uint32_t elem_num;
	uint32_t node_num;
	uint32_t vx_bc_num;
	uint32_t vy_bc_num;

	float* pcl_m;
	PclBodyForce* pcl_bf;
	PclTraction* pcl_t;
	PclPos* pcl_pos;
	MatModel::MaterialModel** pcl_mat_model;

	uint32_t pcl_sorted_var_id;
	PclSortedVarArray pcl_sorted_var_array[2];

	ElemNodeIndex* elem_node_id;
	float* elem_area;
	ElemShapeFuncAB* elem_sf_ab;
	ElemShapeFuncC* elem_sf_c;

	float* elem_density;
	ElemStrainInc* elem_de;
	ElemStress* elem_stress;
	float* elem_am;
	float* elem_am_de_vol;

	ElemNodeVM* elem_node_vm;
	ElemNodeForce* elem_node_force;

	uint32_t* elem_id_array;
	uint32_t* node_elem_id_array;
	uint32_t* node_elem_list;
	float* node_ax;
	float* node_ay;
	float* node_vx;
	float* node_vy;
	float* node_am;
	float* node_de_vol;

	uint32_t* vx_bcs;
	uint32_t* vy_bcs;

	uint32_t thread_num;

	uint32_t* pcl_range;
	uint32_t* elem_range;
	uint32_t* node_elem_range;
	uint32_t* node_range;
	uint32_t* vx_bc_range;
	uint32_t* vy_bc_range;

	uint32_t* new_to_ori_pcl_map0;
	uint32_t* new_to_ori_pcl_map1;
	uint32_t* pcl_in_elem_array0;
	uint32_t* pcl_in_elem_array1;

	uint32_t* elem_bin;
	uint32_t* elem_bin_range;
	uint32_t* elem_bin_offset;
}

uint32_t Step_T2D_ME_mt::get_pcl_num()
{
	return Step_T2D_ME_mt_Internal::pcl_num;
}

uint32_t Step_T2D_ME_mt::get_sorted_var_id()
{
	return Step_T2D_ME_mt_Internal::pcl_sorted_var_id;
}

int Step_T2D_ME_mt::init_calculation()
{
	using Step_T2D_ME_mt_Internal::divide_task_to_thread;
	using Step_T2D_ME_mt_Internal::swap;

	Model_T2D_ME_mt& md = *(Model_T2D_ME_mt*)model;

	omp_set_num_threads(thread_num);

	uint32_t th_num_1 = thread_num + 1;
	pcl_range = new uint32_t[size_t(th_num_1) * 8];
	elem_range = pcl_range + th_num_1;
	node_elem_range = elem_range + th_num_1;
	node_range = node_elem_range + th_num_1;
	vx_bc_range = node_range + th_num_1;
	vy_bc_range = vx_bc_range + th_num_1;
	elem_bin_range = vy_bc_range + th_num_1;
	elem_bin_offset = elem_bin_range + th_num_1;
	elem_range[0] = 0;
	divide_task_to_thread(thread_num, md.node_num, node_range);
	node_elem_range[0] = 0;
	for (uint32_t th_id = 0; th_id < thread_num; ++th_id)
		node_elem_range[th_id + 1] = md.node_elem_list[node_range[th_id]];
	divide_task_to_thread(thread_num, md.vx_bc_num, vx_bc_range);
	divide_task_to_thread(thread_num, md.vy_bc_num, vy_bc_range);
	divide_task_to_thread(thread_num, 0x100, elem_bin_range);
	elem_bin_offset[0] = 0;

	new_to_ori_pcl_map0 = new uint32_t[size_t(md.pcl_num) * 4];
	new_to_ori_pcl_map1 = new_to_ori_pcl_map0 + md.pcl_num;
	pcl_in_elem_array0 = new_to_ori_pcl_map1 + md.pcl_num;
	pcl_in_elem_array1 = pcl_in_elem_array0 + md.pcl_num;

	elem_bin = new uint32_t[size_t(thread_num) << 8];
	
	// accelerate substep function
	Step_T2D_ME_mt_Internal::pcl_num = md.pcl_num;
	Step_T2D_ME_mt_Internal::elem_num = md.elem_num;
	Step_T2D_ME_mt_Internal::node_num = md.node_num;
	Step_T2D_ME_mt_Internal::vx_bc_num = md.vx_bc_num;
	Step_T2D_ME_mt_Internal::vy_bc_num = md.vy_bc_num;

	Step_T2D_ME_mt_Internal::pcl_m = md.pcl_m;
	Step_T2D_ME_mt_Internal::pcl_bf = md.pcl_bf;
	Step_T2D_ME_mt_Internal::pcl_t = md.pcl_t;
	Step_T2D_ME_mt_Internal::pcl_pos = md.pcl_pos;
	Step_T2D_ME_mt_Internal::pcl_mat_model = md.pcl_mat_model;

	Step_T2D_ME_mt_Internal::pcl_sorted_var_id = 0;

	PclSortedVarArray& md_pscv0 = md.pcl_sorted_var_array[0];
	PclSortedVarArray& pscv0 = Step_T2D_ME_mt_Internal::pcl_sorted_var_array[0];
	pscv0.pcl_index = md_pscv0.pcl_index;
	pscv0.pcl_density = md_pscv0.pcl_density;
	pscv0.pcl_disp = md_pscv0.pcl_disp;
	pscv0.pcl_v = md_pscv0.pcl_v;
	pscv0.pcl_N = md_pscv0.pcl_N;
	pscv0.pcl_stress = md_pscv0.pcl_stress;
	pscv0.elem_has_pcl_num = md_pscv0.elem_has_pcl_num;

	PclSortedVarArray& md_pscv1 = md.pcl_sorted_var_array[1];
	PclSortedVarArray& pscv1 = Step_T2D_ME_mt_Internal::pcl_sorted_var_array[1];
	pscv1.pcl_index = md_pscv1.pcl_index;
	pscv1.pcl_density = md_pscv1.pcl_density;
	pscv1.pcl_disp = md_pscv1.pcl_disp;
	pscv1.pcl_v = md_pscv1.pcl_v;
	pscv1.pcl_N = md_pscv1.pcl_N;
	pscv1.pcl_stress = md_pscv1.pcl_stress;
	pscv1.elem_has_pcl_num = md_pscv1.elem_has_pcl_num;

	Step_T2D_ME_mt_Internal::elem_node_id = md.elem_node_id;
	Step_T2D_ME_mt_Internal::elem_area = md.elem_area;
	Step_T2D_ME_mt_Internal::elem_sf_ab = md.elem_sf_ab;
	Step_T2D_ME_mt_Internal::elem_sf_c = md.elem_sf_c;

	Step_T2D_ME_mt_Internal::elem_density = md.elem_density;
	Step_T2D_ME_mt_Internal::elem_de = md.elem_de;
	Step_T2D_ME_mt_Internal::elem_stress = md.elem_stress;
	Step_T2D_ME_mt_Internal::elem_am = md.elem_am;
	Step_T2D_ME_mt_Internal::elem_am_de_vol = md.elem_am_de_vol;

	Step_T2D_ME_mt_Internal::elem_node_vm = md.elem_node_vm;
	Step_T2D_ME_mt_Internal::elem_node_force = md.elem_node_force;

	Step_T2D_ME_mt_Internal::elem_id_array = md.elem_id_array;
	Step_T2D_ME_mt_Internal::node_elem_id_array = md.node_elem_id_array;
	Step_T2D_ME_mt_Internal::node_elem_list = md.node_elem_list;
	Step_T2D_ME_mt_Internal::node_ax = md.node_ax;
	Step_T2D_ME_mt_Internal::node_ay = md.node_ay;
	Step_T2D_ME_mt_Internal::node_vx = md.node_vx;
	Step_T2D_ME_mt_Internal::node_vy = md.node_vy;
	Step_T2D_ME_mt_Internal::node_am = md.node_am;
	Step_T2D_ME_mt_Internal::node_de_vol = md.node_de_vol;

	Step_T2D_ME_mt_Internal::vx_bcs = md.vx_bcs;
	Step_T2D_ME_mt_Internal::vy_bcs = md.vy_bcs;

	Step_T2D_ME_mt_Internal::thread_num = thread_num;

	Step_T2D_ME_mt_Internal::pcl_range = pcl_range;
	Step_T2D_ME_mt_Internal::elem_range = elem_range;
	Step_T2D_ME_mt_Internal::node_elem_range = node_elem_range;
	Step_T2D_ME_mt_Internal::node_range = node_range;
	Step_T2D_ME_mt_Internal::vx_bc_range = vx_bc_range;
	Step_T2D_ME_mt_Internal::vy_bc_range = vy_bc_range;

	Step_T2D_ME_mt_Internal::new_to_ori_pcl_map0 = new_to_ori_pcl_map0;
	Step_T2D_ME_mt_Internal::new_to_ori_pcl_map1 = new_to_ori_pcl_map1;
	Step_T2D_ME_mt_Internal::pcl_in_elem_array0 = pcl_in_elem_array0;
	Step_T2D_ME_mt_Internal::pcl_in_elem_array1 = pcl_in_elem_array1;

	Step_T2D_ME_mt_Internal::elem_bin = elem_bin;
	Step_T2D_ME_mt_Internal::elem_bin_range = elem_bin_range;
	Step_T2D_ME_mt_Internal::elem_bin_offset = elem_bin_offset;

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
	new_to_ori_pcl_map = Step_T2D_ME_mt_Internal::new_to_ori_pcl_map0;
	new_to_ori_pcl_map_tmp = Step_T2D_ME_mt_Internal::new_to_ori_pcl_map1;
	pcl_in_elem_array = Step_T2D_ME_mt_Internal::pcl_in_elem_array0;
	pcl_in_elem_array_tmp = Step_T2D_ME_mt_Internal::pcl_in_elem_array1;
	
	divide_task_to_thread(thread_num, md.pcl_num, pcl_range);
	PclSortedVarArray &psva = md.pcl_sorted_var_array[0];
	uint32_t pcl_num1 = 0;
#pragma omp parallel reduction(+:pcl_num1)
	{
		uint32_t my_th_id = omp_get_thread_num();
		uint32_t p_id0 = pcl_range[my_th_id];
		uint32_t p_id1 = pcl_range[my_th_id + 1];
		uint32_t pcl_in_mesh_num = 0;
		for (uint32_t p_id = p_id0; p_id < p_id1; ++p_id)
		{
			PclDisp& p_d = psva.pcl_disp[p_id];
			p_d.ux = 0.0f;
			p_d.uy = 0.0f;

			PclPos& p_p = md.pcl_pos[p_id];
			PclShapeFunc& p_N = psva.pcl_N[p_id];
			uint32_t pcl_in_elem_id = md.find_pcl_in_which_elem(p_p.x, p_p.y, p_N);
			if (pcl_in_elem_id != UINT32_MAX)
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
	}

	// sort particle variables
	for (uint32_t digit_disp = 0; digit_disp < 32; digit_disp += 8)
	{
#pragma omp parallel
		{
			// cal historgram
			uint32_t my_th_id = omp_get_thread_num();
			uint32_t *my_bin = elem_bin + (size_t(my_th_id) << 8);
			memset(my_bin, 0, 0x100 * sizeof(uint32_t));

			uint32_t p_id0 = pcl_range[my_th_id];
			uint32_t p_id1 = pcl_range[my_th_id + 1];
			uint32_t data_digit;
			for (uint32_t p_id = p_id0; p_id < p_id1; ++p_id)
			{
				data_digit = (pcl_in_elem_array[p_id] >> digit_disp) & 0xFF;
				++my_bin[data_digit];
			}
	#pragma omp barrier

			// cal global offset
			uint32_t bin_id0 = elem_bin_range[my_th_id];
			uint32_t bin_id1 = elem_bin_range[my_th_id + 1];
			uint32_t th_id, bin_id;
			if (bin_id0 != bin_id1)
			{
				for (th_id = 1; th_id < thread_num; ++th_id)
					elem_bin[(th_id << 8) + bin_id0] += elem_bin[((th_id - 1) << 8) + bin_id0];
				for (bin_id = bin_id0 + 1; bin_id < bin_id1; ++bin_id)
				{
					elem_bin[bin_id] += elem_bin[((thread_num - 1) << 8) + bin_id - 1];
					for (th_id = 1; th_id < thread_num; ++th_id)
						elem_bin[(th_id << 8) + bin_id] += elem_bin[((th_id - 1) << 8) + bin_id];
				}
			}
	#pragma omp barrier

			// init offset for all threads
	#pragma omp single
			for (th_id = 1; th_id < thread_num; ++th_id)
				elem_bin_offset[th_id] = elem_bin_offset[th_id - 1] + elem_bin[((thread_num - 1) << 8) + elem_bin_range[th_id] - 1];

			// add offset to each bin
			if (my_th_id != 0)
			{
				uint32_t my_bin_offset = elem_bin_offset[my_th_id];
				for (th_id = 0; th_id < thread_num; ++th_id)
					for (bin_id = bin_id0; bin_id < bin_id1; ++bin_id)
						elem_bin[(th_id << 8) + bin_id] += my_bin_offset;
			}
	#pragma omp barrier

			// reorder memory
			for (uint32_t p_id = p_id1; p_id-- > p_id0;)
			{
				data_digit = (pcl_in_elem_array[p_id] >> digit_disp) & 0xFF;
				pcl_in_elem_array_tmp[--my_bin[data_digit]] = pcl_in_elem_array[p_id];
				new_to_ori_pcl_map_tmp[my_bin[data_digit]] = new_to_ori_pcl_map[p_id];
			}
		}

		swap(new_to_ori_pcl_map_ui, new_to_ori_pcl_map_tmp_ui);
		swap(pcl_in_elem_array_ui, pcl_in_elem_array_tmp_ui);
	}

	// divde pcl and element task
	Step_T2D_ME_mt_Internal::pcl_num = pcl_num1;
	md.pcl_num = pcl_num1;
	divide_task_to_thread(thread_num, md.pcl_num, pcl_range);
	for (uint32_t th_id = 0; th_id < thread_num; ++th_id)
	{
		uint32_t p_id1 = pcl_range[th_id + 1];
		uint32_t pcl_in_elem_id = pcl_in_elem_array[p_id1 - 1];
		while (p_id1 < md.pcl_num && pcl_in_elem_id == pcl_in_elem_array[p_id1])
			++p_id1;
		pcl_range[th_id + 1] = p_id1;
		elem_range[th_id + 1] = pcl_in_elem_id + 1;
	}

#pragma omp parallel
	{
		uint32_t my_th_id = omp_get_thread_num();

		// cal elem_has_pcl_num
		uint32_t e_id0 = elem_range[my_th_id];
		uint32_t e_id1 = elem_range[my_th_id + 1];
		memset(psva.elem_has_pcl_num + e_id0, 0, (e_id1 - e_id0) * sizeof(uint32_t));
		uint32_t p_id0 = pcl_range[my_th_id];
		uint32_t p_id1 = pcl_range[my_th_id + 1];
		PclSortedVarArray& psva1 = md.pcl_sorted_var_array[1];
		for (uint32_t p_id = p_id0; p_id < p_id1; ++p_id)
		{
			++psva.elem_has_pcl_num[pcl_in_elem_array[p_id]];

			psva1.pcl_index[p_id] = psva.pcl_index[p_id];
			psva1.pcl_density[p_id] = psva.pcl_density[p_id];

			PclDisp& p_d1 = psva1.pcl_disp[p_id];
			PclDisp& p_d0 = psva.pcl_disp[p_id];
			p_d1.ux = p_d0.ux;
			p_d1.uy = p_d0.uy;

			PclV& p_v1 = psva1.pcl_v[p_id];
			PclV& p_v0 = psva.pcl_v[p_id];
			p_v1.vx = p_v0.vx;
			p_v1.vy = p_v0.vy;

			PclShapeFunc& p_sf1 = psva1.pcl_N[p_id];
			PclShapeFunc& p_sf0 = psva.pcl_N[p_id];
			p_sf1.N1 = p_sf0.N1;
			p_sf1.N2 = p_sf0.N2;
			p_sf1.N3 = p_sf0.N3;

			PclStress& p_s1 = psva1.pcl_stress[p_id];
			PclStress& p_s0 = psva.pcl_stress[p_id];
			p_s1.s11 = p_s0.s11;
			p_s1.s22 = p_s0.s22;
			p_s1.s12 = p_s0.s12;
		}
#pragma omp barrier

		// reorder particle data
		for (uint32_t p_id = p_id0; p_id < p_id1; ++p_id)
		{
			uint32_t ori_pcl_id = new_to_ori_pcl_map[p_id];

			psva.pcl_index[p_id] = psva1.pcl_index[ori_pcl_id];

			psva.pcl_density[p_id] = psva1.pcl_density[ori_pcl_id];

			PclDisp& p_d0 = psva.pcl_disp[ori_pcl_id];
			PclDisp& p_d1 = psva1.pcl_disp[p_id];
			p_d1.ux = p_d0.ux;
			p_d1.uy = p_d0.uy;

			PclV& p_v0 = psva.pcl_v[ori_pcl_id];
			PclV& p_v1 = psva1.pcl_v[p_id];
			p_v1.vx = p_v0.vx;
			p_v1.vy = p_v0.vy;

			PclShapeFunc& p_sf0 = psva.pcl_N[ori_pcl_id];
			PclShapeFunc& p_sf1 = psva1.pcl_N[p_id];
			p_sf1.N1 = p_sf0.N1;
			p_sf1.N2 = p_sf0.N2;
			p_sf1.N3 = p_sf0.N3;

			PclStress& p_s0 = psva.pcl_stress[ori_pcl_id];
			PclStress& p_s1 = psva1.pcl_stress[p_id];
			p_s1.s11 = p_s0.s11;
			p_s1.s22 = p_s0.s22;
			p_s1.s12 = p_s0.s12;
		}
	}

	return 0;
}

int Step_T2D_ME_mt::finalize_calculation()
{
	Model_T2D_ME_mt &md = *(Model_T2D_ME_mt *)model;
	md.pcl_num = Step_T2D_ME_mt_Internal::pcl_num;
	clear_mem();
	return 0;
}

int solve_substep_T2D_ME_mt(void* _self)
{
	using Step_T2D_ME_mt_Internal::PclBodyForce;
	using Step_T2D_ME_mt_Internal::PclBodyForce;
	using Step_T2D_ME_mt_Internal::PclTraction;
	using Step_T2D_ME_mt_Internal::PclPos;
	using Step_T2D_ME_mt_Internal::PclSortedVarArray;
	using Step_T2D_ME_mt_Internal::PclDisp;
	using Step_T2D_ME_mt_Internal::PclV;
	using Step_T2D_ME_mt_Internal::PclShapeFunc;
	using Step_T2D_ME_mt_Internal::PclStress;
	using Step_T2D_ME_mt_Internal::ElemNodeIndex;
	using Step_T2D_ME_mt_Internal::ElemShapeFuncAB;
	using Step_T2D_ME_mt_Internal::ElemShapeFuncC;
	using Step_T2D_ME_mt_Internal::ElemStrainInc;
	using Step_T2D_ME_mt_Internal::ElemStress;
	using Step_T2D_ME_mt_Internal::ElemNodeVM;
	using Step_T2D_ME_mt_Internal::ElemNodeForce;

	using Step_T2D_ME_mt_Internal::divide_task_to_thread;
	using Step_T2D_ME_mt_Internal::swap;

	using Step_T2D_ME_mt_Internal::pcl_num;
	using Step_T2D_ME_mt_Internal::elem_num;
	using Step_T2D_ME_mt_Internal::node_num;
	using Step_T2D_ME_mt_Internal::vx_bc_num;
	using Step_T2D_ME_mt_Internal::vy_bc_num;

	using Step_T2D_ME_mt_Internal::pcl_m;
	using Step_T2D_ME_mt_Internal::pcl_bf;
	using Step_T2D_ME_mt_Internal::pcl_t;
	using Step_T2D_ME_mt_Internal::pcl_pos;
	using Step_T2D_ME_mt_Internal::pcl_mat_model;

	using Step_T2D_ME_mt_Internal::pcl_sorted_var_id;
	using Step_T2D_ME_mt_Internal::pcl_sorted_var_array;

	using Step_T2D_ME_mt_Internal::elem_node_id;
	using Step_T2D_ME_mt_Internal::elem_area;
	using Step_T2D_ME_mt_Internal::elem_sf_ab;
	using Step_T2D_ME_mt_Internal::elem_sf_c;

	using Step_T2D_ME_mt_Internal::elem_density;
	using Step_T2D_ME_mt_Internal::elem_de;
	using Step_T2D_ME_mt_Internal::elem_stress;
	using Step_T2D_ME_mt_Internal::elem_am;
	using Step_T2D_ME_mt_Internal::elem_am_de_vol;

	using Step_T2D_ME_mt_Internal::elem_node_vm;
	using Step_T2D_ME_mt_Internal::elem_node_force;
	
	using Step_T2D_ME_mt_Internal::elem_id_array;
	using Step_T2D_ME_mt_Internal::node_elem_id_array;
	using Step_T2D_ME_mt_Internal::node_elem_list;
	using Step_T2D_ME_mt_Internal::node_ax;
	using Step_T2D_ME_mt_Internal::node_ay;
	using Step_T2D_ME_mt_Internal::node_vx;
	using Step_T2D_ME_mt_Internal::node_vy;
	using Step_T2D_ME_mt_Internal::node_am;
	using Step_T2D_ME_mt_Internal::node_de_vol;

	using Step_T2D_ME_mt_Internal::vx_bcs;
	using Step_T2D_ME_mt_Internal::vy_bcs;

	using Step_T2D_ME_mt_Internal::thread_num;

	using Step_T2D_ME_mt_Internal::pcl_range;
	using Step_T2D_ME_mt_Internal::elem_range;
	using Step_T2D_ME_mt_Internal::node_elem_range;
	using Step_T2D_ME_mt_Internal::node_range;
	using Step_T2D_ME_mt_Internal::vx_bc_range;
	using Step_T2D_ME_mt_Internal::vy_bc_range;

	using Step_T2D_ME_mt_Internal::elem_bin;
	using Step_T2D_ME_mt_Internal::elem_bin_range;
	using Step_T2D_ME_mt_Internal::elem_bin_offset;

	PclSortedVarArray& pscv0 = pcl_sorted_var_array[pcl_sorted_var_id];
	uint32_t* pcl_index0 = pscv0.pcl_index;
	float* pcl_density0 = pscv0.pcl_density;
	PclDisp* pcl_disp0 = pscv0.pcl_disp; 
	PclV* pcl_v0 = pscv0.pcl_v;
	PclShapeFunc* pcl_N0 = pscv0.pcl_N;
	PclStress* pcl_stress0 = pscv0.pcl_stress;
	uint32_t* elem_has_pcl_num0 = pscv0.elem_has_pcl_num;
	pcl_sorted_var_id ^= 1;
	PclSortedVarArray& pscv1 = pcl_sorted_var_array[pcl_sorted_var_id];
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
	new_to_ori_pcl_map = Step_T2D_ME_mt_Internal::new_to_ori_pcl_map0;
	new_to_ori_pcl_map_tmp = Step_T2D_ME_mt_Internal::new_to_ori_pcl_map1;
	pcl_in_elem_array = Step_T2D_ME_mt_Internal::pcl_in_elem_array0;
	pcl_in_elem_array_tmp = Step_T2D_ME_mt_Internal::pcl_in_elem_array1;
	
	Step_T2D_ME_mt &self = *(Step_T2D_ME_mt *)(_self);
	Model_T2D_ME_mt &md = *(Model_T2D_ME_mt *)(self.model);
	float dtime = float(self.dtime);

	uint32_t pcl_num1 = 0;
#pragma omp parallel reduction(+:pcl_num1)
	{
		uint32_t my_th_id = omp_get_thread_num();

		// velocity mapping
		// nodal force integration
		uint32_t g_p_id = pcl_range[my_th_id];
		uint32_t e_id0 = elem_range[my_th_id];
		uint32_t e_id1 = elem_range[my_th_id + 1];
		for (uint32_t e_id = e_id0; e_id < e_id1; ++e_id)
		{
			uint32_t ne_var_start_id = e_id * 3;
			ElemNodeVM& en_vm1 = elem_node_vm[ne_var_start_id];
			ElemNodeVM& en_vm2 = elem_node_vm[ne_var_start_id + 1];
			ElemNodeVM& en_vm3 = elem_node_vm[ne_var_start_id + 2];
			en_vm1.vm = 0.0f;
			en_vm1.vmx = 0.0f;
			en_vm1.vmy = 0.0f;
			en_vm2.vm = 0.0f;
			en_vm2.vmx = 0.0f;
			en_vm2.vmy = 0.0f;
			en_vm3.vm = 0.0f;
			en_vm3.vmx = 0.0f;
			en_vm3.vmy = 0.0f;
			ElemNodeForce& en_f1 = elem_node_force[ne_var_start_id];
			ElemNodeForce& en_f2 = elem_node_force[ne_var_start_id + 1];
			ElemNodeForce& en_f3 = elem_node_force[ne_var_start_id + 2];
			en_f1.fx = 0.0f;
			en_f1.fy = 0.0f;
			en_f2.fx = 0.0f;
			en_f2.fy = 0.0f;
			en_f3.fx = 0.0f;
			en_f3.fy = 0.0f;
			elem_am[e_id] = 0.0f;
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
			for (uint32_t p_id = 0; p_id < p_num; ++p_id, ++g_p_id)
			{
				uint32_t ori_pcl_id = pcl_index0[g_p_id];
				p_m = pcl_m[ori_pcl_id];
				e_pcl_m += p_m;
				PclBodyForce& p_bf = pcl_bf[ori_pcl_id];
				e_pcl_bfx += p_bf.bfx;
				e_pcl_bfy += p_bf.bfy;
				p_vol = p_m / pcl_density0[g_p_id];
				e_pcl_vol += p_vol;
				PclStress& p_s = pcl_stress0[g_p_id];
				e_s11 += p_s.s11 * p_vol;
				e_s22 += p_s.s22 * p_vol;
				e_s12 += p_s.s12 * p_vol;
				PclShapeFunc& p_N = pcl_N0[g_p_id];
				PclV& p_v = pcl_v0[g_p_id];
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
				PclTraction& p_t = pcl_t[ori_pcl_id];
				en_f1.fx += p_N.N1 * p_t.tx;
				en_f1.fy += p_N.N1 * p_t.ty;
				en_f2.fx += p_N.N2 * p_t.tx;
				en_f2.fy += p_N.N2 * p_t.ty;
				en_f3.fx += p_N.N3 * p_t.tx;
				en_f3.fy += p_N.N3 * p_t.ty;
			}
			elem_density[e_id] = e_pcl_m / e_pcl_vol;
			e_s11 /= e_pcl_vol;
			e_s22 /= e_pcl_vol;
			e_s12 /= e_pcl_vol;
			if (e_pcl_vol > elem_area[e_id])
				e_pcl_vol = elem_area[e_id];
			e_pcl_m *= one_third;
			e_pcl_bfx *= one_third;
			e_pcl_bfy *= one_third;
			elem_am[e_id] = e_pcl_m;
			ElemShapeFuncAB& e_sf = elem_sf_ab[e_id];
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
		uint32_t n_id0 = node_range[my_th_id];
		uint32_t n_id1 = node_range[my_th_id + 1];
		uint32_t n_id;
		uint32_t ne_id = node_elem_range[my_th_id];
		uint32_t ne_id1;
		for (uint32_t n_id = n_id0; n_id < n_id1; ++n_id)
		{
			float n_am = 0.0f;
			float n_fx = 0.0f;
			float n_fy = 0.0f;
			float n_vm = 0.0f;
			float n_vmx = 0.0f;
			float n_vmy = 0.0f;
			ne_id1 = node_elem_list[n_id];
			for (; ne_id < ne_id1; ++ne_id)
			{
				n_am += elem_am[elem_id_array[ne_id]];
				uint32_t node_var_id = node_elem_id_array[ne_id];
				ElemNodeForce& nf = elem_node_force[node_var_id];
				n_fx += nf.fx;
				n_fy += nf.fy;
				ElemNodeVM& nvm = elem_node_vm[node_var_id];
				n_vm += nvm.vm;
				n_vmx += nvm.vmx;
				n_vmy += nvm.vmy;
			}
			if (n_am != 0.0f && n_vm != 0.0f)
			{
				node_am[n_id] = n_am;
				node_ax[n_id] = n_fx / n_am;
				node_vx[n_id] = n_vmx / n_vm + node_ax[n_id] * dtime;
				node_ay[n_id] = n_fy / n_am;
				node_vy[n_id] = n_vmy / n_vm + node_ay[n_id] * dtime;
			}
		}
#pragma omp barrier

		// apply velocity bc
		uint32_t vbc_id, vbc_id0, vbc_id1;
		vbc_id0 = vx_bc_range[my_th_id];
		vbc_id1 = vx_bc_range[my_th_id + 1];
		for (vbc_id = vbc_id0; vbc_id < vbc_id1; ++vbc_id)
		{
			n_id = vx_bcs[vbc_id];
			node_ax[n_id] = 0.0f;
			node_vx[n_id] = 0.0f;
		}
		vbc_id0 = vy_bc_range[my_th_id];
		vbc_id1 = vy_bc_range[my_th_id + 1];
		for (vbc_id = vbc_id0; vbc_id < vbc_id1; ++vbc_id)
		{
			n_id = vy_bcs[vbc_id];
			node_ay[n_id] = 0.0f;
			node_vy[n_id] = 0.0f;
		}
#pragma omp barrier

		// cal element strain and strain enhancement
		e_id0 = elem_range[my_th_id];
		e_id1 = elem_range[my_th_id + 1];
		for (uint32_t e_id = e_id0; e_id < e_id1; ++e_id)
		{
			ElemNodeIndex& e_n_id = elem_node_id[e_id];
			ElemShapeFuncAB& e_sf = elem_sf_ab[e_id];
			ElemStrainInc& e_de = elem_de[e_id];
			e_de.de11 = (e_sf.dN1_dx * node_vx[e_n_id.n1] + e_sf.dN2_dx * node_vx[e_n_id.n2] + e_sf.dN3_dx * node_vx[e_n_id.n3]) * dtime;
			e_de.de22 = (e_sf.dN1_dy * node_vy[e_n_id.n1] + e_sf.dN2_dy * node_vy[e_n_id.n2] + e_sf.dN3_dy * node_vy[e_n_id.n3]) * dtime;
			e_de.de12 = (e_sf.dN1_dx * node_vy[e_n_id.n1] + e_sf.dN2_dx * node_vy[e_n_id.n2] + e_sf.dN3_dx * node_vy[e_n_id.n3]
					   + e_sf.dN1_dy * node_vx[e_n_id.n1] + e_sf.dN2_dy * node_vx[e_n_id.n2] + e_sf.dN3_dy * node_vx[e_n_id.n3]) * dtime * 0.5f;
			float e_de_vol = e_de.de11 + e_de.de22;
			elem_am_de_vol[e_id] = elem_am[e_id] * e_de_vol;
			e_de_vol *= one_third;
			e_de.de11 -= e_de_vol;
			e_de.de22 -= e_de_vol;
		}
#pragma omp barrier

		n_id0 = node_range[my_th_id];
		n_id1 = node_range[my_th_id + 1];
		ne_id = node_elem_range[my_th_id];
		for (n_id = n_id0; n_id < n_id1; ++n_id)
		{
			if (node_am[n_id] != 0.0f)
			{
				float n_am_de_vol = 0.0f;
				ne_id1 = node_elem_list[n_id];
				for (; ne_id < ne_id1; ++ne_id)
					n_am_de_vol += elem_am_de_vol[elem_id_array[ne_id]];
				node_de_vol[n_id] = n_am_de_vol / node_am[n_id];
			}
			else
			{
				ne_id = node_elem_list[n_id];
				node_de_vol[n_id] = 0.0f;
			}
		}
#pragma omp barrier

		// update particle variables
		double dstrain[6] = { 0.0 };
		uint32_t pcl_in_mesh_num = 0;
		g_p_id = pcl_range[my_th_id];
		e_id0 = elem_range[my_th_id];
		e_id1 = elem_range[my_th_id + 1];
		for (uint32_t e_id = e_id0; e_id < e_id1; ++e_id)
		{
			ElemNodeIndex& e_n_id = elem_node_id[e_id];
			float e_de_vol = (node_de_vol[e_n_id.n1]
							+ node_de_vol[e_n_id.n2]
							+ node_de_vol[e_n_id.n3]) * one_third;
			elem_density[e_id] /= (1.0f + e_de_vol);
			e_de_vol *= one_third;
			ElemStrainInc& e_de = elem_de[e_id];
			e_de.de11 += e_de_vol;
			e_de.de22 += e_de_vol;
			dstrain[0] = double(e_de.de11);
			dstrain[1] = double(e_de.de22);
			dstrain[3] = double(e_de.de12);
			uint32_t p_num = elem_has_pcl_num0[e_id];
			for (uint32_t p_id = 0; p_id < p_num; ++p_id, ++g_p_id)
			{
				uint32_t ori_pcl_id = pcl_index0[g_p_id];

				pcl_density0[g_p_id] = elem_density[e_id];

				MatModel::MaterialModel &pcl_mm = *pcl_mat_model[ori_pcl_id];
				int32_t mm_res = pcl_mm.integrate(dstrain);
				const double *dstress = pcl_mm.get_dstress();
				PclStress& p_s = pcl_stress0[g_p_id];
				p_s.s11 += float(dstress[0]);
				p_s.s22 += float(dstress[1]);
				p_s.s12 += float(dstress[3]);

				PclShapeFunc& p_N = pcl_N0[g_p_id];

				PclV& p_v = pcl_v0[g_p_id];
				p_v.vx += (p_N.N1 * node_ax[e_n_id.n1] + p_N.N2 * node_ax[e_n_id.n2] + p_N.N3 * node_ax[e_n_id.n3]) * dtime;
				p_v.vy += (p_N.N1 * node_ay[e_n_id.n1] + p_N.N2 * node_ay[e_n_id.n2] + p_N.N3 * node_ay[e_n_id.n3]) * dtime;

				PclDisp& p_d = pcl_disp0[g_p_id];
				p_d.ux += (p_N.N1 * node_vx[e_n_id.n1] + p_N.N2 * node_vx[e_n_id.n2] + p_N.N3 * node_vx[e_n_id.n3]) * dtime;
				p_d.uy += (p_N.N1 * node_vy[e_n_id.n1] + p_N.N2 * node_vy[e_n_id.n2] + p_N.N3 * node_vy[e_n_id.n3]) * dtime;

				PclPos& p_p = pcl_pos[ori_pcl_id];
				float pcl_x = p_p.x + p_d.ux;
				float pcl_y = p_p.y + p_d.uy;

				uint32_t pcl_in_elem_id = e_id;
				if (!md.is_in_element(pcl_x, pcl_y, e_id, p_N))
				{
					pcl_in_elem_id = md.find_pcl_in_which_elem(pcl_x, pcl_y, p_N);
					if (pcl_in_elem_id != UINT32_MAX)
					{
						if (p_N.N1 < N_min)
							p_N.N1 = N_min;
						if (p_N.N2 < N_min)
							p_N.N2 = N_min;
						if (p_N.N3 < N_min)
							p_N.N3 = N_min;
						++pcl_in_mesh_num;
					}
				}
				new_to_ori_pcl_map[g_p_id] = g_p_id;
				pcl_in_elem_array[g_p_id] = pcl_in_elem_id;
			}
		}
		pcl_num1 += pcl_in_mesh_num;
	}

	// sort particle variables
	for (uint32_t digit_disp = 0; digit_disp < 32; digit_disp += 8)
	{
#pragma omp parallel
		{
			// cal historgram
			uint32_t my_th_id = omp_get_thread_num();
			uint32_t *my_bin = elem_bin + (my_th_id << 8);
			memset(my_bin, 0, 0x100 * sizeof(uint32_t));

			uint32_t p_id0 = pcl_range[my_th_id];
			uint32_t p_id1 = pcl_range[my_th_id + 1];
			uint32_t p_id, data_digit;
			for (p_id = p_id0; p_id < p_id1; ++p_id)
			{
				data_digit = (pcl_in_elem_array[p_id] >> digit_disp) & 0xFF;
				++my_bin[data_digit];
			}
#pragma omp barrier

			// cal global offset
			uint32_t bin_id0 = elem_bin_range[my_th_id];
			uint32_t bin_id1 = elem_bin_range[my_th_id + 1];
			uint32_t th_id, bin_id;
			if (bin_id0 != bin_id1)
			{
				for (th_id = 1; th_id < thread_num; ++th_id)
					elem_bin[(th_id << 8) + bin_id0] += elem_bin[((th_id-1) << 8) + bin_id0];
				for (bin_id = bin_id0 + 1; bin_id < bin_id1; ++bin_id)
				{
					elem_bin[bin_id] += elem_bin[((thread_num-1) << 8) + bin_id-1];
					for (th_id = 1; th_id < thread_num; ++th_id)
						elem_bin[(th_id << 8) + bin_id] += elem_bin[((th_id-1) << 8) + bin_id];
				}
			}
#pragma omp barrier

			// init offset for all threads
#pragma omp single
			for (th_id = 1; th_id < thread_num; ++th_id)
				elem_bin_offset[th_id] = elem_bin_offset[th_id-1] + elem_bin[((thread_num-1) << 8) + elem_bin_range[th_id]-1];

			// add offset to each bin
			if (my_th_id != 0)
			{
				uint32_t my_bin_offset = elem_bin_offset[my_th_id];
				for (th_id = 0; th_id < thread_num; ++th_id)
					for (bin_id = bin_id0; bin_id < bin_id1; ++bin_id)
						elem_bin[(th_id << 8) + bin_id] += my_bin_offset;
			}
#pragma omp barrier

			// reorder memory
			for (p_id = p_id1; p_id-- > p_id0;)
			{
				data_digit = (pcl_in_elem_array[p_id] >> digit_disp) & 0xFF;
				pcl_in_elem_array_tmp[--my_bin[data_digit]] = pcl_in_elem_array[p_id];
				new_to_ori_pcl_map_tmp[my_bin[data_digit]] = new_to_ori_pcl_map[p_id];
			}
		}

		swap(new_to_ori_pcl_map_ui, new_to_ori_pcl_map_tmp_ui);
		swap(pcl_in_elem_array_ui, pcl_in_elem_array_tmp_ui);
	}

	// divde pcl and element task
	pcl_num = pcl_num1;
	divide_task_to_thread(thread_num, pcl_num, pcl_range);
	for (uint32_t th_id = 0; th_id < thread_num; ++th_id)
	{
		uint32_t p_id1 = pcl_range[th_id + 1];
		uint32_t pcl_in_elem_id = pcl_in_elem_array[p_id1 - 1];
		while (p_id1 < pcl_num && pcl_in_elem_id == pcl_in_elem_array[p_id1])
			++p_id1;
		pcl_range[th_id + 1] = p_id1;
		elem_range[th_id + 1] = pcl_in_elem_id + 1;
	}

#pragma omp parallel
	{
		uint32_t my_th_id = omp_get_thread_num();

		// cal elem_has_pcl_num
		uint32_t e_id0 = elem_range[my_th_id];
		uint32_t e_id1 = elem_range[my_th_id + 1];
		memset(elem_has_pcl_num1 + e_id0, 0, (e_id1 - e_id0) * sizeof(uint32_t));
		uint32_t p_id0 = pcl_range[my_th_id];
		uint32_t p_id1 = pcl_range[my_th_id + 1];
		uint32_t p_id;
		for (p_id = p_id0; p_id < p_id1; ++p_id)
			++elem_has_pcl_num1[pcl_in_elem_array[p_id]];

		// reorder particle data
		for (p_id = p_id0; p_id < p_id1; ++p_id)
		{
			uint32_t ori_pcl_id = new_to_ori_pcl_map[p_id];

			pcl_index1[p_id] = pcl_index0[ori_pcl_id];

			pcl_density1[p_id] = pcl_density0[ori_pcl_id];

			PclDisp& p_d0 = pcl_disp0[ori_pcl_id];
			PclDisp& p_d1 = pcl_disp1[p_id];
			p_d1.ux = p_d0.ux;
			p_d1.uy = p_d0.uy;

			PclV& p_v0 = pcl_v0[ori_pcl_id];
			PclV& p_v1 = pcl_v1[p_id];
			p_v1.vx = p_v0.vx;
			p_v1.vy = p_v0.vy;

			PclShapeFunc& p_sf0 = pcl_N0[ori_pcl_id];
			PclShapeFunc& p_sf1 = pcl_N1[p_id];
			p_sf1.N1 = p_sf0.N1;
			p_sf1.N2 = p_sf0.N2;
			p_sf1.N3 = p_sf0.N3;

			PclStress& p_s0 = pcl_stress0[ori_pcl_id];
			PclStress& p_s1 = pcl_stress1[p_id];
			p_s1.s11 = p_s0.s11;
			p_s1.s22 = p_s0.s22;
			p_s1.s12 = p_s0.s12;
		}
	}

	return 0;
}
