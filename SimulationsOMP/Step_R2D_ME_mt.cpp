#include "SimulationsOMP_pcp.h"

#include <fstream>
#include <omp.h>

#include "Step_R2D_ME_mt.h"

#define one_third (1.0/3.0)
#define one_fourth (0.25)
#define N_min (1.0e-8)
#define N_low(xi)  ((1.0 - (xi)) * 0.5)
#define N_high(xi) ((1.0 + (xi)) * 0.5)
#define dN_dxi_low -0.5
#define dN_dxi_high 0.5
#define Block_Low(th_id, th_num, data_num) ((th_id)*(data_num)/(th_num))

static std::fstream res_file_r2d_me_mt;

Step_R2D_ME_mt::Step_R2D_ME_mt(const char* _name) : 
	Step_OMP(_name, "Step_R2D_ME_mt", &substep_func_omp_R2D_ME_mt) {}

Step_R2D_ME_mt::~Step_R2D_ME_mt() {}

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
}

int Step_R2D_ME_mt::init_calculation()
{
	res_file_r2d_me_mt.open("r2d_me_mt_res.txt", std::ios::out | std::ios::binary);
	
	Model_R2D_ME_mt &md = *(Model_R2D_ME_mt *)model;

	omp_set_num_threads(thread_num);

	size_t th_num_1 = thread_num + 1;
	char *mem_range = (char *)task_range_mem.alloc((sizeof(PclRange) + sizeof(size_t) * 2) * th_num_1);
	pcl_range = (PclRange *)mem_range;
	mem_range += sizeof(PclRange) * th_num_1;
	elem_range = (size_t *)mem_range;
	mem_range += sizeof(size_t) * th_num_1;
	node_range = (size_t*)mem_range;

	pcl_range[0].id = 0;
	divide_task_to_thread(thread_num, md.elem_num, elem_range);
	divide_task_to_thread(thread_num, md.node_num, node_range);

	radix_sort_var_id = 0;
	new_to_prev_pcl_maps[0] = (size_t *)radix_sort_var_mem.alloc(sizeof(size_t) * (md.pcl_num * 4 + 2) + Cache_Alignment * 3);
	new_to_prev_pcl_maps[1] = cache_aligned(new_to_prev_pcl_maps[0] + md.pcl_num);
	pcl_in_elem_arrays[0] = cache_aligned(new_to_prev_pcl_maps[1] + md.pcl_num);
	pcl_in_elem_arrays[1] = cache_aligned(pcl_in_elem_arrays[0] + md.pcl_num + 1);
	pcl_in_elem_arrays[0][md.pcl_num] = md.elem_num;
	pcl_in_elem_arrays[1][md.pcl_num] = md.elem_num;

	elem_count_bin = (size_t *)elem_bin_mem.alloc(sizeof(size_t) * thread_num * 0x100 * 2);
	elem_sum_bin = elem_count_bin + thread_num * 0x100;

	pcl_num = md.pcl_num;
	elem_num = md.elem_num;
	node_num = md.node_num;

	pcl_m = md.pcl_m;
	pcl_bf = md.pcl_bf;
	pcl_t = md.pcl_t;
	pcl_pos = md.pcl_pos;
	pcl_vol = md.pcl_vol;
	pcl_mat_model = md.pcl_mat_model;

	pcl_sorted_var_id = 1;

	PclSortedVarArray &md_pscv0 = md.pcl_sorted_var_array[0];
	PclSortedVarArray &pscv0 = pcl_sorted_var_array[0];
	pscv0.pcl_index = md_pscv0.pcl_index;
	pscv0.pcl_density = md_pscv0.pcl_density;
	pscv0.pcl_disp = md_pscv0.pcl_disp;
	pscv0.pcl_v = md_pscv0.pcl_v;
	pscv0.pcl_N = md_pscv0.pcl_N;
	pscv0.pcl_stress = md_pscv0.pcl_stress;
	pscv0.pcl_strain = md_pscv0.pcl_strain;
	pscv0.pcl_estrain = md_pscv0.pcl_estrain;
	pscv0.pcl_pstrain = md_pscv0.pcl_pstrain;

	PclSortedVarArray& md_pscv1 = md.pcl_sorted_var_array[1];
	PclSortedVarArray& pscv1 = pcl_sorted_var_array[1];
	pscv1.pcl_index = md_pscv1.pcl_index;
	pscv1.pcl_density = md_pscv1.pcl_density;
	pscv1.pcl_disp = md_pscv1.pcl_disp;
	pscv1.pcl_v = md_pscv1.pcl_v;
	pscv1.pcl_stress = md_pscv1.pcl_stress;
	pscv1.pcl_strain = md_pscv1.pcl_strain;
	pscv1.pcl_estrain = md_pscv1.pcl_estrain;
	pscv1.pcl_pstrain = md_pscv1.pcl_pstrain;

	elem_sf = md.elem_sf;

	elem_density = md.elem_density;
	elem_pcl_m = md.elem_pcl_m;
	
	elem_node_vm = md.elem_node_vm;
	elem_node_force = md.elem_node_force;

	node_a = md.node_a;
	node_v = md.node_v;
	node_has_vbc = md.node_has_vbc;
	
	if (md.has_rigid_rect())
	{
		RigidRect &rr = md.get_rigid_rect();
		rr_fx_cont = 0.0;
		rr_fy_cont = 0.0;
		rr_m_cont = 0.0;
		rr.reset_f_contact();
	}

	for (size_t th_id = 1; th_id < th_num_1; ++th_id)
		pcl_range[th_id].id = Block_Low(th_id, thread_num, pcl_num);
	PclSortedVarArray &psva = md.pcl_sorted_var_array[0];
	PclSortedVarArray& psva1 = md.pcl_sorted_var_array[1];
	size_t next_pcl_num = 0;
#pragma omp parallel
	{
		size_t my_th_id = size_t(omp_get_thread_num());
		
		// init elem vars
		size_t e_id0 = elem_range[my_th_id];
		size_t e_id1 = elem_range[my_th_id + 1];
		memset(elem_pcl_m + e_id0, 0, (e_id1 - e_id0) * sizeof(double));
		memset(elem_pcl_vol + e_id0, 0, (e_id1 - e_id0) * sizeof(double));
		memset(elem_stress + e_id0, 0, (e_id1 - e_id0) * sizeof(ElemStress));
		memset(elem_node_vm + e_id0 * 3, 0, (e_id1 - e_id0) * 3 * sizeof(ElemNodeVM));
		memset(elem_node_force + e_id0 * 3, 0, (e_id1 - e_id0) * 3 * sizeof(ElemNodeForce));
		memset(elem_m_de_vol + e_id0 * 3, 0, (e_id1 - e_id0) * 3 * sizeof(double));

		size_t sort_var_id = radix_sort_var_id;
		size_t* new_to_prev_pcl_map = new_to_prev_pcl_maps[sort_var_id];
		size_t* new_to_prev_pcl_map_tmp = new_to_prev_pcl_maps[sort_var_id ^ 1];
		size_t* pcl_in_elem_array = pcl_in_elem_arrays[sort_var_id];
		size_t* pcl_in_elem_array_tmp = pcl_in_elem_arrays[sort_var_id ^ 1];

		size_t p_id0 = pcl_range[my_th_id].id;
		size_t p_id1 = pcl_range[my_th_id + 1].id;
		size_t pcl_in_mesh_num = 0;
		size_t e_id, p_id;
		for (p_id = p_id0; p_id < p_id1; ++p_id)
		{
			PclDisp& p_d = psva.pcl_disp[p_id];
			p_d.ux = 0.0;
			p_d.uy = 0.0;
			PclDisp& p_d1 = psva1.pcl_disp[p_id];
			p_d1.ux = 0.0;
			p_d1.uy = 0.0;

			PclPos& p_p = md.pcl_pos[p_id];
			PclShapeFunc& p_N = psva.pcl_N[p_id];
			e_id = md.find_pcl_in_which_elem(p_p.x, p_p.y, p_N);
			if (e_id != elem_num)
			{
				if (p_N.N1 < N_min)
					p_N.N1 = N_min;
				if (p_N.N2 < N_min)
					p_N.N2 = N_min;
				if (p_N.N3 < N_min)
					p_N.N3 = N_min;
				++pcl_in_mesh_num;
			}
			new_to_prev_pcl_map[p_id] = p_id;
			pcl_in_elem_array[p_id] = e_id;
		}

#pragma omp critical
		next_pcl_num += pcl_in_mesh_num;

		if (my_th_id == 0) // master thread
		{
			pcl_range[thread_num].id = next_pcl_num;
			radix_sort_var_id = sort_var_id;
			pcl_num = next_pcl_num;
		}
		else // non master thread
		{
			p_id = Block_Low(my_th_id, thread_num, next_pcl_num);
			if (p_id < next_pcl_num)
			{
				e_id = pcl_in_elem_array[p_id];
				while (e_id == pcl_in_elem_array[p_id + 1])
					++p_id;
				pcl_range[my_th_id].id = p_id + 1;
			}
			else
			{
				pcl_range[my_th_id].id = next_pcl_num;
			}
		}
	}

	K_cont = md.K_cont;
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
		contact_substep_id[p_id] = SIZE_MAX;
	
	return 0;
}

int Step_R2D_ME_mt::finalize_calculation()
{
	Model_R2D_ME_mt &md = *(Model_R2D_ME_mt *)model;
	md.pcl_num = pcl_num;
	return 0;
}

int substep_func_omp_R2D_ME_mt(
	void *_self,
	size_t my_th_id,
	double dt,
	double cur_time,
	size_t substp_id
	)
{
	typedef Model_R2D_ME_mt::ShapeFunc ShapeFunc;
	typedef Model_R2D_ME_mt::DShapeFunc DShapeFunc;
	typedef Model_R2D_ME_mt::Stress Stress;
	typedef Model_R2D_ME_mt::Strain Strain;
	typedef Model_R2D_ME_mt::StrainInc StrainInc;
	typedef Model_R2D_ME_mt::PclBodyForce PclBodyForce;
	typedef Model_R2D_ME_mt::PclTraction PclTraction;
	typedef Model_R2D_ME_mt::PclPos PclPos;
	typedef Model_R2D_ME_mt::PclSortedVarArray PclSortedVarArray;
	typedef Model_R2D_ME_mt::PclDisp PclDisp;
	typedef Model_R2D_ME_mt::PclV PclV;
	typedef Model_R2D_ME_mt::ElemNodeVM ElemNodeVM;
	typedef Model_R2D_ME_mt::ElemNodeForce ElemNodeForce;
	typedef Model_R2D_ME_mt::NodeA NodeA;
	typedef Model_R2D_ME_mt::NodeV NodeV;
	typedef Model_R2D_ME_mt::NodeHasVBC NodeHasVBC;
	typedef Model_R2D_ME_mt::ContPos ContPos;

	Step_R2D_ME_mt &self = *(Step_R2D_ME_mt*)(_self);
	Model_R2D_ME_mt &md = *(Model_R2D_ME_mt*)(self.model);

	PclSortedVarArray& pscv0 = self.pcl_sorted_var_array[self.pcl_sorted_var_id];
	size_t* pcl_index0 = pscv0.pcl_index;
	double *pcl_density0 = pscv0.pcl_density;
	PclDisp *pcl_disp0 = pscv0.pcl_disp;
	PclV *pcl_v0 = pscv0.pcl_v;
	Stress* pcl_stress0 = pscv0.pcl_stress;
	Strain *pcl_strain0 = pscv0.pcl_strain;
	Strain* pcl_estrain0 = pscv0.pcl_estrain;
	Strain* pcl_pstrain0 = pscv0.pcl_pstrain;

	PclSortedVarArray& pscv1 = self.pcl_sorted_var_array[self.pcl_sorted_var_id ^ 1];
	size_t* pcl_index1 = pscv1.pcl_index;
	double* pcl_density1 = pscv1.pcl_density;
	PclDisp* pcl_disp1 = pscv1.pcl_disp;
	PclV *pcl_v1 = pscv1.pcl_v;
	Stress *pcl_stress1 = pscv1.pcl_stress;
	Strain* pcl_strain1 = pscv1.pcl_strain;
	Strain* pcl_estrain1 = pscv1.pcl_estrain;
	Strain* pcl_pstrain1 = pscv1.pcl_pstrain;

	size_t sort_var_id = self.radix_sort_var_id;
	size_t* pcl_in_elem_array = self.pcl_in_elem_arrays[sort_var_id];
	size_t* pcl_in_elem_array_tmp = self.pcl_in_elem_arrays[sort_var_id ^ 1];

	// get p_id0, p_id1...
	size_t p_id0, p_id1, p_id;

	size_t ae_x_id, ae_y_id;
	size_t ori_pcl_id;
	double p_m, p_vol, p_N_m;
	double one_fourth_bfx, one_fourth_bfy;
	double xi, eta;
	double e_pcl_m = 0.0;
	double e_pcl_vol = 0.0;
	double e_s11 = 0.0;
	double e_s22 = 0.0;
	double e_s12 = 0.0;
	double en1_vm = 0.0;
	double en1_vmx = 0.0;
	double en1_vmy = 0.0;
	double en2_vm = 0.0;
	double en2_vmx = 0.0;
	double en2_vmy = 0.0;
	double en3_vm = 0.0;
	double en3_vmx = 0.0;
	double en3_vmy = 0.0;
	double en4_vm = 0.0;
	double en4_vmx = 0.0;
	double en4_vmy = 0.0;
	double en1_fx = 0.0;
	double en1_fy = 0.0;
	double en2_fx = 0.0;
	double en2_fy = 0.0;
	double en3_fx = 0.0;
	double en3_fy = 0.0;
	double en4_fx = 0.0;
	double en4_fy = 0.0;
	DShapeFunc& e_dN = self.elem_dN;
	size_t ae_id = pcl_in_elem_array[p_id0];
	for (p_id = p_id0; p_id < p_id1; ++p_id)
	{
		// elem has pcl
		self.elem_substp_id[ae_id] = substp_id;

		// map pcl mass
		ori_pcl_id = pcl_index0[p_id];
		p_m = self.pcl_m[ori_pcl_id];
		e_pcl_m += p_m;

		// map pcl volume
		p_vol = p_m / pcl_density0[p_id];
		self.pcl_vol[p_id] = p_vol;
		e_pcl_vol += p_vol;

		// map stress
		Stress &p_s = pcl_stress0[p_id];
		e_s11 += p_s.s11 * p_vol;
		e_s22 += p_s.s22 * p_vol;
		e_s12 += p_s.s12 * p_vol;

		// cal shape function
		ae_x_id = ae_id % self.actual_elem_x_num;
		ae_y_id = ae_id / self.actual_elem_x_num;
		PclPos &p_p = self.pcl_pos[ori_pcl_id];
		PclDisp &p_d = pcl_disp0[p_id];
		xi = (p_p.x + p_d.ux - self.mh_xl) * self.inv_elem_hx - double(ae_x_id);
		eta = (p_p.y + p_d.uy - self.mh_yl) * self.inv_elem_hy - double(ae_y_id);
		ShapeFunc& p_N = self.pcl_N[p_id];
		p_N.N1 = N_low(xi) * N_low(eta);
		if (p_N.N1 < N_min)
			p_N.N1 = N_min;
		p_N.N2 = N_high(xi) * N_low(eta);
		if (p_N.N2 < N_min)
			p_N.N2 = N_min;
		p_N.N3 = N_low(xi) * N_high(eta);
		if (p_N.N3 < N_min)
			p_N.N3 = N_min;
		p_N.N4 = N_high(xi) * N_high(eta);
		if (p_N.N4 < N_min)
			p_N.N4 = N_min;
		DShapeFunc& p_dN = self.pcl_dN[p_id];
		p_dN.dN1_dx = dN_dxi_low * N_low(eta) * self.dxi_dx;
		p_dN.dN1_dy = N_low(xi) * dN_dxi_low * self.deta_dy;
		p_dN.dN2_dx = dN_dxi_high * N_low(eta) * self.dxi_dx;
		p_dN.dN2_dy = N_high(xi) * dN_dxi_low * self.deta_dy;
		p_dN.dN3_dx = dN_dxi_low * N_high(eta) * self.dxi_dx;
		p_dN.dN3_dy = N_low(xi) * dN_dxi_high * self.deta_dy;
		p_dN.dN4_dx = dN_dxi_high * N_high(eta) * self.dxi_dx;
		p_dN.dN4_dy = N_high(xi) * dN_dxi_high * self.deta_dy;

		// map velocity
		PclV& p_v = pcl_v0[p_id];
		p_N_m = p_N.N1 * p_m;
		en1_vm += p_N_m;
		en1_vmx += p_N_m * p_v.vx;
		en1_vmy += p_N_m * p_v.vy;
		p_N_m = p_N.N2 * p_m;
		en2_vm += p_N_m;
		en2_vmx += p_N_m * p_v.vx;
		en2_vmy += p_N_m * p_v.vy;
		p_N_m = p_N.N3 * p_m;
		en3_vm += p_N_m;
		en3_vmx += p_N_m * p_v.vx;
		en3_vmy += p_N_m * p_v.vy;
		p_N_m = p_N.N4 * p_m;
		en4_vm += p_N_m;
		en4_vmx += p_N_m * p_v.vx;
		en4_vmy += p_N_m * p_v.vy;

		// external load
		PclBodyForce& p_bf = self.pcl_bf[ori_pcl_id];
		one_fourth_bfx = one_fourth * p_bf.bfx;
		one_fourth_bfy = one_fourth * p_bf.bfy;
		PclTraction& p_t = self.pcl_t[ori_pcl_id];
		en1_fx += one_fourth_bfx + p_N.N1 * p_t.tx;
		en1_fy += one_fourth_bfy + p_N.N1 * p_t.ty;
		en2_fx += one_fourth_bfx + p_N.N2 * p_t.tx;
		en2_fy += one_fourth_bfy + p_N.N2 * p_t.ty;
		en3_fx += one_fourth_bfx + p_N.N3 * p_t.tx;
		en3_fy += one_fourth_bfy + p_N.N3 * p_t.ty;
		en4_fx += one_fourth_bfx + p_N.N4 * p_t.tx;
		en4_fy += one_fourth_bfy + p_N.N4 * p_t.ty;
	
		if (ae_id != pcl_in_elem_array[p_id + 1])
		{
			ae_id = pcl_in_elem_array[p_id + 1];

			self.elem_pcl_m[ae_id] = e_pcl_m;
			self.elem_density[ae_id] = e_pcl_m / e_pcl_vol;

			ElemNodeVM& en1_v = self.elem_node_vm[ae_id * 4];
			en1_v.vm = en1_vm;
			en1_v.vmx = en1_vmx;
			en1_v.vmy = en1_vmy;
			ElemNodeVM& en2_v = self.elem_node_vm[ae_id * 4 + 1];
			en2_v.vm = en2_vm;
			en2_v.vmx = en2_vmx;
			en2_v.vmy = en2_vmy;
			ElemNodeVM& en3_v = self.elem_node_vm[ae_id * 4 + 2];
			en3_v.vm = en3_vm;
			en3_v.vmx = en3_vmx;
			en3_v.vmy = en3_vmy;
			ElemNodeVM& en4_v = self.elem_node_vm[ae_id * 4 + 3];
			en4_v.vm = en4_vm;
			en4_v.vmx = en4_vmx;
			en4_v.vmy = en4_vmy;

			e_s11 /= e_pcl_vol;
			e_s22 /= e_pcl_vol;
			e_s12 /= e_pcl_vol;
			if (e_pcl_vol > self.elem_area)
				e_pcl_vol = self.elem_area;
			// node 1
			en1_fx -= (e_dN.dN1_dx * e_s11 + e_dN.dN1_dy * e_s12) * e_pcl_vol;
			en1_fy -= (e_dN.dN1_dx * e_s12 + e_dN.dN1_dy * e_s22) * e_pcl_vol;
			ElemNodeForce &en1_f = self.elem_node_force[ae_id * 4];
			en1_f.fx = en1_fx;
			en1_f.fy = en1_fy;
			// node 2
			en2_fx -= (e_dN.dN2_dx * e_s11 + e_dN.dN2_dy * e_s12) * e_pcl_vol;
			en2_fy -= (e_dN.dN2_dx * e_s12 + e_dN.dN2_dy * e_s22) * e_pcl_vol;
			ElemNodeForce& en2_f = self.elem_node_force[ae_id * 4 + 1];
			en2_f.fx = en2_fx;
			en2_f.fy = en2_fy;
			// node 3
			en3_fx -= (e_dN.dN3_dx * e_s11 + e_dN.dN3_dy * e_s12) * e_pcl_vol;
			en3_fy -= (e_dN.dN3_dx * e_s12 + e_dN.dN3_dy * e_s22) * e_pcl_vol;
			ElemNodeForce& en3_f = self.elem_node_force[ae_id * 4 + 2];
			en3_f.fx = en3_fx;
			en3_f.fy = en3_fy;
			// node 4
			en4_fx -= (e_dN.dN4_dx * e_s11 + e_dN.dN4_dy * e_s12) * e_pcl_vol;
			en4_fy -= (e_dN.dN4_dx * e_s12 + e_dN.dN4_dy * e_s22) * e_pcl_vol;
			ElemNodeForce& en4_f = self.elem_node_force[ae_id * 4 + 3];
			en4_f.fx = en4_fx;
			en4_f.fy = en4_fy;

			// reset
			e_pcl_m = 0.0;
			e_pcl_vol = 0.0;
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
			en4_vm = 0.0;
			en4_vmx = 0.0;
			en4_vmy = 0.0;
			en1_fx = 0.0;
			en1_fy = 0.0;
			en2_fx = 0.0;
			en2_fy = 0.0;
			en3_fx = 0.0;
			en3_fy = 0.0;
			en4_fx = 0.0;
			en4_fy = 0.0;
		}
	}

	if (md.has_rigid_rect())
	{
		RigidRectForce rr_force;
		rr_force.reset_f_contact();
		self.apply_rigid_rect_avg(
			my_th_id, dt,
			pcl_in_elem_array, pscv0,
			rr_force
			);
		RigidRect& rr = md.get_rigid_rect();
#pragma omp critical
		{
			rr.combine(rr_force);
		}
	}
#pragma omp barrier

	// update node variables
	size_t n_id0 = self.node_range[my_th_id];
	size_t n_id1 = self.node_range[my_th_id + 1];
	size_t n_id, n_x_id, n_y_id;
	double n_vm, n_vmx, n_vmy;
	double n_am, n_fx, n_fy;
	size_t bc_mask;
	for (n_id = n_id0; n_id < n_id1; ++n_id)
	{
		n_x_id = n_id % self.node_x_num;
		n_y_id = n_id / self.node_x_num;
		n_vm = 0.0;
		n_vmx = 0.0;
		n_vmy = 0.0;
		n_am = 0.0;
		n_fx = 0.0;
		n_fy = 0.0;
		ae_id = self.actual_elem_x_num * n_y_id + n_x_id;
		if (self.elem_substp_id[ae_id] == substp_id)
		{
			ElemNodeVM &en_vm = self.elem_node_vm[ae_id * 4 + 3];
			n_vm += en_vm.vm;
			n_vmx += en_vm.vmx;
			n_vmy += en_vm.vmy;
			n_am += self.elem_pcl_m[ae_id];
			ElemNodeForce& en_f = self.elem_node_force[ae_id * 4 + 3];
			n_fx += en_f.fx;
			n_fy += en_f.fy;
		}
		++ae_id;
		if (self.elem_substp_id[ae_id] == substp_id)
		{
			ElemNodeVM& en_vm = self.elem_node_vm[ae_id * 4 + 2];
			n_vm += en_vm.vm;
			n_vmx += en_vm.vmx;
			n_vmy += en_vm.vmy;
			n_am += self.elem_pcl_m[ae_id];
			ElemNodeForce& en_f = self.elem_node_force[ae_id * 4 + 2];
			n_fx += en_f.fx;
			n_fy += en_f.fy;
		}
		ae_id += self.actual_elem_x_num;
		if (self.elem_substp_id[ae_id] == substp_id)
		{
			ElemNodeVM& en_vm = self.elem_node_vm[ae_id * 4];
			n_vm += en_vm.vm;
			n_vmx += en_vm.vmx;
			n_vmy += en_vm.vmy;
			n_am += self.elem_pcl_m[ae_id];
			ElemNodeForce& en_f = self.elem_node_force[ae_id * 4];
			n_fx += en_f.fx;
			n_fy += en_f.fy;
		}
		--ae_id;
		if (self.elem_substp_id[ae_id] == substp_id)
		{
			ElemNodeVM& en_vm = self.elem_node_vm[ae_id * 4 + 1];
			n_vm += en_vm.vm;
			n_vmx += en_vm.vmx;
			n_vmy += en_vm.vmy;
			n_am += self.elem_pcl_m[ae_id];
			ElemNodeForce& en_f = self.elem_node_force[ae_id * 4 + 1];
			n_fx += en_f.fx;
			n_fy += en_f.fy;
		}
		NodeA &node_a = self.node_a[n_id];
		if (n_am != 0.0)
		{
			n_am *= one_fourth;
			node_a.ax = n_fx / n_am;
			node_a.ay = n_fy / n_am;
		}
		if (n_vm != 0.0)
		{
			NodeV& node_v = self.node_v[n_id];
			node_v.vx = n_vmx / n_vm + node_a.ax * dt;
			node_v.vy = n_vmy / n_vm + node_a.ay * dt;
			NodeHasVBC& node_has_vbc = self.node_has_vbc[n_id];
			bc_mask = size_t(node_has_vbc.has_vx_bc) + SIZE_MAX;
			node_a.ax_ui &= bc_mask;
			node_v.vx_ui &= bc_mask;
			bc_mask = size_t(node_has_vbc.has_vy_bc) + SIZE_MAX;
			node_a.ay_ui &= bc_mask;
			node_v.vy_ui &= bc_mask;
		}
	}

#pragma omp master
	{
		if (md.has_rigid_rect())
		{
			RigidRect& rr = md.get_rigid_rect();
			rr.update_motion(dt);
			self.rr_fx_cont = rr.get_fx_contact();
			self.rr_fy_cont = rr.get_fy_contact();
			self.rr_m_cont = rr.get_m_contact();
			rr.reset_f_contact();
		}

		self.pcl_num = 0;
	}
#pragma omp barrier

	// update particle variables
	p_id0 = self.pcl_range[my_th_id].id;
	p_id1 = self.pcl_range[my_th_id + 1].id;
	size_t n1_id, n2_id, n3_id, n4_id;
	NodeV *pn1_v, *pn2_v, *pn3_v, *pn4_v;
	size_t p_ae_x_id, p_ae_y_id;
	double e_de_vol_div_3, p_de_vol_div_3;
	double pcl_x, pcl_y;
	union
	{
		double pcl_dstrain[6];
		struct { double p_de11, p_de22, p_de33, p_de12, p_de23, p_de31; };
	};
	p_de33 = 0.0;
	p_de23 = 0.0;
	p_de31 = 0.0;
	size_t last_ae_id = SIZE_MAX;
	size_t pcl_in_mesh_num = 0;
	for (p_id = p_id0; p_id < p_id1; ++p_id)
	{
		ae_id = pcl_in_elem_array[p_id];
		if (last_ae_id != ae_id)
		{
			ae_x_id = ae_id % self.actual_elem_x_num;
			ae_y_id = ae_id / self.actual_elem_x_num;
			n4_id = self.node_x_num * ae_x_id + ae_y_id;
			pn4_v = self.node_v + n4_id;
			n3_id = n4_id - 1;
			pn3_v = self.node_v + n3_id;
			n2_id = n4_id - self.node_x_num;
			pn2_v = self.node_v + n2_id;
			n1_id = n2_id - 1;
			pn1_v = self.node_v + n1_id;
			e_de_vol_div_3 = (e_dN.dN1_dx * pn1_v->vx + e_dN.dN2_dx * pn2_v->vx
							+ e_dN.dN3_dx * pn3_v->vx + e_dN.dN4_dx * pn4_v->vx
							+ e_dN.dN1_dy * pn1_v->vy + e_dN.dN2_dy * pn2_v->vy
							+ e_dN.dN3_dy * pn3_v->vy + e_dN.dN4_dy * pn4_v->vy) * dt;
			self.elem_density[ae_id] /= 1.0 + e_de_vol_div_3;
			e_de_vol_div_3 *= one_third;
			last_ae_id = ae_id;
		}
		
		// update strain
		DShapeFunc& p_dN = self.pcl_dN[p_id];
		p_de11 = (p_dN.dN1_dx * pn1_v->vx + p_dN.dN2_dx * pn2_v->vx
				+ p_dN.dN3_dx * pn3_v->vx + p_dN.dN4_dx * pn4_v->vx) * dt;
		p_de22 = (p_dN.dN1_dy * pn1_v->vy + p_dN.dN2_dy * pn2_v->vy
				+ p_dN.dN3_dy * pn3_v->vy + p_dN.dN4_dy * pn4_v->vy) * dt;
		p_de12 = (p_dN.dN1_dx * pn1_v->vy + p_dN.dN2_dx * pn2_v->vy
				+ p_dN.dN3_dx * pn3_v->vy + p_dN.dN4_dx * pn4_v->vy
				+ p_dN.dN1_dy * pn1_v->vx + p_dN.dN2_dy * pn2_v->vx
				+ p_dN.dN3_dy * pn3_v->vx + p_dN.dN4_dy * pn4_v->vx) * dt;
		p_de_vol_div_3 = (p_de11 + p_de22) * one_third;
		p_de11 += -p_de_vol_div_3 + e_de_vol_div_3;
		p_de22 += -p_de_vol_div_3 + e_de_vol_div_3;
		Strain &p_e = pcl_strain0[p_id];
		p_e.e11 += p_de11;
		p_e.e22 += p_de22;
		p_e.e12 += p_de12;

		// update stress
		ori_pcl_id = pcl_index0[p_id];
		MatModel::MaterialModel& pcl_mm = *self.pcl_mat_model[ori_pcl_id];
		int mm_res = pcl_mm.integrate(pcl_dstrain);
		const double* dstress = pcl_mm.get_dstress();
		Stress &p_s = pcl_stress0[p_id];
		p_s.s11 += dstress[0];
		p_s.s22 += dstress[1];
		p_s.s12 += dstress[3];
		const double* dee = pcl_mm.get_dstrain_e();
		Strain &pcl_ee = pcl_estrain0[p_id];
		pcl_ee.e11 += dee[0];
		pcl_ee.e22 += dee[1];
		pcl_ee.e12 += dee[3];
		const double* dpe = pcl_mm.get_dstrain_p();
		Strain &pcl_pe = pcl_pstrain0[p_id];
		pcl_pe.e11 += dpe[0];
		pcl_pe.e22 += dpe[1];
		pcl_pe.e12 += dpe[3];

		// update density
		pcl_density0[p_id] = self.elem_density[ae_id];

		// update velocity
		ShapeFunc& p_N = self.pcl_N[p_id];
		NodeA &n1_a = self.node_a[n1_id];
		NodeA &n2_a = self.node_a[n2_id];
		NodeA &n3_a = self.node_a[n3_id];
		NodeA &n4_a = self.node_a[n4_id];
		PclV& p_v = pcl_v0[p_id];
		p_v.vx += (p_N.N1 * n1_a.ax + p_N.N2 * n2_a.ax
				 + p_N.N3 * n3_a.ax + p_N.N4 * n4_a.ax) * dt;
		p_v.vy += (p_N.N1 * n1_a.ay + p_N.N2 * n2_a.ay
				 + p_N.N3 * n3_a.ay + p_N.N4 * n4_a.ay) * dt;

		// update displacement
		NodeV& n1_v = self.node_v[n1_id];
		NodeV& n2_v = self.node_v[n2_id];
		NodeV& n3_v = self.node_v[n3_id];
		NodeV& n4_v = self.node_v[n4_id];
		PclDisp& p_d = pcl_disp0[p_id];
		p_d.ux += (p_N.N1 * n1_v.vx + p_N.N2 * n2_v.vx
				 + p_N.N3 * n3_v.vx + p_N.N4 * n4_v.vx) * dt;
		p_d.uy += (p_N.N1 * n1_v.vy + p_N.N2 * n2_v.vy
				 + p_N.N3 * n3_v.vy + p_N.N4 * n4_v.vy) * dt;
		
		// update location (in which element)
		PclPos& p_p = self.pcl_pos[ori_pcl_id];
		pcl_x = p_p.x + p_d.ux;
		pcl_y = p_p.y + p_d.uy;
		if (pcl_x < self.mh_xl || pcl_x > self.mh_xu ||
			pcl_y < self.mh_yl || pcl_y > self.mh_yu)
		{
			// not in bg mesh
			pcl_in_elem_array[p_id] = self.actual_elem_num;
		}
		else
		{
			// in bg mesh
			p_ae_x_id = (pcl_x - self.mh_xl) * self.inv_elem_hx;
			p_ae_y_id = (pcl_y - self.mh_yl) * self.inv_elem_hy;
			pcl_in_elem_array[p_id] = self.actual_elem_x_num * p_ae_y_id + p_ae_x_id;
			++pcl_in_mesh_num;
		}
	}

#pragma omp critical
	{
		self.pcl_num += pcl_in_mesh_num;
	}
#pragma omp barrier

#pragma omp master
	{
		if (self.pcl_num)
			self.continue_calculation();
		else
			self.exit_calculation();
	}

	// sort particle variables
	size_t* my_cbin = self.elem_count_bin + my_th_id * 0x100;
	size_t* my_sbin = self.elem_sum_bin + my_th_id * 0x100;
	size_t data_digit, bin_id, th_id;
	size_t* other_cbin;
	for (size_t digit_disp = 0, elem_num_tmp = self.actual_elem_num;
		elem_num_tmp; digit_disp += 8, elem_num_tmp >>= 8)
	{
		memset(my_cbin, 0, 0x100 * sizeof(size_t));

		for (p_id = p_id0; p_id < p_id1; ++p_id)
		{
			data_digit = (pcl_in_elem_array[p_id] >> digit_disp) & 0xFF;
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
			other_cbin = self.elem_count_bin + th_id * 0x100;
			for (bin_id = 0; bin_id < 0x100; ++bin_id)
				my_sbin[bin_id] += other_cbin[bin_id];
		}
		for (th_id = my_th_id + 1; th_id < self.thread_num; ++th_id)
		{
			other_cbin = self.elem_count_bin + th_id * 0x100;
			for (bin_id = 1; bin_id < 0x100; ++bin_id)
				my_sbin[bin_id] += other_cbin[bin_id - 1];
		}

		for (p_id = p_id1; p_id-- > p_id0;)
		{
			data_digit = (pcl_in_elem_array[p_id] >> digit_disp) & 0xFF;
			pcl_in_elem_array_tmp[--my_sbin[data_digit]] = pcl_in_elem_array[p_id];
		}

		pcl_in_elem_array_tmp = pcl_in_elem_array;
		sort_var_id ^= 1;
		pcl_in_elem_array = self.pcl_in_elem_arrays[sort_var_id];
#pragma omp barrier
	}

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
	return 0;
}

int Step_R2D_ME_mt::apply_rigid_rect_avg(
	size_t p_id0,
	size_t p_id1,
	size_t* pcl_in_elem,
	PclSortedVarArray& cur_pscv,
	RigidRectForce& rr_force
	)
{
	double p_x, p_y;
	double dist, norm_x, norm_y;
	double f_cont, fx_cont, fy_cont;
	size_t e_id, pcl_ori_id;
	Model_R2D_ME_mt& md = *(Model_R2D_ME_mt*)model;
	RigidRect& rr = md.get_rigid_rect();
	const Point2D &rr_centre = rr.get_centre();
	for (size_t p_id = p_id0; p_id < p_id1; ++p_id)
	{
		pcl_ori_id = cur_pscv.pcl_index[p_id];
		PclPos& p_p = pcl_pos[pcl_ori_id];
		PclDisp& p_d = cur_pscv.pcl_disp[p_id];
		p_x = p_p.x + p_d.ux;
		p_y = p_p.y + p_d.uy;
		if (rr.detect_collision_with_point(
			p_x, p_y, pcl_vol[p_id],
			dist, norm_x, norm_y))
		{
			f_cont = K_cont * dist;
			fx_cont = f_cont * norm_x;
			fy_cont = f_cont * norm_y;
			// apply contact force to rigid body
			rr_force.add_f_contact(
				p_x, p_y,
				-fx_cont, -fy_cont,
				rr_centre.x, rr_centre.y
				);
			// apply contact force to mesh
			ShapeFunc &p_N = pcl_N[p_id];
			e_id = pcl_in_elem[p_id];
			ElemNodeForce& en1_f = elem_node_force[e_id * 4];
			en1_f.fx += p_N.N1 * fx_cont;
			en1_f.fy += p_N.N1 * fy_cont;
			ElemNodeForce& en2_f = elem_node_force[e_id * 4 + 1];
			en2_f.fx += p_N.N2 * fx_cont;
			en2_f.fy += p_N.N2 * fy_cont;
			ElemNodeForce& en3_f = elem_node_force[e_id * 4 + 2];
			en3_f.fx += p_N.N3 * fx_cont;
			en3_f.fy += p_N.N3 * fy_cont;
			ElemNodeForce& en4_f = elem_node_force[e_id * 4 + 3];
			en4_f.fx += p_N.N4 * fx_cont;
			en4_f.fy += p_N.N4 * fy_cont;
		}
	}

	return 0;
}
