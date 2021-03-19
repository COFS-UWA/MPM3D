#include "SimulationsOMP_pcp.h"

#include <assert.h>
#include "tbb/task_arena.h"

#include "ParallelUtils.h"
#include "Step_T2D_CHM_Task.h"
#include "Step_T2D_CHM_TBB.h"

namespace Step_T2D_CHM_Task
{
	constexpr double one_third = 1.0 / 3.0;
	constexpr double one_fourth = 0.25;

	void CalData::set_model(Model_T2D_CHM_mt &md) noexcept
	{
		pmodel = &md;

		pcl_m_s = md.pcl_m_s;
		pcl_density_s = md.pcl_density_s;
		pcl_vol_s = md.pcl_vol_s;
		pcl_bf_s = md.pcl_bf_s;
		pcl_bf_f = md.pcl_bf_f;
		pcl_t = md.pcl_t;
		pcl_pos = md.pcl_pos;
		pcl_vol = md.pcl_vol;
		pcl_mat_model = md.pcl_mat_model;
		
		auto& spva0 = spvas[0];
		const auto& md_spva0 = md.sorted_pcl_var_arrays[0];
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

		auto& spva1 = spvas[1];
		const auto& md_spva1 = md.sorted_pcl_var_arrays[1];
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

		elem_node_id = md.elem_node_id;
		elem_N_ab = md.elem_N_ab;
		elem_N_c = md.elem_N_c;
		elem_area = md.elem_area;
		elem_density_f = md.elem_density_f;
		elem_pcl_n = md.elem_pcl_n;
		elem_pcl_m_s = md.elem_pcl_m_s;
		elem_pcl_m_f = md.elem_pcl_m_f;
		elem_de = md.elem_de;
		elem_p = md.elem_p;
		elem_n2_miu_div_k_vol = md.elem_n2_miu_div_k_vol;
		elem_seep_force = md.elem_seep_force;
		elem_m_de_vol_s = md.elem_m_de_vol_s;
		elem_m_de_vol_f = md.elem_m_de_vol_f;

		// element-node data
		elem_node_vm_s = md.elem_node_vm_s;
		elem_node_vm_f = md.elem_node_vm_f;
		elem_node_force_s = md.elem_node_force_s;
		elem_node_force_f = md.elem_node_force_f;

		// node data
		node_pos = md.node_pos;
		node_a_s = md.node_a_s;
		node_a_f = md.node_a_f;
		node_v_s = md.node_v_s;
		node_v_f = md.node_v_f;
		node_has_vbc_s = md.node_has_vbc_s;
		node_has_vbc_f = md.node_has_vbc_f;
		node_am_s = md.node_am_s;
		node_am_f = md.node_am_f;
		node_de_vol_s = md.node_de_vol_s;
		node_de_vol_f = md.node_de_vol_f;

		Kf = md.Kf;
		miu = md.miu;
		k = md.k;

#ifdef _DEBUG
		ori_pcl_num = md.ori_pcl_num;
		elem_num = md.elem_num;
		node_num = md.node_num;
#endif
	}

	void InitPcl::operator() (size_t wk_id, size_t &pcl_in_mesh_num) const
	{
		const auto& pcl_sort_mem = cd.pcl_sort_mem;
		size_t* const ori_pcl_in_elem = pcl_sort_mem.ori_keys;
		size_t* const ori_cur_to_prev_pcl = pcl_sort_mem.ori_vals;
		Model_T2D_CHM_mt& md = *cd.pmodel;
		Position* const pcl_pos = const_cast<Position* const>(cd.pcl_pos);
		const auto& spva0 = cd.spvas[0];
		const size_t p_id0 = ParallelUtils::block_low(wk_id, task_num, cd.prev_valid_pcl_num);
		const size_t p_id1 = ParallelUtils::block_low(wk_id + 1, task_num, cd.prev_valid_pcl_num);
		size_t my_valid_pcl_num = 0;
		for (size_t p_id = p_id0; p_id < p_id1; ++p_id)
		{
			const size_t ori_p_id = spva0.pcl_index[p_id];
			Position& p_p = pcl_pos[ori_p_id];
			Displacement& p_d = spva0.pcl_u_s[p_id];
			p_p.x += p_d.ux;
			p_p.y += p_d.uy;
			p_d.ux = 0.0;
			p_d.uy = 0.0;
			ShapeFunc& p_N = spva0.pcl_N[p_id];
			size_t e_id = md.find_pcl_in_which_elem(p_p.x, p_p.y, p_N);
			if (e_id == SIZE_MAX)
				e_id = md.find_pcl_in_which_elem_tol(p_p.x, p_p.y, p_N);
			if (e_id != SIZE_MAX)
				++my_valid_pcl_num;
			ori_pcl_in_elem[p_id] = e_id;
			ori_cur_to_prev_pcl[p_id] = p_id;
		}
		pcl_in_mesh_num = my_valid_pcl_num;
	}
	
	void MapPclToBgMesh::operator() (size_t wk_id) const
	{
		size_t e_id;
		size_t p_id0 = block_low(wk_id, task_num, valid_pcl_num);
		e_id = pcl_in_elem[p_id0];
		while (p_id0 != SIZE_MAX && e_id == pcl_in_elem[--p_id0]);
		++p_id0;
		assert(p_id0 <= valid_pcl_num);
		size_t p_id1 = block_low(wk_id + 1, task_num, valid_pcl_num);
		e_id = pcl_in_elem[p_id1];
		while (p_id1 != SIZE_MAX && e_id == pcl_in_elem[--p_id1]);
		++p_id1;
		assert(p_id1 <= valid_pcl_num);
		
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
		double en1_vm_f = 0.0;
		double en1_vmx_f = 0.0;
		double en1_vmy_f = 0.0;
		double en1_vmz_f = 0.0;
		double en2_vm_f = 0.0;
		double en2_vmx_f = 0.0;
		double en2_vmy_f = 0.0;
		double en2_vmz_f = 0.0;
		double en3_vm_f = 0.0;
		double en3_vmx_f = 0.0;
		double en3_vmy_f = 0.0;
		double en3_vmz_f = 0.0;
		double en4_vm_f = 0.0;
		double en4_vmx_f = 0.0;
		double en4_vmy_f = 0.0;
		double en4_vmz_f = 0.0;
		double en1_fx_seep = 0.0;
		double en2_fx_seep = 0.0;
		double en3_fx_seep = 0.0;
		double en4_fx_seep = 0.0;
		double en1_fy_seep = 0.0;
		double en2_fy_seep = 0.0;
		double en3_fy_seep = 0.0;
		double en4_fy_seep = 0.0;
		double en1_fz_seep = 0.0;
		double en2_fz_seep = 0.0;
		double en3_fz_seep = 0.0;
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
		e_id = pcl_in_elem[p_id0];
#ifdef _DEBUG
		assert(e_id < cd.elem_num);
#endif
		double p_N_m;
		for (size_t p_id = p_id0; p_id < p_id1; ++p_id)
		{
			// pcl index
			const size_t prev_p_id = cur_to_prev_pcl[p_id];
			assert(prev_p_id < cd.prev_valid_pcl_num);
			const size_t ori_p_id = pcl_index1[prev_p_id];
#ifdef _DEBUG
			assert(ori_p_id < cd.ori_pcl_num);
#endif
			pcl_index0[p_id] = ori_p_id;

			// m_s
			const double p_m_s = pcl_m_s[ori_p_id];
			e_p_m_s += p_m_s;
			// e_n
			e_n += pcl_vol_s[ori_p_id];
			// vol
			const double p_n = pcl_n1[prev_p_id];
			const double p_vol = pcl_vol_s[ori_p_id] / (1.0 - p_n);
			pcl_vol[p_id] = p_vol;
			e_p_vol += p_vol;
			// vol_f
			const double p_vol_f = p_n * p_vol;
			e_p_vol_f += p_vol_f;
			// m_f
			const double p_m_f = pcl_density_f1[prev_p_id] * p_vol_f;
			e_p_m_f += p_m_f;

			// map pcl stress
			const Stress& p_s1 = pcl_stress1[prev_p_id];
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
			
			// shape function
			const ShapeFunc& p_N1 = pcl_N1[prev_p_id];
			ShapeFunc& p_N0 = pcl_N0[p_id];
			p_N0.N1 = p_N1.N1;
			p_N0.N2 = p_N1.N2;
			p_N0.N3 = p_N1.N3;
			p_N0.N4 = p_N1.N4;
			// solid velocity
			const Velocity& p_v_s1 = pcl_v_s1[prev_p_id];
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
			const Velocity& p_v_f1 = pcl_v_f1[prev_p_id];
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
			const Displacement& p_u_s1 = pcl_u_s1[prev_p_id];
			Displacement& p_u_s0 = pcl_u_s0[p_id];
			p_u_s0.ux = p_u_s1.ux;
			p_u_s0.uy = p_u_s1.uy;
			p_u_s0.uz = p_u_s1.uz;
			const Displacement& p_u_f1 = pcl_u_f1[prev_p_id];
			Displacement& p_u_f0 = pcl_u_f0[p_id];
			p_u_f0.ux = p_u_f1.ux;
			p_u_f0.uy = p_u_f1.uy;
			p_u_f0.uz = p_u_f1.uz;

			// seepage force
			const double f_seep_tmp = cd.miu / cd.k * p_n * p_n * p_vol;
			const double e_fx_seep = (p_v_f0.vx - p_v_s0.vx) * f_seep_tmp;
			en1_fx_seep += p_N0.N1 * e_fx_seep;
			en2_fx_seep += p_N0.N2 * e_fx_seep;
			en3_fx_seep += p_N0.N3 * e_fx_seep;
			en4_fx_seep += p_N0.N4 * e_fx_seep;
			const double e_fy_seep = (p_v_f0.vy - p_v_s0.vy) * f_seep_tmp;
			en1_fy_seep += p_N0.N1 * e_fy_seep;
			en2_fy_seep += p_N0.N2 * e_fy_seep;
			en3_fy_seep += p_N0.N3 * e_fy_seep;
			en4_fy_seep += p_N0.N4 * e_fy_seep;
			const double e_fz_seep = (p_v_f0.vz - p_v_s0.vz) * f_seep_tmp;
			en1_fz_seep += p_N0.N1 * e_fz_seep;
			en2_fz_seep += p_N0.N2 * e_fz_seep;
			en3_fz_seep += p_N0.N3 * e_fz_seep;
			en4_fz_seep += p_N0.N4 * e_fz_seep;
			
			// solid external load
			const Force& p_bf_s = pcl_bf_s[ori_p_id];
			const double one_fourth_bfx_s = one_fourth * p_bf_s.fx;
			const double one_fourth_bfy_s = one_fourth * p_bf_s.fy;
			const double one_fourth_bfz_s = one_fourth * p_bf_s.fz;
			const Force& p_t = pcl_t[ori_p_id];
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
			const Force& p_bf_f = pcl_bf_f[ori_p_id];
			const double one_fourth_bfx_f = one_fourth * p_bf_f.fx;
			const double one_fourth_bfy_f = one_fourth * p_bf_f.fy;
			const double one_fourth_bfz_f = one_fourth * p_bf_f.fz;
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
			
			if (e_id != pcl_in_elem[p_id + 1])
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
				
				e_id = pcl_in_elem[p_id + 1];
#ifdef _DEBUG
				assert(e_id < cd.elem_num || e_id == SIZE_MAX);
#endif

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
				en2_fx_seep = 0.0;
				en3_fx_seep = 0.0;
				en4_fx_seep = 0.0;
				en1_fy_seep = 0.0;
				en2_fy_seep = 0.0;
				en3_fy_seep = 0.0;
				en4_fy_seep = 0.0;
				en1_fz_seep = 0.0;
				en2_fz_seep = 0.0;
				en3_fz_seep = 0.0;
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
	}

	//void ContactRigidRect::operator() (size_t wk_id, Force3D& rr_cf) const
	//{
	//	size_t e_id;
	//	size_t p_id0 = block_low(wk_id, task_num, valid_pcl_num);
	//	e_id = pcl_in_elem[p_id0];
	//	while (p_id0 != SIZE_MAX && e_id == pcl_in_elem[--p_id0]);
	//	++p_id0;
	//	assert(p_id0 <= valid_pcl_num);
	//	size_t p_id1 = block_low(wk_id + 1, task_num, valid_pcl_num);
	//	e_id = pcl_in_elem[p_id1];
	//	while (p_id1 != SIZE_MAX && e_id == pcl_in_elem[--p_id1]);
	//	++p_id1;
	//	assert(p_id1 <= valid_pcl_num);

	//	double dist;
	//	Vector2D lnorm;
	//	Force lcont_f, gcont_f;
	//	Point2D cur_cont_pos;
	//	Force2D rcf;
	//	rcf.reset();
	//	for (size_t p_id = p_id0; p_id < p_id1; ++p_id)
	//	{
	//		const size_t ori_p_id = pcl_index[p_id];
	//		const Position& p_p = pcl_pos[ori_p_id];
	//		const Displacement& p_d = pcl_disp[p_id];
	//		const double p_x = p_p.x + p_d.ux;
	//		const double p_y = p_p.y + p_d.uy;
	//		const double p_r = 0.5 * sqrt(pcl_vol[p_id]);
	//		if (prr->detect_collision_with_point(
	//			p_x, p_y, p_r, dist, lnorm, cur_cont_pos))
	//		{
	//			const double f_cont = K_cont * dist;
	//			lcont_f.fx = f_cont * lnorm.x;
	//			lcont_f.fy = f_cont * lnorm.y;
	//			prr->get_global_vector(lcont_f.vec, gcont_f.vec);
	//			// apply contact force to mesh
	//			const ShapeFunc& p_N = pcl_N[p_id];
	//			const size_t e_id = pcl_in_elem[p_id];
	//			Force &en_f1 = elem_node_force[e_id * 3];
	//			en_f1.fx += p_N.N1 * gcont_f.fx;
	//			en_f1.fy += p_N.N1 * gcont_f.fy;
	//			Force& en_f2 = elem_node_force[e_id * 3 + 1];
	//			en_f2.fx += p_N.N2 * gcont_f.fx;
	//			en_f2.fy += p_N.N2 * gcont_f.fy;
	//			Force& en_f3 = elem_node_force[e_id * 3 + 2];
	//			en_f3.fx += p_N.N3 * gcont_f.fx;
	//			en_f3.fy += p_N.N3 * gcont_f.fy;
	//			// apply contact force to rigid body
	//			const Point2D& rr_cen = prr->get_centre();
	//			rcf.add_force(p_x, p_y,
	//				-gcont_f.fx, -gcont_f.fy,
	//				rr_cen.x, rr_cen.y);
	//		}
	//	}
	//	rr_cf.fx = rcf.fx;
	//	rr_cf.fy = rcf.fy;
	//	rr_cf.m = rcf.m;
	//}

	void UpdateAccelerationAndVelocity::operator() (size_t wk_id) const
	{
		size_t n_id;
		size_t ve_id0 = block_low(wk_id, task_num, four_valid_elem_num);
		n_id = node_has_elem[ve_id0];
		while (ve_id0 != SIZE_MAX && n_id == node_has_elem[--ve_id0]);
		++ve_id0;
		assert(ve_id0 <= four_valid_elem_num);
		size_t ve_id1 = block_low(wk_id + 1, task_num, four_valid_elem_num);
		n_id = node_has_elem[ve_id1];
		while (ve_id1 != SIZE_MAX && n_id == node_has_elem[--ve_id1]);
		++ve_id1;
		assert(ve_id1 <= four_valid_elem_num);
		
		size_t bc_mask, ne_id;
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
		n_id = node_has_elem[ve_id0];
#ifdef _DEBUG
		assert(n_id < cd.node_num);
#endif
		for (size_t ve_id = ve_id0; ve_id < ve_id1; ++ve_id)
		{
			ne_id = node_elem_pair[ve_id];
#ifdef _DEBUG
			assert(ne_id < cd.elem_num * 4);
#endif
			n_am_s += elem_pcl_m_s[ne_id / 4];
			n_am_f += elem_pcl_m_f[ne_id / 4];
			const Force& nf_s = elem_node_force_s[ne_id];
			n_fx_s += nf_s.fx;
			n_fy_s += nf_s.fy;
			n_fz_s += nf_s.fz;
			const Force& nf_f = elem_node_force_f[ne_id];
			n_fx_f += nf_f.fx;
			n_fy_f += nf_f.fy;
			n_fz_f += nf_f.fz;

			const ElemNodeVM& nvm_s = elem_node_vm_s[ne_id];
			n_vm_s += nvm_s.vm;
			n_vmx_s += nvm_s.vmx;
			n_vmy_s += nvm_s.vmy;
			n_vmz_s += nvm_s.vmz;
			const ElemNodeVM& nvm_f = elem_node_vm_f[ne_id];
			n_vm_f += nvm_f.vm;
			n_vmx_f += nvm_f.vmx;
			n_vmy_f += nvm_f.vmy;
			n_vmz_f += nvm_f.vmz;

			if (n_id != node_has_elem[ve_id + 1])
			{
				// solid
				n_am_s *= one_fourth;
				node_am_s[n_id] = n_am_s;
				Acceleration& n_a_s = node_a_s[n_id];
				n_a_s.ax = n_fx_s / n_am_s;
				n_a_s.ay = n_fy_s / n_am_s;
				n_a_s.az = n_fz_s / n_am_s;
				Velocity& n_v_s = node_v_s[n_id];
				n_v_s.vx = n_vmx_s / n_vm_s + n_a_s.ax * cd.dt;
				n_v_s.vy = n_vmy_s / n_vm_s + n_a_s.ay * cd.dt;
				n_v_s.vz = n_vmz_s / n_vm_s + n_a_s.az * cd.dt;
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
				n_v_f.vx = n_vmx_f / n_vm_f + n_a_f.ax * cd.dt;
				n_v_f.vy = n_vmy_f / n_vm_f + n_a_f.ay * cd.dt;
				n_v_f.vz = n_vmz_f / n_vm_f + n_a_f.az * cd.dt;
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

				n_id = node_has_elem[ve_id + 1];
#ifdef _DEBUG
				assert(n_id < cd.node_num || n_id == SIZE_MAX);
#endif // _DEBUG

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
	}

	void CalElemDeAndMapToNode::operator() (size_t wk_id) const
	{
		const size_t ve_id0 = block_low(wk_id, task_num, valid_elem_num);
		const size_t ve_id1 = block_low(wk_id + 1, task_num, valid_elem_num);
		double e_de_vol_s, e_de_vol_f;
		for (size_t ve_id = ve_id0; ve_id < ve_id1; ++ve_id)
		{
			const size_t e_id = valid_elems[ve_id];
#ifdef _DEBUG
			assert(e_id < cd.elem_num);
#endif
			const ElemNodeIndex& eni = elem_node_id[e_id];
			const Velocity& n1_v_s = node_v_s[eni.n1];
			const Velocity& n2_v_s = node_v_s[eni.n2];
			const Velocity& n3_v_s = node_v_s[eni.n3];
			const Velocity& n4_v_s = node_v_s[eni.n4];
			const DShapeFuncABC& e_dN = elem_N_abc[e_id];
			StrainInc& e_de = elem_de[e_id];
			e_de.de11 = (e_dN.dN1_dx * n1_v_s.vx + e_dN.dN2_dx * n2_v_s.vx + e_dN.dN3_dx * n3_v_s.vx + e_dN.dN4_dx * n4_v_s.vx) * cd.dt;
			e_de.de22 = (e_dN.dN1_dy * n1_v_s.vy + e_dN.dN2_dy * n2_v_s.vy + e_dN.dN3_dy * n3_v_s.vy + e_dN.dN4_dy * n4_v_s.vy) * cd.dt;
			e_de.de33 = (e_dN.dN1_dz * n1_v_s.vz + e_dN.dN2_dz * n2_v_s.vz + e_dN.dN3_dz * n3_v_s.vz + e_dN.dN4_dz * n4_v_s.vz) * cd.dt;
			e_de.de12 = (e_dN.dN1_dx * n1_v_s.vy + e_dN.dN2_dx * n2_v_s.vy + e_dN.dN3_dx * n3_v_s.vy + e_dN.dN4_dx * n4_v_s.vy
					   + e_dN.dN1_dy * n1_v_s.vx + e_dN.dN2_dy * n2_v_s.vx + e_dN.dN3_dy * n3_v_s.vx + e_dN.dN4_dy * n4_v_s.vx) * cd.dt * 0.5;
			e_de.de23 = (e_dN.dN1_dy * n1_v_s.vz + e_dN.dN2_dy * n2_v_s.vz + e_dN.dN3_dy * n3_v_s.vz + e_dN.dN4_dy * n4_v_s.vz
					   + e_dN.dN1_dz * n1_v_s.vy + e_dN.dN2_dz * n2_v_s.vy + e_dN.dN3_dz * n3_v_s.vy + e_dN.dN4_dz * n4_v_s.vy) * cd.dt * 0.5;
			e_de.de31 = (e_dN.dN1_dz * n1_v_s.vx + e_dN.dN2_dz * n2_v_s.vx + e_dN.dN3_dz * n3_v_s.vx + e_dN.dN4_dz * n4_v_s.vx
					   + e_dN.dN1_dx * n1_v_s.vz + e_dN.dN2_dx * n2_v_s.vz + e_dN.dN3_dx * n3_v_s.vz + e_dN.dN4_dx * n4_v_s.vz) * cd.dt * 0.5;
			e_de_vol_s = e_de.de11 + e_de.de22 + e_de.de33;
			elem_m_de_vol_s[e_id] = elem_pcl_m_s[e_id] * e_de_vol_s;
			const Velocity& n1_v_f = node_v_f[eni.n1];
			const Velocity& n2_v_f = node_v_f[eni.n2];
			const Velocity& n3_v_f = node_v_f[eni.n3];
			const Velocity& n4_v_f = node_v_f[eni.n4];
			e_de_vol_f = (1.0 - elem_pcl_n[e_id]) / elem_pcl_n[e_id] * -e_de_vol_s
					- (e_dN.dN1_dx * n1_v_f.vx + e_dN.dN2_dx * n2_v_f.vx + e_dN.dN3_dx * n3_v_f.vx + e_dN.dN4_dx * n4_v_f.vx
					+ e_dN.dN1_dy * n1_v_f.vy + e_dN.dN2_dy * n2_v_f.vy + e_dN.dN3_dy * n3_v_f.vy + e_dN.dN4_dy * n4_v_f.vy
					+ e_dN.dN1_dz * n1_v_f.vz + e_dN.dN2_dz * n2_v_f.vz + e_dN.dN3_dz * n3_v_f.vz + e_dN.dN4_dz * n4_v_f.vz) * cd.dt;
			elem_m_de_vol_f[e_id] = elem_pcl_m_f[e_id] * e_de_vol_f;
			e_de_vol_s *= one_third;
			e_de.de11 -= e_de_vol_s;
			e_de.de22 -= e_de_vol_s;
			e_de.de33 -= e_de_vol_s;
		}
	}
	
	void CalNodeDe::operator() (size_t wk_id) const
	{
		size_t n_id;
		size_t ve_id0 = block_low(wk_id, task_num, four_valid_elem_num);
		n_id = node_has_elem[ve_id0];
		while (ve_id0 != SIZE_MAX && n_id == node_has_elem[--ve_id0]);
		++ve_id0;
		assert(ve_id0 <= four_valid_elem_num);
		size_t ve_id1 = block_low(wk_id + 1, task_num, four_valid_elem_num);
		n_id = node_has_elem[ve_id1];
		while (ve_id1 != SIZE_MAX && n_id == node_has_elem[--ve_id1]);
		++ve_id1;
		assert(ve_id1 <= four_valid_elem_num);

		size_t e_id;
		double n_am_de_vol_s = 0.0;
		double n_am_de_vol_f = 0.0;
		n_id = node_has_elem[ve_id0];
#ifdef _DEBUG
		assert(n_id < cd.node_num || n_id == SIZE_MAX);
#endif
		for (size_t ve_id = ve_id0; ve_id < ve_id1; ++ve_id)
		{
			e_id = node_elem_pair[ve_id] / 4;
#ifdef _DEBUG
			assert(e_id < cd.elem_num);
#endif
			n_am_de_vol_s += elem_m_de_vol_s[e_id];
			n_am_de_vol_f += elem_m_de_vol_f[e_id];
			if (n_id != node_has_elem[ve_id + 1])
			{
				node_de_vol_s[n_id] = n_am_de_vol_s * one_fourth / node_am_s[n_id];
				node_de_vol_f[n_id] = n_am_de_vol_f * one_fourth / node_am_f[n_id];
				n_id = node_has_elem[ve_id + 1];
#ifdef _DEBUG
				assert(n_id < cd.node_num || n_id == SIZE_MAX);
#endif
				n_am_de_vol_s = 0.0;
				n_am_de_vol_f = 0.0;
			}
		}
	}

	void MapBgMeshToPcl::operator() (size_t wk_id, size_t& pcl_in_mesh_num) const
	{
		size_t e_id;
		size_t p_id0 = block_low(wk_id, task_num, valid_pcl_num);
		e_id = pcl_in_elem[p_id0];
		while (p_id0 != SIZE_MAX && e_id == pcl_in_elem[--p_id0]);
		++p_id0;
		assert(p_id0 <= valid_pcl_num);
		size_t p_id1 = block_low(wk_id + 1, task_num, valid_pcl_num);
		e_id = pcl_in_elem[p_id1];
		while (p_id1 != SIZE_MAX && e_id == pcl_in_elem[--p_id1]);
		++p_id1;
		assert(p_id1 <= valid_pcl_num);
		
		const Acceleration* pn1_a_s, * pn2_a_s, * pn3_a_s, * pn4_a_s;
		const Acceleration* pn1_a_f, * pn2_a_f, * pn3_a_f, * pn4_a_f;
		const Velocity* pn1_v_s, * pn2_v_s, * pn3_v_s, * pn4_v_s;
		const Velocity* pn1_v_f, * pn2_v_f, * pn3_v_f, * pn4_v_f;
		double e_de_vol_s, e_de_vol_f, e_density_f, e_n, e_p;
		StrainInc* pe_de;
		double dstrain[6];
		Model_T2D_CHM_mt& md = *(cd.pmodel);
		size_t my_valid_pcl_num = 0;
		e_id = SIZE_MAX;
		for (size_t p_id = p_id0; p_id < p_id1; ++p_id)
		{
			if (e_id != pcl_in_elem[p_id])
			{
				e_id = pcl_in_elem[p_id];
#ifdef _DEBUG
				assert(e_id < cd.elem_num);
#endif
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
				e_p = elem_p[e_id] + cd.Kf * e_de_vol_f;

				pe_de = elem_de + e_id;
				e_de_vol_s *= one_third;
				pe_de->de11 += e_de_vol_s;
				pe_de->de22 += e_de_vol_s;
				pe_de->de33 += e_de_vol_s;
			}

			// update velocity
			ShapeFunc& p_N = pcl_N0[p_id];
			Velocity& p_v_s0 = pcl_v_s0[p_id];
			p_v_s0.vx += (p_N.N1 * pn1_a_s->ax + p_N.N2 * pn2_a_s->ax + p_N.N3 * pn3_a_s->ax + p_N.N4 * pn4_a_s->ax) * cd.dt;
			p_v_s0.vy += (p_N.N1 * pn1_a_s->ay + p_N.N2 * pn2_a_s->ay + p_N.N3 * pn3_a_s->ay + p_N.N4 * pn4_a_s->ay) * cd.dt;
			p_v_s0.vz += (p_N.N1 * pn1_a_s->az + p_N.N2 * pn2_a_s->az + p_N.N3 * pn3_a_s->az + p_N.N4 * pn4_a_s->az) * cd.dt;
			Velocity& p_v_f0 = pcl_v_f0[p_id];
			p_v_f0.vx += (p_N.N1 * pn1_a_f->ax + p_N.N2 * pn2_a_f->ax + p_N.N3 * pn3_a_f->ax + p_N.N4 * pn4_a_f->ax) * cd.dt;
			p_v_f0.vy += (p_N.N1 * pn1_a_f->ay + p_N.N2 * pn2_a_f->ay + p_N.N3 * pn3_a_f->ay + p_N.N4 * pn4_a_f->ay) * cd.dt;
			p_v_f0.vz += (p_N.N1 * pn1_a_f->az + p_N.N2 * pn2_a_f->az + p_N.N3 * pn3_a_f->az + p_N.N4 * pn4_a_f->az) * cd.dt;

			// update displacement
			Displacement& p_u_s0 = pcl_u_s0[p_id];
			p_u_s0.ux += (p_N.N1 * pn1_v_s->vx + p_N.N2 * pn2_v_s->vx + p_N.N3 * pn3_v_s->vx + p_N.N4 * pn4_v_s->vx) * cd.dt;
			p_u_s0.uy += (p_N.N1 * pn1_v_s->vy + p_N.N2 * pn2_v_s->vy + p_N.N3 * pn3_v_s->vy + p_N.N4 * pn4_v_s->vy) * cd.dt;
			p_u_s0.uz += (p_N.N1 * pn1_v_s->vz + p_N.N2 * pn2_v_s->vz + p_N.N3 * pn3_v_s->vz + p_N.N4 * pn4_v_s->vz) * cd.dt;
			Displacement& p_u_f0 = pcl_u_f0[p_id];
			p_u_f0.ux += (p_N.N1 * pn1_v_f->vx + p_N.N2 * pn2_v_f->vx + p_N.N3 * pn3_v_f->vx + p_N.N4 * pn4_v_f->vx) * cd.dt;
			p_u_f0.uy += (p_N.N1 * pn1_v_f->vy + p_N.N2 * pn2_v_f->vy + p_N.N3 * pn3_v_f->vy + p_N.N4 * pn4_v_f->vy) * cd.dt;
			p_u_f0.uz += (p_N.N1 * pn1_v_f->vz + p_N.N2 * pn2_v_f->vz + p_N.N3 * pn3_v_f->vz + p_N.N4 * pn4_v_f->vz) * cd.dt;

			// update location (in which element)
			const size_t ori_p_id = pcl_index0[p_id];
#ifdef _DEBUG
			assert(ori_p_id < cd.ori_pcl_num);
#endif

			const Position& p_p = pcl_pos[ori_p_id];
			const double p_x = p_p.x + p_u_s0.ux;
			const double p_y = p_p.y + p_u_s0.uy;
			const double p_z = p_p.z + p_u_s0.uz;
			size_t p_e_id = e_id;
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
				++my_valid_pcl_num;
			new_pcl_in_elem[p_id] = p_e_id;
			new_cur_to_prev_pcl[p_id] = p_id;
#ifdef _DEBUG
			assert(p_e_id < cd.elem_num || p_e_id == SIZE_MAX);
#endif

			// update n
			pcl_n0[p_id] = e_n;
			// update density
			pcl_density_f0[p_id] = e_density_f;
			// update pore pressure
			pcl_p0[p_id] = e_p;

			// update stress
			MatModel::MaterialModel& pcl_mm = *pcl_mat_model[ori_p_id];
			dstrain[0] = pe_de->de11;
			dstrain[1] = pe_de->de22;
			dstrain[2] = pe_de->de33;
			dstrain[3] = pe_de->de12;
			dstrain[4] = pe_de->de23;
			dstrain[5] = pe_de->de31;
			pcl_mm.integrate(dstrain);
			const double* dstress = pcl_mm.get_dstress();
			Stress& p_s = pcl_stress0[p_id];
			p_s.s11 += dstress[0];
			p_s.s22 += dstress[1];
			p_s.s33 += dstress[2];
			p_s.s12 += dstress[3];
			p_s.s23 += dstress[4];
			p_s.s31 += dstress[5];

			const size_t prev_p_id = cur_to_prev_pcl[p_id];
			assert(prev_p_id < cd.prev_valid_pcl_num);

			const Strain& p_e1 = pcl_strain1[prev_p_id];
			Strain& p_e0 = pcl_strain0[p_id];
			p_e0.e11 = p_e1.e11 + pe_de->de11;
			p_e0.e22 = p_e1.e22 + pe_de->de22;
			p_e0.e33 = p_e1.e33 + pe_de->de33;
			p_e0.e12 = p_e1.e12 + pe_de->de12;
			p_e0.e23 = p_e1.e23 + pe_de->de23;
			p_e0.e31 = p_e1.e31 + pe_de->de31;

			const double* estrain = pcl_mm.get_dstrain_e();
			const Strain& p_ee1 = pcl_estrain1[prev_p_id];
			Strain& p_ee0 = pcl_estrain0[p_id];
			p_ee0.e11 = p_ee1.e11 + estrain[0];
			p_ee0.e22 = p_ee1.e22 + estrain[1];
			p_ee0.e33 = p_ee1.e33 + estrain[2];
			p_ee0.e12 = p_ee1.e12 + estrain[3];
			p_ee0.e23 = p_ee1.e23 + estrain[4];
			p_ee0.e31 = p_ee1.e31 + estrain[5];

			const double *pstrain = pcl_mm.get_dstrain_p();
			const Strain& p_pe1 = pcl_pstrain1[prev_p_id];
			Strain& p_pe0 = pcl_pstrain0[p_id];
			p_pe0.e11 = p_pe1.e11 + pstrain[0];
			p_pe0.e22 = p_pe1.e22 + pstrain[1];
			p_pe0.e33 = p_pe1.e33 + pstrain[2];
			p_pe0.e12 = p_pe1.e12 + pstrain[3];
			p_pe0.e23 = p_pe1.e23 + pstrain[4];
			p_pe0.e31 = p_pe1.e31 + pstrain[5];
		}

		pcl_in_mesh_num = my_valid_pcl_num;
	}
}
