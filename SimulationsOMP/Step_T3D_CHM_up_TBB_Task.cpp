#include "SimulationsOMP_pcp.h"

#include <iostream>
#include <assert.h>

#include "ParaUtil.h"
#include "Step_T3D_CHM_up_TBB.h"

#define Block_Low(blk_id, blk_num, data_num) ((blk_id) * (data_num) / (blk_num))

namespace Step_T3D_CHM_up_TBB_Task
{
	constexpr double one_third = 1.0 / 3.0;
	constexpr double one_fourth = 0.25;

	void InitPcl::init(size_t thread_num) noexcept
	{
		Model_T3D_CHM_up_mt& md = *stp.pmodel;
		pcl_vol_s = md.pcl_vol_s;
		const auto& spva0 = stp.spvas[0];
		pcl_n0 = spva0.pcl_n;
		in_pcl_in_elems = stp.in_pcl_in_elems;
		in_prev_pcl_ids = stp.in_prev_pcl_ids;
		task_num = ParaUtil::cal_task_num<
			min_pcl_num_per_task, init_pcl_task_num_per_thread>(
				thread_num, stp.prev_valid_pcl_num);
	}

	void InitPcl::work(size_t tsk_id, InitPclRes &res) const
	{
		Model_T3D_CHM_up_mt& md = *stp.pmodel;
		Position* const pcl_pos = const_cast<Position* const>(md.pcl_pos);
		const auto& spva0 = stp.spvas[0];
		const size_t p_id0 = Block_Low(tsk_id, task_num, stp.prev_valid_pcl_num);
		const size_t p_id1 = Block_Low(tsk_id + 1, task_num, stp.prev_valid_pcl_num);
		size_t valid_pcl_num = 0;
		double max_pcl_vol = 0.0;
		for (size_t p_id = p_id0; p_id < p_id1; ++p_id)
		{
			const size_t ori_p_id = spva0.pcl_index[p_id];
			Position& p_p = pcl_pos[ori_p_id];
			Displacement& p_u = spva0.pcl_u[p_id];
			p_p.x += p_u.ux;
			p_p.y += p_u.uy;
			p_p.z += p_u.uz;
			p_u.ux = 0.0;
			p_u.uy = 0.0;
			p_u.uz = 0.0;
			ShapeFunc& p_N = spva0.pcl_N[p_id];
			size_t e_id = md.find_pcl_in_which_elem(p_p.x, p_p.y, p_p.z, p_N);
			if (e_id == SIZE_MAX)
				e_id = md.find_pcl_in_which_elem_tol(p_p.x, p_p.y, p_p.z, p_N);
			if (e_id != SIZE_MAX)
				++valid_pcl_num;
			const double pcl_vol = pcl_vol_s[ori_p_id] / (1.0 - pcl_n0[ori_p_id]);
			if (max_pcl_vol < pcl_vol)
				max_pcl_vol = pcl_vol;
			in_pcl_in_elems[p_id] = e_id;
			in_prev_pcl_ids[p_id] = p_id;
		}
		res.pcl_num = valid_pcl_num;
		res.max_pcl_vol = max_pcl_vol;
	}
	
	void MapPclToBgMesh::init() noexcept
	{
		Model_T3D_CHM_up_mt& md = *stp.pmodel;
		
		Kf0 = md.Kf;
		k = md.k;
		dyn_viscosity = md.dyn_viscosity;
		m_cav = md.m_cav;
		f_cav_end = md.f_cav_end;
		u_cav_off = pow(1.0 / f_cav_end - 1.0, 1.0 / m_cav) - 1.0;
		u_div_u_cav_lim = pow(Step_T3D_CHM_up_TBB_Task::max_Kf_ratio_divider - 1.0, 1.0 / m_cav);

		// pcl range
		pcl_in_elems = stp.pcl_in_elems;
		prev_pcl_ids = stp.prev_pcl_ids;
		
		// pcl data
		pcl_m_s = md.pcl_m_s;
		pcl_vol_s = md.pcl_vol_s;
		pcl_vol = md.pcl_vol;
		pcl_bf_s = md.pcl_bf_s;
		pcl_bf_f = md.pcl_bf_f;
		pcl_t = md.pcl_t;
		
		// bg mesh data
		elem_dN_abc = md.elem_dN_abc;
		elem_vol = md.elem_vol;
		elem_has_pcls = md.elem_has_pcls;
		elem_pcl_m = md.elem_pcl_m;
		elem_pcl_pm = md.elem_pcl_pm;
		elem_pcl_n = md.elem_pcl_n;
		elem_p = md.elem_p;
		elem_density_f = md.elem_density_f;
		elem_pcl_vol = md.elem_pcl_vol;
		elem_node_vm_s = md.elem_node_vm_s;
		elem_node_p = md.elem_node_p;
		elem_node_force = md.elem_node_force;
		elem_node_p_force = md.elem_node_p_force;
		elem_node_at_surface = md.elem_node_at_surface;
		elem_u_cav = md.elem_u_cav;
	}

	void MapPclToBgMesh::update(size_t tsk_num) noexcept
	{
		const auto& spva0 = stp.spvas[stp.next_spva_id()];
		const auto& spva1 = stp.spvas[stp.prev_spva_id()];
		// pcl_vars0
		pcl_index0 = spva0.pcl_index;
		pcl_n0 = spva0.pcl_n;
		pcl_density_f0 = spva0.pcl_density_f;
		pcl_v_s0 = spva0.pcl_v_s;
		pcl_u0 = spva0.pcl_u;
		pcl_stress0 = spva0.pcl_stress;
		//pcl_p0 = spva0.pcl_p;
		pcl_N0 = spva0.pcl_N;
		// pcl_vars1
		pcl_index1 = spva1.pcl_index;
		pcl_n1 = spva1.pcl_n;
		pcl_density_f1 = spva1.pcl_density_f;
		pcl_v_s1 = spva1.pcl_v_s;
		pcl_u1 = spva1.pcl_u;
		pcl_stress1 = spva1.pcl_stress;
		pcl_p1 = spva1.pcl_p;
		pcl_N1 = spva1.pcl_N;

		substep_index = stp.substep_index;
		dtime = stp.dtime;
		pcl_num = stp.valid_pcl_num;
		task_num = tsk_num;
	}

	void MapPclToBgMesh::work(size_t tsk_id, MapPclToBgMeshRes &res) const
	{
		size_t e_id;
		size_t p_id0 = Block_Low(tsk_id, task_num, pcl_num);
		e_id = pcl_in_elems[p_id0];
		while (p_id0 != SIZE_MAX && e_id == pcl_in_elems[--p_id0]);
		++p_id0;
		assert(p_id0 <= pcl_num);
		size_t p_id1 = Block_Low(tsk_id + 1, task_num, pcl_num);
		e_id = pcl_in_elems[p_id1];
		while (p_id1 != SIZE_MAX && e_id == pcl_in_elems[--p_id1]);
		++p_id1;
		assert(p_id1 <= pcl_num);

		double e_p_m_s = 0.0;
		double e_p_m_f = 0.0;
		double e_p_vol_s = 0.0;
		double e_p_vol_f = 0.0;
		double e_p_vol = 0.0;
		double e_p_pm = 0.0;
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
		double en1_p = 0.0;
		double en2_p = 0.0;
		double en3_p = 0.0;
		double en4_p = 0.0;
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
		double e_pf_bfx = 0.0;
		double e_pf_bfy = 0.0;
		double e_pf_bfz = 0.0;
		const double pf_bf_tmp = k / dyn_viscosity * dtime;
		e_id = pcl_in_elems[p_id0];
#ifdef _DEBUG
		assert(e_id < stp.elem_num);
#endif
		double p_N_m;
		for (size_t p_id = p_id0; p_id < p_id1; ++p_id)
		{
			// pcl index
			const size_t prev_p_id = prev_pcl_ids[p_id];
			assert(prev_p_id < stp.prev_valid_pcl_num);
			const size_t ori_p_id = pcl_index1[prev_p_id];
#ifdef _DEBUG
			assert(ori_p_id < stp.ori_pcl_num);
#endif
			pcl_index0[p_id] = ori_p_id;

			// m_s
			const double p_m_s = pcl_m_s[ori_p_id];
			e_p_m_s += p_m_s;
			const double p_vol_s = pcl_vol_s[ori_p_id];
			e_p_vol_s += p_vol_s;
			// vol
			const double p_n = pcl_n1[prev_p_id];
			const double p_vol = p_vol_s / (1.0 - p_n);
			pcl_vol[p_id] = p_vol;
			e_p_vol += p_vol;
			// vol_f
			const double p_vol_f = p_n * p_vol;
			e_p_vol_f += p_vol_f;
			// m_f
			const double p_den_f = pcl_density_f1[prev_p_id];
			const double p_m_f = p_den_f * p_vol_f;
			e_p_m_f += p_m_f;
			// m_total
			const double p_m = p_m_s + p_m_f;

			// map pcl stress and pore pressure to elements
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
			const double p_p = pcl_p1[prev_p_id];
			e_p += p_p * p_vol;
			//pcl_p0[p_id] = p_p;
			
			// e_p_pm
			e_p_pm += p_n / Kf0 * p_vol;

			// shape function to nodes
			const ShapeFunc& p_N1 = pcl_N1[prev_p_id];
			ShapeFunc& p_N0 = pcl_N0[p_id];
			p_N0.N1 = p_N1.N1;
			p_N0.N2 = p_N1.N2;
			p_N0.N3 = p_N1.N3;
			p_N0.N4 = p_N1.N4;
			// map velocity and pore pressure
			const Velocity& p_v_s1 = pcl_v_s1[prev_p_id];
			Velocity& p_v_s0 = pcl_v_s0[p_id];
			p_v_s0.vx = p_v_s1.vx;
			p_v_s0.vy = p_v_s1.vy;
			p_v_s0.vz = p_v_s1.vz;
			p_N_m = (p_N0.N1 > N_tol ? p_N0.N1 : N_tol) * p_m;
			en1_vm_s += p_N_m;
			en1_vmx_s += p_N_m * p_v_s0.vx;
			en1_vmy_s += p_N_m * p_v_s0.vy;
			en1_vmz_s += p_N_m * p_v_s0.vz;
			en1_p += p_N_m * p_p;
			p_N_m = (p_N0.N2 > N_tol ? p_N0.N2 : N_tol) * p_m;
			en2_vm_s += p_N_m;
			en2_vmx_s += p_N_m * p_v_s0.vx;
			en2_vmy_s += p_N_m * p_v_s0.vy;
			en2_vmz_s += p_N_m * p_v_s0.vz;
			en2_p += p_N_m * p_p;
			p_N_m = (p_N0.N3 > N_tol ? p_N0.N3 : N_tol) * p_m;
			en3_vm_s += p_N_m;
			en3_vmx_s += p_N_m * p_v_s0.vx;
			en3_vmy_s += p_N_m * p_v_s0.vy;
			en3_vmz_s += p_N_m * p_v_s0.vz;
			en3_p += p_N_m * p_p;
			p_N_m = (p_N0.N4 > N_tol ? p_N0.N4 : N_tol) * p_m;
			en4_vm_s += p_N_m;
			en4_vmx_s += p_N_m * p_v_s0.vx;
			en4_vmy_s += p_N_m * p_v_s0.vy;
			en4_vmz_s += p_N_m * p_v_s0.vz;
			en4_p += p_N_m * p_p;

			// displacement (for contact)
			const Displacement& p_u1 = pcl_u1[prev_p_id];
			Displacement& p_u0 = pcl_u0[p_id];
			p_u0.ux = p_u1.ux;
			p_u0.uy = p_u1.uy;
			p_u0.uz = p_u1.uz;

			// cal external load
			const Force& p_bf_s = pcl_bf_s[ori_p_id];
			const Force& p_bf_f = pcl_bf_f[ori_p_id];
			const double one_fourth_bfx = one_fourth * (p_bf_s.fx + p_bf_f.fx * p_m_f);
			const double one_fourth_bfy = one_fourth * (p_bf_s.fy + p_bf_f.fy * p_m_f);
			const double one_fourth_bfz = one_fourth * (p_bf_s.fz + p_bf_f.fz * p_m_f);
			const Force& p_t = pcl_t[ori_p_id];
			en1_fx_s += one_fourth_bfx + p_N0.N1 * p_t.fx;
			en1_fy_s += one_fourth_bfy + p_N0.N1 * p_t.fy;
			en1_fz_s += one_fourth_bfz + p_N0.N1 * p_t.fz;
			en2_fx_s += one_fourth_bfx + p_N0.N2 * p_t.fx;
			en2_fy_s += one_fourth_bfy + p_N0.N2 * p_t.fy;
			en2_fz_s += one_fourth_bfz + p_N0.N2 * p_t.fz;
			en3_fx_s += one_fourth_bfx + p_N0.N3 * p_t.fx;
			en3_fy_s += one_fourth_bfy + p_N0.N3 * p_t.fy;
			en3_fz_s += one_fourth_bfz + p_N0.N3 * p_t.fz;
			en4_fx_s += one_fourth_bfx + p_N0.N4 * p_t.fx;
			en4_fy_s += one_fourth_bfy + p_N0.N4 * p_t.fy;
			en4_fz_s += one_fourth_bfz + p_N0.N4 * p_t.fz;

			// external load of pore pressure
			e_pf_bfx += p_den_f * p_bf_f.fx * p_vol;
			e_pf_bfy += p_den_f * p_bf_f.fy * p_vol;
			e_pf_bfz += p_den_f * p_bf_f.fz * p_vol;

			if (e_id != pcl_in_elems[p_id + 1])
			{
				elem_has_pcls[e_id] = substep_index;
				elem_node_at_surface[e_id * 4] = 0;
				elem_node_at_surface[e_id * 4 + 1] = 0;
				elem_node_at_surface[e_id * 4 + 2] = 0;
				elem_node_at_surface[e_id * 4 + 3] = 0;
				
				elem_pcl_m[e_id] = e_p_m_s + e_p_m_f;
				elem_pcl_n[e_id] = 1.0 - e_p_vol_s / e_p_vol;
				elem_density_f[e_id] = e_p_m_f / e_p_vol_f;

				// map velocity
				ElemNodeVM& en1_v = elem_node_vm_s[e_id * 4];
				en1_v.vm = en1_vm_s;
				en1_v.vmx = en1_vmx_s;
				en1_v.vmy = en1_vmy_s;
				en1_v.vmz = en1_vmz_s;
				ElemNodeVM& en2_v = elem_node_vm_s[e_id * 4 + 1];
				en2_v.vm = en2_vm_s;
				en2_v.vmx = en2_vmx_s;
				en2_v.vmy = en2_vmy_s;
				en2_v.vmz = en2_vmz_s;
				ElemNodeVM& en3_v = elem_node_vm_s[e_id * 4 + 2];
				en3_v.vm = en3_vm_s;
				en3_v.vmx = en3_vmx_s;
				en3_v.vmy = en3_vmy_s;
				en3_v.vmz = en3_vmz_s;
				ElemNodeVM& en4_v = elem_node_vm_s[e_id * 4 + 3];
				en4_v.vm = en4_vm_s;
				en4_v.vmx = en4_vmx_s;
				en4_v.vmy = en4_vmy_s;
				en4_v.vmz = en4_vmz_s;
				// map pore pressure
				elem_node_p[e_id * 4] = en1_p;
				elem_node_p[e_id * 4 + 1] = en2_p;
				elem_node_p[e_id * 4 + 2] = en3_p;
				elem_node_p[e_id * 4 + 3] = en4_p;

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
				elem_pcl_vol[e_id] = e_p_vol;
				const DShapeFuncABC& e_dN = elem_dN_abc[e_id];
				// node 1
				Force& en1_f = elem_node_force[e_id * 4];
				en1_fx_s -= (e_dN.dN1_dx * (e_s11 - e_p) + e_dN.dN1_dy * e_s12 + e_dN.dN1_dz * e_s31) * e_p_vol;
				en1_f.fx = en1_fx_s;
				en1_fy_s -= (e_dN.dN1_dx * e_s12 + e_dN.dN1_dy * (e_s22 - e_p) + e_dN.dN1_dz * e_s23) * e_p_vol;
				en1_f.fy = en1_fy_s;
				en1_fz_s -= (e_dN.dN1_dx * e_s31 + e_dN.dN1_dy * e_s23 + e_dN.dN1_dz * (e_s33 - e_p)) * e_p_vol;
				en1_f.fz = en1_fz_s;
				// node 2
				Force& en2_f = elem_node_force[e_id * 4 + 1];
				en2_fx_s -= (e_dN.dN2_dx * (e_s11 - e_p) + e_dN.dN2_dy * e_s12 + e_dN.dN2_dz * e_s31) * e_p_vol;
				en2_f.fx = en2_fx_s;
				en2_fy_s -= (e_dN.dN2_dx * e_s12 + e_dN.dN2_dy * (e_s22 - e_p) + e_dN.dN2_dz * e_s23) * e_p_vol;
				en2_f.fy = en2_fy_s;
				en2_fz_s -= (e_dN.dN2_dx * e_s31 + e_dN.dN2_dy * e_s23 + e_dN.dN2_dz * (e_s33 - e_p)) * e_p_vol;
				en2_f.fz = en2_fz_s;
				// node 3
				Force& en3_f = elem_node_force[e_id * 4 + 2];
				en3_fx_s -= (e_dN.dN3_dx * (e_s11 - e_p) + e_dN.dN3_dy * e_s12 + e_dN.dN3_dz * e_s31) * e_p_vol;
				en3_f.fx = en3_fx_s;
				en3_fy_s -= (e_dN.dN3_dx * e_s12 + e_dN.dN3_dy * (e_s22 - e_p) + e_dN.dN3_dz * e_s23) * e_p_vol;
				en3_f.fy = en3_fy_s;
				en3_fz_s -= (e_dN.dN3_dx * e_s31 + e_dN.dN3_dy * e_s23 + e_dN.dN3_dz * (e_s33 - e_p)) * e_p_vol;
				en3_f.fz = en3_fz_s;
				// node 4
				Force& en4_f = elem_node_force[e_id * 4 + 3];
				en4_fx_s -= (e_dN.dN4_dx * (e_s11 - e_p) + e_dN.dN4_dy * e_s12 + e_dN.dN4_dz * e_s31) * e_p_vol;
				en4_f.fx = en4_fx_s;
				en4_fy_s -= (e_dN.dN4_dx * e_s12 + e_dN.dN4_dy * (e_s22 - e_p) + e_dN.dN4_dz * e_s23) * e_p_vol;
				en4_f.fy = en4_fy_s;
				en4_fz_s -= (e_dN.dN4_dx * e_s31 + e_dN.dN4_dy * e_s23 + e_dN.dN4_dz * (e_s33 - e_p)) * e_p_vol;
				en4_f.fz = en4_fz_s;

				// cavitation
				double Kf_ratio = 1.0;
				if (m_cav != 0.0) // consider cavitation
				{
					if (e_p < 0.0)
					{
						const double tmp = e_p / elem_u_cav[e_id] + u_cav_off;
						Kf_ratio = tmp < u_div_u_cav_lim ? (1.0 / (1.0 + pow(tmp, m_cav))) : (1.0 / max_Kf_ratio_divider);
					}
				}
				elem_pcl_pm[e_id] = e_p_pm / Kf_ratio;

				elem_node_p_force[e_id * 4    ] = pf_bf_tmp * (e_dN.dN1_dx * e_pf_bfx + e_dN.dN1_dy * e_pf_bfy + e_dN.dN1_dz * e_pf_bfz);
				elem_node_p_force[e_id * 4 + 1] = pf_bf_tmp * (e_dN.dN2_dx * e_pf_bfx + e_dN.dN2_dy * e_pf_bfy + e_dN.dN2_dz * e_pf_bfz);
				elem_node_p_force[e_id * 4 + 2] = pf_bf_tmp * (e_dN.dN3_dx * e_pf_bfx + e_dN.dN3_dy * e_pf_bfy + e_dN.dN3_dz * e_pf_bfz);
				elem_node_p_force[e_id * 4 + 3] = pf_bf_tmp * (e_dN.dN4_dx * e_pf_bfx + e_dN.dN4_dy * e_pf_bfy + e_dN.dN4_dz * e_pf_bfz);

				e_id = pcl_in_elems[p_id + 1];
#ifdef _DEBUG
				assert(e_id < stp.elem_num || e_id == SIZE_MAX);
#endif

				e_p_m_s = 0.0;
				e_p_m_f = 0.0;
				e_p_vol_s = 0.0;
				e_p_vol_f = 0.0;
				e_p_vol = 0.0;
				e_p_pm = 0.0;
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
				en1_p = 0.0;
				en2_p = 0.0;
				en3_p = 0.0;
				en4_p = 0.0;
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
				e_pf_bfx = 0.0;
				e_pf_bfy = 0.0;
				e_pf_bfz = 0.0;
			}
		}

		// contact force calculation
		ContactRigidBody& crb = stp.cont_rigid_body;
		if (crb.has_rigid_mesh())
			crb.apply_t3d_rigid_object(p_id0, p_id1, res.react_force);
	}

	void FindSoilSurface::init() noexcept
	{
		Model_T3D_CHM_up_mt& md = *stp.pmodel;
		
		elem_adj_elems = md.elem_adj_elems;
		elem_has_pcls = md.elem_has_pcls;
		elem_node_at_surface = md.elem_node_at_surface;
		//
		elem_ids = stp.elem_ids;
	}

	void FindSoilSurface::update(size_t tsk_num) noexcept
	{
		substep_index = stp.substep_index;
		elem_num = stp.valid_elem_num;
		task_num = tsk_num;
	}

	void FindSoilSurface::work(size_t tsk_id) const
	{
		const size_t ve_id0 = Block_Low(tsk_id, task_num, elem_num);
		const size_t ve_id1 = Block_Low(tsk_id + 1, task_num, elem_num);
		for (size_t ve_id = ve_id0; ve_id < ve_id1; ++ve_id)
		{
			const size_t e_id = elem_ids[ve_id];
			
			uint16_t en1_at_surface = 0;
			uint16_t en2_at_surface = 0;
			uint16_t en3_at_surface = 0;
			uint16_t en4_at_surface = 0;
			const AdjElemIndex& adj_elem = elem_adj_elems[e_id];
			if (adj_elem.adj_e1 != SIZE_MAX && elem_has_pcls[adj_elem.adj_e1] != substep_index) // adjacent element is empty
			{
				en1_at_surface = 1;
				en2_at_surface = 1;
				en3_at_surface = 1;
			}
			if (adj_elem.adj_e2 != SIZE_MAX && elem_has_pcls[adj_elem.adj_e2] != substep_index) // adjacent elem is empty
			{
				en1_at_surface = 1;
				en4_at_surface = 1;
				en2_at_surface = 1;
			}
			if (adj_elem.adj_e3 != SIZE_MAX && elem_has_pcls[adj_elem.adj_e3] != substep_index) // adjacent elem is empty
			{
				en2_at_surface = 1;
				en4_at_surface = 1;
				en3_at_surface = 1;
			}
			if (adj_elem.adj_e4 != SIZE_MAX && elem_has_pcls[adj_elem.adj_e4] != substep_index) // adjacent elem is empty
			{
				en1_at_surface = 1;
				en3_at_surface = 1;
				en4_at_surface = 1;
			}
			elem_node_at_surface[4 * e_id + 0] |= en1_at_surface;
			elem_node_at_surface[4 * e_id + 1] |= en2_at_surface;
			elem_node_at_surface[4 * e_id + 2] |= en3_at_surface;
			elem_node_at_surface[4 * e_id + 3] |= en4_at_surface;
		}

	}

	void UpdateAccelerationAndVelocity::init() noexcept
	{
		Model_T3D_CHM_up_mt& md = *stp.pmodel;

		elem_pcl_m = md.elem_pcl_m;
		elem_node_force = md.elem_node_force;
		elem_node_vm_s = md.elem_node_vm_s;
		elem_node_p = md.elem_node_p;
		elem_node_at_surface = md.elem_node_at_surface;

		node_has_vbc = md.node_has_vbc;
		node_vbc_vec_s = md.node_vbc_vec_s;
		node_am = md.node_am;
		node_a_s = md.node_a_s;
		node_v_s = md.node_v_s;
		node_p = md.node_p;
		node_at_surface = md.node_at_surface;

		// node ranges
		node_ids = stp.node_ids;
		node_elem_offs = stp.node_elem_offs;
	}

	void UpdateAccelerationAndVelocity::update(size_t tsk_num) noexcept
	{
		dtime = stp.dtime;
		four_elem_num = stp.valid_elem_num * 4;
		task_num = tsk_num;
	}

	void UpdateAccelerationAndVelocity::work(size_t tsk_id) const
	{
		size_t n_id;
		size_t ve_id0 = Block_Low(tsk_id, task_num, four_elem_num);
		n_id = node_ids[ve_id0];
		while (ve_id0 != SIZE_MAX && n_id == node_ids[--ve_id0]);
		++ve_id0;
		assert(ve_id0 <= four_elem_num);
		size_t ve_id1 = Block_Low(tsk_id + 1, task_num, four_elem_num);
		n_id = node_ids[ve_id1];
		while (ve_id1 != SIZE_MAX && n_id == node_ids[--ve_id1]);
		++ve_id1;
		assert(ve_id1 <= four_elem_num);

		size_t bc_mask, ne_id;
		double vbc_len;
		double n_am = 0.0;
		double n_fx = 0.0;
		double n_fy = 0.0;
		double n_fz = 0.0;
		double n_vm = 0.0;
		double n_vmx = 0.0;
		double n_vmy = 0.0;
		double n_vmz = 0.0;
		double n_pm = 0.0;
		uint16_t n_at_surface = 0;
		n_id = node_ids[ve_id0];
#ifdef _DEBUG
		assert(n_id < stp.node_num);
#endif
		for (size_t ve_id = ve_id0; ve_id < ve_id1; ++ve_id)
		{
			ne_id = node_elem_offs[ve_id];
#ifdef _DEBUG
			assert(ne_id < stp.elem_num * 4);
#endif
			n_am += elem_pcl_m[ne_id / 4];
			const Force& nf = elem_node_force[ne_id];
			n_fx += nf.fx;
			n_fy += nf.fy;
			n_fz += nf.fz;
			const ElemNodeVM& nvm = elem_node_vm_s[ne_id];
			n_vm += nvm.vm;
			n_vmx += nvm.vmx;
			n_vmy += nvm.vmy;
			n_vmz += nvm.vmz;
			n_pm += elem_node_p[ne_id];

			n_at_surface |= elem_node_at_surface[ne_id];

			if (n_id != node_ids[ve_id + 1])
			{
				node_at_surface[n_id] = n_at_surface;

				Acceleration& n_a = node_a_s[n_id];
				n_am *= one_fourth;
				node_am[n_id] = n_am;
				n_a.ax = n_fx / n_am;
				n_a.ay = n_fy / n_am;
				n_a.az = n_fz / n_am;
				Velocity& n_v = node_v_s[n_id];
				n_v.vx = n_vmx / n_vm + n_a.ax * dtime;
				n_v.vy = n_vmy / n_vm + n_a.ay * dtime;
				n_v.vz = n_vmz / n_vm + n_a.az * dtime;
				const NodeVBCVec& n_vbc_vec = node_vbc_vec_s[n_id];
				vbc_len = n_a.ax * n_vbc_vec.x + n_a.ay * n_vbc_vec.y + n_a.az * n_vbc_vec.z;
				n_a.ax -= vbc_len * n_vbc_vec.x;
				n_a.ay -= vbc_len * n_vbc_vec.y;
				n_a.az -= vbc_len * n_vbc_vec.z;
				vbc_len = n_v.vx * n_vbc_vec.x + n_v.vy * n_vbc_vec.y + n_v.vz * n_vbc_vec.z;
				n_v.vx -= vbc_len * n_vbc_vec.x;
				n_v.vy -= vbc_len * n_vbc_vec.y;
				n_v.vz -= vbc_len * n_vbc_vec.z;
				const NodeHasVBC& n_has_vbc = node_has_vbc[n_id];
				bc_mask = size_t(n_has_vbc.has_vx_bc) + SIZE_MAX;
				n_a.iax &= bc_mask;
				n_v.ivx &= bc_mask;
				bc_mask = size_t(n_has_vbc.has_vy_bc) + SIZE_MAX;
				n_a.iay &= bc_mask;
				n_v.ivy &= bc_mask;
				bc_mask = size_t(n_has_vbc.has_vz_bc) + SIZE_MAX;
				n_a.iaz &= bc_mask;
				n_v.ivz &= bc_mask;

				double& n_p = node_p[n_id];
				n_p = n_pm / n_vm;
				bc_mask = size_t(n_has_vbc.is_drained) + SIZE_MAX;
				*(reinterpret_cast<size_t *>(&n_p)) &= bc_mask;
				if (n_at_surface == 1) // at surface but in contact
					n_p = 0.0;

				n_id = node_ids[ve_id + 1];
#ifdef _DEBUG
				assert(n_id < stp.node_num || n_id == SIZE_MAX);
#endif

				n_am = 0.0;
				n_fx = 0.0;
				n_fy = 0.0;
				n_fz = 0.0;
				n_vm = 0.0;
				n_vmx = 0.0;
				n_vmy = 0.0;
				n_vmz = 0.0;
				n_pm = 0.0;
				n_at_surface = 0;
			}
		}
	}

	void CalElemDeAndMapToNode::init() noexcept
	{
		Model_T3D_CHM_up_mt& md = *stp.pmodel;

		k = md.k;
		dyn_viscosity = md.dyn_viscosity;

		elem_node_id = md.elem_node_id;
		elem_dN_abc = md.elem_dN_abc;
		elem_pcl_m = md.elem_pcl_m;
		elem_pcl_vol = md.elem_pcl_vol;
		elem_density_f = md.elem_density_f;
		elem_de = md.elem_de;
		elem_m_de_vol_s = md.elem_m_de_vol_s;

		elem_node_p_force = md.elem_node_p_force;
		elem_node_at_surface = md.elem_node_at_surface;

		node_a_s = md.node_a_s;
		node_v_s = md.node_v_s;
		node_p = md.node_p;

		//
		elem_ids = stp.elem_ids;
	}

	void CalElemDeAndMapToNode::update(size_t tsk_num) noexcept
	{
		dtime = stp.dtime;
		elem_num = stp.valid_elem_num;
		task_num = tsk_num;
	}

	void CalElemDeAndMapToNode::work(size_t tsk_id) const
	{
		const size_t ve_id0 = Block_Low(tsk_id, task_num, elem_num);
		const size_t ve_id1 = Block_Low(tsk_id + 1, task_num, elem_num);
		const double pf_coef = k / dyn_viscosity * dtime;
		for (size_t ve_id = ve_id0; ve_id < ve_id1; ++ve_id)
		{
			const size_t e_id = elem_ids[ve_id];
#ifdef _DEBUG
			assert(e_id < stp.elem_num);
#endif
			const ElemNodeIndex& eni = elem_node_id[e_id];
			const Velocity& n_v1 = node_v_s[eni.n1];
			const Velocity& n_v2 = node_v_s[eni.n2];
			const Velocity& n_v3 = node_v_s[eni.n3];
			const Velocity& n_v4 = node_v_s[eni.n4];
			const DShapeFuncABC& e_dN = elem_dN_abc[e_id];
			StrainInc& e_de = elem_de[e_id];
			e_de.de11 = (e_dN.dN1_dx * n_v1.vx + e_dN.dN2_dx * n_v2.vx + e_dN.dN3_dx * n_v3.vx + e_dN.dN4_dx * n_v4.vx) * dtime;
			e_de.de22 = (e_dN.dN1_dy * n_v1.vy + e_dN.dN2_dy * n_v2.vy + e_dN.dN3_dy * n_v3.vy + e_dN.dN4_dy * n_v4.vy) * dtime;
			e_de.de33 = (e_dN.dN1_dz * n_v1.vz + e_dN.dN2_dz * n_v2.vz + e_dN.dN3_dz * n_v3.vz + e_dN.dN4_dz * n_v4.vz) * dtime;
			e_de.de12 = (e_dN.dN1_dx * n_v1.vy + e_dN.dN2_dx * n_v2.vy + e_dN.dN3_dx * n_v3.vy + e_dN.dN4_dx * n_v4.vy
					   + e_dN.dN1_dy * n_v1.vx + e_dN.dN2_dy * n_v2.vx + e_dN.dN3_dy * n_v3.vx + e_dN.dN4_dy * n_v4.vx) * dtime * 0.5;
			e_de.de23 = (e_dN.dN1_dy * n_v1.vz + e_dN.dN2_dy * n_v2.vz + e_dN.dN3_dy * n_v3.vz + e_dN.dN4_dy * n_v4.vz
					   + e_dN.dN1_dz * n_v1.vy + e_dN.dN2_dz * n_v2.vy + e_dN.dN3_dz * n_v3.vy + e_dN.dN4_dz * n_v4.vy) * dtime * 0.5;
			e_de.de31 = (e_dN.dN1_dz * n_v1.vx + e_dN.dN2_dz * n_v2.vx + e_dN.dN3_dz * n_v3.vx + e_dN.dN4_dz * n_v4.vx
					   + e_dN.dN1_dx * n_v1.vz + e_dN.dN2_dx * n_v2.vz + e_dN.dN3_dx * n_v3.vz + e_dN.dN4_dx * n_v4.vz) * dtime * 0.5;
			const double e_de_vol = e_de.de11 + e_de.de22 + e_de.de33;

			const double e_pcl_vol = elem_pcl_vol[e_id];
			double n1_pf = -one_fourth * e_de_vol * e_pcl_vol;
			double n2_pf = n1_pf;
			double n3_pf = n1_pf;
			double n4_pf = n1_pf;

			// p_force
			const double n1_p = node_p[eni.n1];
			const double n2_p = node_p[eni.n2];
			const double n3_p = node_p[eni.n3];
			const double n4_p = node_p[eni.n4];
			const double dp_dx = e_dN.dN1_dx * n1_p + e_dN.dN2_dx * n2_p + e_dN.dN3_dx * n3_p + e_dN.dN4_dx * n4_p;
			const double dp_dy = e_dN.dN1_dy * n1_p + e_dN.dN2_dy * n2_p + e_dN.dN3_dy * n3_p + e_dN.dN4_dy * n4_p;
			const double dp_dz = e_dN.dN1_dz * n1_p + e_dN.dN2_dz * n2_p + e_dN.dN3_dz * n3_p + e_dN.dN4_dz * n4_p;
			const double pf_tmp = pf_coef * e_pcl_vol;
			n1_pf -= pf_tmp * (e_dN.dN1_dx * dp_dx + e_dN.dN1_dy * dp_dy + e_dN.dN1_dz * dp_dz);
			n2_pf -= pf_tmp * (e_dN.dN2_dx * dp_dx + e_dN.dN2_dy * dp_dy + e_dN.dN2_dz * dp_dz);
			n3_pf -= pf_tmp * (e_dN.dN3_dx * dp_dx + e_dN.dN3_dy * dp_dy + e_dN.dN3_dz * dp_dz);
			n4_pf -= pf_tmp * (e_dN.dN4_dx * dp_dx + e_dN.dN4_dy * dp_dy + e_dN.dN4_dz * dp_dz);

			const Acceleration &n_a1 = node_a_s[eni.n1];
			const Acceleration& n_a2 = node_a_s[eni.n2];
			const Acceleration& n_a3 = node_a_s[eni.n3];
			const Acceleration& n_a4 = node_a_s[eni.n4];
			const double e_ax_s = (n_a1.ax + n_a2.ax + n_a3.ax + n_a4.ax) * one_fourth;
			const double e_ay_s = (n_a1.ay + n_a2.ay + n_a3.ay + n_a4.ay) * one_fourth;
			const double e_az_s = (n_a1.az + n_a2.az + n_a3.az + n_a4.az) * one_fourth;
			const double pf_a_tmp = pf_tmp * elem_density_f[e_id];
			// this causes instability
			//n1_pf -= pf_a_tmp * (e_dN.dN1_dx * e_ax_s + e_dN.dN1_dy * e_ay_s + e_dN.dN1_dz * e_az_s);
			//n2_pf -= pf_a_tmp * (e_dN.dN2_dx * e_ax_s + e_dN.dN2_dy * e_ay_s + e_dN.dN2_dz * e_az_s);
			//n3_pf -= pf_a_tmp * (e_dN.dN3_dx * e_ax_s + e_dN.dN3_dy * e_ay_s + e_dN.dN3_dz * e_az_s);
			//n4_pf -= pf_a_tmp * (e_dN.dN4_dx * e_ax_s + e_dN.dN4_dy * e_ay_s + e_dN.dN4_dz * e_az_s);
			
			elem_node_p_force[e_id * 4] += n1_pf;
			elem_node_p_force[e_id * 4 + 1] += n2_pf;
			elem_node_p_force[e_id * 4 + 2] += n3_pf;
			elem_node_p_force[e_id * 4 + 3] += n4_pf;

			// strain enhancement
			//elem_m_de_vol[e_id] = elem_pcl_m[e_id] * e_de_vol;
			//e_de_vol *= one_third;
			//e_de.de11 -= e_de_vol;
			//e_de.de22 -= e_de_vol;
			//e_de.de33 -= e_de_vol;
		}
	}
	
	void CalNodeDe::init() noexcept
	{
		Model_T3D_CHM_up_mt& md = *stp.pmodel;

		elem_pcl_pm = md.elem_pcl_pm;
		elem_node_p_force = md.elem_node_p_force;
		node_has_vbc = md.node_has_vbc;
		node_dp = md.node_dp;
		node_at_surface = md.node_at_surface;

		elem_m_de_vol_s = md.elem_m_de_vol_s;
		node_am = md.node_am;
		node_de_vol_s = md.node_de_vol_s;

		// node ranges
		node_ids = stp.node_ids;
		node_elem_offs = stp.node_elem_offs;
	}

	void CalNodeDe::update(size_t tsk_num) noexcept
	{
		four_elem_num = stp.valid_elem_num * 4;
		task_num = tsk_num;
	}

	void CalNodeDe::work(size_t tsk_id) const
	{
		size_t n_id;
		size_t ve_id0 = Block_Low(tsk_id, task_num, four_elem_num);
		n_id = node_ids[ve_id0];
		while (ve_id0 != SIZE_MAX && n_id == node_ids[--ve_id0]);
		++ve_id0;
		assert(ve_id0 <= four_elem_num);
		size_t ve_id1 = Block_Low(tsk_id + 1, task_num, four_elem_num);
		n_id = node_ids[ve_id1];
		while (ve_id1 != SIZE_MAX && n_id == node_ids[--ve_id1]);
		++ve_id1;
		assert(ve_id1 <= four_elem_num);

		size_t ne_id;
		//double n_am_de_vol_s = 0.0;
		double n_pm = 0.0;
		double n_pf = 0.0;
		n_id = node_ids[ve_id0];
#ifdef _DEBUG
		assert(n_id < stp.node_num || n_id == SIZE_MAX);
#endif
		for (size_t ve_id = ve_id0; ve_id < ve_id1; ++ve_id)
		{
			ne_id = node_elem_offs[ve_id];
#ifdef _DEBUG
			assert(ne_id < stp.elem_num * 4);
#endif
			n_pm += elem_pcl_pm[ne_id / 4];
			n_pf += elem_node_p_force[ne_id];
			// strain enhancement
			//n_am_de_vol_s += elem_m_de_vol_s[e_id];
			if (n_id != node_ids[ve_id + 1])
			{
				n_pm *= one_fourth;
				node_dp[n_id] = n_pf / n_pm;
				// the node is at surface without contact
				if (node_has_vbc[n_id].is_drained || node_at_surface[n_id] == 1)
					node_dp[n_id] = 0.0;

				// strain enhancement
				//node_de_vol[n_id] = n_am_de_vol * one_fourth / node_am[n_id];

				n_id = node_ids[ve_id + 1];
#ifdef _DEBUG
				assert(n_id < stp.node_num || n_id == SIZE_MAX);
#endif
				// strain enhancement
				n_pm = 0.0;
				n_pf = 0.0;
				//n_am_de_vol_s = 0.0;
			}
		}
	}

	void MapBgMeshToPcl::init() noexcept
	{
		Model_T3D_CHM_up_mt& md = *stp.pmodel;

		Kf0 = md.Kf;
		m_cav = md.m_cav;
		f_cav_end = md.f_cav_end;
		u_cav_off = pow(1.0 / f_cav_end - 1.0, 1.0 / m_cav) - 1.0;
		u_div_u_cav_lim = pow(Step_T3D_CHM_up_TBB_Task::max_Kf_ratio_divider - 1.0, 1.0 / m_cav);
		
		pcl_pos = md.pcl_pos;
		pcl_mat_model = md.pcl_mat_model;
		pcl_is_cavitated = md.pcl_is_cavitated;
		//
		elem_node_id = md.elem_node_id;
		node_a_s = md.node_a_s;
		node_v_s = md.node_v_s;
		node_p = md.node_p;
		node_dp = md.node_dp;
		elem_p = md.elem_p;
		elem_pcl_n = md.elem_pcl_n;
		elem_density_f = md.elem_density_f;
		elem_de = md.elem_de;
		elem_u_cav = md.elem_u_cav;
		node_de_vol_s = md.node_de_vol_s;
		//
		auto& pcl_sort = stp.pcl_sort;
		pcl_in_elems = pcl_sort.out_pcl_in_elems();
		prev_pcl_ids = pcl_sort.out_prev_pcl_ids();
		in_pcl_in_elems = pcl_sort.in_pcl_in_elems();
		in_prev_pcl_ids = pcl_sort.in_prev_pcl_ids();
	}

	void MapBgMeshToPcl::update(size_t tsk_num) noexcept
	{
		const auto& spva0 = stp.spvas[stp.next_spva_id()];
		pcl_index0 = spva0.pcl_index;
		pcl_n0 = spva0.pcl_n;
		pcl_density_f0 = spva0.pcl_density_f;
		pcl_v_s0 = spva0.pcl_v_s;
		pcl_u0 = spva0.pcl_u;
		pcl_N0 = spva0.pcl_N;
		pcl_stress0 = spva0.pcl_stress;
		pcl_p0 = spva0.pcl_p;
		pcl_strain0 = spva0.pcl_strain;
		pcl_estrain0 = spva0.pcl_estrain;
		pcl_pstrain0 = spva0.pcl_pstrain;
		const auto& spva1 = stp.spvas[stp.prev_spva_id()];
		pcl_strain1 = spva1.pcl_strain;
		pcl_estrain1 = spva1.pcl_estrain;
		pcl_pstrain1 = spva1.pcl_pstrain;
		
		dtime = stp.dtime;
		pcl_num = stp.prev_valid_pcl_num;
		task_num = tsk_num;
	}

	void MapBgMeshToPcl::work(size_t tsk_id, MapBgMeshToPclRes &res) const
	{
		size_t e_id;
		size_t p_id0 = Block_Low(tsk_id, task_num, pcl_num);
		e_id = pcl_in_elems[p_id0];
		while (p_id0 != SIZE_MAX && e_id == pcl_in_elems[--p_id0]);
		++p_id0;
		assert(p_id0 <= stp.prev_valid_pcl_num);
		size_t p_id1 = Block_Low(tsk_id + 1, task_num, pcl_num);
		e_id = pcl_in_elems[p_id1];
		while (p_id1 != SIZE_MAX && e_id == pcl_in_elems[--p_id1]);
		++p_id1;
		assert(p_id1 <= stp.prev_valid_pcl_num);

		const Acceleration *pn_a1, *pn_a2, *pn_a3, *pn_a4;
		const Velocity *pn_v1, *pn_v2, *pn_v3, *pn_v4;
		StrainInc dstrain;
		double e_p, e_n, e_density_f, Kf_ratio;
		size_t valid_pcl_num = 0;
		e_id = SIZE_MAX;
		Model_T3D_CHM_up_mt& md = *stp.pmodel;
		for (size_t p_id = p_id0; p_id < p_id1; ++p_id)
		{
			if (e_id != pcl_in_elems[p_id])
			{
				e_id = pcl_in_elems[p_id];
#ifdef _DEBUG
				assert(e_id < stp.elem_num);
#endif
				const ElemNodeIndex& eni = elem_node_id[e_id];
				pn_a1 = node_a_s + eni.n1;
				pn_a2 = node_a_s + eni.n2;
				pn_a3 = node_a_s + eni.n3;
				pn_a4 = node_a_s + eni.n4;
				pn_v1 = node_v_s + eni.n1;
				pn_v2 = node_v_s + eni.n2;
				pn_v3 = node_v_s + eni.n3;
				pn_v4 = node_v_s + eni.n4;

				//double e_de_vol = (node_de_vol_s[eni.n1]
				//	+ node_de_vol_s[eni.n2] + node_de_vol_s[eni.n3]
				//	+ node_de_vol_s[eni.n4]) * one_fourth;
				//e_density = elem_density[e_id] / (1.0 + e_de_vol);
				
				//pe_de = elem_de + e_id;
				//e_de_vol *= one_third;
				//dstrain.de11 = pe_de->de11 + e_de_vol;
				//dstrain.de22 = pe_de->de22 + e_de_vol;
				//dstrain.de33 = pe_de->de33 + e_de_vol;
				//dstrain.de12 = pe_de->de12;
				//dstrain.de23 = pe_de->de23;
				//dstrain.de31 = pe_de->de31;

				StrainInc &e_de = elem_de[e_id];
				dstrain.de11 = e_de.de11;
				dstrain.de22 = e_de.de22;
				dstrain.de33 = e_de.de33;
				dstrain.de12 = e_de.de12;
				dstrain.de23 = e_de.de23;
				dstrain.de31 = e_de.de31;

				const double de_vol_s = e_de.de11 + e_de.de22 + e_de.de33;
				e_n = (elem_pcl_n[e_id] + de_vol_s) / (1.0 + de_vol_s);

				const double e_dp = (node_dp[eni.n1] + node_dp[eni.n2]
					+ node_dp[eni.n3] + node_dp[eni.n4]) * one_fourth;
				e_p = elem_p[e_id] + e_dp;
				// cavitation
				Kf_ratio = 1.0;
				if (m_cav != 0.0) // consider cavitation
				{
					if (e_p < 0.0)
					{
						const double tmp = e_p / elem_u_cav[e_id] + u_cav_off;
						Kf_ratio = tmp < u_div_u_cav_lim ? (1.0 / (1.0 + pow(tmp, m_cav))) : (1.0 / max_Kf_ratio_divider);
					}
				}
				const double de_vol_f = -e_dp / (Kf_ratio * Kf0);
				e_density_f = elem_density_f[e_id] / (1.0 + de_vol_f);
			}

			// update velocity
			ShapeFunc& p_N = pcl_N0[p_id];
			Velocity& p_v0 = pcl_v_s0[p_id];
			p_v0.vx += (p_N.N1 * pn_a1->ax + p_N.N2 * pn_a2->ax + p_N.N3 * pn_a3->ax + p_N.N4 * pn_a4->ax) * dtime;
			p_v0.vy += (p_N.N1 * pn_a1->ay + p_N.N2 * pn_a2->ay + p_N.N3 * pn_a3->ay + p_N.N4 * pn_a4->ay) * dtime;
			p_v0.vz += (p_N.N1 * pn_a1->az + p_N.N2 * pn_a2->az + p_N.N3 * pn_a3->az + p_N.N4 * pn_a4->az) * dtime;

			// update displacement
			Displacement& p_d0 = pcl_u0[p_id];
			p_d0.ux += (p_N.N1 * pn_v1->vx + p_N.N2 * pn_v2->vx + p_N.N3 * pn_v3->vx + p_N.N4 * pn_v4->vx) * dtime;
			p_d0.uy += (p_N.N1 * pn_v1->vy + p_N.N2 * pn_v2->vy + p_N.N3 * pn_v3->vy + p_N.N4 * pn_v4->vy) * dtime;
			p_d0.uz += (p_N.N1 * pn_v1->vz + p_N.N2 * pn_v2->vz + p_N.N3 * pn_v3->vz + p_N.N4 * pn_v4->vz) * dtime;

			// update location (in which element)
			const size_t ori_p_id = pcl_index0[p_id];
#ifdef _DEBUG
			assert(ori_p_id < stp.ori_pcl_num);
#endif
			const Position& p_p = pcl_pos[ori_p_id];
			const double p_x = p_p.x + p_d0.ux;
			const double p_y = p_p.y + p_d0.uy;
			const double p_z = p_p.z + p_d0.uz;
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
			in_pcl_in_elems[p_id] = p_e_id;
			in_prev_pcl_ids[p_id] = p_id;
			if (p_e_id != SIZE_MAX)
				++valid_pcl_num;
#ifdef _DEBUG
			assert(p_e_id < stp.elem_num || p_e_id == SIZE_MAX);
#endif

			// update density
			pcl_n0[p_id] = e_n;
			pcl_p0[p_id] = e_p;
			pcl_density_f0[p_id] = e_density_f;
			pcl_is_cavitated[p_id] = Kf_ratio;

			// update stress
			MatModel::MaterialModel& pcl_mm = *pcl_mat_model[ori_p_id];
			pcl_mm.integrate(dstrain.de);
			const double* dstress = pcl_mm.get_dstress();
			Stress& p_s = pcl_stress0[p_id];
			p_s.s11 += dstress[0];
			p_s.s22 += dstress[1];
			p_s.s33 += dstress[2];
			p_s.s12 += dstress[3];
			p_s.s23 += dstress[4];
			p_s.s31 += dstress[5];

			const size_t prev_p_id = prev_pcl_ids[p_id];
#ifdef _DEBUG
			assert(prev_p_id < stp.prev_valid_pcl_num_tmp);
#endif
			const Strain& p_e1 = pcl_strain1[prev_p_id];
			Strain& p_e0 = pcl_strain0[p_id];
			p_e0.e11 = p_e1.e11 + dstrain.de11;
			p_e0.e22 = p_e1.e22 + dstrain.de22;
			p_e0.e33 = p_e1.e33 + dstrain.de33;
			p_e0.e12 = p_e1.e12 + dstrain.de12;
			p_e0.e23 = p_e1.e23 + dstrain.de23;
			p_e0.e31 = p_e1.e31 + dstrain.de31;

			const double* estrain = pcl_mm.get_dstrain_e();
			const Strain& p_ee1 = pcl_estrain1[prev_p_id];
			Strain& p_ee0 = pcl_estrain0[p_id];
			p_ee0.e11 = p_ee1.e11 + estrain[0];
			p_ee0.e22 = p_ee1.e22 + estrain[1];
			p_ee0.e33 = p_ee1.e33 + estrain[2];
			p_ee0.e12 = p_ee1.e12 + estrain[3];
			p_ee0.e23 = p_ee1.e23 + estrain[4];
			p_ee0.e31 = p_ee1.e31 + estrain[5];

			const double* pstrain = pcl_mm.get_dstrain_p();
			const Strain& p_pe1 = pcl_pstrain1[prev_p_id];
			Strain& p_pe0 = pcl_pstrain0[p_id];
			p_pe0.e11 = p_pe1.e11 + pstrain[0];
			p_pe0.e22 = p_pe1.e22 + pstrain[1];
			p_pe0.e33 = p_pe1.e33 + pstrain[2];
			p_pe0.e12 = p_pe1.e12 + pstrain[3];
			p_pe0.e23 = p_pe1.e23 + pstrain[4];
			p_pe0.e31 = p_pe1.e31 + pstrain[5];
		}

		res.pcl_num = valid_pcl_num;
	}

// ================== contact funct ================
	void ContactRigidBody::init(
		double max_pcl_vol,
		bool is_first_step) noexcept
	{
		Model_T3D_CHM_up_mt& md = *stp.pmodel;
		pcm = md.get_contact_model();
		prmesh = nullptr;
		if (md.has_t3d_rigid_mesh())
		{
			prmesh = &md.get_t3d_rigid_mesh();
			// set max dist for efficiency
			prmesh->init_max_dist(0.5 * pow(max_pcl_vol, one_third) * 4.0);
			if (is_first_step)
				prmesh->reset_cont_force();
		}
		//
		pcl_pos = md.pcl_pos;
		pcl_vol = md.pcl_vol;
		elem_node_force = md.elem_node_force;
		elem_node_at_surface = md.elem_node_at_surface;
		//
		pcl_in_elems = stp.pcl_in_elems;
	}

	void ContactRigidBody::update() noexcept
	{
		const auto& spva0 = stp.spvas[stp.next_spva_id()];
		pcl_index = spva0.pcl_index;
		pcl_u = spva0.pcl_u;
		pcl_N = spva0.pcl_N;
		substep_index = stp.substep_index;
	}

	void ContactRigidBody::apply_t3d_rigid_object(
		size_t p_id0,
		size_t p_id1,
		Force3D& rc_cf)
		const noexcept
	{
		double dist;
		Vector3D lnorm, gnorm;
		Force lcont_f, gcont_f;
		Point3D cur_cont_pos;
		ParticleVariablesGetter pv_getter;
		Force3D rc_cf_tmp;
		rc_cf_tmp.reset();
		for (size_t p_id = p_id0; p_id < p_id1; ++p_id)
		{
			const size_t ori_p_id = pcl_index[p_id];
			const Position& p_p = pcl_pos[ori_p_id];
			const Displacement& p_u = pcl_u[p_id];
			const double p_x = p_p.x + p_u.ux;
			const double p_y = p_p.y + p_u.uy;
			const double p_z = p_p.z + p_u.uz;
			const double p_r = 0.5 * pow(pcl_vol[p_id], one_third);
			if (prmesh->detect_collision_with_point(
				p_x, p_y, p_z, p_r, dist, lnorm, cur_cont_pos))
			{
				prmesh->get_global_vector(lnorm, gnorm);
				pcm->cal_contact_force(
					substep_index,
					ori_p_id,
					dist,
					lnorm,
					cur_cont_pos,
					p_r + p_r,
					pv_getter,
					lcont_f.vec);
				prmesh->get_global_vector(lcont_f.vec, gcont_f.vec);
				// apply contact force to mesh
				// denoted the element as in contact
				const size_t e_id = pcl_in_elems[p_id];
				elem_node_at_surface[e_id * 4] |= 2;
				elem_node_at_surface[e_id * 4 + 1] |= 2;
				elem_node_at_surface[e_id * 4 + 2] |= 2;
				elem_node_at_surface[e_id * 4 + 3] |= 2;
				// contact force
				const ShapeFunc& p_N = pcl_N[p_id];
				Force& en_f1 = elem_node_force[e_id * 4];
				en_f1.fx += p_N.N1 * gcont_f.fx;
				en_f1.fy += p_N.N1 * gcont_f.fy;
				en_f1.fz += p_N.N1 * gcont_f.fz;
				Force& en_f2 = elem_node_force[e_id * 4 + 1];
				en_f2.fx += p_N.N2 * gcont_f.fx;
				en_f2.fy += p_N.N2 * gcont_f.fy;
				en_f2.fz += p_N.N2 * gcont_f.fz;
				Force& en_f3 = elem_node_force[e_id * 4 + 2];
				en_f3.fx += p_N.N3 * gcont_f.fx;
				en_f3.fy += p_N.N3 * gcont_f.fy;
				en_f3.fz += p_N.N3 * gcont_f.fz;
				Force& en_f4 = elem_node_force[e_id * 4 + 3];
				en_f4.fx += p_N.N4 * gcont_f.fx;
				en_f4.fy += p_N.N4 * gcont_f.fy;
				en_f4.fz += p_N.N4 * gcont_f.fz;
				// apply contact force to rigid body
				const Point3D& rm_cen = prmesh->get_pos();
				rc_cf_tmp.add_force(p_x, p_y, p_z,
					-gcont_f.fx, -gcont_f.fy, -gcont_f.fz,
					rm_cen.x, rm_cen.y, rm_cen.z);
			}
		}
		rc_cf = rc_cf_tmp;
	}
}
