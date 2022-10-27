#include "SimulationsOMP_pcp.h"

#include <assert.h>

#include "ParaUtil.h"
#include "Step_T3D_CHM_TBB.h"

#define Block_Low(blk_id, blk_num, data_num) ((blk_id) * (data_num) / (blk_num))

namespace Step_T3D_CHM_TBB_Task
{
	constexpr double one_third = 1.0 / 3.0;
	constexpr double one_fourth = 0.25;
	
	void InitPcl::init(size_t thread_num) noexcept
	{
		pcl_m_s = stp.pcl_m_s;
		pcl_density_s = stp.pcl_density_s;
		in_pcl_in_elems = stp.in_pcl_in_elems;
		in_prev_pcl_ids = stp.in_prev_pcl_ids;
		task_num = ParaUtil::cal_task_num<
			min_pcl_num_per_task, init_pcl_task_num_per_thread>(
				thread_num, stp.prev_valid_pcl_num);
	}

	void InitPcl::work(size_t tsk_id, InitPclRes& res) const
	{
		Model_T3D_CHM_mt& md = *stp.pmodel;
		Position* const pcl_pos = const_cast<Position* const>(stp.pcl_pos);
		const auto& spva0 = stp.spvas[0];
		const size_t p_id0 = Block_Low(tsk_id, task_num, stp.prev_valid_pcl_num);
		const size_t p_id1 = Block_Low(tsk_id + 1, task_num, stp.prev_valid_pcl_num);
		size_t valid_pcl_num = 0;
		double max_pcl_vol = 0.0;
		for (size_t p_id = p_id0; p_id < p_id1; ++p_id)
		{
			const size_t ori_p_id = spva0.pcl_index[p_id];
			Position& p_p = pcl_pos[ori_p_id];
			Displacement& p_d = spva0.pcl_u_s[p_id];
			p_p.x += p_d.ux;
			p_p.y += p_d.uy;
			p_p.z += p_d.uz;
			p_d.ux = 0.0;
			p_d.uy = 0.0;
			p_d.uz = 0.0;
			ShapeFunc& p_N = spva0.pcl_N[p_id];
			size_t e_id = md.find_pcl_in_which_elem(p_p.x, p_p.y, p_p.z, p_N);
			if (e_id == SIZE_MAX)
				e_id = md.find_pcl_in_which_elem_tol(p_p.x, p_p.y, p_p.z, p_N);
			if (e_id != SIZE_MAX)
				++valid_pcl_num;
			const double p_vol = pcl_m_s[ori_p_id]
				/ (pcl_density_s[ori_p_id] * (1.0 - spva0.pcl_n[p_id]));
			if (max_pcl_vol < p_vol)
				max_pcl_vol = p_vol;
			in_pcl_in_elems[p_id] = e_id;
			in_prev_pcl_ids[p_id] = p_id;
		}
		res.pcl_num = valid_pcl_num;
		res.max_pcl_vol = max_pcl_vol;
	}
	
	void MapPclToBgMesh::init() noexcept
	{
		// pcl data
		pcl_m_s = stp.pcl_m_s;
		pcl_vol_s = stp.pcl_vol_s;
		pcl_bf_s = stp.pcl_bf_s;
		pcl_bf_f = stp.pcl_bf_f;
		pcl_t = stp.pcl_t;
		pcl_vol = stp.pcl_vol;
		// bg mesh data
		elem_N_abc = stp.elem_N_abc;
		elem_vol = stp.elem_vol;
		elem_pcl_m_s = stp.elem_pcl_m_s;
		elem_pcl_m_f = stp.elem_pcl_m_f;
		elem_pcl_n = stp.elem_pcl_n;
		elem_density_f = stp.elem_density_f;
		elem_p = stp.elem_p;
		elem_node_vm_s = stp.elem_node_vm_s;
		elem_node_vm_f = stp.elem_node_vm_f;
		elem_node_force_s = stp.elem_node_force_s;
		elem_node_force_f = stp.elem_node_force_f;
		// pcl range
		pcl_in_elems = stp.pcl_in_elems;
		prev_pcl_ids = stp.prev_pcl_ids;
	}

	void MapPclToBgMesh::update(size_t tsk_num) noexcept
	{
		const auto& spva0 = stp.spvas[stp.next_spva_id()];
		const auto& spva1 = stp.spvas[stp.prev_spva_id()];
		pcl_index0 = spva0.pcl_index;
		pcl_n0 = spva0.pcl_n;
		pcl_density_f0 = spva0.pcl_density_f;
		pcl_v_s0 = spva0.pcl_v_s;
		pcl_v_f0 = spva0.pcl_v_f;
		pcl_u_s0 = spva0.pcl_u_s;
		pcl_u_f0 = spva0.pcl_u_f;
		pcl_stress0 = spva0.pcl_stress;
		pcl_N0 = spva0.pcl_N;
		pcl_index1 = spva1.pcl_index;
		pcl_n1 = spva1.pcl_n;
		pcl_density_f1 = spva1.pcl_density_f;
		pcl_v_s1 = spva1.pcl_v_s;
		pcl_v_f1 = spva1.pcl_v_f;
		pcl_u_s1 = spva1.pcl_u_s;
		pcl_u_f1 = spva1.pcl_u_f;
		pcl_stress1 = spva1.pcl_stress;
		pcl_p1 = spva1.pcl_p;
		pcl_N1 = spva1.pcl_N;

		task_num = tsk_num;
	}
	
	void MapPclToBgMesh::work(size_t tsk_id, MapPclToBgMeshRes& res) const
	{
		size_t e_id;
		size_t p_id0 = Block_Low(tsk_id, task_num, stp.valid_pcl_num);
		e_id = pcl_in_elems[p_id0];
		while (p_id0 != SIZE_MAX && e_id == pcl_in_elems[--p_id0]);
		++p_id0;
		assert(p_id0 <= stp.valid_pcl_num);
		size_t p_id1 = Block_Low(tsk_id + 1, task_num, stp.valid_pcl_num);
		e_id = pcl_in_elems[p_id1];
		while (p_id1 != SIZE_MAX && e_id == pcl_in_elems[--p_id1]);
		++p_id1;
		assert(p_id1 <= stp.valid_pcl_num);

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
			const double f_seep_tmp = stp.miu / stp.k * p_n * p_n * p_vol;
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

			if (e_id != pcl_in_elems[p_id + 1])
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

				e_id = pcl_in_elems[p_id + 1];
#ifdef _DEBUG
				assert(e_id < stp.elem_num || e_id == SIZE_MAX);
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

		// contact force calculation
		Force3D rcy_cf;
		ContactRigidBody& crb = stp.cont_rigid_body;
		if (crb.has_rigid_cylinder())
		{
			crb.apply_rigid_cylinder(p_id0, p_id1, rcy_cf);
			res.react_force = rcy_cf;
		}
		if (crb.has_rigid_mesh())
		{
			crb.apply_t3d_rigid_object(p_id0, p_id1, rcy_cf);
			res.react_force = rcy_cf;
		}
	}
	
	void UpdateAccelerationAndVelocity::init() noexcept
	{
		elem_pcl_m_s = stp.elem_pcl_m_s;
		elem_pcl_m_f = stp.elem_pcl_m_f;
		elem_node_force_s = stp.elem_node_force_s;
		elem_node_force_f = stp.elem_node_force_f;
		elem_node_vm_s = stp.elem_node_vm_s;
		elem_node_vm_f = stp.elem_node_vm_f;
		node_am_s = stp.node_am_s;
		node_am_f = stp.node_am_f;
		node_a_s = stp.node_a_s;
		node_a_f = stp.node_a_f;
		node_v_s = stp.node_v_s;
		node_v_f = stp.node_v_f;
		node_has_vbc_s = stp.node_has_vbc_s;
		node_has_vbc_f = stp.node_has_vbc_f;
		node_vbc_vec_s = stp.node_vbc_vec_s;
		node_vbc_vec_f = stp.node_vbc_vec_f;
		// node ranges
		node_ids = stp.node_ids;
		node_elem_offs = stp.node_elem_offs;
	}

	void UpdateAccelerationAndVelocity::update(size_t tsk_num)
	{
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

			if (n_id != node_ids[ve_id + 1])
			{
				// solid
				n_am_s *= one_fourth;
				node_am_s[n_id] = n_am_s;
				Acceleration& n_a_s = node_a_s[n_id];
				n_a_s.ax = n_fx_s / n_am_s;
				n_a_s.ay = n_fy_s / n_am_s;
				n_a_s.az = n_fz_s / n_am_s;
				Velocity& n_v_s = node_v_s[n_id];
				n_v_s.vx = n_vmx_s / n_vm_s + n_a_s.ax * stp.dtime;
				n_v_s.vy = n_vmy_s / n_vm_s + n_a_s.ay * stp.dtime;
				n_v_s.vz = n_vmz_s / n_vm_s + n_a_s.az * stp.dtime;
				const NodeVBCVec& n_vbc_vec_s = node_vbc_vec_s[n_id];
				vbc_len = n_a_s.ax * n_vbc_vec_s.x + n_a_s.ay * n_vbc_vec_s.y + n_a_s.az * n_vbc_vec_s.z;
				n_a_s.ax -= vbc_len * n_vbc_vec_s.x;
				n_a_s.ay -= vbc_len * n_vbc_vec_s.y;
				n_a_s.az -= vbc_len * n_vbc_vec_s.z;
				vbc_len = n_v_s.vx * n_vbc_vec_s.x + n_v_s.vy * n_vbc_vec_s.y + n_v_s.vz * n_vbc_vec_s.z;
				n_v_s.vx -= vbc_len * n_vbc_vec_s.x;
				n_v_s.vy -= vbc_len * n_vbc_vec_s.y;
				n_v_s.vz -= vbc_len * n_vbc_vec_s.z;
				const NodeHasVBC& n_has_vbc_s = node_has_vbc_s[n_id];
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
				n_v_f.vx = n_vmx_f / n_vm_f + n_a_f.ax * stp.dtime;
				n_v_f.vy = n_vmy_f / n_vm_f + n_a_f.ay * stp.dtime;
				n_v_f.vz = n_vmz_f / n_vm_f + n_a_f.az * stp.dtime;
				const NodeVBCVec& n_vbc_vec_f = node_vbc_vec_f[n_id];
				vbc_len = n_a_f.ax * n_vbc_vec_f.x + n_a_f.ay * n_vbc_vec_f.y + n_a_f.az * n_vbc_vec_f.z;
				n_a_f.ax -= vbc_len * n_vbc_vec_f.x;
				n_a_f.ay -= vbc_len * n_vbc_vec_f.y;
				n_a_f.az -= vbc_len * n_vbc_vec_f.z;
				vbc_len = n_v_f.vx * n_vbc_vec_f.x + n_v_f.vy * n_vbc_vec_f.y + n_v_f.vz * n_vbc_vec_f.z;
				n_v_f.vx -= vbc_len * n_vbc_vec_f.x;
				n_v_f.vy -= vbc_len * n_vbc_vec_f.y;
				n_v_f.vz -= vbc_len * n_vbc_vec_f.z;
				const NodeHasVBC& n_has_vbc_f = node_has_vbc_f[n_id];
				bc_mask = SIZE_MAX + size_t(n_has_vbc_f.has_vx_bc);
				n_a_f.iax &= bc_mask;
				n_v_f.ivx &= bc_mask;
				bc_mask = SIZE_MAX + size_t(n_has_vbc_f.has_vy_bc);
				n_a_f.iay &= bc_mask;
				n_v_f.ivy &= bc_mask;
				bc_mask = SIZE_MAX + size_t(n_has_vbc_f.has_vz_bc);
				n_a_f.iaz &= bc_mask;
				n_v_f.ivz &= bc_mask;

				n_id = node_ids[ve_id + 1];
#ifdef _DEBUG
				assert(n_id < stp.node_num || n_id == SIZE_MAX);
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

	void CalElemDeAndMapToNode::init() noexcept
	{
		elem_node_id = stp.elem_node_id;
		elem_N_abc = stp.elem_N_abc;
		elem_pcl_m_s = stp.elem_pcl_m_s;
		elem_pcl_m_f = stp.elem_pcl_m_f;
		elem_pcl_n = stp.elem_pcl_n;
		node_v_s = stp.node_v_s;
		node_v_f = stp.node_v_f;
		elem_de = stp.elem_de;
		elem_m_de_vol_s = stp.elem_m_de_vol_s;
		elem_m_de_vol_f = stp.elem_m_de_vol_f;
		elem_ids = stp.elem_ids;
	}

	void CalElemDeAndMapToNode::update(size_t tsk_num)
	{
		elem_num = stp.valid_elem_num;
		task_num = tsk_num;
	}

	void CalElemDeAndMapToNode::work(size_t tsk_id) const
	{
		const size_t ve_id0 = Block_Low(tsk_id, task_num, elem_num);
		const size_t ve_id1 = Block_Low(tsk_id + 1, task_num, elem_num);
		double e_de_vol_s, e_de_vol_f;
		for (size_t ve_id = ve_id0; ve_id < ve_id1; ++ve_id)
		{
			const size_t e_id = elem_ids[ve_id];
#ifdef _DEBUG
			assert(e_id < stp.elem_num);
#endif
			const ElemNodeIndex& eni = elem_node_id[e_id];
			const Velocity& n1_v_s = node_v_s[eni.n1];
			const Velocity& n2_v_s = node_v_s[eni.n2];
			const Velocity& n3_v_s = node_v_s[eni.n3];
			const Velocity& n4_v_s = node_v_s[eni.n4];
			const DShapeFuncABC& e_dN = elem_N_abc[e_id];
			StrainInc& e_de = elem_de[e_id];
			e_de.de11 = (e_dN.dN1_dx * n1_v_s.vx + e_dN.dN2_dx * n2_v_s.vx + e_dN.dN3_dx * n3_v_s.vx + e_dN.dN4_dx * n4_v_s.vx) * stp.dtime;
			e_de.de22 = (e_dN.dN1_dy * n1_v_s.vy + e_dN.dN2_dy * n2_v_s.vy + e_dN.dN3_dy * n3_v_s.vy + e_dN.dN4_dy * n4_v_s.vy) * stp.dtime;
			e_de.de33 = (e_dN.dN1_dz * n1_v_s.vz + e_dN.dN2_dz * n2_v_s.vz + e_dN.dN3_dz * n3_v_s.vz + e_dN.dN4_dz * n4_v_s.vz) * stp.dtime;
			e_de.de12 = (e_dN.dN1_dx * n1_v_s.vy + e_dN.dN2_dx * n2_v_s.vy + e_dN.dN3_dx * n3_v_s.vy + e_dN.dN4_dx * n4_v_s.vy
					   + e_dN.dN1_dy * n1_v_s.vx + e_dN.dN2_dy * n2_v_s.vx + e_dN.dN3_dy * n3_v_s.vx + e_dN.dN4_dy * n4_v_s.vx) * stp.dtime * 0.5;
			e_de.de23 = (e_dN.dN1_dy * n1_v_s.vz + e_dN.dN2_dy * n2_v_s.vz + e_dN.dN3_dy * n3_v_s.vz + e_dN.dN4_dy * n4_v_s.vz
					   + e_dN.dN1_dz * n1_v_s.vy + e_dN.dN2_dz * n2_v_s.vy + e_dN.dN3_dz * n3_v_s.vy + e_dN.dN4_dz * n4_v_s.vy) * stp.dtime * 0.5;
			e_de.de31 = (e_dN.dN1_dz * n1_v_s.vx + e_dN.dN2_dz * n2_v_s.vx + e_dN.dN3_dz * n3_v_s.vx + e_dN.dN4_dz * n4_v_s.vx
					   + e_dN.dN1_dx * n1_v_s.vz + e_dN.dN2_dx * n2_v_s.vz + e_dN.dN3_dx * n3_v_s.vz + e_dN.dN4_dx * n4_v_s.vz) * stp.dtime * 0.5;
			e_de_vol_s = e_de.de11 + e_de.de22 + e_de.de33;
			elem_m_de_vol_s[e_id] = elem_pcl_m_s[e_id] * e_de_vol_s;
			const Velocity& n1_v_f = node_v_f[eni.n1];
			const Velocity& n2_v_f = node_v_f[eni.n2];
			const Velocity& n3_v_f = node_v_f[eni.n3];
			const Velocity& n4_v_f = node_v_f[eni.n4];
			e_de_vol_f = (1.0 - elem_pcl_n[e_id]) / elem_pcl_n[e_id] * -e_de_vol_s
					- (e_dN.dN1_dx * n1_v_f.vx + e_dN.dN2_dx * n2_v_f.vx + e_dN.dN3_dx * n3_v_f.vx + e_dN.dN4_dx * n4_v_f.vx
					+ e_dN.dN1_dy * n1_v_f.vy + e_dN.dN2_dy * n2_v_f.vy + e_dN.dN3_dy * n3_v_f.vy + e_dN.dN4_dy * n4_v_f.vy
					+ e_dN.dN1_dz * n1_v_f.vz + e_dN.dN2_dz * n2_v_f.vz + e_dN.dN3_dz * n3_v_f.vz + e_dN.dN4_dz * n4_v_f.vz) * stp.dtime;
			elem_m_de_vol_f[e_id] = elem_pcl_m_f[e_id] * e_de_vol_f;
			e_de_vol_s *= one_third;
			e_de.de11 -= e_de_vol_s;
			e_de.de22 -= e_de_vol_s;
			e_de.de33 -= e_de_vol_s;
		}
	}
	
	void CalNodeDe::init() noexcept
	{
		elem_m_de_vol_s = stp.elem_m_de_vol_s;
		elem_m_de_vol_f = stp.elem_m_de_vol_f;
		node_am_s = stp.node_am_s;
		node_am_f = stp.node_am_f;
		node_de_vol_s = stp.node_de_vol_s;
		node_de_vol_f = stp.node_de_vol_f;
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

		size_t e_id;
		double n_am_de_vol_s = 0.0;
		double n_am_de_vol_f = 0.0;
		n_id = node_ids[ve_id0];
#ifdef _DEBUG
		assert(n_id < stp.node_num || n_id == SIZE_MAX);
#endif
		for (size_t ve_id = ve_id0; ve_id < ve_id1; ++ve_id)
		{
			e_id = node_elem_offs[ve_id] / 4;
#ifdef _DEBUG
			assert(e_id < stp.elem_num);
#endif
			n_am_de_vol_s += elem_m_de_vol_s[e_id];
			n_am_de_vol_f += elem_m_de_vol_f[e_id];
			if (n_id != node_ids[ve_id + 1])
			{
				node_de_vol_s[n_id] = n_am_de_vol_s * one_fourth / node_am_s[n_id];
				node_de_vol_f[n_id] = n_am_de_vol_f * one_fourth / node_am_f[n_id];
				n_id = node_ids[ve_id + 1];
#ifdef _DEBUG
				assert(n_id < stp.node_num || n_id == SIZE_MAX);
#endif
				n_am_de_vol_s = 0.0;
				n_am_de_vol_f = 0.0;
			}
		}
	}

	void MapBgMeshToPcl::init() noexcept
	{
		elem_node_id = stp.elem_node_id;
		node_a_s = stp.node_a_s;
		node_a_f = stp.node_a_f;
		node_v_s = stp.node_v_s;
		node_v_f = stp.node_v_f;
		elem_density_f = stp.elem_density_f;
		elem_pcl_n = stp.elem_pcl_n;
		elem_p = stp.elem_p;
		elem_de = stp.elem_de;

		node_de_vol_s = stp.node_de_vol_s;
		node_de_vol_f = stp.node_de_vol_f;
		pcl_pos = stp.pcl_pos;
		pcl_mat_model = stp.pcl_mat_model;
		// range
		pcl_in_elems = stp.pcl_in_elems;
		prev_pcl_ids = stp.prev_pcl_ids;
		in_pcl_in_elems = stp.in_pcl_in_elems;
		in_prev_pcl_ids = stp.in_prev_pcl_ids;

		// cavitation
		pcl_u_cav = stp.pcl_u_cav;
		pcl_is_cavitated = stp.pcl_is_cavitated;
		elem_u_cav = stp.elem_u_cav;
	}

	void MapBgMeshToPcl::update(size_t tsk_num) noexcept
	{
		const auto& spva0 = stp.spvas[stp.next_spva_id()];
		const auto& spva1 = stp.spvas[stp.prev_spva_id()];
		pcl_index0 = spva0.pcl_index;
		pcl_density_f0 = spva0.pcl_density_f;
		pcl_n0 = spva0.pcl_n;
		pcl_p0 = spva0.pcl_p;
		pcl_N0 = spva0.pcl_N;
		pcl_v_s0 = spva0.pcl_v_s;
		pcl_v_f0 = spva0.pcl_v_f;
		pcl_u_s0 = spva0.pcl_u_s;
		pcl_u_f0 = spva0.pcl_u_f;
		pcl_stress0 = spva0.pcl_stress;
		pcl_strain0 = spva0.pcl_strain;
		pcl_estrain0 = spva0.pcl_estrain;
		pcl_pstrain0 = spva0.pcl_pstrain;
		pcl_strain1 = spva1.pcl_strain;
		pcl_estrain1 = spva1.pcl_estrain;
		pcl_pstrain1 = spva1.pcl_pstrain;

		pcl_num = stp.prev_valid_pcl_num;
		task_num = tsk_num;
	}
	
	void MapBgMeshToPcl::work(size_t tsk_id, MapBgMeshToPclRes &res) const
	{
		size_t e_id;
		size_t p_id0 = Block_Low(tsk_id, task_num, stp.prev_valid_pcl_num);
		e_id = pcl_in_elems[p_id0];
		while (p_id0 != SIZE_MAX && e_id == pcl_in_elems[--p_id0]);
		++p_id0;
		assert(p_id0 <= stp.prev_valid_pcl_num);
		size_t p_id1 = Block_Low(tsk_id + 1, task_num, stp.prev_valid_pcl_num);
		e_id = pcl_in_elems[p_id1];
		while (p_id1 != SIZE_MAX && e_id == pcl_in_elems[--p_id1]);
		++p_id1;
		assert(p_id1 <= stp.prev_valid_pcl_num);

		const Acceleration* pn1_a_s, * pn2_a_s, * pn3_a_s, * pn4_a_s;
		const Acceleration* pn1_a_f, * pn2_a_f, * pn3_a_f, * pn4_a_f;
		const Velocity* pn1_v_s, * pn2_v_s, * pn3_v_s, * pn4_v_s;
		const Velocity* pn1_v_f, * pn2_v_f, * pn3_v_f, * pn4_v_f;
		double e_de_vol_s, e_de_vol_f, e_density_f, e_n, e_p;
		StrainInc* pe_de;
		double dstrain[6];
		Model_T3D_CHM_mt& md = *(stp.pmodel);
		size_t valid_pcl_num = 0;
		e_id = SIZE_MAX;
		// cavitation
		double e_u_cav, Kf_ratio;
		for (size_t p_id = p_id0; p_id < p_id1; ++p_id)
		{
			if (e_id != pcl_in_elems[p_id])
			{
				e_id = pcl_in_elems[p_id];
#ifdef _DEBUG
				assert(e_id < stp.elem_num);
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
				//e_p = elem_p[e_id] + stp.Kf * e_de_vol_f;
				// cavitation
				Kf_ratio = 1.0;
				if (stp.m_cav != 0.0) // consider cavitation
				{
					e_u_cav = elem_u_cav[e_id];
					if (elem_p[e_id] < 0.0)
					{
						const double tmp = elem_p[e_id] / elem_u_cav[e_id] + stp.u_cav_off;
						Kf_ratio = tmp < stp.u_div_u_cav_lim ? (1.0 / (1.0 + pow(tmp, stp.m_cav))) : (1.0 / max_Kf_ratio_divider);
					}
				}
				e_p = elem_p[e_id] + Kf_ratio * stp.Kf * e_de_vol_f;
				
				pe_de = elem_de + e_id;
				e_de_vol_s *= one_third;
				pe_de->de11 += e_de_vol_s;
				pe_de->de22 += e_de_vol_s;
				pe_de->de33 += e_de_vol_s;
			}

			// update velocity
			ShapeFunc& p_N = pcl_N0[p_id];
			Velocity& p_v_s0 = pcl_v_s0[p_id];
			p_v_s0.vx += (p_N.N1 * pn1_a_s->ax + p_N.N2 * pn2_a_s->ax + p_N.N3 * pn3_a_s->ax + p_N.N4 * pn4_a_s->ax) * stp.dtime;
			p_v_s0.vy += (p_N.N1 * pn1_a_s->ay + p_N.N2 * pn2_a_s->ay + p_N.N3 * pn3_a_s->ay + p_N.N4 * pn4_a_s->ay) * stp.dtime;
			p_v_s0.vz += (p_N.N1 * pn1_a_s->az + p_N.N2 * pn2_a_s->az + p_N.N3 * pn3_a_s->az + p_N.N4 * pn4_a_s->az) * stp.dtime;
			Velocity& p_v_f0 = pcl_v_f0[p_id];
			p_v_f0.vx += (p_N.N1 * pn1_a_f->ax + p_N.N2 * pn2_a_f->ax + p_N.N3 * pn3_a_f->ax + p_N.N4 * pn4_a_f->ax) * stp.dtime;
			p_v_f0.vy += (p_N.N1 * pn1_a_f->ay + p_N.N2 * pn2_a_f->ay + p_N.N3 * pn3_a_f->ay + p_N.N4 * pn4_a_f->ay) * stp.dtime;
			p_v_f0.vz += (p_N.N1 * pn1_a_f->az + p_N.N2 * pn2_a_f->az + p_N.N3 * pn3_a_f->az + p_N.N4 * pn4_a_f->az) * stp.dtime;

			// update displacement
			Displacement& p_u_s0 = pcl_u_s0[p_id];
			p_u_s0.ux += (p_N.N1 * pn1_v_s->vx + p_N.N2 * pn2_v_s->vx + p_N.N3 * pn3_v_s->vx + p_N.N4 * pn4_v_s->vx) * stp.dtime;
			p_u_s0.uy += (p_N.N1 * pn1_v_s->vy + p_N.N2 * pn2_v_s->vy + p_N.N3 * pn3_v_s->vy + p_N.N4 * pn4_v_s->vy) * stp.dtime;
			p_u_s0.uz += (p_N.N1 * pn1_v_s->vz + p_N.N2 * pn2_v_s->vz + p_N.N3 * pn3_v_s->vz + p_N.N4 * pn4_v_s->vz) * stp.dtime;
			Displacement& p_u_f0 = pcl_u_f0[p_id];
			p_u_f0.ux += (p_N.N1 * pn1_v_f->vx + p_N.N2 * pn2_v_f->vx + p_N.N3 * pn3_v_f->vx + p_N.N4 * pn4_v_f->vx) * stp.dtime;
			p_u_f0.uy += (p_N.N1 * pn1_v_f->vy + p_N.N2 * pn2_v_f->vy + p_N.N3 * pn3_v_f->vy + p_N.N4 * pn4_v_f->vy) * stp.dtime;
			p_u_f0.uz += (p_N.N1 * pn1_v_f->vz + p_N.N2 * pn2_v_f->vz + p_N.N3 * pn3_v_f->vz + p_N.N4 * pn4_v_f->vz) * stp.dtime;

			// update location (in which element)
			const size_t ori_p_id = pcl_index0[p_id];
#ifdef _DEBUG
			assert(ori_p_id < stp.ori_pcl_num);
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
				++valid_pcl_num;
			in_pcl_in_elems[p_id] = p_e_id;
			in_prev_pcl_ids[p_id] = p_id;
#ifdef _DEBUG
			assert(p_e_id < stp.elem_num || p_e_id == SIZE_MAX);
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

			const size_t prev_p_id = prev_pcl_ids[p_id];
#ifdef _DEBUG
			assert(prev_p_id < stp.prev_valid_pcl_num_tmp);
#endif
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

			if (stp.m_cav != 0.0)
			{
				pcl_u_cav[p_id] = e_u_cav;
				pcl_is_cavitated[p_id] = Kf_ratio;
			}
		}

		res.pcl_num = valid_pcl_num;
	}

	void ContactRigidBody::init(
		double max_pcl_vol,
		bool is_first_step) noexcept
	{
		Model_T3D_CHM_mt& md = *stp.pmodel;
		pcm_s = md.get_contact_model_s();
		pcm_f = md.get_contact_model_f();
		prcy = nullptr;
		prmesh = nullptr;
		if (md.has_rigid_cylinder())
		{
			prcy = &md.get_rigid_cylinder();
			if (is_first_step)
				prcy->reset_cont_force();
		}
		if (md.has_t3d_rigid_mesh())
		{
			prmesh = &md.get_t3d_rigid_mesh();
			// set max dist for efficiency
			prmesh->init_max_dist(0.5 * pow(max_pcl_vol, one_third) * 4.0);
			if (is_first_step)
				prmesh->reset_cont_force();
		}
		//
		pcl_pos = stp.pcl_pos;
		pcl_vol = stp.pcl_vol;
		elem_node_force_s = stp.elem_node_force_s;
		elem_node_force_f = stp.elem_node_force_f;
		// pcl range
		pcl_in_elems = stp.pcl_in_elems;
	}

	void ContactRigidBody::update() noexcept
	{
		const auto& spva0 = stp.spvas[stp.next_spva_id()];
		pcl_index = spva0.pcl_index;
		pcl_u_s = spva0.pcl_u_s;
		pcl_u_f = spva0.pcl_u_f;
		pcl_N = spva0.pcl_N;
	}

	void ContactRigidBody::apply_rigid_cylinder(
		size_t p_id0,
		size_t p_id1,
		Force3D& rc_cf)
		const noexcept
	{
		double dist;
		Vector3D lnorm, gnorm;
		Point3D cur_cont_pos;
		Force lcont_fs, gcont_fs;
		Force lcont_ff, gcont_ff;
		rc_cf.reset();
		const size_t substp_id = stp.substep_index;
		for (size_t p_id = p_id0; p_id < p_id1; ++p_id)
		{
			const size_t ori_p_id = pcl_index[p_id];
			const Position& p_p = pcl_pos[ori_p_id];
			const Displacement& p_u = pcl_u_s[p_id];
			const double p_x = p_p.x + p_u.ux;
			const double p_y = p_p.y + p_u.uy;
			const double p_z = p_p.z + p_u.uz;
			const double p_r = 0.5 * pow(pcl_vol[p_id], one_third);
			const ShapeFunc& p_N = pcl_N[p_id];
			const size_t e_id = pcl_in_elems[p_id];
			ParticleVariablesGetter pv_place_holder;
			if (prcy->detect_collision_with_point(
				p_x, p_y, p_z, p_r, dist, lnorm, cur_cont_pos))
			{
				prcy->get_global_vector(lnorm, gnorm);
				// solid pcl
				pcm_s->cal_contact_force(
					substp_id,
					ori_p_id,
					dist,
					lnorm,
					cur_cont_pos,
					p_r + p_r,
					pv_place_holder,
					lcont_fs.vec);
				prcy->get_global_vector(lcont_fs.vec, gcont_fs.vec);
				Force& en_f_s1 = elem_node_force_s[e_id * 4];
				en_f_s1.fx += p_N.N1 * gcont_fs.fx;
				en_f_s1.fy += p_N.N1 * gcont_fs.fy;
				en_f_s1.fz += p_N.N1 * gcont_fs.fz;
				Force& en_f_s2 = elem_node_force_s[e_id * 4 + 1];
				en_f_s2.fx += p_N.N2 * gcont_fs.fx;
				en_f_s2.fy += p_N.N2 * gcont_fs.fy;
				en_f_s2.fz += p_N.N2 * gcont_fs.fz;
				Force& en_f_s3 = elem_node_force_s[e_id * 4 + 2];
				en_f_s3.fx += p_N.N3 * gcont_fs.fx;
				en_f_s3.fy += p_N.N3 * gcont_fs.fy;
				en_f_s3.fz += p_N.N3 * gcont_fs.fz;
				Force& en_f_s4 = elem_node_force_s[e_id * 4 + 3];
				en_f_s4.fx += p_N.N4 * gcont_fs.fx;
				en_f_s4.fy += p_N.N4 * gcont_fs.fy;
				en_f_s4.fz += p_N.N4 * gcont_fs.fz;
				// fluid pcl
				pcm_f->cal_contact_force(
					substp_id,
					ori_p_id,
					dist,
					lnorm,
					cur_cont_pos,
					p_r + p_r,
					pv_place_holder,
					lcont_ff.vec);
				prcy->get_global_vector(lcont_ff.vec, gcont_ff.vec);
				Force& en_f_f1 = elem_node_force_f[e_id * 4];
				en_f_f1.fx += p_N.N1 * gcont_ff.fx;
				en_f_f1.fy += p_N.N1 * gcont_ff.fy;
				en_f_f1.fz += p_N.N1 * gcont_ff.fz;
				Force& en_f_f2 = elem_node_force_f[e_id * 4 + 1];
				en_f_f2.fx += p_N.N2 * gcont_ff.fx;
				en_f_f2.fy += p_N.N2 * gcont_ff.fy;
				en_f_f2.fz += p_N.N2 * gcont_ff.fz;
				Force& en_f_f3 = elem_node_force_f[e_id * 4 + 2];
				en_f_f3.fx += p_N.N3 * gcont_ff.fx;
				en_f_f3.fy += p_N.N3 * gcont_ff.fy;
				en_f_f3.fz += p_N.N3 * gcont_ff.fz;
				Force& en_f_f4 = elem_node_force_f[e_id * 4 + 3];
				en_f_f4.fx += p_N.N4 * gcont_ff.fx;
				en_f_f4.fy += p_N.N4 * gcont_ff.fy;
				en_f_f4.fz += p_N.N4 * gcont_ff.fz;
				// apply contact force to rigid body
				const Point3D& rc_cen = prcy->get_centre();
				rc_cf.add_force(p_x, p_y, p_z,
					-(gcont_fs.fx + gcont_ff.fx),
					-(gcont_fs.fy + gcont_ff.fy),
					-(gcont_fs.fz + gcont_ff.fz),
					rc_cen.x, rc_cen.y, rc_cen.z);
			}
		}
	}

	void ContactRigidBody::apply_t3d_rigid_object(
		size_t p_id0,
		size_t p_id1,
		Force3D& rc_cf)
		const noexcept
	{
		double dist;
		Vector3D lnorm, gnorm;
		Point3D cur_cont_pos;
		Force lcont_fs, gcont_fs;
		Force lcont_ff, gcont_ff;
		rc_cf.reset();
		const size_t substp_id = stp.substep_index;
		for (size_t p_id = p_id0; p_id < p_id1; ++p_id)
		{
			const size_t ori_p_id = pcl_index[p_id];
			const Position& p_p = pcl_pos[ori_p_id];
			const Displacement& p_u = pcl_u_s[p_id];
			const double p_x = p_p.x + p_u.ux;
			const double p_y = p_p.y + p_u.uy;
			const double p_z = p_p.z + p_u.uz;
			const double p_r = 0.5 * pow(pcl_vol[p_id], one_third);
			const ShapeFunc& p_N = pcl_N[p_id];
			const size_t e_id = pcl_in_elems[p_id];
			ParticleVariablesGetter pv_place_holder;
			if (prmesh->detect_collision_with_point(
				p_x, p_y, p_z, p_r, dist, lnorm, cur_cont_pos))
			{
				prmesh->get_global_vector(lnorm, gnorm);
				// solid pcl
				pcm_s->cal_contact_force(
					substp_id,
					ori_p_id,
					dist,
					lnorm,
					cur_cont_pos,
					p_r + p_r,
					pv_place_holder,
					lcont_fs.vec);
				prmesh->get_global_vector(lcont_fs.vec, gcont_fs.vec);
				Force& en_f_s1 = elem_node_force_s[e_id * 4];
				en_f_s1.fx += p_N.N1 * gcont_fs.fx;
				en_f_s1.fy += p_N.N1 * gcont_fs.fy;
				en_f_s1.fz += p_N.N1 * gcont_fs.fz;
				Force& en_f_s2 = elem_node_force_s[e_id * 4 + 1];
				en_f_s2.fx += p_N.N2 * gcont_fs.fx;
				en_f_s2.fy += p_N.N2 * gcont_fs.fy;
				en_f_s2.fz += p_N.N2 * gcont_fs.fz;
				Force& en_f_s3 = elem_node_force_s[e_id * 4 + 2];
				en_f_s3.fx += p_N.N3 * gcont_fs.fx;
				en_f_s3.fy += p_N.N3 * gcont_fs.fy;
				en_f_s3.fz += p_N.N3 * gcont_fs.fz;
				Force& en_f_s4 = elem_node_force_s[e_id * 4 + 3];
				en_f_s4.fx += p_N.N4 * gcont_fs.fx;
				en_f_s4.fy += p_N.N4 * gcont_fs.fy;
				en_f_s4.fz += p_N.N4 * gcont_fs.fz;
				// fluid pcl
				pcm_f->cal_contact_force(
					substp_id,
					ori_p_id,
					dist,
					lnorm,
					cur_cont_pos,
					p_r + p_r,
					pv_place_holder,
					lcont_ff.vec);
				prmesh->get_global_vector(lcont_ff.vec, gcont_ff.vec);
				Force& en_f_f1 = elem_node_force_f[e_id * 4];
				en_f_f1.fx += p_N.N1 * gcont_ff.fx;
				en_f_f1.fy += p_N.N1 * gcont_ff.fy;
				en_f_f1.fz += p_N.N1 * gcont_ff.fz;
				Force& en_f_f2 = elem_node_force_f[e_id * 4 + 1];
				en_f_f2.fx += p_N.N2 * gcont_ff.fx;
				en_f_f2.fy += p_N.N2 * gcont_ff.fy;
				en_f_f2.fz += p_N.N2 * gcont_ff.fz;
				Force& en_f_f3 = elem_node_force_f[e_id * 4 + 2];
				en_f_f3.fx += p_N.N3 * gcont_ff.fx;
				en_f_f3.fy += p_N.N3 * gcont_ff.fy;
				en_f_f3.fz += p_N.N3 * gcont_ff.fz;
				Force& en_f_f4 = elem_node_force_f[e_id * 4 + 3];
				en_f_f4.fx += p_N.N4 * gcont_ff.fx;
				en_f_f4.fy += p_N.N4 * gcont_ff.fy;
				en_f_f4.fz += p_N.N4 * gcont_ff.fz;
				// apply contact force to rigid body
				const Point3D& rc_cen = prmesh->get_pos();
				rc_cf.add_force(p_x, p_y, p_z,
					-(gcont_fs.fx + gcont_ff.fx),
					-(gcont_fs.fy + gcont_ff.fy),
					-(gcont_fs.fz + gcont_ff.fz),
					rc_cen.x, rc_cen.y, rc_cen.z);
			}
		}
	}
}
