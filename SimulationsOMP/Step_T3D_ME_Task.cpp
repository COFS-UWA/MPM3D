#include "SimulationsOMP_pcp.h"

#include <assert.h>
#include "tbb/task_arena.h"

#include "ParallelUtils.h"
#include "Step_T3D_ME_Task.h"
#include "Step_T3D_ME_TBB.h"

namespace Step_T3D_ME_Task
{
	constexpr double one_third = 1.0 / 3.0;
	constexpr double one_fourth = 0.25;

	void CalData::set_model(Model_T3D_ME_mt &md) noexcept
	{
		pmodel = &md;
		pcl_m = md.pcl_m;
		pcl_bf = md.pcl_bf;
		pcl_t = md.pcl_t;
		pcl_pos = md.pcl_pos;
		pcl_vol = md.pcl_vol;
		pcl_mat_model = md.pcl_mat_model;
		
		auto& spva0 = spvas[0];
		const auto& md_spva0 = md.sorted_pcl_var_arrays[0];
		spva0.pcl_index = md_spva0.pcl_index;
		spva0.pcl_density = md_spva0.pcl_density;
		spva0.pcl_v = md_spva0.pcl_v;
		spva0.pcl_disp = md_spva0.pcl_disp;
		spva0.pcl_stress = md_spva0.pcl_stress;
		spva0.pcl_strain = md_spva0.pcl_strain;
		spva0.pcl_estrain = md_spva0.pcl_estrain;
		spva0.pcl_pstrain = md_spva0.pcl_pstrain;
		spva0.pcl_N = md_spva0.pcl_N;

		auto& spva1 = spvas[1];
		const auto& md_spva1 = md.sorted_pcl_var_arrays[1];
		spva1.pcl_index = md_spva1.pcl_index;
		spva1.pcl_density = md_spva1.pcl_density;
		spva1.pcl_v = md_spva1.pcl_v;
		spva1.pcl_disp = md_spva1.pcl_disp;
		spva1.pcl_stress = md_spva1.pcl_stress;
		spva1.pcl_strain = md_spva1.pcl_strain;
		spva1.pcl_estrain = md_spva1.pcl_estrain;
		spva1.pcl_pstrain = md_spva1.pcl_pstrain;
		spva1.pcl_N = md_spva1.pcl_N;

		elem_node_id = md.elem_node_id;
		elem_dN_abc = md.elem_dN_abc;
		elem_dN_d = md.elem_dN_d;
		elem_vol = md.elem_vol;
		elem_pcl_m = md.elem_pcl_m;
		elem_density = md.elem_density;
		elem_de = md.elem_de;
		elem_m_de_vol = md.elem_m_de_vol;
		elem_node_vm = md.elem_node_vm;
		elem_node_force = md.elem_node_force;
		node_a = md.node_a;
		node_v = md.node_v;
		node_has_vbc = md.node_has_vbc;
		node_am = md.node_am;
		node_de_vol = md.node_de_vol;

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
		Model_T3D_ME_mt& md = *cd.pmodel;
		Position* const pcl_pos = const_cast<Position* const>(cd.pcl_pos);
		const auto& spva0 = cd.spvas[0];
		const size_t p_id0 = ParallelUtils::block_low(wk_id, task_num, cd.prev_valid_pcl_num);
		const size_t p_id1 = ParallelUtils::block_low(wk_id + 1, task_num, cd.prev_valid_pcl_num);
		size_t valid_pcl_num = 0;
		for (size_t p_id = p_id0; p_id < p_id1; ++p_id)
		{
			const size_t ori_p_id = spva0.pcl_index[p_id];
			Position& p_p = pcl_pos[ori_p_id];
			Displacement& p_d = spva0.pcl_disp[p_id];
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
			ori_pcl_in_elem[p_id] = e_id;
			ori_cur_to_prev_pcl[p_id] = p_id;
		}
		pcl_in_mesh_num = valid_pcl_num;
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
		e_id = pcl_in_elem[p_id0];
#ifdef _DEBUG
		assert(e_id < cd.elem_num);
#endif
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

			// map pcl mass
			const double p_m = pcl_m[ori_p_id];
			e_p_m += p_m;

			// map pcl volume
			const double p_vol = p_m / pcl_density1[prev_p_id];
			pcl_vol[p_id] = p_vol;
			e_p_vol += p_vol;

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

			// shape function
			const ShapeFunc& p_N1 = pcl_N1[prev_p_id];
			ShapeFunc& p_N0 = pcl_N0[p_id];
			p_N0.N1 = p_N1.N1;
			p_N0.N2 = p_N1.N2;
			p_N0.N3 = p_N1.N3;
			p_N0.N4 = p_N1.N4;
			// map velocity
			double p_N_m;
			const Velocity& p_v1 = pcl_v1[prev_p_id];
			Velocity& p_v0 = pcl_v0[p_id];
			p_v0.vx = p_v1.vx;
			p_v0.vy = p_v1.vy;
			p_v0.vz = p_v1.vz;
			p_N_m = (p_N0.N1 > N_tol ? p_N0.N1 : N_tol) * p_m;
			en1_vm += p_N_m;
			en1_vmx += p_N_m * p_v0.vx;
			en1_vmy += p_N_m * p_v0.vy;
			en1_vmz += p_N_m * p_v0.vz;
			p_N_m = (p_N0.N2 > N_tol ? p_N0.N2 : N_tol) * p_m;
			en2_vm += p_N_m;
			en2_vmx += p_N_m * p_v0.vx;
			en2_vmy += p_N_m * p_v0.vy;
			en2_vmz += p_N_m * p_v0.vz;
			p_N_m = (p_N0.N3 > N_tol ? p_N0.N3 : N_tol) * p_m;
			en3_vm += p_N_m;
			en3_vmx += p_N_m * p_v0.vx;
			en3_vmy += p_N_m * p_v0.vy;
			en3_vmz += p_N_m * p_v0.vz;
			p_N_m = (p_N0.N4 > N_tol ? p_N0.N4 : N_tol) * p_m;
			en4_vm += p_N_m;
			en4_vmx += p_N_m * p_v0.vx;
			en4_vmy += p_N_m * p_v0.vy;
			en4_vmz += p_N_m * p_v0.vz;

			// displacement (for contact)
			const Displacement& p_d1 = pcl_disp1[prev_p_id];
			Displacement& p_d0 = pcl_disp0[p_id];
			p_d0.ux = p_d1.ux;
			p_d0.uy = p_d1.uy;
			p_d0.uz = p_d1.uz;

			// cal external load
			const Force& p_bf = pcl_bf[ori_p_id];
			const double one_fourth_bfx = one_fourth * p_bf.fx;
			const double one_fourth_bfy = one_fourth * p_bf.fy;
			const double one_fourth_bfz = one_fourth * p_bf.fz;
			const Force& p_t = pcl_t[ori_p_id];
			en1_fx += one_fourth_bfx + p_N0.N1 * p_t.fx;
			en1_fy += one_fourth_bfy + p_N0.N1 * p_t.fy;
			en1_fz += one_fourth_bfz + p_N0.N1 * p_t.fz;
			en2_fx += one_fourth_bfx + p_N0.N2 * p_t.fx;
			en2_fy += one_fourth_bfy + p_N0.N2 * p_t.fy;
			en2_fz += one_fourth_bfz + p_N0.N2 * p_t.fz;
			en3_fx += one_fourth_bfx + p_N0.N3 * p_t.fx;
			en3_fy += one_fourth_bfy + p_N0.N3 * p_t.fy;
			en3_fz += one_fourth_bfz + p_N0.N3 * p_t.fz;
			en4_fx += one_fourth_bfx + p_N0.N4 * p_t.fx;
			en4_fy += one_fourth_bfy + p_N0.N4 * p_t.fy;
			en4_fz += one_fourth_bfz + p_N0.N4 * p_t.fz;

			if (e_id != pcl_in_elem[p_id + 1])
			{
				elem_pcl_m[e_id] = e_p_m;
				elem_density[e_id] = e_p_m / e_p_vol;

				ElemNodeVM& en1_v = elem_node_vm[e_id * 4];
				en1_v.vm = en1_vm;
				en1_v.vmx = en1_vmx;
				en1_v.vmy = en1_vmy;
				en1_v.vmz = en1_vmz;
				ElemNodeVM& en2_v = elem_node_vm[e_id * 4 + 1];
				en2_v.vm = en2_vm;
				en2_v.vmx = en2_vmx;
				en2_v.vmy = en2_vmy;
				en2_v.vmz = en2_vmz;
				ElemNodeVM& en3_v = elem_node_vm[e_id * 4 + 2];
				en3_v.vm = en3_vm;
				en3_v.vmx = en3_vmx;
				en3_v.vmy = en3_vmy;
				en3_v.vmz = en3_vmz;
				ElemNodeVM& en4_v = elem_node_vm[e_id * 4 + 3];
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
				const DShapeFuncABC& e_dN = elem_dN_abc[e_id];
				// node 1
				Force& en1_f = elem_node_force[e_id * 4];
				en1_fx -= (e_dN.dN1_dx * e_s11 + e_dN.dN1_dy * e_s12 + e_dN.dN1_dz * e_s31) * e_p_vol;
				en1_f.fx = en1_fx;
				en1_fy -= (e_dN.dN1_dx * e_s12 + e_dN.dN1_dy * e_s22 + e_dN.dN1_dz * e_s23) * e_p_vol;
				en1_f.fy = en1_fy;
				en1_fz -= (e_dN.dN1_dx * e_s31 + e_dN.dN1_dy * e_s23 + e_dN.dN1_dz * e_s33) * e_p_vol;
				en1_f.fz = en1_fz;
				// node 2
				Force& en2_f = elem_node_force[e_id * 4 + 1];
				en2_fx -= (e_dN.dN2_dx * e_s11 + e_dN.dN2_dy * e_s12 + e_dN.dN2_dz * e_s31) * e_p_vol;
				en2_f.fx = en2_fx;
				en2_fy -= (e_dN.dN2_dx * e_s12 + e_dN.dN2_dy * e_s22 + e_dN.dN2_dz * e_s23) * e_p_vol;
				en2_f.fy = en2_fy;
				en2_fz -= (e_dN.dN2_dx * e_s31 + e_dN.dN2_dy * e_s23 + e_dN.dN2_dz * e_s33) * e_p_vol;
				en2_f.fz = en2_fz;
				// node 3
				Force& en3_f = elem_node_force[e_id * 4 + 2];
				en3_fx -= (e_dN.dN3_dx * e_s11 + e_dN.dN3_dy * e_s12 + e_dN.dN3_dz * e_s31) * e_p_vol;
				en3_f.fx = en3_fx;
				en3_fy -= (e_dN.dN3_dx * e_s12 + e_dN.dN3_dy * e_s22 + e_dN.dN3_dz * e_s23) * e_p_vol;
				en3_f.fy = en3_fy;
				en3_fz -= (e_dN.dN3_dx * e_s31 + e_dN.dN3_dy * e_s23 + e_dN.dN3_dz * e_s33) * e_p_vol;
				en3_f.fz = en3_fz;
				// node 4
				Force& en4_f = elem_node_force[e_id * 4 + 3];
				en4_fx -= (e_dN.dN4_dx * e_s11 + e_dN.dN4_dy * e_s12 + e_dN.dN4_dz * e_s31) * e_p_vol;
				en4_f.fx = en4_fx;
				en4_fy -= (e_dN.dN4_dx * e_s12 + e_dN.dN4_dy * e_s22 + e_dN.dN4_dz * e_s23) * e_p_vol;
				en4_f.fy = en4_fy;
				en4_fz -= (e_dN.dN4_dx * e_s31 + e_dN.dN4_dy * e_s23 + e_dN.dN4_dz * e_s33) * e_p_vol;
				en4_f.fz = en4_fz;

				e_id = pcl_in_elem[p_id + 1];
#ifdef _DEBUG
				assert(e_id < cd.elem_num || e_id == SIZE_MAX);
#endif

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
		double n_am = 0.0;
		double n_fx = 0.0;
		double n_fy = 0.0;
		double n_fz = 0.0;
		double n_vm = 0.0;
		double n_vmx = 0.0;
		double n_vmy = 0.0;
		double n_vmz = 0.0;
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
			n_am += elem_pcl_m[ne_id / 4];
			const Force& nf = elem_node_force[ne_id];
			n_fx += nf.fx;
			n_fy += nf.fy;
			n_fz += nf.fz;
			const ElemNodeVM& nvm = elem_node_vm[ne_id];
			n_vm += nvm.vm;
			n_vmx += nvm.vmx;
			n_vmy += nvm.vmy;
			n_vmz += nvm.vmz;

			if (n_id != node_has_elem[ve_id + 1])
			{
				Acceleration& n_a = node_a[n_id];
				n_am *= one_third;
				node_am[n_id] = n_am;
				n_a.ax = n_fx / n_am;
				n_a.ay = n_fy / n_am;
				n_a.az = n_fz / n_am;
				Velocity& n_v = node_v[n_id];
				n_v.vx = n_vmx / n_vm + n_a.ax * cd.dt;
				n_v.vy = n_vmy / n_vm + n_a.ay * cd.dt;
				n_v.vz = n_vmz / n_vm + n_a.az * cd.dt;
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

				n_id = node_has_elem[ve_id + 1];
#ifdef _DEBUG
				assert(n_id < cd.node_num || n_id == SIZE_MAX);
#endif

				n_am = 0.0;
				n_fx = 0.0;
				n_fy = 0.0;
				n_fz = 0.0;
				n_vm = 0.0;
				n_vmx = 0.0;
				n_vmy = 0.0;
				n_vmz = 0.0;
			}
		}
	}

	void CalElemDeAndMapToNode::operator() (size_t wk_id) const
	{
		const size_t ve_id0 = block_low(wk_id, task_num, valid_elem_num);
		const size_t ve_id1 = block_low(wk_id + 1, task_num, valid_elem_num);
		for (size_t ve_id = ve_id0; ve_id < ve_id1; ++ve_id)
		{
			const size_t e_id = valid_elems[ve_id];
#ifdef _DEBUG
			assert(e_id < cd.elem_num);
#endif

			const ElemNodeIndex& eni = elem_node_id[e_id];
			const Velocity& n_v1 = node_v[eni.n1];
			const Velocity& n_v2 = node_v[eni.n2];
			const Velocity& n_v3 = node_v[eni.n3];
			const Velocity& n_v4 = node_v[eni.n4];
			const DShapeFuncABC& e_dN = elem_dN_abc[e_id];
			StrainInc& e_de = elem_de[e_id];
			e_de.de11 = (e_dN.dN1_dx * n_v1.vx + e_dN.dN2_dx * n_v2.vx + e_dN.dN3_dx * n_v3.vx + e_dN.dN4_dx * n_v4.vx) * cd.dt;
			e_de.de22 = (e_dN.dN1_dy * n_v1.vy + e_dN.dN2_dy * n_v2.vy + e_dN.dN3_dy * n_v3.vy + e_dN.dN4_dy * n_v4.vy) * cd.dt;
			e_de.de33 = (e_dN.dN1_dz * n_v1.vz + e_dN.dN2_dz * n_v2.vz + e_dN.dN3_dz * n_v3.vz + e_dN.dN4_dz * n_v4.vz) * cd.dt;
			e_de.de12 = (e_dN.dN1_dx * n_v1.vy + e_dN.dN2_dx * n_v2.vy + e_dN.dN3_dx * n_v3.vy + e_dN.dN4_dx * n_v4.vy
					   + e_dN.dN1_dy * n_v1.vx + e_dN.dN2_dy * n_v2.vx + e_dN.dN3_dy * n_v3.vx + e_dN.dN4_dy * n_v4.vx) * cd.dt * 0.5;
			e_de.de23 = (e_dN.dN1_dy * n_v1.vz + e_dN.dN2_dy * n_v2.vz + e_dN.dN3_dy * n_v3.vz + e_dN.dN4_dy * n_v4.vz
					   + e_dN.dN1_dz * n_v1.vy + e_dN.dN2_dz * n_v2.vy + e_dN.dN3_dz * n_v3.vy + e_dN.dN4_dz * n_v4.vy) * cd.dt * 0.5;
			e_de.de31 = (e_dN.dN1_dz * n_v1.vx + e_dN.dN2_dz * n_v2.vx + e_dN.dN3_dz * n_v3.vx + e_dN.dN4_dz * n_v4.vx
					   + e_dN.dN1_dx * n_v1.vz + e_dN.dN2_dx * n_v2.vz + e_dN.dN3_dx * n_v3.vz + e_dN.dN4_dx * n_v4.vz) * cd.dt * 0.5;
			double e_de_vol = e_de.de11 + e_de.de22 + e_de.de33;
			elem_m_de_vol[e_id] = elem_pcl_m[e_id] * e_de_vol;
			e_de_vol *= one_third;
			e_de.de11 -= e_de_vol;
			e_de.de22 -= e_de_vol;
			e_de.de33 -= e_de_vol;
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

		double n_am_de_vol = 0.0;
		n_id = node_has_elem[ve_id0];
#ifdef _DEBUG
		assert(n_id < cd.node_num || n_id == SIZE_MAX);
#endif
		for (size_t ve_id = ve_id0; ve_id < ve_id1; ++ve_id)
		{
#ifdef _DEBUG
			assert(node_elem_pair[ve_id] / 4 < cd.elem_num);
#endif
			n_am_de_vol += elem_m_de_vol[node_elem_pair[ve_id] / 4];
			if (n_id != node_has_elem[ve_id + 1])
			{
				node_de_vol[n_id] = n_am_de_vol * one_fourth / node_am[n_id];
				n_id = node_has_elem[ve_id + 1];
#ifdef _DEBUG
				assert(n_id < cd.node_num || n_id == SIZE_MAX);
#endif
				n_am_de_vol = 0.0;
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
		
		const Acceleration *pn_a1, *pn_a2, *pn_a3, *pn_a4;
		const Velocity *pn_v1, *pn_v2, *pn_v3, *pn_v4;
		StrainInc* pe_de;
		double dstrain[6];
		const double *estrain, *pstrain, *dstress;
		size_t valid_pcl_num = 0;
		for (size_t p_id = p_id0; p_id < p_id1; ++p_id)
		{
			if (e_id != pcl_in_elem[p_id])
			{
				e_id = pcl_in_elem[p_id];
#ifdef _DEBUG
				assert(e_id < cd.elem_num);
#endif
				const ElemNodeIndex& eni = elem_node_id[e_id];
				pn_a1 = node_a + eni.n1;
				pn_a2 = node_a + eni.n2;
				pn_a3 = node_a + eni.n3;
				pn_a4 = node_a + eni.n4;
				pn_v1 = node_v + eni.n1;
				pn_v2 = node_v + eni.n2;
				pn_v3 = node_v + eni.n3;
				pn_v4 = node_v + eni.n4;

				double e_de_vol = (node_de_vol[eni.n1]
					+ node_de_vol[eni.n2] + node_de_vol[eni.n3]
					+ node_de_vol[eni.n4]) * one_fourth;
				elem_density[e_id] /= 1.0 + e_de_vol;
				
				pe_de = elem_de + e_id;
				e_de_vol *= one_third;
				pe_de->de11 += e_de_vol;
				pe_de->de22 += e_de_vol;
				pe_de->de33 += e_de_vol;
			}

			// update velocity
			ShapeFunc& p_N = pcl_N0[p_id];
			Velocity& p_v0 = pcl_v0[p_id];
			p_v0.vx += (p_N.N1 * pn_a1->ax + p_N.N2 * pn_a2->ax + p_N.N3 * pn_a3->ax + p_N.N4 * pn_a4->ax) * cd.dt;
			p_v0.vy += (p_N.N1 * pn_a1->ay + p_N.N2 * pn_a2->ay + p_N.N3 * pn_a3->ay + p_N.N4 * pn_a4->ay) * cd.dt;
			p_v0.vz += (p_N.N1 * pn_a1->az + p_N.N2 * pn_a2->az + p_N.N3 * pn_a3->az + p_N.N4 * pn_a4->az) * cd.dt;

			// update displacement
			Displacement& p_d0 = pcl_disp0[p_id];
			p_d0.ux += (p_N.N1 * pn_v1->vx + p_N.N2 * pn_v2->vx + p_N.N3 * pn_v3->vx + p_N.N4 * pn_v4->vx) * cd.dt;
			p_d0.uy += (p_N.N1 * pn_v1->vy + p_N.N2 * pn_v2->vy + p_N.N3 * pn_v3->vy + p_N.N4 * pn_v4->vy) * cd.dt;
			p_d0.uz += (p_N.N1 * pn_v1->vz + p_N.N2 * pn_v2->vz + p_N.N3 * pn_v3->vz + p_N.N4 * pn_v4->vz) * cd.dt;

			// update location (in which element)
			const size_t ori_p_id = pcl_index0[p_id];
#ifdef _DEBUG
			assert(ori_p_id < cd.ori_pcl_num);
#endif
			const Position& p_p = pcl_pos[ori_p_id];
			const double p_x = p_p.x + p_d0.ux;
			const double p_y = p_p.y + p_d0.uy;
			const double p_z = p_p.z + p_d0.uz;
			size_t p_e_id = e_id;
			Model_T3D_ME_mt& md = *cd.pmodel;
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
			if (p_e_id != SIZE_MAX)
				++valid_pcl_num;
			new_pcl_in_elem[p_id] = p_e_id;
			new_cur_to_prev_pcl[p_id] = p_id;
#ifdef _DEBUG
			assert(p_e_id < cd.elem_num || p_e_id == SIZE_MAX);
#endif

			// update density
			pcl_density0[p_id] = elem_density[e_id];

			// update stress
			MatModel::MaterialModel& pcl_mm = *pcl_mat_model[ori_p_id];
			dstrain[0] = pe_de->de11;
			dstrain[1] = pe_de->de22;
			dstrain[2] = pe_de->de33;
			dstrain[3] = pe_de->de12;
			dstrain[4] = pe_de->de23;
			dstrain[5] = pe_de->de31;
			pcl_mm.integrate(dstrain);
			dstress = pcl_mm.get_dstress();
			Stress& p_s = pcl_stress0[p_id];
			p_s.s11 += dstress[0];
			p_s.s22 += dstress[1];
			p_s.s33 += dstress[2];
			p_s.s12 += dstress[3];
			p_s.s23 += dstress[4];
			p_s.s31 += dstress[5];

			const size_t prev_p_id = cur_to_prev_pcl[p_id];
#ifdef _DEBUG
			assert(prev_p_id < cd.prev_valid_pcl_num);
#endif
			const Strain& p_e1 = pcl_strain1[prev_p_id];
			Strain& p_e0 = pcl_strain0[p_id];
			p_e0.e11 = p_e1.e11 + pe_de->de11;
			p_e0.e22 = p_e1.e22 + pe_de->de22;
			p_e0.e33 = p_e1.e33 + pe_de->de33;
			p_e0.e12 = p_e1.e12 + pe_de->de12;
			p_e0.e23 = p_e1.e23 + pe_de->de23;
			p_e0.e31 = p_e1.e31 + pe_de->de31;

			estrain = pcl_mm.get_dstrain_e();
			const Strain& p_ee1 = pcl_estrain1[prev_p_id];
			Strain& p_ee0 = pcl_estrain0[p_id];
			p_ee0.e11 = p_ee1.e11 + estrain[0];
			p_ee0.e22 = p_ee1.e22 + estrain[1];
			p_ee0.e33 = p_ee1.e33 + estrain[2];
			p_ee0.e12 = p_ee1.e12 + estrain[3];
			p_ee0.e23 = p_ee1.e23 + estrain[4];
			p_ee0.e31 = p_ee1.e31 + estrain[5];

			pstrain = pcl_mm.get_dstrain_p();
			const Strain& p_pe1 = pcl_pstrain1[prev_p_id];
			Strain& p_pe0 = pcl_pstrain0[p_id];
			p_pe0.e11 = p_pe1.e11 + pstrain[0];
			p_pe0.e22 = p_pe1.e22 + pstrain[1];
			p_pe0.e33 = p_pe1.e33 + pstrain[2];
			p_pe0.e12 = p_pe1.e12 + pstrain[3];
			p_pe0.e23 = p_pe1.e23 + pstrain[4];
			p_pe0.e31 = p_pe1.e31 + pstrain[5];
		}

		pcl_in_mesh_num = valid_pcl_num;
	}
}
