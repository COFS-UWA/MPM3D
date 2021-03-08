#include "SimulationsOMP_pcp.h"

#include <assert.h>
#include "tbb/task_arena.h"

#include "ParallelUtils.h"
#include "Step_T2D_ME_Task.h"
#include "Step_T2D_ME_TBB.h"

namespace Step_T2D_ME_Task
{
	constexpr double one_third = 1.0 / 3.0;

	void CalData::set_model(Model_T2D_ME_mt &md) noexcept
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
		elem_dN_ab = md.elem_dN_ab;
		elem_dN_c = md.elem_dN_c;
		elem_area = md.elem_area;
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
		Model_T2D_ME_mt& md = *cd.pmodel;
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
			p_d.ux = 0.0;
			p_d.uy = 0.0;
			ShapeFunc& p_N = spva0.pcl_N[p_id];
			size_t e_id = md.find_pcl_in_which_elem(p_p.x, p_p.y, p_N);
			if (e_id == SIZE_MAX)
				e_id = md.find_pcl_in_which_elem_tol(p_p.x, p_p.y, p_N);
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
		double en1_fx = 0.0;
		double en1_fy = 0.0;
		double en2_fx = 0.0;
		double en2_fy = 0.0;
		double en3_fx = 0.0;
		double en3_fy = 0.0;
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
			p_s0.s12 = p_s1.s12;
			e_s11 += p_s0.s11 * p_vol;
			e_s22 += p_s0.s22 * p_vol;
			e_s12 += p_s0.s12 * p_vol;

			// shape function
			const ShapeFunc& p_N1 = pcl_N1[prev_p_id];
			ShapeFunc& p_N0 = pcl_N0[p_id];
			p_N0.N1 = p_N1.N1;
			p_N0.N2 = p_N1.N2;
			p_N0.N3 = p_N1.N3;
			// map velocity
			double p_N_m;
			const Velocity& p_v1 = pcl_v1[prev_p_id];
			Velocity& p_v0 = pcl_v0[p_id];
			p_v0.vx = p_v1.vx;
			p_v0.vy = p_v1.vy;
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

			// displacement (for contact)
			const Displacement& p_d1 = pcl_disp1[prev_p_id];
			Displacement& p_d0 = pcl_disp0[p_id];
			p_d0.ux = p_d1.ux;
			p_d0.uy = p_d1.uy;

			// cal external load
			const Force& p_bf = pcl_bf[ori_p_id];
			const double one_third_bfx = one_third * p_bf.fx;
			const double one_third_bfy = one_third * p_bf.fy;
			const Force& p_t = pcl_t[ori_p_id];
			en1_fx += one_third_bfx + p_N0.N1 * p_t.fx;
			en1_fy += one_third_bfy + p_N0.N1 * p_t.fy;
			en2_fx += one_third_bfx + p_N0.N2 * p_t.fx;
			en2_fy += one_third_bfy + p_N0.N2 * p_t.fy;
			en3_fx += one_third_bfx + p_N0.N3 * p_t.fx;
			en3_fy += one_third_bfy + p_N0.N3 * p_t.fy;

			if (e_id != pcl_in_elem[p_id + 1])
			{
				elem_pcl_m[e_id] = e_p_m;
				elem_density[e_id] = e_p_m / e_p_vol;

				ElemNodeVM& en1_v = elem_node_vm[e_id * 3];
				en1_v.vm = en1_vm;
				en1_v.vmx = en1_vmx;
				en1_v.vmy = en1_vmy;
				ElemNodeVM& en2_v = elem_node_vm[e_id * 3 + 1];
				en2_v.vm = en2_vm;
				en2_v.vmx = en2_vmx;
				en2_v.vmy = en2_vmy;
				ElemNodeVM& en3_v = elem_node_vm[e_id * 3 + 2];
				en3_v.vm = en3_vm;
				en3_v.vmx = en3_vmx;
				en3_v.vmy = en3_vmy;

				e_s11 /= e_p_vol;
				e_s22 /= e_p_vol;
				e_s12 /= e_p_vol;
				if (e_p_vol > elem_area[e_id])
					e_p_vol = elem_area[e_id];
				const ShapeFuncAB& e_dN = elem_dN_ab[e_id];
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

				e_id = pcl_in_elem[p_id + 1];
#ifdef _DEBUG
				assert(e_id < cd.elem_num || e_id == SIZE_MAX);
#endif

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
	}

	void UpdateAccelerationAndVelocity::operator() (size_t wk_id) const
	{
		size_t n_id;
		size_t ve_id0 = block_low(wk_id, task_num, three_valid_elem_num);
		n_id = node_has_elem[ve_id0];
		while (ve_id0 != SIZE_MAX && n_id == node_has_elem[--ve_id0]);
		++ve_id0;
		assert(ve_id0 <= three_valid_elem_num);
		size_t ve_id1 = block_low(wk_id + 1, task_num, three_valid_elem_num);
		n_id = node_has_elem[ve_id1];
		while (ve_id1 != SIZE_MAX && n_id == node_has_elem[--ve_id1]);
		++ve_id1;
		assert(ve_id1 <= three_valid_elem_num);
		
		size_t bc_mask, ne_id;
		double n_am = 0.0;
		double n_fx = 0.0;
		double n_fy = 0.0;
		double n_vm = 0.0;
		double n_vmx = 0.0;
		double n_vmy = 0.0;
		n_id = node_has_elem[ve_id0];
#ifdef _DEBUG
		assert(n_id < cd.node_num);
#endif
		for (size_t ve_id = ve_id0; ve_id < ve_id1; ++ve_id)
		{
			ne_id = node_elem_pair[ve_id];
#ifdef _DEBUG
			assert(ne_id < cd.elem_num * 3);
#endif
			n_am += elem_pcl_m[ne_id / 3];
			const Force& nf = elem_node_force[ne_id];
			n_fx += nf.fx;
			n_fy += nf.fy;
			const ElemNodeVM& nvm = elem_node_vm[ne_id];
			n_vm += nvm.vm;
			n_vmx += nvm.vmx;
			n_vmy += nvm.vmy;

			if (n_id != node_has_elem[ve_id + 1])
			{
				Acceleration& n_a = node_a[n_id];
				n_am *= one_third;
				node_am[n_id] = n_am;
				n_a.ax = n_fx / n_am;
				n_a.ay = n_fy / n_am;
				Velocity& n_v = node_v[n_id];
				n_v.vx = n_vmx / n_vm + n_a.ax * cd.dt;
				n_v.vy = n_vmy / n_vm + n_a.ay * cd.dt;
				NodeHasVBC& n_has_vbc = node_has_vbc[n_id];
				bc_mask = size_t(n_has_vbc.has_vx_bc) + SIZE_MAX;
				n_a.iax &= bc_mask;
				n_v.ivx &= bc_mask;
				bc_mask = size_t(n_has_vbc.has_vy_bc) + SIZE_MAX;
				n_a.iay &= bc_mask;
				n_v.ivy &= bc_mask;

				n_id = node_has_elem[ve_id + 1];
#ifdef _DEBUG
				assert(n_id < cd.node_num || n_id == SIZE_MAX);
#endif

				n_am = 0.0;
				n_fx = 0.0;
				n_fy = 0.0;
				n_vm = 0.0;
				n_vmx = 0.0;
				n_vmy = 0.0;
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
			const ShapeFuncAB& e_dN = elem_dN_ab[e_id];
			StrainInc& e_de = elem_de[e_id];
			e_de.de11 = (e_dN.dN1_dx * n_v1.vx + e_dN.dN2_dx * n_v2.vx + e_dN.dN3_dx * n_v3.vx) * cd.dt;
			e_de.de22 = (e_dN.dN1_dy * n_v1.vy + e_dN.dN2_dy * n_v2.vy + e_dN.dN3_dy * n_v3.vy) * cd.dt;
			e_de.de12 = (e_dN.dN1_dx * n_v1.vy + e_dN.dN2_dx * n_v2.vy + e_dN.dN3_dx * n_v3.vy
					   + e_dN.dN1_dy * n_v1.vx + e_dN.dN2_dy * n_v2.vx + e_dN.dN3_dy * n_v3.vx) * cd.dt * 0.5;
			double e_de_vol = e_de.de11 + e_de.de22;
			elem_m_de_vol[e_id] = elem_pcl_m[e_id] * e_de_vol;
			e_de_vol *= one_third;
			e_de.de11 -= e_de_vol;
			e_de.de22 -= e_de_vol;
		}
	}
	
	void CalNodeDe::operator() (size_t wk_id) const
	{
		size_t n_id;
		size_t ve_id0 = block_low(wk_id, task_num, three_valid_elem_num);
		n_id = node_has_elem[ve_id0];
		while (ve_id0 != SIZE_MAX && n_id == node_has_elem[--ve_id0]);
		++ve_id0;
		assert(ve_id0 <= three_valid_elem_num);
		size_t ve_id1 = block_low(wk_id + 1, task_num, three_valid_elem_num);
		n_id = node_has_elem[ve_id1];
		while (ve_id1 != SIZE_MAX && n_id == node_has_elem[--ve_id1]);
		++ve_id1;
		assert(ve_id1 <= three_valid_elem_num);

		double n_am_de_vol = 0.0;
		n_id = node_has_elem[ve_id0];
#ifdef _DEBUG
		assert(n_id < cd.node_num || n_id == SIZE_MAX);
#endif
		for (size_t ve_id = ve_id0; ve_id < ve_id1; ++ve_id)
		{
#ifdef _DEBUG
			assert(node_elem_pair[ve_id] / 3 < cd.elem_num);
#endif
			n_am_de_vol += elem_m_de_vol[node_elem_pair[ve_id] / 3];
			if (n_id != node_has_elem[ve_id + 1])
			{
				node_de_vol[n_id] = n_am_de_vol * one_third / node_am[n_id];
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
		
		const Acceleration* pn_a1, * pn_a2, * pn_a3;
		const Velocity* pn_v1, * pn_v2, * pn_v3;
		StrainInc* pe_de;
		double dstrain[6];
		dstrain[2] = 0.0;
		dstrain[4] = 0.0;
		dstrain[5] = 0.0;
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
				pn_v1 = node_v + eni.n1;
				pn_v2 = node_v + eni.n2;
				pn_v3 = node_v + eni.n3;

				double e_de_vol = one_third * (node_de_vol[eni.n1]
					+ node_de_vol[eni.n2] + node_de_vol[eni.n3]);
				elem_density[e_id] /= 1.0 + e_de_vol;
				
				pe_de = elem_de + e_id;
				e_de_vol *= one_third;
				pe_de->de11 += e_de_vol;
				pe_de->de22 += e_de_vol;
			}

			// update velocity
			ShapeFunc& p_N = pcl_N0[p_id];
			Velocity& p_v0 = pcl_v0[p_id];
			p_v0.vx += (p_N.N1 * pn_a1->ax + p_N.N2 * pn_a2->ax + p_N.N3 * pn_a3->ax) * cd.dt;
			p_v0.vy += (p_N.N1 * pn_a1->ay + p_N.N2 * pn_a2->ay + p_N.N3 * pn_a3->ay) * cd.dt;

			// update displacement
			Displacement& p_d0 = pcl_disp0[p_id];
			p_d0.ux += (p_N.N1 * pn_v1->vx + p_N.N2 * pn_v2->vx + p_N.N3 * pn_v3->vx) * cd.dt;
			p_d0.uy += (p_N.N1 * pn_v1->vy + p_N.N2 * pn_v2->vy + p_N.N3 * pn_v3->vy) * cd.dt;

			// update location (in which element)
			const size_t ori_p_id = pcl_index0[p_id];
#ifdef _DEBUG
			assert(ori_p_id < cd.ori_pcl_num);
#endif
			const Position& p_p = pcl_pos[ori_p_id];
			const double p_x = p_p.x + p_d0.ux;
			const double p_y = p_p.y + p_d0.uy;
			size_t p_e_id = e_id;
			Model_T2D_ME_mt& md = *cd.pmodel;
			if (!md.is_in_element(p_x, p_y, e_id, p_N))
			{
				p_e_id = md.find_pcl_in_which_elem(p_x, p_y, p_N);
				if (p_e_id == SIZE_MAX)
				{
					if (md.is_in_element_tol(p_x, p_y, e_id, p_N))
						p_e_id = e_id;
					else
						p_e_id = md.find_pcl_in_which_elem_tol(p_x, p_y, p_N);
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
			dstrain[3] = pe_de->de12;
			pcl_mm.integrate(dstrain);
			dstress = pcl_mm.get_dstress();
			Stress& p_s = pcl_stress0[p_id];
			p_s.s11 += dstress[0];
			p_s.s22 += dstress[1];
			p_s.s12 += dstress[3];

			const size_t prev_p_id = cur_to_prev_pcl[p_id];
#ifdef _DEBUG
			assert(prev_p_id < cd.prev_valid_pcl_num);
#endif
			const Strain& p_e1 = pcl_strain1[prev_p_id];
			Strain& p_e0 = pcl_strain0[p_id];
			p_e0.e11 = p_e1.e11 + pe_de->de11;
			p_e0.e22 = p_e1.e22 + pe_de->de22;
			p_e0.e12 = p_e1.e12 + pe_de->de12;

			estrain = pcl_mm.get_dstrain_e();
			const Strain& p_ee1 = pcl_estrain1[prev_p_id];
			Strain& p_ee0 = pcl_estrain0[p_id];
			p_ee0.e11 = p_ee1.e11 + estrain[0];
			p_ee0.e22 = p_ee1.e22 + estrain[1];
			p_ee0.e12 = p_ee1.e12 + estrain[3];

			pstrain = pcl_mm.get_dstrain_p();
			const Strain& p_pe1 = pcl_pstrain1[prev_p_id];
			Strain& p_pe0 = pcl_pstrain0[p_id];
			p_pe0.e11 = p_pe1.e11 + pstrain[0];
			p_pe0.e22 = p_pe1.e22 + pstrain[1];
			p_pe0.e12 = p_pe1.e12 + pstrain[3];
		}

		pcl_in_mesh_num = valid_pcl_num;
	}
}
