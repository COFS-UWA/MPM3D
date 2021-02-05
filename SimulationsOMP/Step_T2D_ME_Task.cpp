#include "SimulationsOMP_pcp.h"

#include <assert.h>
#include "tbb/task_arena.h"

#include "SortParticleTask.h"
#include "SortTriMeshNodeTask.hpp"
#include "Step_T2D_ME_Task.h"
#include "Step_T2D_ME_TBB.h"

#define one_third (1.0/3.0)
#define Block_Low(block_id, block_num, data_num) ((block_id) * (data_num) / (block_num))

namespace Step_T2D_ME_Task
{
	TaskData::TaskData(Step_T2D_ME_TBB& _stp,
		size_t _pcl_num_per_map_pcl_to_mesh_task,
		size_t _node_num_per_update_a_and_v_task,
		size_t _elem_num_per_cal_elem_de_task,
		size_t _node_num_per_cal_node_de_task,
		size_t _pcl_num_per_task_map_mesh_to_pcl) :
		stp(_stp),
		md(*(Model_T2D_ME_mt*)(stp.model)),
		pcl_sort_mem(stp.pcl_sort_mem),
		node_sort_mem(stp.node_sort_mem),
		pcl_m(md.pcl_m),
		pcl_bf(md.pcl_bf),
		pcl_t(md.pcl_t),
		pcl_pos(md.pcl_pos),
		pcl_vol(md.pcl_vol),
		pcl_mat_model(md.pcl_mat_model),
		elem_node_id(md.elem_node_id),
		elem_dN_ab(md.elem_dN_ab),
		elem_dN_c(md.elem_dN_c),
		elem_area(md.elem_area),
		elem_pcl_m(md.elem_pcl_m),
		elem_density(md.elem_density),
		elem_de(md.elem_de),
		elem_m_de_vol(md.elem_m_de_vol),
		valid_elems(stp.valid_elems),
		tmp_valid_elems(stp.tmp_valid_elems),
		elem_node_vm(md.elem_node_vm),
		elem_node_force(md.elem_node_force),
		node_a(md.node_a),
		node_v(md.node_v),
		node_has_vbc(md.node_has_vbc),
		node_am(md.node_am),
		node_de_vol(md.node_de_vol),
		pcl_digit_num(SortUtils::max_digit_num(md.ori_pcl_num)),
		node_digit_num(SortUtils::max_digit_num(md.node_num)),
#ifdef _DEBUG
		ori_pcl_num(md.ori_pcl_num),
		elem_num(md.elem_num),
		node_num(md.node_num),
#endif
		pcl_num_per_map_pcl_to_mesh_task(_pcl_num_per_map_pcl_to_mesh_task),
		node_num_per_update_a_and_v_task(_node_num_per_update_a_and_v_task),
		elem_num_per_cal_elem_de_task(_elem_num_per_cal_elem_de_task),
		node_num_per_cal_node_de_task(_node_num_per_cal_node_de_task),
		pcl_num_per_task_map_mesh_to_pcl(_pcl_num_per_task_map_mesh_to_pcl)
	{
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
	}
	
	tbb::task* MapPclToBgMeshTask::execute()
	{
		size_t* const pcl_in_elem = td.pcl_sort_mem.out_keys;
		size_t e_id;
		e_id = pcl_in_elem[p_id0];
		while (p_id0 != SIZE_MAX && e_id == pcl_in_elem[--p_id0]);
		++p_id0;
		assert(p_id0 <= td.valid_pcl_num);
		e_id = pcl_in_elem[p_id1];
		while (p_id1 != SIZE_MAX && e_id == pcl_in_elem[--p_id1]);
		++p_id1;
		assert(p_id1 <= td.valid_pcl_num);
		
		size_t* const prev_pcl_id = td.pcl_sort_mem.out_vals;
		const double* const pcl_m = td.pcl_m;
		const Force* pcl_bf = td.pcl_bf;
		const Force* pcl_t = td.pcl_t;
		double* const pcl_vol = td.pcl_vol;
		const auto& spva0 = td.spvas[td.sorted_pcl_var_id];
		const auto& spva1 = td.spvas[td.sorted_pcl_var_id ^ 1];
		size_t* const pcl_index0 = spva0.pcl_index;
		double* const pcl_density0 = spva0.pcl_density;
		Velocity* const pcl_v0 = spva0.pcl_v;
		Displacement* const pcl_disp0 = spva0.pcl_disp;
		Stress* const pcl_stress0 = spva0.pcl_stress;
		ShapeFunc* const pcl_N0 = spva0.pcl_N;
		const size_t* const pcl_index1 = spva1.pcl_index;
		const double* const pcl_density1 = spva1.pcl_density;
		const Velocity* const pcl_v1 = spva1.pcl_v;
		const Displacement* const pcl_disp1 = spva1.pcl_disp;
		const Stress* const pcl_stress1 = spva1.pcl_stress;
		const ShapeFunc* const pcl_N1 = spva1.pcl_N;
		const ShapeFuncAB* const elem_dN_ab = td.elem_dN_ab;
		const double* const elem_area = td.elem_area;
		double* const elem_pcl_m = td.elem_pcl_m;
		double* const elem_density = td.elem_density;
		ElemNodeVM *const elem_node_vm = td.elem_node_vm;
		Force *const elem_node_force = td.elem_node_force;
		size_t* const valid_elems = td.valid_elems;
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
		size_t valid_elem_num = 0;
		e_id = pcl_in_elem[p_id0];
		assert(e_id < td.elem_num);
		for (size_t p_id = p_id0; p_id < p_id1; ++p_id)
		{
			// pcl index
			const size_t prev_p_id = prev_pcl_id[p_id];
			assert(prev_p_id < td.prev_valid_pcl_num);
			const size_t ori_p_id = pcl_index1[prev_p_id];
			assert(ori_p_id < td.ori_pcl_num);
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
#define N_tol (1.0e-10)
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

				valid_elems[valid_elem_num++] = e_id;
				
				e_id = pcl_in_elem[p_id + 1];
				assert(e_id < td.elem_num || e_id == SIZE_MAX);

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
		return nullptr;
	}

	tbb::task* UpdateAccelerationAndVelocityTask::execute()
	{
		const size_t* const node_has_elem = td.node_sort_mem.out_keys;
		size_t n_id;
		n_id = node_has_elem[ve_id0];
		while (ve_id0 != SIZE_MAX && n_id == node_has_elem[--ve_id0]);
		++ve_id0;
		assert(ve_id0 <= td.valid_elem_num * 3);
		n_id = node_has_elem[ve_id1];
		while (ve_id1 != SIZE_MAX && n_id == node_has_elem[--ve_id1]);
		++ve_id1;
		assert(ve_id1 <= td.valid_elem_num * 3);
		
		const size_t* const node_elem_pair = td.node_sort_mem.out_vals;
		const double *const elem_pcl_m = td.elem_pcl_m;
		const Force *const elem_node_force = td.elem_node_force;
		const ElemNodeVM *const elem_node_vm = td.elem_node_vm;
		Acceleration *const node_a = td.node_a;
		double *const node_am = td.node_am;
		Velocity *const node_v = td.node_v;
		NodeHasVBC* const node_has_vbc = td.node_has_vbc;
		size_t bc_mask, ne_id, e_id;
		double n_am = 0.0;
		double n_fx = 0.0;
		double n_fy = 0.0;
		double n_vm = 0.0;
		double n_vmx = 0.0;
		double n_vmy = 0.0;
		n_id = node_has_elem[ve_id0];
		assert(n_id < td.node_num);
		for (size_t ve_id = ve_id0; ve_id < ve_id1; ++ve_id)
		{
			ne_id = node_elem_pair[ve_id];
			assert(ne_id < td.elem_num * 3);
			e_id = ne_id / 3;
			n_am += elem_pcl_m[e_id];
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
				n_v.vx = n_vmx / n_vm + n_a.ax * td.dt;
				n_v.vy = n_vmy / n_vm + n_a.ay *td. dt;
				NodeHasVBC& n_has_vbc = node_has_vbc[n_id];
				bc_mask = size_t(n_has_vbc.has_vx_bc) + SIZE_MAX;
				n_a.iax &= bc_mask;
				n_v.ivx &= bc_mask;
				bc_mask = size_t(n_has_vbc.has_vy_bc) + SIZE_MAX;
				n_a.iay &= bc_mask;
				n_v.ivy &= bc_mask;

				n_id = node_has_elem[ve_id + 1];
				assert(n_id < td.node_num || n_id == SIZE_MAX);

				n_am = 0.0;
				n_fx = 0.0;
				n_fy = 0.0;
				n_vm = 0.0;
				n_vmx = 0.0;
				n_vmy = 0.0;
			}
		}
		return nullptr;
	}

	tbb::task* CalElemDeAndMapToNode::execute()
	{
		const size_t* const valid_elems = td.valid_elems;
		const ElemNodeIndex* const elem_node_id = td.elem_node_id;
		const Velocity* const node_v = td.node_v;
		const ShapeFuncAB* const elem_dN_ab = td.elem_dN_ab;
		StrainInc* const elem_de = td.elem_de;
		const double *const elem_pcl_m = td.elem_pcl_m;
		double *const elem_m_de_vol = td.elem_m_de_vol;
		for (size_t ve_id = ve_id0; ve_id < ve_id1; ++ve_id)
		{
			const size_t e_id = valid_elems[ve_id];
			assert(e_id < td.elem_num);

			const ElemNodeIndex& eni = elem_node_id[e_id];
			const Velocity& n_v1 = node_v[eni.n1];
			const Velocity& n_v2 = node_v[eni.n2];
			const Velocity& n_v3 = node_v[eni.n3];
			const ShapeFuncAB& e_dN = elem_dN_ab[e_id];
			StrainInc& e_de = elem_de[e_id];
			e_de.de11 = (e_dN.dN1_dx * n_v1.vx + e_dN.dN2_dx * n_v2.vx + e_dN.dN3_dx * n_v3.vx) * td.dt;
			e_de.de22 = (e_dN.dN1_dy * n_v1.vy + e_dN.dN2_dy * n_v2.vy + e_dN.dN3_dy * n_v3.vy) * td.dt;
			e_de.de12 = (e_dN.dN1_dx * n_v1.vy + e_dN.dN2_dx * n_v2.vy + e_dN.dN3_dx * n_v3.vy
					   + e_dN.dN1_dy * n_v1.vx + e_dN.dN2_dy * n_v2.vx + e_dN.dN3_dy * n_v3.vx) * td.dt * 0.5;
			double e_de_vol = e_de.de11 + e_de.de22;
			elem_m_de_vol[e_id] = elem_pcl_m[e_id] * e_de_vol;
			e_de_vol *= one_third;
			e_de.de11 -= e_de_vol;
			e_de.de22 -= e_de_vol;
		}
		return nullptr;
	}
	
	tbb::task* CalNodeDe::execute()
	{
		const size_t* const node_has_elem = td.node_sort_mem.out_keys;
		size_t n_id;
		n_id = node_has_elem[ve_id0];
		while (ve_id0 != SIZE_MAX && n_id == node_has_elem[--ve_id0]);
		++ve_id0;
		assert(ve_id0 <= td.valid_elem_num * 3);
		n_id = node_has_elem[ve_id1];
		while (ve_id1 != SIZE_MAX && n_id == node_has_elem[--ve_id1]);
		++ve_id1;
		assert(ve_id1 <= td.valid_elem_num * 3);

		const size_t* const node_elem_pair = td.node_sort_mem.out_vals;
		const double* const elem_m_de_vol = td.elem_m_de_vol;
		const double* const node_am = td.node_am;
		double* const node_de_vol = td.node_de_vol;
		double n_am_de_vol = 0.0;
		n_id = node_has_elem[ve_id0];
		assert(n_id < td.node_num);
		for (size_t ve_id = ve_id0; ve_id < ve_id1; ++ve_id)
		{
			const size_t e_id = node_elem_pair[ve_id] / 3;
			assert(e_id < td.elem_num);
			n_am_de_vol += elem_m_de_vol[e_id];
			if (n_id != node_has_elem[ve_id + 1])
			{
				node_de_vol[n_id] = n_am_de_vol * one_third / node_am[n_id];
				n_id = node_has_elem[ve_id + 1];
				assert(n_id < td.node_num || n_id == SIZE_MAX);
				n_am_de_vol = 0.0;
			}
		}
		return nullptr;
	}

	tbb::task* MapBgMeshToPclTask::execute()
	{
		size_t* const pcl_in_elem = td.pcl_sort_mem.out_keys;
		size_t e_id;
		e_id = pcl_in_elem[p_id0];
		while (p_id0 != SIZE_MAX && e_id == pcl_in_elem[--p_id0]);
		++p_id0;
		assert(p_id0 <= td.valid_pcl_num);
		e_id = pcl_in_elem[p_id1];
		while (p_id1 != SIZE_MAX && e_id == pcl_in_elem[--p_id1]);
		++p_id1;
		assert(p_id1 <= td.valid_pcl_num);

		const size_t* const prev_pcl_id = td.pcl_sort_mem.out_vals;
		const ElemNodeIndex* const elem_node_id = td.elem_node_id;
		const Acceleration *const node_a = td.node_a;
		const Velocity* const node_v = td.node_v;
		double* const elem_density = td.elem_density;
		StrainInc* const elem_de = td.elem_de;
		const double* const node_de_vol = td.node_de_vol;
		MatModel::MaterialModel** pcl_mat_model = td.pcl_mat_model;
		const auto &spva0 = td.spvas[td.sorted_pcl_var_id];
		const auto& spva1 = td.spvas[td.sorted_pcl_var_id ^ 1];
		const size_t* const pcl_index0 = spva0.pcl_index;
		const ShapeFunc *const pcl_N0 = spva0.pcl_N;
		double* const pcl_density0 = spva0.pcl_density;
		Velocity *const pcl_v0 = spva0.pcl_v;
		Displacement* const pcl_disp0 = spva0.pcl_disp;
		Stress* const pcl_stress0 = spva0.pcl_stress;
		Strain* const pcl_strain0 = spva0.pcl_strain;
		Strain* const pcl_estrain0 = spva0.pcl_estrain;
		Strain* const pcl_pstrain0 = spva0.pcl_pstrain;
		Strain* const pcl_strain1 = spva1.pcl_strain;
		Strain* const pcl_estrain1 = spva1.pcl_estrain;
		Strain* const pcl_pstrain1 = spva1.pcl_pstrain;
		const Acceleration* pn_a1, * pn_a2, * pn_a3;
		const Velocity* pn_v1, * pn_v2, * pn_v3;
		StrainInc* pe_de;
		double dstrain[6];
		dstrain[2] = 0.0;
		dstrain[4] = 0.0;
		dstrain[5] = 0.0;
		const double *estrain, *pstrain, *dstress;
		for (size_t p_id = p_id0; p_id < p_id1; ++p_id)
		{
			if (e_id != pcl_in_elem[p_id])
			{
				e_id = pcl_in_elem[p_id];
#ifdef _DEBUG
				assert(e_id < td.elem_num);
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
			const ShapeFunc& p_N = pcl_N0[p_id];
			Velocity& p_v0 = pcl_v0[p_id];
			p_v0.vx += (p_N.N1 * pn_a1->ax + p_N.N2 * pn_a2->ax + p_N.N3 * pn_a3->ax) * td.dt;
			p_v0.vy += (p_N.N1 * pn_a1->ay + p_N.N2 * pn_a2->ay + p_N.N3 * pn_a3->ay) * td.dt;

			// update displacement
			Displacement& p_d0 = pcl_disp0[p_id];
			p_d0.ux += (p_N.N1 * pn_v1->vx + p_N.N2 * pn_v2->vx + p_N.N3 * pn_v3->vx) * td.dt;
			p_d0.uy += (p_N.N1 * pn_v1->vy + p_N.N2 * pn_v2->vy + p_N.N3 * pn_v3->vy) * td.dt;

			// update density
			pcl_density0[p_id] = elem_density[e_id];

			// update stress
			const size_t ori_p_id = pcl_index0[p_id];
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

			const size_t prev_p_id = prev_pcl_id[p_id];
#ifdef _DEBUG
			assert(prev_p_id < td.prev_valid_pcl_num);
#endif
			const Strain& p_e1 = pcl_strain1[prev_p_id];
			Strain& p_e0 = pcl_strain0[p_id];
			p_e0.e11 = p_e1.e11 + pe_de->de11;
			p_e0.e22 = p_e1.e22 + pe_de->de22;
			p_e0.e12 = p_e1.e12 + pe_de->de12;

			estrain = pcl_mm.get_dstrain_e();
			Strain& p_ee1 = pcl_estrain1[prev_p_id];
			Strain& p_ee0 = pcl_estrain0[p_id];
			p_ee0.e11 = p_ee1.e11 + estrain[0];
			p_ee0.e22 = p_ee1.e22 + estrain[1];
			p_ee0.e12 = p_ee1.e12 + estrain[3];

			pstrain = pcl_mm.get_dstrain_p();
			Strain& p_pe1 = pcl_pstrain1[prev_p_id];
			Strain& p_pe0 = pcl_pstrain0[p_id];
			p_pe0.e11 = p_pe1.e11 + pstrain[0];
			p_pe0.e22 = p_pe1.e22 + pstrain[1];
			p_pe0.e12 = p_pe1.e12 + pstrain[3];
		}
		return nullptr;
	}

	tbb::task *Step_T2D_ME_Task::execute()
	{
		// sort pcl id
		td.prev_valid_pcl_num = td.valid_pcl_num;
		set_ref_count(2);
		spawn_and_wait_for_all(*new(allocate_child())
			SortUtils::SortParticleTask(
				td.pcl_sort_mem,
				td.valid_pcl_num,
				td.pcl_digit_num));
		if (td.valid_pcl_num == 0)
			return nullptr;

		size_t task_num, task_id, start_id, end_id;
		// sort node id
		task_num = (td.valid_pcl_num + td.pcl_num_per_map_pcl_to_mesh_task - 1)
				 / td.pcl_num_per_map_pcl_to_mesh_task;
		set_ref_count(2 + task_num);
		spawn(*new(allocate_child())
			SortUtils::SortTriMeshNodeTask<Model_T2D_ME_mt::ElemNodeIndex>(
				td.elem_num,
				td.pcl_sort_mem.out_keys,
				td.valid_pcl_num,
				td.node_digit_num,
				td.pcl_sort_mem,
				td.tmp_valid_elems,
				td.valid_elems));
		// map pcl to bg mesh
		start_id = 0;
		for (task_id = 1; task_id < task_num; ++task_id)
		{
			end_id = Block_Low(task_id, task_num, td.valid_pcl_num);
			spawn(*new(allocate_child()) MapPclToBgMeshTask(start_id, end_id ,td));
			start_id = end_id;
		}
		spawn_and_wait_for_all(*new(allocate_child()) MapPclToBgMeshTask(start_id, td.valid_pcl_num, td));

		// cal node acceleration and velocity
		const size_t three_valid_elem_num = td.valid_elem_num * 3;
		task_num = (three_valid_elem_num + td.node_num_per_update_a_and_v_task - 1)
				  / td.node_num_per_update_a_and_v_task;
		set_ref_count(1 + task_num);
		start_id = 0;
		for (task_id = 1; task_id < task_num; ++task_id)
		{
			end_id = Block_Low(task_id, task_num, three_valid_elem_num);
			spawn(*new(allocate_child()) UpdateAccelerationAndVelocityTask(start_id, end_id, td));
			start_id = end_id;
		}
		spawn_and_wait_for_all(*new(allocate_child()) MapPclToBgMeshTask(start_id, three_valid_elem_num, td));

		// cal element strain increment and map to node
		task_num = (td.valid_elem_num + td.elem_num_per_cal_elem_de_task - 1)
				/ td.elem_num_per_cal_elem_de_task;
		set_ref_count(1 + task_num);
		start_id = 0;
		for (task_id = 1; task_id < task_num; ++task_id)
		{
			end_id = Block_Low(task_id, task_num, td.valid_elem_num);
			spawn(*new(allocate_child()) CalElemDeAndMapToNode(start_id, end_id, td));
			start_id = end_id;
		}
		spawn_and_wait_for_all(*new(allocate_child()) CalElemDeAndMapToNode(start_id, td.valid_elem_num, td));

		// cal strain increment at node
		task_num = (three_valid_elem_num + td.node_num_per_cal_node_de_task - 1)
				/ td.node_num_per_cal_node_de_task;
		set_ref_count(1 + task_num);
		start_id = 0;
		for (task_id = 1; task_id < task_num; ++task_id)
		{
			end_id = Block_Low(task_id, task_num, three_valid_elem_num);
			spawn(*new(allocate_child()) CalNodeDe(start_id, end_id, td));
			start_id = end_id;
		}
		spawn_and_wait_for_all(*new(allocate_child()) CalNodeDe(start_id, three_valid_elem_num, td));

		// map bg mesh back to pcl
		task_num = (td.valid_pcl_num + td.pcl_num_per_task_map_mesh_to_pcl - 1)
				/ td.pcl_num_per_task_map_mesh_to_pcl;
		set_ref_count(1 + task_num);
		start_id = 0;
		for (task_id = 1; task_id < task_num; ++task_id)
		{
			end_id = Block_Low(task_id, task_num, td.valid_pcl_num);
			spawn(*new(allocate_child()) MapBgMeshToPclTask(start_id, end_id, td));
			start_id = end_id;
		}
		spawn_and_wait_for_all(*new(allocate_child()) MapBgMeshToPclTask(start_id, td.valid_pcl_num, td));
		
		return nullptr;
	}
}
