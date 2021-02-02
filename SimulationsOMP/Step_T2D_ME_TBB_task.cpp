#include "SimulationsOMP_pcp.h"

#include <assert.h>
#include "tbb/task_arena.h"

#include "Step_T2D_ME_TBB_task.h"
#include "Step_T2D_ME_TBB.h"

#define one_third (1.0/3.0)
#define Block_Low(block_id, block_num, data_num) ((block_id) * (data_num) / (block_num))

namespace Step_T2D_ME_TBB_task
{
	CalData_T2D_ME_TBB::CalData_T2D_ME_TBB(Step_T2D_ME_TBB& _stp) :
		stp(_stp), md(*(Model_T2D_ME_mt*)(stp.model)),
		spva0(md.sorted_pcl_var_arrays[stp.sorted_pcl_var_id]),
		spva1(md.sorted_pcl_var_arrays[stp.sorted_pcl_var_id ^ 1]),
		pcl_sort_mem(stp.pcl_sort_mem),
		pcl_m(md.pcl_m),
		pcl_bf(md.pcl_bf),
		pcl_t(md.pcl_t),
		pcl_pos(md.pcl_pos),
		pcl_vol(md.pcl_vol),
		pcl_mat_model(md.pcl_mat_model),
		pcl_index0(spva0.pcl_index),
		pcl_density0(spva0.pcl_density),
		pcl_v0(spva0.pcl_v),
		pcl_disp0(spva0.pcl_disp),
		pcl_N0(spva0.pcl_N),
		pcl_stress0(spva0.pcl_stress),
		pcl_strain0(spva0.pcl_strain),
		pcl_estrain0(spva0.pcl_estrain),
		pcl_pstrain0(spva0.pcl_pstrain),
		pcl_index1(spva1.pcl_index),
		pcl_density1(spva1.pcl_density),
		pcl_v1(spva1.pcl_v),
		pcl_disp1(spva1.pcl_disp),
		pcl_N1(spva1.pcl_N),
		pcl_stress1(spva1.pcl_stress),
		pcl_strain1(spva1.pcl_strain),
		pcl_estrain1(spva1.pcl_estrain),
		pcl_pstrain1(spva1.pcl_pstrain),
		elem_node_id(md.elem_node_id),
		elem_dN_ab(md.elem_dN_ab),
		elem_dN_c(md.elem_dN_c),
		elem_area(md.elem_area),
		elem_pcl_m(md.elem_pcl_m),
		elem_density(md.elem_density),
		elem_de(md.elem_de),
		elem_m_de_vol(md.elem_m_de_vol),
		elem_node_vm(md.elem_node_vm),
		elem_node_force(md.elem_node_force),
		node_a(md.node_a),
		node_v(md.node_v),
		node_has_vbc(md.node_has_vbc),
		node_am(md.node_am),
		node_de_vol(md.node_de_vol) {}
	
	MapPclToBgMeshTask::MapPclToBgMeshTask(
		size_t _block_id,
		size_t _block_num,
		CalData_T2D_ME_TBB& _cd) :
		block_id(_block_id),
		block_num(_block_num),
		cd(_cd) {}
	
	tbb::task* MapPclToBgMeshTask::execute()
	{
		Step_T2D_ME_TBB& stp = cd.stp;
		size_t* const pcl_in_elem = cd.pcl_sort_mem.out_key;
		size_t* const prev_pcl_id = cd.pcl_sort_mem.out_val;
		size_t e_id;
		size_t p_id0 = Block_Low(block_id, block_num, stp.valid_pcl_num);
		e_id = pcl_in_elem[p_id0];
		while (p_id0 != SIZE_MAX && e_id == pcl_in_elem[--p_id0]);
		++p_id0;
		assert(p_id0 <= stp.valid_pcl_num);
		size_t p_id1 = Block_Low(block_id + 1, block_num, stp.valid_pcl_num);
		e_id = pcl_in_elem[p_id1];
		while (p_id1 != SIZE_MAX && e_id == pcl_in_elem[--p_id1]);
		++p_id1;
		assert(p_id1 <= stp.valid_pcl_num);
		
		Model_T2D_ME_mt& md = cd.md;
		const double* const pcl_m = cd.pcl_m;
		const Force* pcl_bf = cd.pcl_bf;
		const Force* pcl_t = cd.pcl_t;
		double* const pcl_vol = cd.pcl_vol;
		size_t* const pcl_index0 = cd.pcl_index0;
		double* const pcl_density0 = cd.pcl_density0;
		ShapeFunc* const pcl_N0 = cd.pcl_N0;
		Velocity* const pcl_v0 = cd.pcl_v0;
		Displacement* const pcl_disp0 = cd.pcl_disp0;
		Stress* const pcl_stress0 = cd.pcl_stress0;
		const size_t* const pcl_index1 = cd.pcl_index1;
		const double* const pcl_density1 = cd.pcl_density1;
		const ShapeFunc* const pcl_N1 = cd.pcl_N1;
		const Velocity* const pcl_v1 = cd.pcl_v1;
		const Displacement* const pcl_disp1 = cd.pcl_disp1;
		const Stress* const pcl_stress1 = cd.pcl_stress1;
		const ShapeFuncAB* const elem_dN_ab = cd.elem_dN_ab;
		const double* const elem_area = cd.elem_area;
		double* const elem_pcl_m = cd.elem_pcl_m;
		double* const elem_density = cd.elem_density;
		ElemNodeVM *const elem_node_vm = cd.elem_node_vm;
		Force *const elem_node_force = cd.elem_node_force;
		const size_t my_th_id = tbb::task_arena::current_thread_index();
		size_t* const valid_elems = cd.valid_elem_arrays[my_th_id];
		size_t valid_elem_num = 0;
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
		assert(e_id < md.elem_num);
		for (size_t p_id = p_id0; p_id < p_id1; ++p_id)
		{
			// pcl index
			const size_t prev_p_id = prev_pcl_id[p_id];
			assert(prev_p_id < stp.prev_valid_pcl_num);
			const size_t ori_p_id = pcl_index1[prev_p_id];
			assert(ori_p_id < md.ori_pcl_num);
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
				assert(e_id < md.elem_num || e_id == SIZE_MAX);

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
	
}
