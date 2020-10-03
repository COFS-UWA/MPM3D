#include "Simulations_pcp.h"

#include "Step_T2D_ME_mt.h"

#define one_third (1.0f/3.0f)
#define N_min (1.0e-10f)

int Step_T2D_ME_mt::init_calculation()
{
	Model_T2D_ME_mt& md = *(Model_T2D_ME_mt*)model;

	md.cur_pcl_sorted_var_id = 0;

	uint32_t p_id, pcl_in_elem_id;
	uint32_t count_sort_id = 0;
	PclSortedVarArray &psva = md.pcl_sorted_var_array[0];
	for (p_id = 0; p_id < md.pcl_num; ++p_id)
	{
		PclDisp &p_d = psva.pcl_disp[p_id];
		p_d.ux = 0.0f;
		p_d.uy = 0.0f;
	
		PclPos& p_p = md.pcl_pos[p_id];
		PclShapeFunc& p_N = psva.pcl_N[p_id];
		pcl_in_elem_id = md.find_pcl_in_which_elem(p_p.x, p_p.y, p_N);
		if (pcl_in_elem_id == UINT32_MAX)
		{
			--md.pcl_num;
			continue;
		}
		md.pcl_unsorted_id_array[count_sort_id] = p_id;
		md.pcl_in_elem_id_array[count_sort_id] = pcl_in_elem_id;
		++count_sort_id;
		if (p_N.N1 < N_min)
			p_N.N1 = N_min;
		if (p_N.N2 < N_min)
			p_N.N2 = N_min;
		if (p_N.N3 < N_min)
			p_N.N3 = N_min;
	}

	memset(md.elem_has_pcl_num_array, 0, sizeof(uint32_t) * md.elem_num);
	for (p_id = 0; p_id < md.pcl_num; ++p_id)
		++md.elem_has_pcl_num_array[md.pcl_in_elem_id_array[p_id]];
	ElemPclList &ep_list = psva.elem_pcl_list[0];
	ep_list.end_id = md.elem_has_pcl_num_array[0];
	for (uint32_t e_id = 1; e_id < md.elem_num; ++e_id)
	{
		md.elem_has_pcl_num_array[e_id] += md.elem_has_pcl_num_array[e_id - 1];
		ep_list.end_id = md.elem_has_pcl_num_array[e_id];
	}
	for (p_id = md.pcl_num - 1; p_id < md.pcl_num; --p_id)
	{
		pcl_in_elem_id = md.pcl_in_elem_id_array[p_id];
		--md.elem_has_pcl_num_array[pcl_in_elem_id];
		md.pcl_new_to_cur_map[md.elem_has_pcl_num_array[pcl_in_elem_id]]
			= md.pcl_unsorted_id_array[p_id];
	}

	// reorder particle variables
	PclSortedVarArray& psva0 = md.pcl_sorted_var_array[0];
	PclSortedVarArray& psva1 = md.pcl_sorted_var_array[1];
	for (p_id = 0; p_id < md.pcl_num; ++p_id)
	{
		psva1.pcl_index[p_id].id = psva0.pcl_index[p_id].id;
		psva1.pcl_density[p_id].density = psva0.pcl_density[p_id].density;
		PclV& p_v0 = psva0.pcl_v[p_id];
		PclV& p_v1 = psva1.pcl_v[p_id];
		p_v1.vx = p_v0.vx;
		p_v1.vy = p_v0.vy;
		PclStress& p_s0 = psva0.pcl_stress[p_id];
		PclStress& p_s1 = psva1.pcl_stress[p_id];
		p_s1.s11 = p_s0.s11;
		p_s1.s22 = p_s0.s22;
		p_s1.s12 = p_s0.s12;
	}

	uint32_t cur_pcl_id;
	for (p_id = 0; p_id < md.pcl_num; ++p_id)
	{
		cur_pcl_id = md.pcl_new_to_cur_map[p_id];

		psva0.pcl_index[p_id].id = psva1.pcl_index[cur_pcl_id].id;

		psva0.pcl_density[p_id].density = psva1.pcl_density[cur_pcl_id].density;

		PclDisp& p_d0 = psva0.pcl_disp[p_id];
		PclDisp& p_d1 = psva1.pcl_disp[cur_pcl_id];
		p_d0.ux = p_d1.ux;
		p_d0.uy = p_d1.uy;

		PclV& p_v0 = psva0.pcl_v[p_id];
		PclV& p_v1 = psva1.pcl_v[cur_pcl_id];
		p_v0.vx = p_v1.vx;
		p_v0.vy = p_v1.vy;

		PclShapeFunc& p_sf0 = psva0.pcl_N[p_id];
		PclShapeFunc& p_sf1 = psva1.pcl_N[cur_pcl_id];
		p_sf0.N1 = p_sf1.N1;
		p_sf0.N2 = p_sf1.N2;
		p_sf0.N3 = p_sf1.N3;

		PclStress& p_s0 = psva0.pcl_stress[p_id];
		PclStress& p_s1 = psva1.pcl_stress[cur_pcl_id];
		p_s0.s11 = p_s1.s11;
		p_s0.s22 = p_s1.s22;
		p_s0.s12 = p_s1.s12;
	}

	return 0;
}

int Step_T2D_ME_mt::finalize_calculation() { return 0; }

int solve_substep_T2D_ME_mt(void* _self)
{
	typedef Model_T2D_ME_mt::PclMass PclMass;
	typedef Model_T2D_ME_mt::PclBodyForce PclBodyForce;
	typedef Model_T2D_ME_mt::PclTraction PclTraction;
	typedef Model_T2D_ME_mt::PclPos PclPos;
	typedef Model_T2D_ME_mt::PclSortedVarArray PclSortedVarArray;
	typedef Model_T2D_ME_mt::PclIndex PclIndex;
	typedef Model_T2D_ME_mt::PclDensity PclDensity;
	typedef Model_T2D_ME_mt::PclDisp PclDisp;
	typedef Model_T2D_ME_mt::PclV PclV;
	typedef Model_T2D_ME_mt::PclShapeFunc PclShapeFunc;
	typedef Model_T2D_ME_mt::PclStress PclStress;
	typedef Model_T2D_ME_mt::ElemPclList ElemPclList;
	typedef Model_T2D_ME_mt::ElemNodeIndex ElemNodeIndex;
	typedef Model_T2D_ME_mt::ElemArea ElemArea;
	typedef Model_T2D_ME_mt::ElemShapeFuncAB ElemShapeFuncAB;
	typedef Model_T2D_ME_mt::ElemShapeFuncC ElemShapeFuncC;
	typedef Model_T2D_ME_mt::ElemDensity ElemDensity;
	typedef Model_T2D_ME_mt::ElemStrainInc ElemStrainInc;
	typedef Model_T2D_ME_mt::ElemStress ElemStress;
	typedef Model_T2D_ME_mt::ElemAm ElemAm;
	typedef Model_T2D_ME_mt::ElemAmDeVol ElemAmDeVol;
	typedef Model_T2D_ME_mt::ElemNodeVM ElemNodeVM;
	typedef Model_T2D_ME_mt::ElemNodeForce ElemNodeForce;
	typedef Model_T2D_ME_mt::NodeElemList NodeElemList;
	typedef Model_T2D_ME_mt::NodeA NodeA;
	typedef Model_T2D_ME_mt::NodeV NodeV;
	typedef Model_T2D_ME_mt::NodeAm NodeAm;
	typedef Model_T2D_ME_mt::NodeDeVol NodeDeVol;
	
	Step_T2D_ME_mt& self = *(Step_T2D_ME_mt *)(_self);
	Model_T2D_ME_mt& md = *(Model_T2D_ME_mt*)(self.model);

	uint32_t pcl_num = md.pcl_num;
	uint32_t p_id;
	uint32_t elem_num = md.elem_num;
	uint32_t e_id;
	uint32_t node_num = md.node_num;
	uint32_t n_id;

	uint32_t vx_bc_n_num = md.vx_bc_n_num;
	uint32_t* vx_bc_n_ids = md.vx_bc_n_ids;
	uint32_t vy_bc_n_num = md.vy_bc_n_num;
	uint32_t* vy_bc_n_ids = md.vy_bc_n_ids;
	
	PclMass *pcl_m = md.pcl_m;
	PclBodyForce* pcl_bf = md.pcl_bf;
	PclTraction* pcl_t = md.pcl_t;
	PclPos* pcl_pos = md.pcl_pos;
	MatModel::MaterialModel **pcl_mat_model = md.pcl_mat_model;

	PclSortedVarArray &pscv0 = md.pcl_sorted_var_array[md.cur_pcl_sorted_var_id];
	PclIndex* pcl_index0 = pscv0.pcl_index;
	PclDensity *pcl_density0 = pscv0.pcl_density;
	PclDisp* pcl_disp0 = pscv0.pcl_disp;
	PclV* pcl_v0 = pscv0.pcl_v;
	PclShapeFunc* pcl_N0 = pscv0.pcl_N;
	PclStress* pcl_stress0 = pscv0.pcl_stress;
	ElemPclList* elem_pcl_list0 = pscv0.elem_pcl_list;
	md.cur_pcl_sorted_var_id ^= 1;
	PclSortedVarArray& pscv1 = md.pcl_sorted_var_array[md.cur_pcl_sorted_var_id];
	PclIndex* pcl_index1 = pscv1.pcl_index;
	PclDensity* pcl_density1 = pscv1.pcl_density;
	PclDisp* pcl_disp1 = pscv1.pcl_disp;
	PclV* pcl_v1 = pscv1.pcl_v;
	PclShapeFunc* pcl_N1 = pscv1.pcl_N;
	PclStress* pcl_stress1 = pscv1.pcl_stress;
	ElemPclList* elem_pcl_list1 = pscv1.elem_pcl_list;

	ElemNodeIndex *elem_node_id = md.elem_node_id;
	ElemArea* elem_area = md.elem_area;
	ElemShapeFuncAB* elem_sf_ab = md.elem_sf_ab;
	ElemShapeFuncC* elem_sf_c = md.elem_sf_c;

	ElemDensity *elem_density = md.elem_density;
	ElemStrainInc *elem_de = md.elem_de;
	ElemStress* elem_stress = md.elem_stress;
	ElemAm* elem_am = md.elem_am;
	ElemAmDeVol* elem_am_de_vol = md.elem_am_de_vol;

	ElemNodeVM *elem_node_vm = md.elem_node_vm;
	ElemNodeForce* elem_node_force = md.elem_node_force;
	
	uint32_t *elem_id_array = md.elem_id_array;
	uint32_t *node_elem_id_array = md.node_elem_id_array;
	NodeElemList *node_elem_list = md.node_elem_list;
	NodeA *node_a = md.node_a;
	NodeV* node_v = md.node_v;
	NodeAm* node_am = md.node_am;
	NodeDeVol* node_de_vol = md.node_de_vol;

	uint32_t* pcl_unsorted_id_array = md.pcl_unsorted_id_array;
	uint32_t* pcl_in_elem_id_array = md.pcl_in_elem_id_array;
	uint32_t* elem_has_pcl_num_array = md.elem_has_pcl_num_array;
	uint32_t* pcl_new_to_cur_map = md.pcl_new_to_cur_map;

	float dtime = float(self.dtime);
	
	union // local variables
	{
		// velocity mapping
		// nodal force integration
		struct
		{
			uint32_t ne_var_start_id;
			uint32_t ori_pcl_id;
			float p_vol, p_m, p_N_m;
			float e_pcl_m, e_pcl_bfx, e_pcl_bfy;
			float e_pcl_vol, e_s11, e_s22, e_s12;
		};
		// update nodal variables
		struct
		{
			uint32_t ne_id;
			uint32_t node_var_id;
			float n_am, n_fx, n_fy;
			float n_vm, n_vmx, n_vmy;
			uint32_t vbc_id;
		};
		// element strain
		// strain enhancement method
		struct
		{
			float e_de_vol, pcl_x, pcl_y;
			float n_am, n_am_de_vol;
		};
		// update particle variables
		struct
		{
			uint32_t count_sort_id;
			uint32_t pcl_in_elem_id;
			int32_t mm_int_res;
		};
		// sort particle variables
		struct
		{
			uint32_t cur_pcl_id;
		};
	};
	
	// velocity mapping
	// nodal force integration
	ne_var_start_id = 0;
	p_id = 0;
	for (e_id = 0; e_id < elem_num; ++e_id)
	{
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
		ElemAm& e_am = elem_am[e_id];
		e_am.am = 0.0f;
		ElemPclList& e_pcl_list = elem_pcl_list0[e_id];
		if (p_id == e_pcl_list.end_id)
			continue; // no particle
		e_pcl_m = 0.0f;
		e_pcl_bfx = 0.0f;
		e_pcl_bfy = 0.0f;
		e_pcl_vol = 0.0f;
		e_s11 = 0.0f;
		e_s22 = 0.0f;
		e_s12 = 0.0f;
		for (; p_id < e_pcl_list.end_id; ++p_id)
		{
			ori_pcl_id = pcl_index0[p_id].id;
			p_m = pcl_m[ori_pcl_id].m;
			e_pcl_m += p_m;
			PclBodyForce& p_bf = pcl_bf[ori_pcl_id];
			e_pcl_bfx += p_bf.bfx;
			e_pcl_bfy += p_bf.bfy;
			p_vol = p_m / pcl_density0[p_id].density;
			e_pcl_vol += p_vol;
			PclStress& p_s = pcl_stress0[p_id];
			e_s11 += p_s.s11 * p_vol;
			e_s22 += p_s.s22 * p_vol;
			e_s12 += p_s.s12 * p_vol;
			PclShapeFunc& p_N = pcl_N0[p_id];
			PclV& p_v = pcl_v0[p_id];
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
		elem_density[e_id].density = e_pcl_m / e_pcl_vol;
		e_s11 /= e_pcl_vol;
		e_s22 /= e_pcl_vol;
		e_s12 /= e_pcl_vol;
		if (e_pcl_vol > elem_area[e_id].area)
			e_pcl_vol = elem_area[e_id].area;
		e_pcl_m *= one_third;
		e_pcl_bfx *= one_third;
		e_pcl_bfy *= one_third;
		e_am.am = e_pcl_m;
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
		ne_var_start_id += 3;
	}

	// update nodal variables
	ne_id = 0;
	for (n_id = 0; n_id < node_num; ++n_id)
	{
		n_am = 0.0f;
		n_fx = 0.0f;
		n_fy = 0.0f;
		n_vm = 0.0f;
		n_vmx = 0.0f;
		n_vmy = 0.0f;
		NodeElemList& ne_ids = node_elem_list[n_id];
		for (; ne_id < ne_ids.end_id; ++ne_id)
		{
			node_var_id = elem_id_array[ne_id];
			ElemAm &eam = elem_am[node_var_id];
			n_am += eam.am;
			node_var_id = node_elem_id_array[ne_id];
			ElemNodeForce &nf = elem_node_force[node_var_id];
			n_fx += nf.fx;
			n_fy += nf.fy;
			ElemNodeVM &nvm = elem_node_vm[node_var_id];
			n_vm += nvm.vm;
			n_vmx += nvm.vmx;
			n_vmy += nvm.vmy;
		}
		if (n_am != 0.0f)
		{
			node_am[n_id].am = n_am;
			NodeA& n_a = node_a[n_id];
			n_a.ax = n_fx / n_am;
			n_a.ay = n_fy / n_am;
			NodeV& n_v = node_v[n_id];
			n_v.vx = n_vmx / n_vm + n_a.ax * dtime;
			n_v.vy = n_vmy / n_vm + n_a.ay * dtime;
		}
	}

	// apply velocity bc
	for (vbc_id = 0; vbc_id < vx_bc_n_num; ++vbc_id)
	{
		n_id = vx_bc_n_ids[vbc_id];
		node_a[n_id].ax = 0.0f;
		node_v[n_id].vx = 0.0f;
	}
	for (vbc_id = 0; vbc_id < vy_bc_n_num; ++vbc_id)
	{
		n_id = vy_bc_n_ids[vbc_id];
		node_a[n_id].ax = 0.0f;
		node_v[n_id].vx = 0.0f;
	}

	// cal element strain and strain enhancement
	for (e_id = 0; e_id < elem_num; ++e_id)
	{
		ElemNodeIndex& e_n_id = elem_node_id[e_id];
		NodeV& n_v1 = node_v[e_n_id.n1];
		NodeV& n_v2 = node_v[e_n_id.n2];
		NodeV& n_v3 = node_v[e_n_id.n3];
		ElemShapeFuncAB &e_sf = elem_sf_ab[e_id];
		ElemStrainInc &e_de = elem_de[e_id];
		e_de.de11 = (e_sf.dN1_dx * n_v1.vx + e_sf.dN2_dx * n_v2.vx + e_sf.dN3_dx * n_v3.vx) * dtime;
		e_de.de22 = (e_sf.dN1_dy * n_v1.vy + e_sf.dN2_dy * n_v2.vy + e_sf.dN3_dy * n_v3.vy) * dtime;
		e_de.de12 = (e_sf.dN1_dx * n_v1.vy + e_sf.dN2_dx * n_v2.vy + e_sf.dN3_dx * n_v3.vy
				   + e_sf.dN1_dy * n_v1.vx + e_sf.dN2_dy * n_v2.vx + e_sf.dN3_dy * n_v3.vx) * dtime * 0.5f;
		e_de_vol = e_de.de11 + e_de.de22;
		ElemAmDeVol& e_am_de_vol = elem_am_de_vol[e_id];
		e_am_de_vol.am_de_vol = elem_am[e_id].am * e_de_vol;
		e_de_vol *= one_third;
		e_de.de11 -= e_de_vol;
		e_de.de22 -= e_de_vol;
	}

	ne_id = 0;
	for (n_id = 0; n_id < node_num; ++n_id)
	{
		NodeElemList& ne_ids = node_elem_list[n_id];
		n_am_de_vol = 0.0f;
		for (; ne_id < ne_ids.end_id; ++ne_id)
		{
			node_var_id = elem_id_array[ne_id];
			n_am_de_vol += elem_am_de_vol[node_var_id].am_de_vol;
		}
		n_am = node_am[n_id].am;
		if (n_am != 0.0f)
			node_de_vol[n_id].de_vol = n_am_de_vol / n_am;
	}

	// update particle variables
	double dstrain[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	const double *dstress;
	count_sort_id = 0;
	p_id = 0;
	for (e_id = 0; e_id < elem_num; ++e_id)
	{
		ElemNodeIndex &e_n_id = elem_node_id[e_id];
		NodeA& n_a1 = node_a[e_n_id.n1];
		NodeA& n_a2 = node_a[e_n_id.n2];
		NodeA& n_a3 = node_a[e_n_id.n3];
		NodeV& n_v1 = node_v[e_n_id.n1];
		NodeV& n_v2 = node_v[e_n_id.n2];
		NodeV& n_v3 = node_v[e_n_id.n3];
		e_de_vol = (node_de_vol[e_n_id.n1].de_vol
				  + node_de_vol[e_n_id.n2].de_vol
				  + node_de_vol[e_n_id.n3].de_vol) * one_third;
		ElemDensity& e_d = elem_density[e_id];
		e_d.density /= (1.0f + e_de_vol);
		e_de_vol *= one_third;
		ElemStrainInc& e_de = elem_de[e_id];
		e_de.de11 += e_de_vol;
		e_de.de12 += e_de_vol;
		dstrain[0] = double(e_de.de11);
		dstrain[1] = double(e_de.de22);
		dstrain[3] = double(e_de.de12);
		ElemPclList& e_pcl_list = elem_pcl_list0[e_id];
		for (; p_id < e_pcl_list.end_id; ++p_id)
		{
			PclShapeFunc& p_N = pcl_N0[p_id];
			
			PclV& p_v = pcl_v0[p_id];
			p_v.vx += (p_N.N1 * n_a1.ax + p_N.N2 * n_a2.ax + p_N.N3 * n_a3.ax) * dtime;
			p_v.vy += (p_N.N1 * n_a1.ay + p_N.N2 * n_a2.ay + p_N.N3 * n_a3.ay) * dtime;

			PclDisp& p_d = pcl_disp0[p_id];
			p_d.ux += (p_N.N1 * n_v1.vx + p_N.N2 * n_v2.vx + p_N.N3 * n_v3.vx) * dtime;
			p_d.uy += (p_N.N1 * n_v1.vy + p_N.N2 * n_v2.vy + p_N.N3 * n_v3.vy) * dtime;
			
			ori_pcl_id = pcl_index0[p_id].id;
			PclPos& p_p = pcl_pos[ori_pcl_id];
			pcl_x = p_p.x + p_d.ux;
			pcl_y = p_p.y + p_d.uy;

			pcl_in_elem_id = e_id;
			if (!md.is_in_element(pcl_x, pcl_y, e_id, p_N))
			{
				pcl_in_elem_id = md.find_pcl_in_which_elem(pcl_x, pcl_y, p_N);
				if (pcl_in_elem_id == UINT32_MAX)
				{
					--pcl_num;
					continue;
				}
				// add to sort list
				md.pcl_unsorted_id_array[count_sort_id] = p_id;
				md.pcl_in_elem_id_array[count_sort_id] = pcl_in_elem_id;
				++count_sort_id;
			}
			if (p_N.N1 < N_min)
				p_N.N1 = N_min;
			if (p_N.N2 < N_min)
				p_N.N2 = N_min;
			if (p_N.N3 < N_min)
				p_N.N3 = N_min;

			pcl_density0[p_id].density = e_d.density;

			MatModel::MaterialModel& pcl_mm = *pcl_mat_model[ori_pcl_id];
			mm_int_res = pcl_mm.integrate(dstrain);
			dstress = pcl_mm.get_dstress();
			PclStress& p_s0 = pcl_stress0[p_id];
			p_s0.s11 += float(dstress[0]);
			p_s0.s22 += float(dstress[1]);
			p_s0.s12 += float(dstress[3]);
		}
	}

	// sort particle variables
	memset(elem_has_pcl_num_array, 0, sizeof(uint32_t) * elem_num);
	for (p_id = 0; p_id < pcl_num; ++p_id)
		++elem_has_pcl_num_array[pcl_in_elem_id_array[p_id]];
	elem_pcl_list1[0].end_id = elem_has_pcl_num_array[0];
	for (e_id = 1; e_id < elem_num; ++e_id)
	{
		elem_has_pcl_num_array[e_id] += elem_has_pcl_num_array[e_id-1];
		elem_pcl_list1[e_id].end_id = elem_has_pcl_num_array[e_id-1];
	}
	for (p_id = pcl_num-1; p_id < pcl_num; --p_id)
	{
		pcl_in_elem_id = pcl_in_elem_id_array[p_id];
		--elem_has_pcl_num_array[pcl_in_elem_id];
		pcl_new_to_cur_map[elem_has_pcl_num_array[pcl_in_elem_id]] = pcl_unsorted_id_array[p_id];
	}

	// reorder particle variables
	for (p_id = 0; p_id < pcl_num; ++p_id)
	{
		cur_pcl_id = pcl_new_to_cur_map[p_id];

		pcl_index1[p_id].id = pcl_index0[cur_pcl_id].id;

		pcl_density1[p_id].density = pcl_density0[cur_pcl_id].density;
		
		PclDisp& p_d0 = pcl_disp0[cur_pcl_id];
		PclDisp& p_d1 = pcl_disp1[p_id];
		p_d1.ux = p_d0.ux;
		p_d1.uy = p_d0.uy;

		PclV &p_v0 = pcl_v0[cur_pcl_id];
		PclV &p_v1 = pcl_v1[p_id];
		p_v1.vx = p_v0.vx;
		p_v1.vy = p_v0.vy;

		PclShapeFunc &p_sf0 = pcl_N0[cur_pcl_id];
		PclShapeFunc& p_sf1 = pcl_N1[p_id];
		p_sf1.N1 = p_sf0.N1;
		p_sf1.N2 = p_sf0.N2;
		p_sf1.N3 = p_sf0.N3;

		PclStress& p_s0 = pcl_stress0[cur_pcl_id];
		PclStress& p_s1 = pcl_stress1[p_id];
		p_s1.s11 = p_s0.s11;
		p_s1.s22 = p_s0.s22;
		p_s1.s12 = p_s0.s12;
	}

	md.pcl_num = pcl_num;

	return 0;
}
