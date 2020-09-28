#include "Simulations_pcp.h"

#include "Step_T2D_ME_mt.h"

#define one_third (1.0/3.0)

int Step_T2D_ME_mt::init_calculation()
{
	return 0;
}

int Step_T2D_ME_mt::finalize_calculation()
{

	return 0;
}

int solve_substep_T2D_ME_mt(void* _self)
{
	typedef Model_T2D_ME_mt::PclIndex PclIndex;
	typedef Model_T2D_ME_mt::PclMass PclMass;
	typedef Model_T2D_ME_mt::PclDensity PclDensity;
	typedef Model_T2D_ME_mt::PclStress PclStress;
	typedef Model_T2D_ME_mt::PclShapeFunc PclShapeFunc;
	typedef Model_T2D_ME_mt::PclV PclV;
	typedef Model_T2D_ME_mt::PclBodyForce PclBodyForce;
	typedef Model_T2D_ME_mt::PclTraction PclTraction;
	typedef Model_T2D_ME_mt::ElemPclList ElemPclList;
	typedef Model_T2D_ME_mt::ElemNodeVarOffset ElemNodeVarOffset;
	typedef Model_T2D_ME_mt::ElemDensity ElemDensity;
	typedef Model_T2D_ME_mt::ElemShapeFuncAB ElemShapeFuncAB;
	typedef Model_T2D_ME_mt::ElemShapeFuncC ElemShapeFuncC;
	typedef Model_T2D_ME_mt::NodeElemDeVol NodeElemDeVol;
	typedef Model_T2D_ME_mt::NodeElemVM NodeElemVM;
	typedef Model_T2D_ME_mt::NodeElemAF NodeElemAF;
	typedef Model_T2D_ME_mt::NodeElemList NodeElemList;
	typedef Model_T2D_ME_mt::NodeMotion NodeMotion;

	Step_T2D_ME_mt& self = *(Step_T2D_ME_mt *)(_self);
	Model_T2D_ME_mt& md = *(Model_T2D_ME_mt*)(self.model);

	PclIndex* pcl_index = md.pcl_index;
	PclMass *pcl_m = md.pcl_m;
	PclDensity *pcl_density = md.pcl_density;
	PclStress* pcl_stress = md.pcl_stress;
	PclShapeFunc* pcl_N = md.pcl_N;
	PclV *pcl_v = md.pcl_v;
	PclBodyForce *pcl_bf = md.pcl_bf;
	PclTraction *pcl_t = md.pcl_t;

	ElemPclList* elem_pcl_list = md.elem_pcl_list;
	ElemNodeVarOffset *elem_node_var_offset = md.elem_node_var_offset;
	ElemDensity *elem_density = md.elem_density;
	ElemShapeFuncAB* elem_sf_ab = md.elem_sf_ab;
	ElemShapeFuncC* elem_sf_c = md.elem_sf_c;

	NodeElemAF *ne_af = md.ne_af;
	NodeElemVM* ne_vm = md.ne_vm;
	NodeElemDeVol* ne_de_vol = md.ne_de_vol;
	
	uint32_t *node_elem_offset = md.node_elem_offset;
	NodeElemList *node_elem_list = md.node_elem_list;
	NodeMotion *node_motion = md.node_motion;

	float p_vol, p_m, p_N_m;
	float e_pcl_m, e_pcl_bfx, e_pcl_bfy;
	float e_pcl_vol, e_s11, e_s22, e_s12;
	for (uint32_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		e_pcl_m = 0.0f;
		e_pcl_bfx = 0.0f;
		e_pcl_bfy = 0.0f;
		e_pcl_vol = 0.0f;
		e_s11 = 0.0f;
		e_s22 = 0.0f;
		e_s12 = 0.0f;
		ElemNodeVarOffset& e_nv_off = elem_node_var_offset[e_id];
		NodeElemVM& ne_vm1 = ne_vm[e_nv_off.n1];
		NodeElemVM& ne_vm2 = ne_vm[e_nv_off.n2];
		NodeElemVM& ne_vm3 = ne_vm[e_nv_off.n3];
		ne_vm1.vm = 0.0f;
		ne_vm1.vmx = 0.0f;
		ne_vm1.vmy = 0.0f;
		ne_vm2.vm = 0.0f;
		ne_vm2.vmx = 0.0f;
		ne_vm2.vmy = 0.0f;
		ne_vm3.vm = 0.0f;
		ne_vm3.vmx = 0.0f;
		ne_vm3.vmy = 0.0f;
		NodeElemAF& ne_af1 = ne_af[e_nv_off.n1];
		NodeElemAF& ne_af2 = ne_af[e_nv_off.n2];
		NodeElemAF& ne_af3 = ne_af[e_nv_off.n3];
		ne_af1.am = 0.0f;
		ne_af1.fx = 0.0f;
		ne_af1.fy = 0.0f;
		ne_af2.am = 0.0f;
		ne_af2.fx = 0.0f;
		ne_af2.fy = 0.0f;
		ne_af3.am = 0.0f;
		ne_af3.fx = 0.0f;
		ne_af3.fy = 0.0f;
		ElemPclList& e_pcl_list = elem_pcl_list[e_id];
		for (uint32_t pcl_id = e_pcl_list.start_id; pcl_id < e_pcl_list.end_id; ++pcl_id)
		{
			p_m = pcl_m[pcl_id].m;
			e_pcl_m += p_m;
			e_pcl_bfx += pcl_bf[pcl_id].bfx;
			e_pcl_bfy += pcl_bf[pcl_id].bfy;
			p_vol = p_m / pcl_density[pcl_id].density;
			e_pcl_vol += p_vol;
			PclStress& p_s = pcl_stress[pcl_id];
			e_s11 += p_s.s11 * p_vol;
			e_s22 += p_s.s22 * p_vol;
			e_s12 += p_s.s12 * p_vol;
			PclV& p_v = pcl_v[pcl_id];
			PclShapeFunc &p_N = pcl_N[pcl_id];
			p_N_m = p_N.N1 * p_m;
			ne_vm1.vm += p_N_m;
			ne_vm1.vmx += p_N_m * p_v.vx;
			ne_vm1.vmy += p_N_m * p_v.vy;
			p_N_m = p_N.N2 * p_m;
			ne_vm2.vm += p_N_m;
			ne_vm2.vmx += p_N_m * p_v.vx;
			ne_vm2.vmy += p_N_m * p_v.vy;
			p_N_m = p_N.N3 * p_m;
			ne_vm3.vm += p_N_m;
			ne_vm3.vmx += p_N_m * p_v.vx;
			ne_vm3.vmy += p_N_m * p_v.vy;
			PclTraction& p_t = pcl_t[pcl_id];
			ne_af1.fx += p_N.N1 * p_t.tx;
			ne_af1.fy += p_N.N1 * p_t.ty;
			ne_af2.fx += p_N.N2 * p_t.tx;
			ne_af2.fy += p_N.N2 * p_t.ty;
			ne_af3.fx += p_N.N3 * p_t.tx;
			ne_af3.fy += p_N.N3 * p_t.ty;
		}
		if (e_pcl_m != 0.0)
		{
			elem_density[e_id].density = e_pcl_m / e_pcl_vol;
			e_s11 /= e_pcl_vol;
			e_s22 /= e_pcl_vol;
			e_s12 /= e_pcl_vol;
			e_pcl_m *= one_third;
			e_pcl_bfx *= one_third;
			e_pcl_bfy *= one_third;
			ElemShapeFuncAB& e_sf = elem_sf_ab[e_id];
			ne_af1.am = e_pcl_m;
			ne_af1.fx += e_pcl_bfx;
			ne_af1.fx -= (e_sf.dN1_dx * e_s11 + e_sf.dN1_dy * e_s12) * e_pcl_vol;
			ne_af1.fy += e_pcl_bfy;
			ne_af1.fy -= (e_sf.dN1_dx * e_s12 + e_sf.dN1_dy * e_s22) * e_pcl_vol;
			ne_af2.am = e_pcl_m;
			ne_af2.fx += e_pcl_bfx;
			ne_af2.fx -= (e_sf.dN2_dx * e_s11 + e_sf.dN2_dy * e_s12) * e_pcl_vol;
			ne_af2.fy += e_pcl_bfy;
			ne_af2.fy -= (e_sf.dN2_dx * e_s12 + e_sf.dN2_dy * e_s22) * e_pcl_vol;
			ne_af3.am = e_pcl_m;
			ne_af3.fx += e_pcl_bfx;
			ne_af3.fx -= (e_sf.dN3_dx * e_s11 + e_sf.dN3_dy * e_s12) * e_pcl_vol;
			ne_af3.fy += e_pcl_bfy;
			ne_af3.fy -= (e_sf.dN3_dx * e_s12 + e_sf.dN3_dy * e_s22) * e_pcl_vol;
		}
	}

	uint32_t ne_off;
	float n_am, n_fx, n_fy;
	float n_vm, n_vmx, n_vmy;
	for (uint32_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		NodeElemList &ne_ids = node_elem_list[n_id];
		if (ne_ids.start_id == ne_ids.end_id)
			continue; // has boundary velocity
		n_am = 0.0f;
		n_fx = 0.0f;
		n_fy = 0.0f;
		n_vm = 0.0f;
		n_vmx = 0.0f;
		n_vmy = 0.0f;
		for (uint32_t ne_id = ne_ids.start_id; ne_id < ne_ids.end_id; ++ne_id)
		{
			ne_off = node_elem_offset[ne_id];
			NodeElemAF &naf = ne_af[ne_off];
			n_am += naf.am;
			n_fx += naf.fx;
			n_fy += naf.fy;
			NodeElemVM &nvm = ne_vm[ne_off];
			n_vm += nvm.vm;
			n_vmx += nvm.vmx;
			n_vmy += nvm.vmy;
		}
		if (n_am != 0.0f)
		{
			NodeMotion& n_m = node_motion[n_id];
			n_m.ax = n_fx / n_am;
			n_m.ay = n_fy / n_am;
			n_m.vx = n_vmx / n_vm + n_m.ax * self.dtime;
			n_m.vy = n_vmy / n_vm + n_m.ay * self.dtime;
		}
	}

	float de_vol;
	for (uint32_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		//ElemStrainInc;
	}
	return 0;
}
