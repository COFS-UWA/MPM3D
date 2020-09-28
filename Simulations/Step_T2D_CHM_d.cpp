#include "Simulations_pcp.h"

#include "Step_T2D_CHM_d.h"

#define one_third (1.0/3.0)

Step_T2D_CHM_d::Step_T2D_CHM_d(const char* _name) :
	Step(_name, "Step_T2D_CHM_d", &solve_substep_T2D_CHM_d),
	model(nullptr), damping_ratio(0.0) {}

Step_T2D_CHM_d::~Step_T2D_CHM_d() {}

int Step_T2D_CHM_d::init_calculation()
{
	Model_T2D_CHM_d &md = *model;

	if (is_first_step) {}

	for (size_t pcl_id = 0; pcl_id < md.spcl_num; ++pcl_id)
	{
		SolidParticle &pcl = md.spcls[pcl_id];
		pcl.pe = md.elems;
		pcl.vol_s = pcl.m / pcl.density;
		pcl.x_ori = pcl.x;
		pcl.y_ori = pcl.y;
		pcl.ux = 0.0;
		pcl.uy = 0.0;
	}

	for (size_t pcl_id = 0; pcl_id < md.fpcl_num; ++pcl_id)
	{
		FluidParticle& pcl = md.fpcls[pcl_id];
		pcl.pe = md.elems;
		pcl.x_ori = pcl.x;
		pcl.y_ori = pcl.y;
		pcl.ux = 0.0;
		pcl.uy = 0.0;
	}
	
	return 0;
}

int Step_T2D_CHM_d::finalize_calculation() { return 0; }

int solve_substep_T2D_CHM_d(void *_self)
{
	typedef Model_T2D_CHM_d::Node Node;
	typedef Model_T2D_CHM_d::Element Element;
	typedef Model_T2D_CHM_d::SolidParticle SolidParticle;
	typedef Model_T2D_CHM_d::FluidParticle FluidParticle;
	Step_T2D_CHM_d &self = *(Step_T2D_CHM_d *)(_self);
	Model_T2D_CHM_d &md = *self.model;

	// init nodes
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];

		// solid phase
		n.am_s = 0.0;
		n.vm_s = 0.0;
		n.vx_s = 0.0;
		n.vy_s = 0.0;
		n.fx_ext_s = 0.0;
		n.fy_ext_s = 0.0;
		n.fx_int_s = 0.0;
		n.fy_int_s = 0.0;
		n.de_vol_s = 0.0;

		// fluid phase
		n.am_fs = 0.0;
		n.am_fp = 0.0;
		n.vm_f = 0.0;
		n.vx_f = 0.0;
		n.vy_f = 0.0;
		n.fx_ext_f = 0.0;
		n.fy_ext_f = 0.0;
		n.fx_int_f = 0.0;
		n.fy_int_f = 0.0;
		n.de_vol_fs = 0.0;
		n.de_vol_fp = 0.0;

		// solid - fluid interaction
		n.fx_drag = 0.0;
		n.fy_drag = 0.0;
	}

	// init elements
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];

		// solid phase
		e.pcl_m_s = 0.0;
		e.pcl_n = 0.0;
		e.pcl_vol_s = 0.0;
		e.s11 = 0.0;
		e.s22 = 0.0;
		e.s12 = 0.0;

		// fluid phase
		// seepage fluid
		e.pcl_m_fs = 0.0;
		e.pcl_vol_fs = 0.0;
		e.ps = 0.0;
		// pure fluid
		e.pcl_m_fp = 0.0;
		e.pcl_vol_fp = 0.0;
		e.pp = 0.0;
	}

	// init particles
	double N_m_s;
	for (size_t pcl_id = 0; pcl_id < md.spcl_num; ++pcl_id)
	{
		SolidParticle &pcl = md.spcls[pcl_id];
		if (pcl.pe)
		{
			if (!(pcl.pe = md.find_in_which_element(pcl)))
				continue;

			pcl.vol = pcl.vol_s / (1.0 - pcl.n);

			Element &e = *pcl.pe;
			e.pcl_m_s += pcl.m;
			e.pcl_n += pcl.vol_s;
			e.pcl_vol_s += pcl.vol;
			e.s11 += pcl.s11 * pcl.vol;
			e.s22 += pcl.s22 * pcl.vol;
			e.s12 += pcl.s12 * pcl.vol;

			Node &n1 = md.nodes[e.n1];
			N_m_s = pcl.N1 * pcl.m;
			n1.vm_s += N_m_s;
			n1.vx_s += N_m_s * pcl.vx;
			n1.vy_s += N_m_s * pcl.vy;

			Node &n2 = md.nodes[e.n2];
			N_m_s = pcl.N2 * pcl.m;
			n2.vm_s += N_m_s;
			n2.vx_s += N_m_s * pcl.vx;
			n2.vy_s += N_m_s * pcl.vy;

			Node &n3 = md.nodes[e.n3];
			N_m_s = pcl.N3 * pcl.m;
			n3.vm_s += N_m_s;
			n3.vx_s += N_m_s * pcl.vx;
			n3.vy_s += N_m_s * pcl.vy;
		}
	}

	md.reset_spcl_grids();

	double N_m_f;
	for (size_t pcl_id = 0; pcl_id < md.fpcl_num; ++pcl_id)
	{
		FluidParticle &pcl = md.fpcls[pcl_id];
		if (pcl.pe)
		{
			if (!(pcl.pe = md.find_in_which_element(pcl)))
				continue;

			pcl.vol = pcl.m / pcl.density;
			pcl.is_seepage = (pcl.pe->pcl_m_s != 0.0
							&& md.is_in_solid(pcl.x, pcl.y)) ? true : false;

			Element &e = *pcl.pe;
			if (pcl.is_seepage)
			{
				e.pcl_m_fs += pcl.m;
				e.pcl_vol_fs += pcl.vol;
				e.ps += pcl.p * pcl.vol;
			}
			else
			{
				e.pcl_m_fp += pcl.m;
				e.pcl_vol_fp += pcl.vol;
				e.pp += pcl.p * pcl.vol;
			}

			Node &n1 = md.nodes[e.n1];
			Node &n2 = md.nodes[e.n2];
			Node &n3 = md.nodes[e.n3];

			N_m_f = pcl.N1 * pcl.m;
			n1.vm_f += N_m_f;
			n1.vx_f += N_m_f * pcl.vx;
			n1.vy_f += N_m_f * pcl.vy;

			N_m_f = pcl.N2 * pcl.m;
			n2.vm_f += N_m_f;
			n2.vx_f += N_m_f * pcl.vx;
			n2.vy_f += N_m_f * pcl.vy;

			N_m_f = pcl.N3 * pcl.m;
			n3.vm_f += N_m_f;
			n3.vx_f += N_m_f * pcl.vx;
			n3.vy_f += N_m_f * pcl.vy;
		}
	}
	
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.vm_s != 0.0)
		{
			n.vx_s /= n.vm_s;
			n.vy_s /= n.vm_s;
		}
		if (n.vm_f != 0.0)
		{
			n.vx_f /= n.vm_f;
			n.vy_f /= n.vm_f;
		}
	}
	
	for (size_t v_id = 0; v_id < md.vsx_num; ++v_id)
	{
		Node& n = md.nodes[md.vsxs[v_id].node_id];
		n.vx_s = md.vsxs[v_id].v;
	}
	for (size_t v_id = 0; v_id < md.vsy_num; ++v_id)
	{
		Node& n = md.nodes[md.vsys[v_id].node_id];
		n.vy_s = md.vsys[v_id].v;
	}
	for (size_t v_id = 0; v_id < md.vfx_num; ++v_id)
	{
		Node& n = md.nodes[md.vfxs[v_id].node_id];
		n.vx_f = md.vfxs[v_id].v;
	}
	for (size_t v_id = 0; v_id < md.vfy_num; ++v_id)
	{
		Node& n = md.nodes[md.vfys[v_id].node_id];
		n.vy_f = md.vfys[v_id].v;
	}

	double e_dN1_dx, e_dN1_dy;
	double e_dN2_dx, e_dN2_dy;
	double e_dN3_dx, e_dN3_dy;
	double n_am, e_pcl_f_vol_max;
	double vis_s11, vis_s22, vis_s12;
	double n2_miu_div_k_vol, pcl_v_f, pcl_v_s;
	double n2_miu_div_k_vrx_vol, n2_miu_div_k_vry_vol;
	double e_one_n, p_vol_f;
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];

		PointInTriangle &epit = e.pt_in_tri;
		e_dN1_dx = epit.dN1_dx();
		e_dN1_dy = epit.dN1_dy();
		e_dN2_dx = epit.dN2_dx();
		e_dN2_dy = epit.dN2_dy();
		e_dN3_dx = epit.dN3_dx();
		e_dN3_dy = epit.dN3_dy();
		
		Node& n1 = md.nodes[e.n1];
		Node& n2 = md.nodes[e.n2];
		Node& n3 = md.nodes[e.n3];

		if (e.pcl_m_s != 0.0)
		{
			n_am = one_third * e.pcl_m_s;
			n1.am_s += n_am;
			n2.am_s += n_am;
			n3.am_s += n_am;

			e.pcl_n = 1.0 - e.pcl_n / e.pcl_vol_s; // 1.0 - Vs / V
			e.s11 /= e.pcl_vol_s;
			e.s22 /= e.pcl_vol_s;
			e.s12 /= e.pcl_vol_s;
			if (e.pcl_vol_s > e.area)
				e.pcl_vol_s = e.area;

			n1.fx_int_s += (e_dN1_dx * e.s11 + e_dN1_dy * e.s12) * e.pcl_vol_s;
			n1.fy_int_s += (e_dN1_dx * e.s12 + e_dN1_dy * e.s22) * e.pcl_vol_s;
			n2.fx_int_s += (e_dN2_dx * e.s11 + e_dN2_dy * e.s12) * e.pcl_vol_s;
			n2.fy_int_s += (e_dN2_dx * e.s12 + e_dN2_dy * e.s22) * e.pcl_vol_s;
			n3.fx_int_s += (e_dN3_dx * e.s11 + e_dN3_dy * e.s12) * e.pcl_vol_s;
			n3.fy_int_s += (e_dN3_dx * e.s12 + e_dN3_dy * e.s22) * e.pcl_vol_s;
		}

		if (e.pcl_m_fp != 0.0)
		{
			n_am = one_third * e.pcl_m_fp;
			n1.am_fp += n_am;
			n2.am_fp += n_am;
			n3.am_fp += n_am;

			e.pp /= e.pcl_vol_fp;
			e.pcl_density_fp = e.pcl_m_fp / e.pcl_vol_fp;
			e_pcl_f_vol_max = e.area - e.pcl_vol_s;
			if (e.pcl_vol_fp > e_pcl_f_vol_max)
				e.pcl_vol_fp = e_pcl_f_vol_max;

			vis_s11 = md.miu * 2.0 * (e_dN1_dx * n1.vx_f + e_dN2_dx * n2.vx_f + e_dN3_dx * n3.vx_f) - e.pp;
			vis_s22 = md.miu * 2.0 * (e_dN1_dy * n1.vy_f + e_dN2_dy * n2.vy_f + e_dN3_dy * n3.vy_f) - e.pp;
			vis_s12 = md.miu * (e_dN1_dx * n1.vy_f + e_dN2_dx * n2.vy_f + e_dN3_dx * n3.vy_f
							  + e_dN1_dy * n1.vx_f + e_dN2_dy * n2.vx_f + e_dN3_dy * n3.vx_f);

			n1.fx_int_f += (e_dN1_dx * vis_s11 + e_dN1_dy * vis_s12) * e.pcl_vol_fp;
			n1.fy_int_f += (e_dN1_dx * vis_s12 + e_dN1_dy * vis_s22) * e.pcl_vol_fp;

			n2.fx_int_f += (e_dN2_dx * vis_s11 + e_dN2_dy * vis_s12) * e.pcl_vol_fp;
			n2.fy_int_f += (e_dN2_dx * vis_s12 + e_dN2_dy * vis_s22) * e.pcl_vol_fp;

			n3.fx_int_f += (e_dN3_dx * vis_s11 + e_dN3_dy * vis_s12) * e.pcl_vol_fp;
			n3.fy_int_f += (e_dN3_dx * vis_s12 + e_dN3_dy * vis_s22) * e.pcl_vol_fp;
		}

		if (e.pcl_m_fs != 0.0)
		{
			n_am = one_third * e.pcl_m_fs;
			n1.am_fs += n_am;
			n2.am_fs += n_am;
			n3.am_fs += n_am;
			
			e.ps /= e.pcl_vol_fs;
			e.pcl_density_fs = e.pcl_m_fs / e.pcl_vol_fs;
			e_pcl_f_vol_max = e.pcl_n * e.pcl_vol_s;
			if (e.pcl_vol_fs > e_pcl_f_vol_max)
				e.pcl_vol_fs = e_pcl_f_vol_max;

			n2_miu_div_k_vol = one_third * e.pcl_n * e.pcl_n * md.miu / md.k * e.pcl_vol_fs;
			pcl_v_f = one_third * (n1.vx_f + n2.vx_f + n3.vx_f);
			pcl_v_s = one_third * (n1.vx_s + n2.vx_s + n3.vx_s);
			n2_miu_div_k_vrx_vol = n2_miu_div_k_vol * (pcl_v_f - pcl_v_s);
			pcl_v_f = one_third * (n1.vy_f + n2.vy_f + n3.vy_f);
			pcl_v_s = one_third * (n1.vy_s + n2.vy_s + n3.vy_s);
			n2_miu_div_k_vry_vol = n2_miu_div_k_vol * (pcl_v_f - pcl_v_s);
			
			n1.fx_drag += n2_miu_div_k_vrx_vol;
			n1.fy_drag += n2_miu_div_k_vry_vol;
			n2.fx_drag += n2_miu_div_k_vrx_vol;
			n2.fy_drag += n2_miu_div_k_vry_vol;
			n3.fx_drag += n2_miu_div_k_vrx_vol;
			n3.fy_drag += n2_miu_div_k_vry_vol;

			e_one_n = 1.0 - e.pcl_n;
			p_vol_f = -e.ps * e.pcl_vol_fs;

			n1.fx_int_s += e_dN1_dx * e_one_n * p_vol_f;
			n1.fy_int_s += e_dN1_dy * e_one_n * p_vol_f;
			n1.fx_int_f += e_dN1_dx * e.pcl_n * p_vol_f;
			n1.fy_int_f += e_dN1_dy * e.pcl_n * p_vol_f;

			n2.fx_int_s += e_dN2_dx * e_one_n * p_vol_f;
			n2.fy_int_s += e_dN2_dy * e_one_n * p_vol_f;
			n2.fx_int_f += e_dN2_dx * e.pcl_n * p_vol_f;
			n2.fy_int_f += e_dN2_dy * e.pcl_n * p_vol_f;

			n3.fx_int_s += e_dN3_dx * e_one_n * p_vol_f;
			n3.fy_int_s += e_dN3_dy * e_one_n * p_vol_f;
			n3.fx_int_f += e_dN3_dx * e.pcl_n * p_vol_f;
			n3.fy_int_f += e_dN3_dy * e.pcl_n * p_vol_f;
		}
	}

	// body force
	double bf_mag;
	for (size_t bf_id = 0; bf_id < md.bfsx_num; ++bf_id)
	{
		BodyForceAtPcl &bf = md.bfsxs[bf_id];
		SolidParticle &pcl = md.spcls[bf.pcl_id];
		if (pcl.pe)
		{
			bf_mag = one_third * pcl.m * bf.bf;
			Element &e = *pcl.pe;
			Node &n1 = md.nodes[e.n1];
			Node &n2 = md.nodes[e.n2];
			Node &n3 = md.nodes[e.n3];
			n1.fx_ext_s += bf_mag;
			n2.fx_ext_s += bf_mag;
			n3.fx_ext_s += bf_mag;
		}
	}
	for (size_t bf_id = 0; bf_id < md.bfsy_num; ++bf_id)
	{
		BodyForceAtPcl &bf = md.bfsys[bf_id];
		SolidParticle &pcl = md.spcls[bf.pcl_id];
		if (pcl.pe)
		{
			bf_mag = one_third * pcl.m * bf.bf;
			Element &e = *pcl.pe;
			Node &n1 = md.nodes[e.n1];
			Node &n2 = md.nodes[e.n2];
			Node &n3 = md.nodes[e.n3];
			n1.fy_ext_s += bf_mag;
			n2.fy_ext_s += bf_mag;
			n3.fy_ext_s += bf_mag;
		}
	}
	for (size_t bf_id = 0; bf_id < md.bffx_num; ++bf_id)
	{
		BodyForceAtPcl& bf = md.bffxs[bf_id];
		FluidParticle& pcl = md.fpcls[bf.pcl_id];
		if (pcl.pe)
		{
			bf_mag = one_third * pcl.m * bf.bf;
			Element& e = *pcl.pe;
			Node& n1 = md.nodes[e.n1];
			Node& n2 = md.nodes[e.n2];
			Node& n3 = md.nodes[e.n3];
			n1.fx_ext_f += bf_mag;
			n2.fx_ext_f += bf_mag;
			n3.fx_ext_f += bf_mag;
		}
	}
	for (size_t bf_id = 0; bf_id < md.bffy_num; ++bf_id)
	{
		BodyForceAtPcl& bf = md.bffys[bf_id];
		FluidParticle& pcl = md.fpcls[bf.pcl_id];
		if (pcl.pe)
		{
			bf_mag = one_third * pcl.m * bf.bf;
			Element& e = *pcl.pe;
			Node& n1 = md.nodes[e.n1];
			Node& n2 = md.nodes[e.n2];
			Node& n3 = md.nodes[e.n3];
			n1.fy_ext_f += bf_mag;
			n2.fy_ext_f += bf_mag;
			n3.fy_ext_f += bf_mag;
		}
	}

	// surface force
	for (size_t tf_id = 0; tf_id < md.tx_num; ++tf_id)
	{
		TractionBCAtPcl &tf = md.txs[tf_id];
		SolidParticle &pcl = md.spcls[tf.pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			Node &n1 = md.nodes[e.n1];
			Node& n2 = md.nodes[e.n2];
			Node& n3 = md.nodes[e.n3];
			n1.fx_ext_s += pcl.N1 * tf.t;
			n2.fx_ext_s += pcl.N2 * tf.t;
			n3.fx_ext_s += pcl.N3 * tf.t;
		}
	}
	for (size_t tf_id = 0; tf_id < md.ty_num; ++tf_id)
	{
		TractionBCAtPcl &tf = md.tys[tf_id];
		SolidParticle &pcl = md.spcls[tf.pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			Node &n1 = md.nodes[e.n1];
			Node& n2 = md.nodes[e.n2];
			Node& n3 = md.nodes[e.n3];
			n1.fy_ext_s += pcl.N1 * tf.t;
			n2.fy_ext_s += pcl.N2 * tf.t;
			n3.fy_ext_s += pcl.N3 * tf.t;
		}
	}

	// contact with rigid circle
	if (md.rigid_circle_is_init)
		self.apply_rigid_circle();

	// update nodal acceleration
	double nf, v_sign;
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.am_s != 0.0)
		{
			nf = n.fx_ext_s - n.fx_int_s;
			n.ax_s = (nf + n.fx_drag - self.damping_ratio * abs(nf) * get_sign(n.vx_s)) / n.am_s;
			nf = n.fy_ext_s - n.fy_int_s;
			n.ay_s = (nf + n.fy_drag - self.damping_ratio * abs(nf) * get_sign(n.vy_s)) / n.am_s;
		}
		if (n.am_fs != 0.0 || n.am_fp != 0.0)
		{
			nf = n.fx_ext_f - n.fx_int_f;
			n.ax_f = (nf - n.fx_drag - self.damping_ratio * abs(nf) * get_sign(n.vx_f)) / (n.am_fs + n.am_fp);
			nf = n.fy_ext_f - n.fy_int_f;
			n.ay_f = (nf - n.fy_drag - self.damping_ratio * abs(nf) * get_sign(n.vy_f)) / (n.am_fs + n.am_fp);
		}
	}

	// apply acceleration bc
	for (size_t a_id = 0; a_id < md.asx_num; ++a_id)
	{
		Node &n = md.nodes[md.asxs[a_id].node_id];
		n.ax_s = md.asxs[a_id].a;
	}
	for (size_t a_id = 0; a_id < md.asy_num; ++a_id)
	{
		Node &n = md.nodes[md.asys[a_id].node_id];
		n.ay_s = md.asys[a_id].a;
	}
	for (size_t a_id = 0; a_id < md.afx_num; ++a_id)
	{
		Node &n = md.nodes[md.afxs[a_id].node_id];
		n.ax_f = md.afxs[a_id].a;
	}
	for (size_t a_id = 0; a_id < md.afy_num; ++a_id)
	{
		Node &n = md.nodes[md.afys[a_id].node_id];
		n.ay_f = md.afys[a_id].a;
	}

	// update nodal velocity
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.am_s != 0.0)
		{
			n.vx_s += n.ax_s * self.dtime;
			n.vy_s += n.ay_s * self.dtime;
		}
		if (n.am_fs != 0.0 || n.am_fp != 0.0)
		{
			n.vx_f += n.ax_f * self.dtime;
			n.vy_f += n.ay_f * self.dtime;
		}
	}

	// apply velocity bc
	for (size_t v_id = 0; v_id < md.vsx_num; ++v_id)
	{
		Node &n = md.nodes[md.vsxs[v_id].node_id];
		n.vx_s = md.vsxs[v_id].v;
		n.ax_s = 0.0;
	}
	for (size_t v_id = 0; v_id < md.vsy_num; ++v_id)
	{
		Node &n = md.nodes[md.vsys[v_id].node_id];
		n.vy_s = md.vsys[v_id].v;
		n.ay_s = 0.0;
	}
	for (size_t v_id = 0; v_id < md.vfx_num; ++v_id)
	{
		Node &n =  md.nodes[md.vfxs[v_id].node_id];
		n.vx_f = md.vfxs[v_id].v;
		n.ax_f = 0.0;
	}
	for (size_t v_id = 0; v_id < md.vfy_num; ++v_id)
	{
		Node &n = md.nodes[md.vfys[v_id].node_id];
		n.vy_f = md.vfys[v_id].v;
		n.ay_f = 0.0;
	}

	// update displacement increment
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.am_s != 0.0)
		{
			n.dux_s = n.vx_s * self.dtime;
			n.duy_s = n.vy_s * self.dtime;
		}
		if (n.am_fs != 0.0 || n.am_fp != 0.0)
		{
			n.dux_f = n.vx_f * self.dtime;
			n.duy_f = n.vy_f * self.dtime;
		}
	}

	double de11, de22, de_vol_s;
	double de_vol_fs, de_vol_fp;
	double pcl_m_de_vol;
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
		PointInTriangle& epit = e.pt_in_tri;
		e_dN1_dx = epit.dN1_dx();
		e_dN1_dy = epit.dN1_dy();
		e_dN2_dx = epit.dN2_dx();
		e_dN2_dy = epit.dN2_dy();
		e_dN3_dx = epit.dN3_dx();
		e_dN3_dy = epit.dN3_dy();

		Node& n1 = md.nodes[e.n1];
		Node& n2 = md.nodes[e.n2];
		Node& n3 = md.nodes[e.n3];
		if (e.pcl_m_s != 0.0)
		{
			de11 = n1.dux_s * e_dN1_dx + n2.dux_s * e_dN2_dx + n3.dux_s * e_dN3_dx;
			de22 = n1.duy_s * e_dN1_dy + n2.duy_s * e_dN2_dy + n3.duy_s * e_dN3_dy;
			e.de12 = (n1.dux_s * e_dN1_dy + n2.dux_s * e_dN2_dy + n3.dux_s * e_dN3_dy
					+ n1.duy_s * e_dN1_dx + n2.duy_s * e_dN2_dx + n3.duy_s * e_dN3_dx) * 0.5;
			// volumetric strain of solid phase
			de_vol_s = one_third * (de11 + de22);
			pcl_m_de_vol = de_vol_s * e.pcl_m_s;
			n1.de_vol_s += pcl_m_de_vol;
			n2.de_vol_s += pcl_m_de_vol;
			n3.de_vol_s += pcl_m_de_vol;
			e.dde11 = de11 - de_vol_s;
			e.dde22 = de22 - de_vol_s;
		}
		if (e.pcl_m_fp != 0.0)
		{
			de_vol_fp = -(n1.dux_f * e_dN1_dx + n2.dux_f * e_dN2_dx + n3.dux_f * e_dN3_dx
						+ n1.duy_f * e_dN1_dy + n2.duy_f * e_dN2_dy + n3.duy_f * e_dN3_dy);
			pcl_m_de_vol = one_third * de_vol_fp * e.pcl_m_fp;
			n1.de_vol_fp += pcl_m_de_vol;
			n2.de_vol_fp += pcl_m_de_vol;
			n3.de_vol_fp += pcl_m_de_vol;
		}
		if (e.pcl_m_fs != 0.0)
		{
			// volumetric strain of fluid phase, compression as positive
			de_vol_fs = (1.0 - e.pcl_n) / e.pcl_n * -de_vol_s
					   - (n1.dux_f * e_dN1_dx + n2.dux_f * e_dN2_dx + n3.dux_f * e_dN3_dx
						+ n1.duy_f * e_dN1_dy + n2.duy_f * e_dN2_dy + n3.duy_f * e_dN3_dy);
			pcl_m_de_vol = one_third * de_vol_fs * e.pcl_m_fs;
			n1.de_vol_fs += pcl_m_de_vol;
			n2.de_vol_fs += pcl_m_de_vol;
			n3.de_vol_fs += pcl_m_de_vol;
		}
	}

	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.am_s != 0.0)
			n.de_vol_s /= n.am_s;
		if (n.am_fs != 0.0)
			n.de_vol_fs /= n.am_fs;
		if (n.am_fp != 0.0)
			n.de_vol_fp /= n.am_fp;
	}

	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
		Node& n1 = md.nodes[e.n1];
		Node& n2 = md.nodes[e.n2];
		Node& n3 = md.nodes[e.n3];
		if (e.pcl_m_s != 0.0)
		{
			e.de_vol_s = one_third * (n1.de_vol_s + n2.de_vol_s + n3.de_vol_s);
			e.pcl_n = (e.de_vol_s + e.pcl_n) / (1.0 + e.de_vol_s);
		}
		if (e.pcl_m_fp != 0.0)
		{
			e.de_vol_fp = one_third * (n1.de_vol_fp + n2.de_vol_fp + n2.de_vol_fp);
			e.pcl_density_fp /= 1.0 - e.de_vol_fp;
			e.pp += md.Kf * e.de_vol_fp;
		}
		if (e.pcl_m_fs != 0.0)
		{
			e.de_vol_fs = one_third * (n1.de_vol_fs + n2.de_vol_fs + n3.de_vol_fs);
			e.pcl_density_fs /= 1.0 - e.de_vol_fs;
			e.ps += md.Kf * e.de_vol_fs;
		}
	}
	
	// map variables back to pcl
	double de12;
	double dstrain[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	int mm_res;
	const double* dstress;
	for (size_t pcl_id = 0; pcl_id < md.spcl_num; ++pcl_id)
	{
		SolidParticle &pcl = md.spcls[pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			Node &n1 = md.nodes[e.n1];
			Node &n2 = md.nodes[e.n2];
			Node &n3 = md.nodes[e.n3];

			pcl.vx += (n1.ax_s * pcl.N1 + n2.ax_s * pcl.N2 + n3.ax_s * pcl.N3) * self.dtime;
			pcl.vy += (n1.ay_s * pcl.N1 + n2.ay_s * pcl.N2 + n3.ay_s * pcl.N3) * self.dtime;
			pcl.ux += n1.dux_s * pcl.N1 + n2.dux_s * pcl.N2 + n3.dux_s * pcl.N3;
			pcl.uy += n1.duy_s * pcl.N1 + n2.duy_s * pcl.N2 + n3.duy_s * pcl.N3;
			pcl.x = pcl.x_ori + pcl.ux;
			pcl.y = pcl.y_ori + pcl.uy;

			de_vol_s = one_third * e.de_vol_s;
			de11 = e.dde11 + de_vol_s;
			de22 = e.dde22 + de_vol_s;
			de12 = e.de12;
			pcl.e11 += de11;
			pcl.e22 += de22;
			pcl.e12 += de12;

			dstrain[0] = de11;
			dstrain[1] = de22;
			dstrain[3] = de12;
			mm_res = pcl.mm->integrate(dstrain);
			dstress = pcl.mm->get_dstress();
			pcl.s11 += dstress[0];
			pcl.s22 += dstress[1];
			pcl.s12 += dstress[3];

			pcl.n = e.pcl_n;
		}
	}
	
	for (size_t pcl_id = 0; pcl_id < md.fpcl_num; ++pcl_id)
	{
		FluidParticle& pcl = md.fpcls[pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			Node &n1 = md.nodes[e.n1];
			Node &n2 = md.nodes[e.n2];
			Node &n3 = md.nodes[e.n3];

			pcl.vx += (n1.ax_f * pcl.N1 + n2.ax_f * pcl.N2 + n3.ax_f * pcl.N3) * self.dtime;
			pcl.vy += (n1.ay_f * pcl.N1 + n2.ay_f * pcl.N2 + n3.ay_f * pcl.N3) * self.dtime;
			pcl.ux += n1.dux_f * pcl.N1 + n2.dux_f * pcl.N2 + n3.dux_f * pcl.N3;
			pcl.uy += n1.duy_f * pcl.N1 + n2.duy_f * pcl.N2 + n3.duy_f * pcl.N3;
			pcl.x = pcl.x_ori + pcl.ux;
			pcl.y = pcl.y_ori + pcl.uy;

			if (pcl.is_seepage)
			{
				pcl.p = e.ps;
				pcl.density = e.pcl_density_fs;
			}
			else
			{
				pcl.p = e.pp;
				pcl.density = e.pcl_density_fp;
			}
		}
	}

	return 0;
}

int Step_T2D_CHM_d::apply_rigid_circle()
{
	Model_T2D_CHM_d& md = *model;
	RigidCircle& rc = md.rigid_circle;

	rc.reset_rf();

	double dist, norm_x, norm_y;
	double f_cont, fx_cont, fy_cont;
	for (size_t p_id = 0; p_id < md.spcl_num; ++p_id)
	{
		SolidParticle& pcl = md.spcls[p_id];
		if (pcl.pe && rc.detect_collision_with_point(
						pcl.x, pcl.y, pcl.vol, dist, norm_x, norm_y))
		{
			f_cont = md.Ks_cont * dist;
			fx_cont = f_cont * norm_x;
			fy_cont = f_cont * norm_y;
			rc.add_rf(pcl.x, pcl.y,	-fx_cont, -fy_cont);
			Element& e = *pcl.pe;
			Node& n1 = md.nodes[e.n1];
			n1.fx_ext_s += pcl.N1 * fx_cont;
			n1.fy_ext_s += pcl.N1 * fy_cont;
			Node& n2 = md.nodes[e.n2];
			n2.fx_ext_s += pcl.N2 * fx_cont;
			n2.fy_ext_s += pcl.N2 * fy_cont;
			Node& n3 = md.nodes[e.n3];
			n3.fx_ext_s += pcl.N3 * fx_cont;
			n3.fy_ext_s += pcl.N3 * fy_cont;
		}
	}
	for (size_t p_id = 0; p_id < md.fpcl_num; ++p_id)
	{
		FluidParticle& pcl = md.fpcls[p_id];
		if (pcl.pe && rc.detect_collision_with_point(
						pcl.x, pcl.y, pcl.vol, dist, norm_x, norm_y))
		{
			if (pcl.is_seepage)
				f_cont = md.Kfs_cont * dist;
			else
				f_cont = md.Kfp_cont * dist;
			fx_cont = f_cont * norm_x;
			fy_cont = f_cont * norm_y;
			rc.add_rf(pcl.x, pcl.y, -fx_cont, -fy_cont);
			Element& e = *pcl.pe;
			Node& n1 = md.nodes[e.n1];
			n1.fx_ext_f += pcl.N1 * fx_cont;
			n1.fy_ext_f += pcl.N1 * fy_cont;
			Node& n2 = md.nodes[e.n2];
			n2.fx_ext_f += pcl.N2 * fx_cont;
			n2.fy_ext_f += pcl.N2 * fy_cont;
			Node& n3 = md.nodes[e.n3];
			n3.fx_ext_f += pcl.N3 * fx_cont;
			n3.fy_ext_f += pcl.N3 * fy_cont;
		}
	}

	rc.update_motion(dtime);
	return 0;
}
