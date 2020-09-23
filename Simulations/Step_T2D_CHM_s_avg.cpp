#include "Simulations_pcp.h"

#include "Step_T2D_CHM_s.h"

#define one_third (1.0/3.0)

int solve_substep_T2D_CHM_s_avg(void *_self)
{
	typedef Model_T2D_CHM_s::Node Node;
	typedef Model_T2D_CHM_s::Element Element;
	typedef Model_T2D_CHM_s::Particle Particle;
	Step_T2D_CHM_s &self = *(Step_T2D_CHM_s *)(_self);
	Model_T2D_CHM_s &md = *self.model;

	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		n.has_mp = false;
		// solid phase
		n.am_s = 0.0;
		n.vm_s = 0.0;
		n.vx_s = 0.0;
		n.vy_s = 0.0;
		n.fx_ext_s = 0.0;
		n.fy_ext_s = 0.0;
		n.fx_int_s = 0.0;
		n.fy_int_s = 0.0;
		// fluid phase
		n.am_f = 0.0;
		n.vm_f = 0.0;
		n.vx_f = 0.0;
		n.vy_f = 0.0;
		n.fx_ext_f = 0.0;
		n.fy_ext_f = 0.0;
		n.fx_int_f = 0.0;
		n.fy_int_f = 0.0;
		// solid - fluid interaction
		n.fx_drag = 0.0;
		n.fy_drag = 0.0;
		// strain enhancement
		n.de_vol_s = 0.0;
		n.de_vol_f = 0.0;
	}

	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
		e.has_mp = false;
		e.pcl_m_s = 0.0;
		e.pcl_m_f = 0.0;
		e.pcl_n = 0.0;
		e.pcl_vol_f = 0.0;
		e.pcl_vol = 0.0;
		e.s11 = 0.0;
		e.s22 = 0.0;
		e.s12 = 0.0;
		e.p = 0.0;
	}

	double pcl_vol_f, N_m_s, N_m_f;
	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle &pcl = md.pcls[pcl_id];
		if (pcl.pe)
		{
			if (!(pcl.pe = md.find_in_which_element(pcl)))
				continue;

			pcl.vol = pcl.vol_s / (1.0 - pcl.n);
			pcl_vol_f = pcl.n * pcl.vol;
			pcl.m_f = pcl.density_f * pcl_vol_f;
			
			Element &e = *pcl.pe;
			e.has_mp = true;
			e.pcl_m_s += pcl.m_s;
			e.pcl_m_f += pcl.m_f;
			e.pcl_n += pcl.vol_s;
			e.pcl_vol_f += pcl_vol_f;
			e.pcl_vol += pcl.vol;
			e.s11 += pcl.s11 * pcl.vol;
			e.s22 += pcl.s22 * pcl.vol;
			e.s12 += pcl.s12 * pcl.vol;
			e.p += pcl.p * pcl.vol;

			Node &n1 = md.nodes[e.n1];
			n1.has_mp = true;
			// solid phase
			N_m_s = pcl.N1 * pcl.m_s;
			n1.vm_s += N_m_s;
			n1.vx_s += N_m_s * pcl.vx_s;
			n1.vy_s += N_m_s * pcl.vy_s;
			// fluid phase
			N_m_f = pcl.N1 * pcl.m_f;
			n1.vm_f += N_m_f;
			n1.vx_f += N_m_f * pcl.vx_f;
			n1.vy_f += N_m_f * pcl.vy_f;

			Node &n2 = md.nodes[e.n2];
			n2.has_mp = true;
			// mixture phase
			N_m_s = pcl.N2 * pcl.m_s;
			n2.vm_s += N_m_s;
			n2.vx_s += N_m_s * pcl.vx_s;
			n2.vy_s += N_m_s * pcl.vy_s;
			// fluid phase
			N_m_f = pcl.N2 * pcl.m_f;
			n2.vm_f += N_m_f;
			n2.vx_f += N_m_f * pcl.vx_f;
			n2.vy_f += N_m_f * pcl.vy_f;

			Node &n3 = md.nodes[e.n3];
			n3.has_mp = true;
			// mixture phase
			N_m_s = pcl.N3 * pcl.m_s;
			n3.vm_s += N_m_s;
			n3.vx_s += N_m_s * pcl.vx_s;
			n3.vy_s += N_m_s * pcl.vy_s;
			// fluid phase
			N_m_f = pcl.N3 * pcl.m_f;
			n3.vm_f += N_m_f;
			n3.vx_f += N_m_f * pcl.vx_f;
			n3.vy_f += N_m_f * pcl.vy_f;
		}
	}

	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node& n = md.nodes[n_id];
		if (n.has_mp)
		{
			n.vx_s /= n.vm_s;
			n.vy_s /= n.vm_s;
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

	double pcl_m_s_by_3, pcl_m_f_by_3;
	double n2_miu_div_k;
	double n2_miu_div_k_vrx_vol, n2_miu_div_k_vry_vol;
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
		if (e.has_mp)
		{
			e.pcl_n = 1.0 - e.pcl_n / e.pcl_vol;
			e.pcl_density_f = e.pcl_m_f / e.pcl_vol_f;
			e.s11 /= e.pcl_vol;
			e.s22 /= e.pcl_vol;
			e.s12 /= e.pcl_vol;
			e.p /= e.pcl_vol;
			if (e.pcl_vol > e.area)
				e.pcl_vol = e.area;

			Node &n1 = md.nodes[e.n1];
			Node &n2 = md.nodes[e.n2];
			Node &n3 = md.nodes[e.n3];

			pcl_m_s_by_3 = one_third * e.pcl_m_s;
			n1.am_s += pcl_m_s_by_3;
			n2.am_s += pcl_m_s_by_3;
			n3.am_s += pcl_m_s_by_3;
			pcl_m_f_by_3 = one_third * e.pcl_m_f;
			n1.am_f += pcl_m_f_by_3;
			n2.am_f += pcl_m_f_by_3;
			n3.am_f += pcl_m_f_by_3;
			
			// solid - fluid interaction
			n2_miu_div_k = e.pcl_n * e.pcl_n * md.miu / md.k;
			n2_miu_div_k_vrx_vol = one_third * one_third * n2_miu_div_k
				* (n1.vx_f + n2.vx_f + n3.vx_f - n1.vx_s - n2.vx_s - n3.vx_s) * e.pcl_vol;
			n1.fx_drag += n2_miu_div_k_vrx_vol;
			n2.fx_drag += n2_miu_div_k_vrx_vol;
			n3.fx_drag += n2_miu_div_k_vrx_vol;
			n2_miu_div_k_vry_vol = one_third * one_third * n2_miu_div_k
				* (n1.vy_f + n2.vy_f + n3.vy_f - n1.vy_s - n2.vy_s - n3.vy_s) * e.pcl_vol;
			n1.fy_drag += n2_miu_div_k_vry_vol;
			n2.fy_drag += n2_miu_div_k_vry_vol;
			n3.fy_drag += n2_miu_div_k_vry_vol;

			// internal force
			// node1			
			n1.fx_int_s += (e.dN1_dx * (e.s11 - (1.0 - e.pcl_n) * e.p) + e.dN1_dy * e.s12) * e.pcl_vol;
			n1.fy_int_s += (e.dN1_dx * e.s12 + e.dN1_dy * (e.s22 - (1.0 - e.pcl_n) * e.p)) * e.pcl_vol;
			n1.fx_int_f += (e.dN1_dx * e.pcl_n * -e.p) * e.pcl_vol;
			n1.fy_int_f += (e.dN1_dy * e.pcl_n * -e.p) * e.pcl_vol;
			// node2
			n2.fx_int_s += (e.dN2_dx * (e.s11 - (1.0 - e.pcl_n) * e.p) + e.dN2_dy * e.s12) * e.pcl_vol;
			n2.fy_int_s += (e.dN2_dx * e.s12 + e.dN2_dy * (e.s22 - (1.0 - e.pcl_n) * e.p)) * e.pcl_vol;
			n2.fx_int_f += (e.dN2_dx * e.pcl_n * -e.p) * e.pcl_vol;
			n2.fy_int_f += (e.dN2_dy * e.pcl_n * -e.p) * e.pcl_vol;
			// node3
			n3.fx_int_s += (e.dN3_dx * (e.s11 - (1.0 - e.pcl_n) * e.p) + e.dN3_dy * e.s12) * e.pcl_vol;
			n3.fy_int_s += (e.dN3_dx * e.s12 + e.dN3_dy * (e.s22 - (1.0 - e.pcl_n) * e.p)) * e.pcl_vol;
			n3.fx_int_f += (e.dN3_dx * e.pcl_n * -e.p) * e.pcl_vol;
			n3.fy_int_f += (e.dN3_dy * e.pcl_n * -e.p) * e.pcl_vol;
		}
	}

	// body force
	double bf_s, bf_f;
	for (size_t bf_id = 0; bf_id < md.bfx_num; ++bf_id)
	{
		BodyForceAtPcl &bf = md.bfxs[bf_id];
		Particle &pcl = md.pcls[bf.pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			//bf_s = one_third * pcl.m_s * bf.bf;
			//bf_f = one_third * pcl.m_f * bf.bf;
			bf_s = one_third * (pcl.density_s - pcl.density_f) * (1.0 - pcl.n) * pcl.vol * bf.bf;
			bf_f = one_third * 0.0;
			Node &n1 = md.nodes[e.n1];
			n1.fx_ext_s += bf_s;
			n1.fx_ext_f += bf_f;
			Node &n2 = md.nodes[e.n2];
			n2.fx_ext_s += bf_s;
			n2.fx_ext_f += bf_f;
			Node &n3 = md.nodes[e.n3];
			n3.fx_ext_s += bf_s;
			n3.fx_ext_f += bf_f;
		}
	}
	for (size_t bf_id = 0; bf_id < md.bfy_num; ++bf_id)
	{
		BodyForceAtPcl &bf = md.bfys[bf_id];
		Particle &pcl = md.pcls[bf.pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			//bf_s = one_third * pcl.m_s * bf.bf;
			//bf_f = one_third * pcl.m_f * bf.bf;
			bf_s = one_third * (pcl.density_s - pcl.density_f) * (1.0 - pcl.n) * pcl.vol * bf.bf;
			bf_f = one_third * 0.0;
			Node &n1 = md.nodes[e.n1];
			n1.fy_ext_s += bf_s;
			n1.fy_ext_f += bf_f;
			Node &n2 = md.nodes[e.n2];
			n2.fy_ext_s += bf_s;
			n2.fy_ext_f += bf_f;
			Node &n3 = md.nodes[e.n3];
			n3.fy_ext_s += bf_s;
			n3.fy_ext_f += bf_f;
		}
	}

	// surface force
	for (size_t tf_id = 0; tf_id < md.tx_num; ++tf_id)
	{
		TractionBCAtPcl &tf = md.txs[tf_id];
		Particle &pcl = md.pcls[tf.pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			Node &n1 = md.nodes[e.n1];
			Node &n2 = md.nodes[e.n2];
			Node &n3 = md.nodes[e.n3];
			n1.fx_ext_s += pcl.N1 * tf.t;
			n2.fx_ext_s += pcl.N2 * tf.t;
			n3.fx_ext_s += pcl.N3 * tf.t;
		}
	}
	for (size_t tf_id = 0; tf_id < md.ty_num; ++tf_id)
	{
		TractionBCAtPcl &tf = md.tys[tf_id];
		Particle &pcl = md.pcls[tf.pcl_id];
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

	// update acceleration
	double nf, v_sign;
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.has_mp)
		{
			nf = n.fx_ext_s - n.fx_int_s;
			n.ax_s = (nf + n.fx_drag - self.damping_ratio * abs(nf) * get_sign(n.vx_s)) / n.am_s;
			nf = n.fy_ext_s - n.fy_int_s;
			n.ay_s = (nf + n.fy_drag - self.damping_ratio * abs(nf) * get_sign(n.vy_s)) / n.am_s;
			nf = n.fx_ext_f - n.fx_int_f;
			n.ax_f = (nf - n.fx_drag - self.damping_ratio * abs(nf) * get_sign(n.vx_f)) / n.am_f;
			nf = n.fy_ext_f - n.fy_int_f;
			n.ay_f = (nf - n.fy_drag - self.damping_ratio * abs(nf) * get_sign(n.vy_f)) / n.am_f;
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

	// update velocity
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.has_mp)
		{
			n.vx_s += n.ax_s * self.dtime;
			n.vy_s += n.ay_s * self.dtime;
			n.vx_f += n.ax_f * self.dtime;
			n.vy_f += n.ay_f * self.dtime;
		}
	}

	// contact with rigid circle
	if (md.rigid_circle_is_init)
		self.apply_rigid_circle_avg(self.dtime);

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
		if (n.has_mp)
		{
			n.dux_s = n.vx_s * self.dtime;
			n.duy_s = n.vy_s * self.dtime;
			n.dux_f = n.vx_f * self.dtime;
			n.duy_f = n.vy_f * self.dtime;
		}
	}

	double de11, de22, de_vol_s_by_3;
	double vol_de_vol_s, vol_de_vol_f;
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
		if (e.has_mp)
		{
			Node &n1 = md.nodes[e.n1];
			Node &n2 = md.nodes[e.n2];
			Node &n3 = md.nodes[e.n3];

			de11 = n1.dux_s * e.dN1_dx + n2.dux_s * e.dN2_dx + n3.dux_s * e.dN3_dx;
			de22 = n1.duy_s * e.dN1_dy + n2.duy_s * e.dN2_dy + n3.duy_s * e.dN3_dy;
			e.de12 = (n1.dux_s * e.dN1_dy + n2.dux_s * e.dN2_dy + n3.dux_s * e.dN3_dy
			 		+ n1.duy_s * e.dN1_dx + n2.duy_s * e.dN2_dx + n3.duy_s * e.dN3_dx) * 0.5;
			e.de_vol_s = de11 + de22;
			de_vol_s_by_3 = one_third * e.de_vol_s;
			e.dde11 = de11 - de_vol_s_by_3;
			e.dde22 = de22 - de_vol_s_by_3;

			e.de_vol_f = (1.0 - e.pcl_n) / e.pcl_n * -e.de_vol_s
					   - (n1.dux_f * e.dN1_dx + n2.dux_f * e.dN2_dx + n3.dux_f * e.dN3_dx
						+ n1.duy_f * e.dN1_dy + n2.duy_f * e.dN2_dy + n3.duy_f * e.dN3_dy);
		
			vol_de_vol_s = one_third * e.de_vol_s * e.pcl_m_s;
			n1.de_vol_s += vol_de_vol_s;
			n2.de_vol_s += vol_de_vol_s;
			n3.de_vol_s += vol_de_vol_s;

			vol_de_vol_f = one_third * e.de_vol_f * e.pcl_m_f;
			n1.de_vol_f += vol_de_vol_f;
			n2.de_vol_f += vol_de_vol_f;
			n3.de_vol_f += vol_de_vol_f;
		}
	}

	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.has_mp)
		{
			n.de_vol_s /= n.am_s;
			n.de_vol_f /= n.am_f;
		}
	}

	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
		if (e.has_mp)
		{
			Node &n1 = md.nodes[e.n1];
			Node &n2 = md.nodes[e.n2];
			Node &n3 = md.nodes[e.n3];
			e.de_vol_s = one_third * (n1.de_vol_s + n2.de_vol_s + n3.de_vol_s);
			e.pcl_n = (e.de_vol_s + e.pcl_n) / (1.0 + e.de_vol_s);
			e.de_vol_f = one_third * (n1.de_vol_f + n2.de_vol_f + n3.de_vol_f);
			e.pcl_density_f /= 1.0 - e.de_vol_f;
			e.p += md.Kf * e.de_vol_f;
		}
	}

	double de12;
	double dstrain[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	int mm_res;
	const double* dstress;
	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle &pcl = md.pcls[pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			Node &n1 = md.nodes[e.n1];
			Node &n2 = md.nodes[e.n2];
			Node &n3 = md.nodes[e.n3];

			pcl.vx_s += (n1.ax_s * pcl.N1 + n2.ax_s * pcl.N2 + n3.ax_s * pcl.N3) * self.dtime;
			pcl.vy_s += (n1.ay_s * pcl.N1 + n2.ay_s * pcl.N2 + n3.ay_s * pcl.N3) * self.dtime;
			pcl.vx_f += (n1.ax_f * pcl.N1 + n2.ax_f * pcl.N2 + n3.ax_f * pcl.N3) * self.dtime;
			pcl.vy_f += (n1.ay_f * pcl.N1 + n2.ay_f * pcl.N2 + n3.ay_f * pcl.N3) * self.dtime;

			pcl.ux_s += n1.dux_s * pcl.N1 + n2.dux_s * pcl.N2 + n3.dux_s * pcl.N3;
			pcl.uy_s += n1.duy_s * pcl.N1 + n2.duy_s * pcl.N2 + n3.duy_s * pcl.N3;
			pcl.ux_f += n1.dux_f * pcl.N1 + n2.dux_f * pcl.N2 + n3.dux_f * pcl.N3;
			pcl.uy_f += n1.duy_f * pcl.N1 + n2.duy_f * pcl.N2 + n3.duy_f * pcl.N3;

			pcl.x = pcl.x_ori + pcl.ux_s;
			pcl.y = pcl.y_ori + pcl.uy_s;

			de_vol_s_by_3 = one_third * e.de_vol_s;
			de11 = e.dde11 + de_vol_s_by_3;
			de22 = e.dde22 + de_vol_s_by_3;
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
			pcl.p = e.p;
			pcl.density_f = e.pcl_density_f;
		}
	}
	
	return 0;
}

int Step_T2D_CHM_s::apply_rigid_circle_avg(double dt)
{
	Model_T2D_CHM_s& md = *model;
	RigidCircle &rc = md.rigid_circle;

	rc.reset_rf();

	double dist, norm_x, norm_y;
	double fs_cont, fsx_cont, fsy_cont;
	double ff_cont, ffx_cont, ffy_cont;
	double nfsx_cont, nfsy_cont, ndasx, ndasy;
	double nffx_cont, nffy_cont, ndafx, ndafy;
	for (size_t p_id = 0; p_id < md.pcl_num; ++p_id)
	{
		Particle& pcl = md.pcls[p_id];
		if (pcl.pe && rc.detect_collision_with_point(
						pcl.x, pcl.y, pcl.vol, dist, norm_x, norm_y))
		{
			fs_cont = md.Ks_cont * dist;
			fsx_cont = fs_cont * norm_x;
			fsy_cont = fs_cont * norm_y;
			ff_cont = md.Kf_cont * dist;
			ffx_cont = ff_cont * norm_x;
			ffy_cont = ff_cont * norm_y;
			// reaction force by the rigid object
			rc.add_rf(pcl.x, pcl.y,
				-(fsx_cont + ffx_cont),
				-(fsy_cont + ffy_cont)
				);
			// adjust velocity at nodes
			Element& e = *pcl.pe;
			// node 1
			Node& n1 = md.nodes[e.n1];
			nfsx_cont = pcl.N1 * fsx_cont;
			nfsy_cont = pcl.N1 * fsy_cont;
			ndasx = nfsx_cont / n1.am_s;
			ndasy = nfsy_cont / n1.am_s;
			n1.ax_s += ndasx;
			n1.ay_s += ndasy;
			n1.vx_s += ndasx * dt;
			n1.vy_s += ndasy * dt;
			nffx_cont = pcl.N1 * ffx_cont;
			nffy_cont = pcl.N1 * ffy_cont;
			ndafx = nffx_cont / n1.am_f;
			ndafy = nffy_cont / n1.am_f;
			n1.ax_f += ndafx;
			n1.ay_f += ndafy;
			n1.vx_f += ndafx * dt;
			n1.vy_f += ndafy * dt;
			// node 2
			Node& n2 = md.nodes[e.n2];
			nfsx_cont = pcl.N2 * fsx_cont;
			nfsy_cont = pcl.N2 * fsy_cont;
			ndasx = nfsx_cont / n2.am_s;
			ndasy = nfsy_cont / n2.am_s;
			n2.ax_s += ndasx;
			n2.ay_s += ndasy;
			n2.vx_s += ndasx * dt;
			n2.vy_s += ndasy * dt;
			nffx_cont = pcl.N2 * ffx_cont;
			nffy_cont = pcl.N2 * ffy_cont;
			ndafx = nffx_cont / n2.am_f;
			ndafy = nffy_cont / n2.am_f;
			n2.ax_f += ndafx;
			n2.ay_f += ndafy;
			n2.vx_f += ndafx * dt;
			n2.vy_f += ndafy * dt;
			// node 3
			Node& n3 = md.nodes[e.n3];
			nfsx_cont = pcl.N3 * fsx_cont;
			nfsy_cont = pcl.N3 * fsy_cont;
			ndasx = nfsx_cont / n3.am_s;
			ndasy = nfsy_cont / n3.am_s;
			n3.ax_s += ndasx;
			n3.ay_s += ndasy;
			n3.vx_s += ndasx * dt;
			n3.vy_s += ndasy * dt;
			nffx_cont = pcl.N3 * ffx_cont;
			nffy_cont = pcl.N3 * ffy_cont;
			ndafx = nffx_cont / n3.am_f;
			ndafy = nffy_cont / n3.am_f;
			n3.ax_f += ndafx;
			n3.ay_f += ndafy;
			n3.vx_f += ndafx * dt;
			n3.vy_f += ndafy * dt;
		}
	}

	rc.update_motion(dt);

	return 0;
}
