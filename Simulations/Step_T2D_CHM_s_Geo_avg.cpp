#include "Simulations_pcp.h"

#include "Step_T2D_CHM_s_Geo.h"

#define one_third (1.0f/3.0f)

int solve_substep_T2D_CHM_s_Geo_avg(void *_self)
{
	typedef Model_T2D_CHM_s::Node Node;
	typedef Model_T2D_CHM_s::Element Element;
	typedef Model_T2D_CHM_s::Particle Particle;
	Step_T2D_CHM_s_Geo &self = *(Step_T2D_CHM_s_Geo *)(_self);
	Model_T2D_CHM_s &md = *self.model;

	// init nodes
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		n.am_s = 0.0;
		n.vm_s = 0.0;
		n.vx_s = 0.0;
		n.vy_s = 0.0;
		n.fx_ext_s = 0.0;
		n.fy_ext_s = 0.0;
		n.fx_int_s = 0.0;
		n.fy_int_s = 0.0;
		n.de_vol_s = 0.0;
	}

	// init elements
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
		e.pcl_m_s = 0.0;
		e.pcl_n = 0.0;
		e.pcl_vol = 0.0;
		e.s11 = 0.0;
		e.s22 = 0.0;
		e.s12 = 0.0;
	}

	double N_m_s;
	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle &pcl = md.pcls[pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			e.s11 += pcl.s11 * pcl.vol;
			e.s22 += pcl.s22 * pcl.vol;
			e.s12 += pcl.s12 * pcl.vol;
			e.pcl_vol += pcl.vol;

			// node 1
			Node &n1 = md.nodes[e.n1];
			n1.am_s += one_third * pcl.m_s;
			N_m_s = pcl.N1 * pcl.m_s;
			n1.vm_s += N_m_s;
			n1.vx_s += N_m_s * pcl.vx_s;
			n1.vy_s += N_m_s * pcl.vy_s;

			// node 2
			Node &n2 = md.nodes[e.n2];
			n2.am_s += one_third * pcl.m_s;
			N_m_s = pcl.N2 * pcl.m_s;
			n2.vm_s += N_m_s;
			n2.vx_s += N_m_s * pcl.vx_s;
			n2.vy_s += N_m_s * pcl.vy_s;

			// node 3
			Node &n3 = md.nodes[e.n3];
			n3.am_s += one_third * pcl.m_s;
			N_m_s = pcl.N3 * pcl.m_s;
			n3.vm_s += N_m_s;
			n3.vx_s += N_m_s * pcl.vx_s;
			n3.vy_s += N_m_s * pcl.vy_s;
		}
	}

	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
		if (e.has_mp)
		{
			e.s11 /= e.pcl_vol;
			e.s22 /= e.pcl_vol;
			e.s12 /= e.pcl_vol;
			if (e.pcl_vol > e.area)
				e.pcl_vol = e.area;

			Node &n1 = md.nodes[e.n1];
			Node &n2 = md.nodes[e.n2];
			Node &n3 = md.nodes[e.n3];

			// node 1
			n1.fx_int_s += (e.dN1_dx * e.s11 + e.dN1_dy * e.s12) * e.pcl_vol;
			n1.fy_int_s += (e.dN1_dx * e.s12 + e.dN1_dy * e.s22) * e.pcl_vol;
			// node 2
			n2.fx_int_s += (e.dN2_dx * e.s11 + e.dN2_dy * e.s12) * e.pcl_vol;
			n2.fy_int_s += (e.dN2_dx * e.s12 + e.dN2_dy * e.s22) * e.pcl_vol;
			// node 3
			n3.fx_int_s += (e.dN3_dx * e.s11 + e.dN3_dy * e.s12) * e.pcl_vol;
			n3.fy_int_s += (e.dN3_dx * e.s12 + e.dN3_dy * e.s22) * e.pcl_vol;
		}
	}

	// body force
	double bf_s;
	for (size_t bf_id = 0; bf_id < md.bfx_num; ++bf_id)
	{
		BodyForceAtPcl &bf = md.bfxs[bf_id];
		Particle &pcl = md.pcls[bf.pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			// body force on particle
			bf_s = (pcl.density_s - pcl.density_f) * (1.0 - pcl.n) * pcl.vol * bf.bf;
			// node 1
			Node &n1 = md.nodes[e.n1];
			n1.fx_ext_s += one_third * bf_s;
			// node 2
			Node &n2 = md.nodes[e.n2];
			n2.fx_ext_s += one_third * bf_s;
			// node 3
			Node &n3 = md.nodes[e.n3];
			n3.fx_ext_s += one_third * bf_s;
		}
	}
	for (size_t bf_id = 0; bf_id < md.bfy_num; ++bf_id)
	{
		BodyForceAtPcl &bf = md.bfys[bf_id];
		Particle &pcl = md.pcls[bf.pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			// body force on particle
			bf_s = (pcl.density_s - pcl.density_f) * (1.0 - pcl.n) * pcl.vol* bf.bf;
			// node 1
			Node &n1 = md.nodes[e.n1];
			n1.fy_ext_s += one_third * bf_s;
			// node 2
			Node &n2 = md.nodes[e.n2];
			n2.fy_ext_s += one_third * bf_s;
			// node 3
			Node &n3 = md.nodes[e.n3];
			n3.fy_ext_s += one_third * bf_s;
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
			// node 1
			Node &n1 = md.nodes[e.n1];
			n1.fx_ext_s += pcl.N1 * tf.t;
			// node 2
			Node &n2 = md.nodes[e.n2];
			n2.fx_ext_s += pcl.N2 * tf.t;
			// node 3
			Node &n3 = md.nodes[e.n3];
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
			// node 1
			Node &n1 = md.nodes[e.n1];
			n1.fy_ext_s += pcl.N1 * tf.t;
			// node 2
			Node &n2 = md.nodes[e.n2];
			n2.fy_ext_s += pcl.N2 * tf.t;
			// node 3
			Node &n3 = md.nodes[e.n3];
			n3.fy_ext_s += pcl.N3 * tf.t;
		}
	}

	// update nodal acceleration
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.has_mp)
		{
			n.ax_s = (n.fx_ext_s - n.fx_int_s) / n.am_s;
			n.ay_s = (n.fy_ext_s - n.fy_int_s) / n.am_s;
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

	// update nodal velocity
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.has_mp)
		{
			n.vx_s /= n.vm_s;
			n.vx_s += n.ax_s * self.dtime;
			n.vy_s /= n.vm_s;
			n.vy_s += n.ay_s * self.dtime;
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

	double f_ub = 0.0;
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.has_mp)
			f_ub += n.am_s * n.am_s * (n.ax_s * n.ax_s + n.ay_s * n.ay_s);
	}
	// initial unbalanced force
	if (!self.init_f_ub_is_init)
	{
		self.init_f_ub_is_init = true;
		self.init_f_ub = f_ub;
	}

	if (self.init_f_ub_is_init)
		self.f_ub_ratio = sqrt(f_ub / self.init_f_ub);

	// cal kinetic energy
	double e_kin = 0.0;
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.has_mp)
			e_kin += n.am_s * (n.vx_s * n.vx_s + n.vy_s * n.vy_s);
	}

	// kinetic damping
	if (e_kin < self.e_kin_prev)
	{
		if (!self.e_kin_max_is_init)
		{
			self.e_kin_max_is_init = true;
			self.e_kin_max = self.e_kin_prev;
		}

		// reset pcl velocity
		for (size_t p_id = 0; p_id < md.pcl_num; ++p_id)
		{
			Particle &pcl = md.pcls[p_id];
			pcl.vx_s = 0.0;
			pcl.vy_s = 0.0;
		}
		self.e_kin_prev = 0.0;
	}
	else
	{
		self.e_kin_prev = e_kin;
	}

	if (self.e_kin_max_is_init)
		self.e_kin_ratio = sqrt(self.e_kin_prev / self.e_kin_max);

	// output to file
	//if (self.substep_num % 100 == 1)
	//{
	//	const char *rcd_str = "%lf, %lf, %lf, %lf, %lf, %lf, %lf\n";
	//	char rcd[100];
	//	snprintf(rcd, sizeof(rcd) / sizeof(rcd[0]), rcd_str,
	//		self.current_time, self.f_ub_ratio, self.e_kin_ratio, f_ub, e_kin,
	//		md.pcls[161].y, md.pcls[161].s22);
	//	self.out_file.write(rcd, strlen(rcd));
	//}

	if (e_kin < self.e_kin_prev)
		if (self.e_kin_ratio < self.e_kin_ratio_bound &&
			self.f_ub_ratio  < self.f_ub_ratio_bound)
			return 2;
		else
			return 1;

	// update displacement increment
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.has_mp)
		{
			n.dux_s = n.vx_s * self.dtime;
			n.duy_s = n.vy_s * self.dtime;
		}
	}

	double de11, de22, de_vol_by_3, vol_de_vol;
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
		if (e.has_mp)
		{
			Node &n1 = md.nodes[e.n1];
			Node &n2 = md.nodes[e.n2];
			Node &n3 = md.nodes[e.n3];
			// strain increment
			de11 = n1.dux_s * e.dN1_dx + n2.dux_s * e.dN2_dx + n3.dux_s * e.dN3_dx;
			de22 = n1.duy_s * e.dN1_dy + n2.duy_s * e.dN2_dy + n3.duy_s * e.dN3_dy;
			e.de12 = (n1.dux_s * e.dN1_dy + n2.dux_s * e.dN2_dy + n3.dux_s * e.dN3_dy
			 		+ n1.duy_s * e.dN1_dx + n2.duy_s * e.dN2_dx + n3.duy_s * e.dN3_dx) * 0.5;
			// volumetric strain of solid phase
			e.de_vol_s = de11 + de22;
			de_vol_by_3 = e.de_vol_s / 3.0;
			e.dde11 = de11 - de_vol_by_3;
			e.dde22 = de22 - de_vol_by_3;

			vol_de_vol = one_third * e.de_vol_s * e.pcl_m_s;
			n1.de_vol_s += vol_de_vol;
			n2.de_vol_s += vol_de_vol;
			n3.de_vol_s += vol_de_vol;
		}
	}

	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.has_mp)
			n.de_vol_s /= n.am_s;
	}
	
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element& e = md.elems[e_id];
		if (e.has_mp)
		{
			Node& n1 = md.nodes[e.n1];
			Node& n2 = md.nodes[e.n2];
			Node& n3 = md.nodes[e.n3];
			e.de_vol_s = (n1.de_vol_s + n2.de_vol_s + n3.de_vol_s) * one_third;
		}
	}
	
	double de12, de_vol_s_by_3;
	double dstrain[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	const double* dstress;
	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle &pcl = md.pcls[pcl_id];
		if (pcl.pe)
		{
			Element& e = *pcl.pe;
			Node& n1 = md.nodes[e.n1];
			Node& n2 = md.nodes[e.n2];
			Node& n3 = md.nodes[e.n3];

			// velocity
			pcl.vx_s += (n1.ax_s * pcl.N1 + n2.ax_s * pcl.N2 + n3.ax_s * pcl.N3) * self.dtime;
			pcl.vy_s += (n1.ay_s * pcl.N1 + n2.ay_s * pcl.N2 + n3.ay_s * pcl.N3) * self.dtime;

			// strain
			de_vol_s_by_3 = e.de_vol_s / 3.0;
			de11 = e.dde11 + de_vol_s_by_3;
			de22 = e.dde22 + de_vol_s_by_3;
			de12 = e.de12;

			// update stress
			dstrain[0] = de11;
			dstrain[1] = de22;
			dstrain[3] = de12;
			pcl.mm->integrate(dstrain);
			dstress = pcl.mm->get_dstress();
			pcl.s11 += dstress[0];
			pcl.s22 += dstress[1];
			pcl.s12 += dstress[3];
		}
	}
	
	return 0;
}
