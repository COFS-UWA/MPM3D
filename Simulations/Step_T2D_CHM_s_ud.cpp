#include "Simulations_pcp.h"

#include "Step_T2D_CHM_s_ud.h"

Step_T2D_CHM_s_ud::Step_T2D_CHM_s_ud(const char* _name) :
	Step(_name, "Step_T2D_CHM_s_ud", &solve_substep_T2D_CHM_s_ud_avg),
	model(nullptr), damping_ratio(0.0) {}

Step_T2D_CHM_s_ud::~Step_T2D_CHM_s_ud() {}

int Step_T2D_CHM_s_ud::init_calculation()
{
	Model_T2D_CHM_s &md = *model;

	if (is_first_step) {}

	double pcl_vol;
	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle &pcl = md.pcls[pcl_id];
		pcl.pe = md.elems;
		pcl.vol_s = pcl.m_s / pcl.density_s;
		pcl_vol = pcl.vol_s / (1.0 - pcl.n);
		// pcl.m_f is total mass of both phase actually
		pcl.m_f = (pcl.density_f * pcl.n + pcl.density_s * (1.0 - pcl.n)) * pcl_vol;
		pcl.x_ori = pcl.x;
		pcl.y_ori = pcl.y;
		pcl.ux_s = 0.0;
		pcl.uy_s = 0.0;
	}

	return 0;
}

int Step_T2D_CHM_s_ud::finalize_calculation() { return 0; }

int solve_substep_T2D_CHM_s_ud(void *_self)
{
	typedef Model_T2D_CHM_s::Particle Particle;
	typedef Model_T2D_CHM_s::Element Element;
	typedef Model_T2D_CHM_s::Node Node;
	
	Step_T2D_CHM_s_ud &self = *(Step_T2D_CHM_s_ud *)(_self);
	Model_T2D_CHM_s &md = *self.model;

	// init nodes
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		// material point
		n.has_mp = false;
		// only cal mixture face
		n.m_s = 0.0;
		n.vx_s = 0.0;
		n.vy_s = 0.0;
		n.fx_ext_s = 0.0;
		n.fy_ext_s = 0.0;
		n.fx_int_s = 0.0;
		n.fy_int_s = 0.0;
		// strain enhancement approach
		n.pcl_vol = 0.0;
		n.de_vol_s = 0.0;
	}

	// init elements
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
		e.pcls = nullptr;
		e.pcl_vol = 0.0;
		e.pcl_n = 0.0;
		e.s11 = 0.0;
		e.s22 = 0.0;
		e.s12 = 0.0;
		e.p = 0.0;
	}

	// init particles
	double N_m_s, N_m_f;
	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle &pcl = md.pcls[pcl_id];
		if (pcl.pe)
		{
			if (!(pcl.pe = const_cast<Element *>(md.find_in_which_element(pcl))))
				continue;
			pcl.pe->add_pcl(pcl);

			pcl.vol = pcl.vol_s / (1.0 - pcl.n);

			Element &e = *pcl.pe;
			e.pcl_vol += pcl.vol;
			e.pcl_n += pcl.vol_s;
			e.s11 += pcl.vol * pcl.s11;
			e.s22 += pcl.vol * pcl.s22;
			e.s12 += pcl.vol * pcl.s12;
			e.p += pcl.vol * pcl.p;

			// node 1
			Node &n1 = md.nodes[e.n1];
			n1.has_mp = true;
			// mixture phase
			N_m_s = pcl.N1 * pcl.m_f;
			n1.m_s += N_m_s;
			n1.vx_s += N_m_s * pcl.vx_s;
			n1.vy_s += N_m_s * pcl.vy_s;

			// node 2
			Node &n2 = md.nodes[e.n2];
			n2.has_mp = true;
			// mixture phase
			N_m_s = pcl.N2 * pcl.m_f;
			n2.m_s += N_m_s;
			n2.vx_s += N_m_s * pcl.vx_s;
			n2.vy_s += N_m_s * pcl.vy_s;

			// node 3
			Node &n3 = md.nodes[e.n3];
			n3.has_mp = true;
			// mixture phase
			N_m_s = pcl.N3 * pcl.m_f;
			n3.m_s += N_m_s;
			n3.vx_s += N_m_s * pcl.vx_s;
			n3.vy_s += N_m_s * pcl.vy_s;
		}
	}

	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
		if (e.pcls)
		{
			e.pcl_n = 1.0 - e.pcl_n / e.pcl_vol;
			e.s11 /= e.pcl_vol;
			e.s22 /= e.pcl_vol;
			e.s12 /= e.pcl_vol;
			e.p /= e.pcl_vol;
			if (e.pcl_vol > e.area)
				e.pcl_vol = e.area;

			Node &n1 = md.nodes[e.n1];
			Node &n2 = md.nodes[e.n2];
			Node &n3 = md.nodes[e.n3];
						
			// node 1
			n1.fx_int_s += (e.dN1_dx * (e.s11 - e.p) + e.dN1_dy * e.s12) * e.pcl_vol;
			n1.fy_int_s += (e.dN1_dx * e.s12 + e.dN1_dy * (e.s22 - e.p)) * e.pcl_vol;
			// node 2
			n2.fx_int_s += (e.dN2_dx * (e.s11 - e.p) + e.dN2_dy * e.s12) * e.pcl_vol;
			n2.fy_int_s += (e.dN2_dx * e.s12 + e.dN2_dy * (e.s22 - e.p)) * e.pcl_vol;
			// node 3
			n3.fx_int_s += (e.dN3_dx * (e.s11 - e.p) + e.dN3_dy * e.s12) * e.pcl_vol;
			n3.fy_int_s += (e.dN3_dx * e.s12 + e.dN3_dy * (e.s22 - e.p)) * e.pcl_vol;
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
			bf_s = pcl.m_f * bf.bf;
			// node 1
			Node &n1 = md.nodes[e.n1];
			n1.fx_ext_s += pcl.N1 * bf_s;
			// node 2
			Node &n2 = md.nodes[e.n2];
			n2.fx_ext_s += pcl.N2 * bf_s;
			// node 3
			Node &n3 = md.nodes[e.n3];
			n3.fx_ext_s += pcl.N3 * bf_s;
		}
	}
	for (size_t bf_id = 0; bf_id < md.bfy_num; ++bf_id)
	{
		BodyForceAtPcl &bf = md.bfys[bf_id];
		Particle &pcl = md.pcls[bf.pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			bf_s = pcl.m_f * bf.bf;
			// node 1
			Node &n1 = md.nodes[e.n1];
			n1.fy_ext_s += pcl.N1 * bf_s;
			// node 2
			Node &n2 = md.nodes[e.n2];
			n2.fy_ext_s += pcl.N2 * bf_s;
			// node 3
			Node &n3 = md.nodes[e.n3];
			n3.fy_ext_s += pcl.N3 * bf_s;
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
	double nf;
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.has_mp)
		{
			// fx_s
			nf = n.fx_ext_s - n.fx_int_s;
			n.ax_s = (nf - self.damping_ratio * abs(nf) * get_sign(n.vx_s)) / n.m_s;
			// fy_s
			nf = n.fy_ext_s - n.fy_int_s;
			n.ay_s = (nf - self.damping_ratio * abs(nf) * get_sign(n.vy_s)) / n.m_s;
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

	// update nodal momentum
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.has_mp)
		{
			n.vx_s /= n.m_s;
			n.vy_s /= n.m_s;
			n.vx_s += n.ax_s * self.dtime;
			n.vy_s += n.ay_s * self.dtime;
		}
	}

	// contact with rigid circle
	if (md.rigid_circle_is_init)
		self.apply_rigid_circle(self.dtime);

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

	// update displacement increment of both phases
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.has_mp)
		{
			// solid phase
			n.dux_s = n.vx_s * self.dtime;
			n.duy_s = n.vy_s * self.dtime;
		}
	}

	// map variables back to particles and update their variables
	double de11, de22, de_mean;
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
		if (e.pcls)
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
			de_mean = e.de_vol_s / 3.0;
			e.dde11 = de11 - de_mean;
			e.dde22 = de22 - de_mean;
		}
	}

	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle &pcl = md.pcls[pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			double vol_de_vol_s = pcl.vol * e.de_vol_s;
			// node 1
			Node &n1 = md.nodes[e.n1];
			n1.pcl_vol += pcl.N1 * pcl.vol;
			n1.de_vol_s += pcl.N1 * vol_de_vol_s;
			// node 2
			Node &n2 = md.nodes[e.n2];
			n2.pcl_vol += pcl.N2 * pcl.vol;
			n2.de_vol_s += pcl.N2 * vol_de_vol_s;
			// node 3
			Node &n3 = md.nodes[e.n3];
			n3.pcl_vol += pcl.N3 * pcl.vol;
			n3.de_vol_s += pcl.N3 * vol_de_vol_s;
		}
	}

	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.has_mp)
			n.de_vol_s /= n.pcl_vol;
	}

	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
		if (e.pcls)
		{
			Node &n1 = md.nodes[e.n1];
			Node &n2 = md.nodes[e.n2];
			Node &n3 = md.nodes[e.n3];
			e.de_vol_s = (n1.de_vol_s + n2.de_vol_s + n3.de_vol_s) / 3.0;
			e.p += md.Kf * -e.de_vol_s / e.pcl_n;
		}
	}

	double de12;
	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle &pcl = md.pcls[pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			Node &n1 = md.nodes[e.n1];
			Node &n2 = md.nodes[e.n2];
			Node &n3 = md.nodes[e.n3];

			// velocity
			pcl.vx_s += (n1.ax_s * pcl.N1 + n2.ax_s * pcl.N2 + n3.ax_s * pcl.N3) * self.dtime;
			pcl.vy_s += (n1.ay_s * pcl.N1 + n2.ay_s * pcl.N2 + n3.ay_s * pcl.N3) * self.dtime;

			// displacement
			pcl.ux_s += n1.dux_s * pcl.N1 + n2.dux_s * pcl.N2 + n3.dux_s * pcl.N3;
			pcl.uy_s += n1.duy_s * pcl.N1 + n2.duy_s * pcl.N2 + n3.duy_s * pcl.N3;

			// update position
			pcl.x = pcl.x_ori + pcl.ux_s;
			pcl.y = pcl.y_ori + pcl.uy_s;

			// strain
			de_mean = e.de_vol_s / 3.0;
			de11 = e.dde11 + de_mean;
			de22 = e.dde22 + de_mean;
			de12 = e.de12;
			// strain
			pcl.e11 += de11;
			pcl.e22 += de22;
			pcl.e12 += de12;

			// update stress using constitutive model
			double dstrain[6] = { de11, de22, 0.0, de12, 0.0, 0.0 };
			pcl.mm->integrate(dstrain);
			const double *dstress = pcl.mm->get_dstress();
			pcl.s11 += dstress[0];
			pcl.s22 += dstress[1];
			pcl.s12 += dstress[3];

			// pore pressure
			pcl.p = e.p;

			// fluid density
			pcl.density_f /= 1.0 - e.de_vol_s / e.pcl_n;
			// porosity
			pcl.n = (e.de_vol_s + e.pcl_n) / (1.0 + e.de_vol_s);
		}
	}
	
	return 0;
}

void Step_T2D_CHM_s_ud::apply_rigid_circle(double dtime)
{
	Model_T2D_CHM_s& md = *model;
	RigidCircle& rc = md.get_rigid_circle();

	rc.reset_rf();
	
	double dist, norm_x, norm_y;
	double fs_cont, fsx_cont, fsy_cont;
	double nfsx_cont, nfsy_cont, ndasx, ndasy;
	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle& pcl = md.pcls[pcl_id];
		if (pcl.pe && rc.detect_collision_with_point(
			pcl.x, pcl.y, pcl.vol, dist, norm_x, norm_y))
		{
			fs_cont = md.Ks_cont * dist;
			fsx_cont = fs_cont * norm_x;
			fsy_cont = fs_cont * norm_y;
			// reaction force by the rigid object
			rc.add_rf(pcl.x, pcl.y,	-fsx_cont, -fsy_cont);
			// adjust velocity at nodes
			Element& e = *pcl.pe;
			// node 1
			Node& n1 = md.nodes[e.n1];
			nfsx_cont = pcl.N1 * fsx_cont;
			nfsy_cont = pcl.N1 * fsy_cont;
			ndasx = nfsx_cont / n1.m_s;
			ndasy = nfsy_cont / n1.m_s;
			n1.ax_s += ndasx;
			n1.ay_s += ndasy;
			n1.vx_s += ndasx * dtime;
			n1.vy_s += ndasy * dtime;
			// node 2
			Node& n2 = md.nodes[e.n2];
			nfsx_cont = pcl.N2 * fsx_cont;
			nfsy_cont = pcl.N2 * fsy_cont;
			ndasx = nfsx_cont / n2.m_s;
			ndasy = nfsy_cont / n2.m_s;
			n2.ax_s += ndasx;
			n2.ay_s += ndasy;
			n2.vx_s += ndasx * dtime;
			n2.vy_s += ndasy * dtime;
			// node 3
			Node& n3 = md.nodes[e.n3];
			nfsx_cont = pcl.N3 * fsx_cont;
			nfsy_cont = pcl.N3 * fsy_cont;
			ndasx = nfsx_cont / n3.m_s;
			ndasy = nfsy_cont / n3.m_s;
			n3.ax_s += ndasx;
			n3.ay_s += ndasy;
			n3.vx_s += ndasx * dtime;
			n3.vy_s += ndasy * dtime;
		}
	}

	rc.update_motion(dtime);
}
