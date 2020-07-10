#include "Simulations_pcp.h"

#include <cmath>

#include "MaterialModel.h"
#include "Step_T2D_CHM_s_Geo.h"

Step_T2D_CHM_s_Geo::Step_T2D_CHM_s_Geo(const char* _name) :
	Step(_name, "Step_T2D_CHM_s_Geo", &solve_substep_T2D_CHM_s_Geo),
	model(nullptr), damping_ratio(0.0) {}

Step_T2D_CHM_s_Geo::~Step_T2D_CHM_s_Geo() {}

int Step_T2D_CHM_s_Geo::init_calculation()
{
	Model_T2D_CHM_s &md = *model;

	if (is_first_step) {}

	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle &pcl = md.pcls[pcl_id];
		pcl.pe = md.find_in_which_element(pcl);
		pcl.vol_s = pcl.m_s / pcl.density_s;
		pcl.vol = pcl.vol_s / (1.0 - pcl.n);
		pcl.m_f = pcl.vol * pcl.n * pcl.density_f;
	}

	// convergence criteria
	// unbalanced force
	init_f_ub = 0.0;
	init_f_ub_is_init = false;
	f_ub_ratio = 1.0;
	// kinematic energy
	e_kin_max = 0.0;
	e_kin_max_is_init = false;
	e_kin_prev = 0.0;
	e_kin_ratio = 1.0;

	//out_file.open("ratio_res.txt", std::ios::out | std::ios::binary);

	return 0;
}

int Step_T2D_CHM_s_Geo::finalize_calculation()
{
	Model_T2D_CHM_s &md = *model;

	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle &pcl = md.pcls[pcl_id];
		pcl.vx_s = 0.0;
		pcl.vy_s = 0.0;
		pcl.vx_f = 0.0;
		pcl.vy_f = 0.0;
		pcl.p = 0.0;
		pcl.e11 = 0.0;
		pcl.e22 = 0.0;
		pcl.e12 = 0.0;
	}

	//out_file.close();

	return 0;
}

int solve_substep_T2D_CHM_s_Geo(void *_self)
{
	typedef Model_T2D_CHM_s::Particle Particle;
	typedef Model_T2D_CHM_s::Element Element;
	typedef Model_T2D_CHM_s::Node Node;
	Step_T2D_CHM_s_Geo &self = *(Step_T2D_CHM_s_Geo *)(_self);
	Model_T2D_CHM_s &md = *self.model;

	// init nodes
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		// material point
		n.has_mp = false;
		// solid phase
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
		e.pcl_vol = 0.0;
		e.pcl_n = 0.0;
		e.s11 = 0.0;
		e.s22 = 0.0;
		e.s12 = 0.0;
		e.pcls = nullptr;
	}

	// init particles
	double N_m_s;
	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle &pcl = md.pcls[pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			e.add_pcl(pcl);
			e.pcl_vol += pcl.vol;
			e.s11 += pcl.vol * pcl.s11;
			e.s22 += pcl.vol * pcl.s22;
			e.s12 += pcl.vol * pcl.s12;

			// node 1
			Node &n1 = md.nodes[e.n1];
			n1.has_mp = true;
			// solid phase
			N_m_s = pcl.N1 * pcl.m_s;
			n1.m_s += N_m_s;
			n1.vx_s += N_m_s * pcl.vx_s;
			n1.vy_s += N_m_s * pcl.vy_s;

			// node 2
			Node &n2 = md.nodes[e.n2];
			n2.has_mp = true;
			// solid phase
			N_m_s = pcl.N2 * pcl.m_s;
			n2.m_s += N_m_s;
			n2.vx_s += N_m_s * pcl.vx_s;
			n2.vy_s += N_m_s * pcl.vy_s;

			// node 3
			Node &n3 = md.nodes[e.n3];
			n3.has_mp = true;
			// solid phase
			N_m_s = pcl.N3 * pcl.m_s;
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
			// body force on particle
			bf_s = (pcl.density_s - pcl.density_f) * (1.0 - pcl.n) * pcl.vol* bf.bf;
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

	// update nodal acceleration of fluid pahse
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.has_mp)
		{
			// fx_s
			n.ax_s = (n.fx_ext_s - n.fx_int_s) / n.m_s;
			// fy_s
			n.ay_s = (n.fy_ext_s - n.fy_int_s) / n.m_s;
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
			n.vx_s += n.ax_s * self.dtime;
			n.vy_s /= n.m_s;
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
			f_ub += n.m_s * n.m_s * (n.ax_s * n.ax_s + n.ay_s * n.ay_s);
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
			e_kin += n.m_s * (n.vx_s * n.vx_s + n.vy_s * n.vy_s);
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
	double ds11, ds22, ds12;
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
	
	double de12;
	double de_vol_s, de_vol_f;
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

			// strain enhancement appraoch
			de_vol_s = n1.de_vol_s * pcl.N1 + n2.de_vol_s * pcl.N2 + n3.de_vol_s * pcl.N3;

			// strain
			de_vol_s /= 3.0;
			de11 = e.dde11 + de_vol_s;
			de22 = e.dde22 + de_vol_s;
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
		}
	}
	
	return 0;
}
