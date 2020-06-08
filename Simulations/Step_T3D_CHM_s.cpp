#include "Simulations_pcp.h"

#include <cmath>
#include "MaterialModel.h"

#include "Step_T3D_CHM_s.h"

Step_T3D_CHM_s::Step_T3D_CHM_s(const char *_name) :
	Step(_name, "Step_T3D_CHM_s", &solve_substep_T3D_CHM_s),
	model(nullptr), damping_ratio(0.0) {}

Step_T3D_CHM_s::~Step_T3D_CHM_s() {}

int Step_T3D_CHM_s::init_calculation()
{
	Model_T3D_CHM_s &md = *model;

	if (is_first_step) {}

	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle &pcl = md.pcls[pcl_id];
		pcl.pe = md.elems;
		pcl.vol_s = pcl.m_s / pcl.density_s;
		pcl.x_ori = pcl.x;
		pcl.y_ori = pcl.y;
		pcl.z_ori = pcl.z;
		pcl.ux_s = 0.0;
		pcl.uy_s = 0.0;
		pcl.uz_s = 0.0;
		pcl.ux_f = 0.0;
		pcl.uy_f = 0.0;
		pcl.uz_f = 0.0;
	}

	return 0;
}

int Step_T3D_CHM_s::finalize_calculation() { return 0; }

namespace
{

inline double get_sign(double var)
{
	if (var > 0.0)
		return 1.0;
	else if (var < 0.0)
		return -1.0;
	return 0.0;
}

}

int solve_substep_T3D_CHM_s(void *_self)
{
	typedef Model_T3D_CHM_s::Node Node;
	typedef Model_T3D_CHM_s::Element Element;
	typedef Model_T3D_CHM_s::Particle Particle;
	
	Step_T3D_CHM_s &self = *(Step_T3D_CHM_s *)(_self);
	Model_T3D_CHM_s &md = *self.model;

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
		n.vz_s = 0.0;
		n.fx_ext_s = 0.0;
		n.fy_ext_s = 0.0;
		n.fz_ext_s = 0.0;
		n.fx_int_s = 0.0;
		n.fy_int_s = 0.0;
		n.fz_int_s = 0.0;
		// fluid phase
		n.m_f = 0.0;
		n.vx_f = 0.0;
		n.vy_f = 0.0;
		n.vz_f = 0.0;
		n.fx_ext_f = 0.0;
		n.fy_ext_f = 0.0;
		n.fz_ext_f = 0.0;
		n.fx_int_f = 0.0;
		n.fy_int_f = 0.0;
		n.fz_int_f = 0.0;
		// solid - fluid interaction
		n.fx_drag = 0.0;
		n.fy_drag = 0.0;
		n.fz_drag = 0.0;
		// strain enhancement approach
		n.pcl_vol = 0.0;
		n.de_vol_s = 0.0;
		n.de_vol_f = 0.0;
	}

	// init elements
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
		e.pcls = nullptr;
		e.pcl_vol = 0.0;
		e.n = 0.0;
		e.s11 = 0.0;
		e.s22 = 0.0;
		e.s33 = 0.0;
		e.s12 = 0.0;
		e.s23 = 0.0;
		e.s31 = 0.0;
		e.p = 0.0;
	}

	// init particles
	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle &pcl = md.pcls[pcl_id];
		if (pcl.pe)
		{
			if (!(pcl.pe = md.find_in_which_element(pcl)))
				continue;
			pcl.pe->add_pcl(pcl);

			pcl.vol = pcl.vol_s / (1.0 - pcl.n);
			pcl.m_f = pcl.n * pcl.density_f * pcl.vol;
			double n2_miu_div_k = pcl.n * pcl.n * md.miu / md.k;
			double n2_miu_div_k_vrx_vol = n2_miu_div_k * (pcl.vx_f - pcl.vx_s) * pcl.vol;
			double n2_miu_div_k_vry_vol = n2_miu_div_k * (pcl.vy_f - pcl.vy_s) * pcl.vol;
			double n2_miu_div_k_vrz_vol = n2_miu_div_k * (pcl.vz_f - pcl.vz_s) * pcl.vol;

			Element &e = *pcl.pe;
			e.pcl_vol += pcl.vol;
			e.n += pcl.vol_s;
			e.s11 += pcl.vol * pcl.s11;
			e.s22 += pcl.vol * pcl.s22;
			e.s33 += pcl.vol * pcl.s33;
			e.s12 += pcl.vol * pcl.s12;
			e.s23 += pcl.vol * pcl.s23;
			e.s31 += pcl.vol * pcl.s31;
			e.p += pcl.vol * pcl.p;

			double mvx_s = pcl.m_s * pcl.vx_s;
			double mvy_s = pcl.m_s * pcl.vy_s;
			double mvz_s = pcl.m_s * pcl.vz_s;
			double mvx_f = pcl.m_f * pcl.vx_f;
			double mvy_f = pcl.m_f * pcl.vy_f;
			double mvz_f = pcl.m_f * pcl.vz_f;
			
			// node 1
			Node &n1 = md.nodes[e.n1];
			n1.has_mp = true;
			// solid phase
			n1.m_s += pcl.N1 * pcl.m_s;
			n1.vx_s += pcl.N1 * mvx_s;
			n1.vy_s += pcl.N1 * mvy_s;
			n1.vz_s += pcl.N1 * mvz_s;
			// fluid phase
			n1.m_f += pcl.N1 * pcl.m_f;
			n1.vx_f += pcl.N1 * mvx_f;
			n1.vy_f += pcl.N1 * mvy_f;
			n1.vz_f += pcl.N1 * mvz_f;
			// solid - fluid interaction
			n1.fx_drag += pcl.N1 * n2_miu_div_k_vrx_vol;
			n1.fy_drag += pcl.N1 * n2_miu_div_k_vry_vol;
			n1.fz_drag += pcl.N1 * n2_miu_div_k_vrz_vol;

			// node 2
			Node &n2 = md.nodes[e.n2];
			n2.has_mp = true;
			// mixture phase
			n2.m_s += pcl.N2 * pcl.m_s;
			n2.vx_s += pcl.N2 * mvx_s;
			n2.vy_s += pcl.N2 * mvy_s;
			n2.vz_s += pcl.N2 * mvz_s;
			// fluid phase
			n2.m_f += pcl.N2 * pcl.m_f;
			n2.vx_f += pcl.N2 * mvx_f;
			n2.vy_f += pcl.N2 * mvy_f;
			n2.vz_f += pcl.N2 * mvz_f;
			// solid - fluid interaction
			n2.fx_drag += pcl.N2 * n2_miu_div_k_vrx_vol;
			n2.fy_drag += pcl.N2 * n2_miu_div_k_vry_vol;
			n2.fz_drag += pcl.N2 * n2_miu_div_k_vrz_vol;

			// node 3
			Node &n3 = md.nodes[e.n3];
			n3.has_mp = true;
			// mixture phase
			n3.m_s += pcl.N3 * pcl.m_s;
			n3.vx_s += pcl.N3 * mvx_s;
			n3.vy_s += pcl.N3 * mvy_s;
			n3.vz_s += pcl.N3 * mvz_s;
			// fluid phase
			n3.m_f += pcl.N3 * pcl.m_f;
			n3.vx_f += pcl.N3 * mvx_f;
			n3.vy_f += pcl.N3 * mvy_f;
			n3.vz_f += pcl.N3 * mvz_f;
			// solid - fluid interaction
			n3.fx_drag += pcl.N3 * n2_miu_div_k_vrx_vol;
			n3.fy_drag += pcl.N3 * n2_miu_div_k_vry_vol;
			n3.fz_drag += pcl.N3 * n2_miu_div_k_vrz_vol;

			// node 3
			Node &n4 = md.nodes[e.n4];
			n4.has_mp = true;
			// mixture phase
			n4.m_s += pcl.N4 * pcl.m_s;
			n4.vx_s += pcl.N4 * mvx_s;
			n4.vy_s += pcl.N4 * mvy_s;
			n4.vz_s += pcl.N4 * mvz_s;
			// fluid phase
			n4.m_f += pcl.N4 * pcl.m_f;
			n4.vx_f += pcl.N4 * mvx_f;
			n4.vy_f += pcl.N4 * mvy_f;
			n4.vz_f += pcl.N4 * mvz_f;
			// solid - fluid interaction
			n4.fx_drag += pcl.N4 * n2_miu_div_k_vrx_vol;
			n4.fy_drag += pcl.N4 * n2_miu_div_k_vry_vol;
			n4.fz_drag += pcl.N4 * n2_miu_div_k_vrz_vol;
		}
	}

	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
		if (e.pcls)
		{
			e.n = 1.0 - e.n / e.pcl_vol; // 1.0 - Vs / V
			e.s11 /= e.pcl_vol;
			e.s22 /= e.pcl_vol;
			e.s33 /= e.pcl_vol;
			e.s12 /= e.pcl_vol;
			e.s23 /= e.pcl_vol;
			e.s31 /= e.pcl_vol;
			e.p /= e.pcl_vol;
			if (e.pcl_vol > e.vol)
				e.pcl_vol = e.vol;

			Node &n1 = md.nodes[e.n1];
			Node &n2 = md.nodes[e.n2];
			Node &n3 = md.nodes[e.n3];
			Node &n4 = md.nodes[e.n4];

			// node 1
			n1.fx_int_s += (e.dN1_dx * (e.s11 - (1.0 - e.n) * e.p) + e.dN1_dy * e.s12 + e.dN1_dz * e.s31) * e.pcl_vol;
			n1.fy_int_s += (e.dN1_dx * e.s12 + e.dN1_dy * (e.s22 - (1.0 - e.n) * e.p) + e.dN1_dz * e.s23) * e.pcl_vol;
			n1.fz_int_s += (e.dN1_dx * e.s31 + e.dN1_dy * e.s23 + e.dN1_dz * (e.s33 - (1.0 - e.n) * e.p)) * e.pcl_vol;
			n1.fx_int_f += (e.dN1_dx * e.n * -e.p) * e.pcl_vol;
			n1.fy_int_f += (e.dN1_dy * e.n * -e.p) * e.pcl_vol;
			n1.fz_int_f += (e.dN1_dz * e.n * -e.p) * e.pcl_vol;
			
			// node 2
			n2.fx_int_s += (e.dN2_dx * (e.s11 - (1.0 - e.n) * e.p) + e.dN2_dy * e.s12 + e.dN2_dz * e.s31) * e.pcl_vol;
			n2.fy_int_s += (e.dN2_dx * e.s12 + e.dN2_dy * (e.s22 - (1.0 - e.n) * e.p) + e.dN2_dz * e.s23) * e.pcl_vol;
			n2.fz_int_s += (e.dN2_dx * e.s31 + e.dN2_dy * e.s23 + e.dN2_dz * (e.s33 - (1.0 - e.n) * e.p)) * e.pcl_vol;
			n2.fx_int_f += (e.dN2_dx * e.n * -e.p) * e.pcl_vol;
			n2.fy_int_f += (e.dN2_dy * e.n * -e.p) * e.pcl_vol;
			n2.fz_int_f += (e.dN2_dz * e.n * -e.p) * e.pcl_vol;

			// node 3
			n3.fx_int_s += (e.dN3_dx * (e.s11 - (1.0 - e.n) * e.p) + e.dN3_dy * e.s12 + e.dN3_dz * e.s31) * e.pcl_vol;
			n3.fy_int_s += (e.dN3_dx * e.s12 + e.dN3_dy * (e.s22 - (1.0 - e.n) * e.p) + e.dN3_dz * e.s23) * e.pcl_vol;
			n3.fz_int_s += (e.dN3_dx * e.s31 + e.dN3_dy * e.s23 + e.dN3_dz * (e.s33 - (1.0 - e.n) * e.p)) * e.pcl_vol;
			n3.fx_int_f += (e.dN3_dx * e.n * -e.p) * e.pcl_vol;
			n3.fy_int_f += (e.dN3_dy * e.n * -e.p) * e.pcl_vol;
			n3.fz_int_f += (e.dN3_dz * e.n * -e.p) * e.pcl_vol;

			// node 4
			n4.fx_int_s += (e.dN4_dx * (e.s11 - (1.0 - e.n) * e.p) + e.dN4_dy * e.s12 + e.dN4_dz * e.s31) * e.pcl_vol;
			n4.fy_int_s += (e.dN4_dx * e.s12 + e.dN4_dy * (e.s22 - (1.0 - e.n) * e.p) + e.dN4_dz * e.s23) * e.pcl_vol;
			n4.fz_int_s += (e.dN4_dx * e.s31 + e.dN4_dy * e.s23 + e.dN4_dz * (e.s33 - (1.0 - e.n) * e.p)) * e.pcl_vol;
			n4.fx_int_f += (e.dN4_dx * e.n * -e.p) * e.pcl_vol;
			n4.fy_int_f += (e.dN4_dy * e.n * -e.p) * e.pcl_vol;
			n4.fz_int_f += (e.dN4_dz * e.n * -e.p) * e.pcl_vol;
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
			// body force on particle
			//bf_s = pcl.m_s * bf.bf;
			//bf_f = pcl.m_f * bf.bf;
			bf_s = (pcl.density_s - pcl.density_f) * (1.0 - pcl.n) * pcl.vol * bf.bf;
			bf_f = 0.0;
			// node1
			Node &n1 = md.nodes[e.n1];
			n1.fx_ext_s += pcl.N1 * bf_s;
			n1.fx_ext_f += pcl.N1 * bf_f;
			// node2
			Node &n2 = md.nodes[e.n2];
			n2.fx_ext_s += pcl.N2 * bf_s;
			n2.fx_ext_f += pcl.N2 * bf_f;
			// node3
			Node &n3 = md.nodes[e.n3];
			n3.fx_ext_s += pcl.N3 * bf_s;
			n3.fx_ext_f += pcl.N3 * bf_f;
			// node4
			Node &n4 = md.nodes[e.n4];
			n4.fx_ext_s += pcl.N4 * bf_s;
			n4.fx_ext_f += pcl.N4 * bf_f;
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
			//bf_s = pcl.m_s * bf.bf;
			//bf_f = pcl.m_f * bf.bf;
			bf_s = (pcl.density_s - pcl.density_f) * (1.0 - pcl.n) * pcl.vol * bf.bf;
			bf_f = 0.0;
			// node 1
			Node &n1 = md.nodes[e.n1];
			n1.fy_ext_s += pcl.N1 * bf_s;
			n1.fy_ext_f += pcl.N1 * bf_f;
			// node 2
			Node &n2 = md.nodes[e.n2];
			n2.fy_ext_s += pcl.N2 * bf_s;
			n2.fy_ext_f += pcl.N2 * bf_f;
			// node 3
			Node &n3 = md.nodes[e.n3];
			n3.fy_ext_s += pcl.N3 * bf_s;
			n3.fy_ext_f += pcl.N3 * bf_f;
			// node4
			Node &n4 = md.nodes[e.n4];
			n4.fy_ext_s += pcl.N4 * bf_s;
			n4.fy_ext_f += pcl.N4 * bf_f;
		}
	}
	for (size_t bf_id = 0; bf_id < md.bfz_num; ++bf_id)
	{
		BodyForceAtPcl &bf = md.bfzs[bf_id];
		Particle &pcl = md.pcls[bf.pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			// body force on particle
			//bf_s = pcl.m_s * bf.bf;
			//bf_f = pcl.m_f * bf.bf;
			bf_s = (pcl.density_s - pcl.density_f) * (1.0 - pcl.n) * pcl.vol * bf.bf;
			bf_f = 0.0;
			// node 1
			Node &n1 = md.nodes[e.n1];
			n1.fz_ext_s += pcl.N1 * bf_s;
			n1.fz_ext_f += pcl.N1 * bf_f;
			// node 2
			Node &n2 = md.nodes[e.n2];
			n2.fz_ext_s += pcl.N2 * bf_s;
			n2.fz_ext_f += pcl.N2 * bf_f;
			// node 3
			Node &n3 = md.nodes[e.n3];
			n3.fz_ext_s += pcl.N3 * bf_s;
			n3.fz_ext_f += pcl.N3 * bf_f;
			// node4
			Node &n4 = md.nodes[e.n4];
			n4.fz_ext_s += pcl.N4 * bf_s;
			n4.fz_ext_f += pcl.N4 * bf_f;
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
			// node 4
			Node &n4 = md.nodes[e.n4];
			n4.fx_ext_s += pcl.N4 * tf.t;
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
			// node 4
			Node &n4 = md.nodes[e.n4];
			n4.fy_ext_s += pcl.N4 * tf.t;
		}
	}
	for (size_t tf_id = 0; tf_id < md.tz_num; ++tf_id)
	{
		TractionBCAtPcl &tf = md.tzs[tf_id];
		Particle &pcl = md.pcls[tf.pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			// node 1
			Node &n1 = md.nodes[e.n1];
			n1.fz_ext_s += pcl.N1 * tf.t;
			// node 2
			Node &n2 = md.nodes[e.n2];
			n2.fz_ext_s += pcl.N2 * tf.t;
			// node 3
			Node &n3 = md.nodes[e.n3];
			n3.fz_ext_s += pcl.N3 * tf.t;
			// node 4
			Node &n4 = md.nodes[e.n4];
			n4.fz_ext_s += pcl.N4 * tf.t;
		}
	}
	
	// update nodal acceleration of fluid pahse
	double nf;
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.has_mp)
		{
			// fx_s
			nf = n.fx_ext_s - n.fx_int_s;
			n.ax_s = (nf + n.fx_drag - self.damping_ratio * abs(nf) * get_sign(n.vx_s)) / n.m_s;
			//n.ax_s = (n.fx_ext_s - n.fx_int_s) / n.m_s;
			// fy_s
			nf = n.fy_ext_s - n.fy_int_s;
			n.ay_s = (nf + n.fy_drag - self.damping_ratio * abs(nf) * get_sign(n.vy_s)) / n.m_s;
			//n.ay_s = (n.fy_ext_s - n.fy_int_s) / n.m_s;
			// fz_s
			nf = n.fz_ext_s - n.fz_int_s;
			n.az_s = (nf + n.fz_drag - self.damping_ratio * abs(nf) * get_sign(n.vz_s)) / n.m_s;
			//n.az_s = (n.fz_ext_s - n.fz_int_s) / n.m_s;
			// fx_f
			nf = n.fx_ext_f - n.fx_int_f;
			n.ax_f = (nf - n.fx_drag - self.damping_ratio * abs(nf) * get_sign(n.vx_f)) / n.m_f;
			//n.ax_f = (n.fx_ext_f - n.fx_int_f) / n.m_f;
			// fy_f
			nf = n.fy_ext_f - n.fy_int_f;
			n.ay_f = (nf - n.fy_drag - self.damping_ratio * abs(nf) * get_sign(n.vy_f)) / n.m_f;
			//n.ay_f = (n.fy_ext_f - n.fy_int_f) / n.m_f;
			// fz_f
			nf = n.fz_ext_f - n.fz_int_f;
			n.az_f = (nf - n.fz_drag - self.damping_ratio * abs(nf) * get_sign(n.vz_f)) / n.m_f;
			//n.az_f = (n.fz_ext_f - n.fz_int_f) / n.m_f;
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
	for (size_t a_id = 0; a_id < md.asz_num; ++a_id)
	{
		Node &n = md.nodes[md.aszs[a_id].node_id];
		n.az_s = md.aszs[a_id].a;
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
	for (size_t a_id = 0; a_id < md.afz_num; ++a_id)
	{
		Node &n = md.nodes[md.afzs[a_id].node_id];
		n.az_f = md.afzs[a_id].a;
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
			n.vz_s /= n.m_s;
			n.vz_s += n.az_s * self.dtime;
			n.vx_f /= n.m_f;
			n.vx_f += n.ax_f * self.dtime;
			n.vy_f /= n.m_f;
			n.vy_f += n.ay_f * self.dtime;
			n.vz_f /= n.m_f;
			n.vz_f += n.az_f * self.dtime;
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
	for (size_t v_id = 0; v_id < md.vsz_num; ++v_id)
	{
		Node &n = md.nodes[md.vszs[v_id].node_id];
		n.vz_s = md.vszs[v_id].v;
		n.az_s = 0.0;
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
	for (size_t v_id = 0; v_id < md.vfz_num; ++v_id)
	{
		Node &n = md.nodes[md.vfzs[v_id].node_id];
		n.vz_f = md.vfzs[v_id].v;
		n.az_f = 0.0;
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
			n.duz_s = n.vz_s * self.dtime;
			// fluid phase
			n.dux_f = n.vx_f * self.dtime;
			n.duy_f = n.vy_f * self.dtime;
			n.duz_f = n.vz_f * self.dtime;
		}
	}

	// map variables back to particles and update their variables
	double de11, de22, de33, de_mean;
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
		if (e.pcls)
		{
			Node &n1 = md.nodes[e.n1];
			Node &n2 = md.nodes[e.n2];
			Node &n3 = md.nodes[e.n3];
			Node &n4 = md.nodes[e.n4];
			// strain increment
			de11 = n1.dux_s * e.dN1_dx + n2.dux_s * e.dN2_dx + n3.dux_s * e.dN3_dx + n4.dux_s * e.dN4_dx;
			de22 = n1.duy_s * e.dN1_dy + n2.duy_s * e.dN2_dy + n3.duy_s * e.dN3_dy + n4.duy_s * e.dN4_dy;
			de33 = n1.duz_s * e.dN1_dz + n2.duz_s * e.dN2_dz + n3.duz_s * e.dN3_dz + n4.duz_s * e.dN4_dz;
			e.de12 = (n1.dux_s * e.dN1_dy + n2.dux_s * e.dN2_dy + n3.dux_s * e.dN3_dy + n4.dux_s * e.dN4_dy
					+ n1.duy_s * e.dN1_dx + n2.duy_s * e.dN2_dx + n3.duy_s * e.dN3_dx + n4.duy_s * e.dN4_dx) * 0.5;
			e.de23 = (n1.duy_s * e.dN1_dz + n2.duy_s * e.dN2_dz + n3.duy_s * e.dN3_dz + n4.duy_s * e.dN4_dz
					+ n1.duz_s * e.dN1_dy + n2.duz_s * e.dN2_dy + n3.duz_s * e.dN3_dy + n4.duz_s * e.dN4_dy) * 0.5;
			e.de31 = (n1.duz_s * e.dN1_dx + n2.duz_s * e.dN2_dx + n3.duz_s * e.dN3_dx + n4.duz_s * e.dN4_dx
					+ n1.dux_s * e.dN1_dz + n2.dux_s * e.dN2_dz + n3.dux_s * e.dN3_dz + n4.dux_s * e.dN4_dz) * 0.5;
			// volumetric strain of solid phase
			e.de_vol_s = de11 + de22 + de33;
			de_mean = e.de_vol_s / 3.0;
			e.dde11 = de11 - de_mean;
			e.dde22 = de22 - de_mean;
			e.dde33 = de33 - de_mean;
			// "volumetric strain" of fluid phase, take compression as positive
			e.de_vol_f = -(1.0 - e.n) / e.n * e.de_vol_s
						-(n1.dux_f * e.dN1_dx + n2.dux_f * e.dN2_dx + n3.dux_f * e.dN3_dx + n4.dux_f * e.dN4_dx
						+ n1.duy_f * e.dN1_dy + n2.duy_f * e.dN2_dy + n3.duy_f * e.dN3_dy + n4.duy_f * e.dN4_dy
						+ n1.duz_f * e.dN1_dz + n2.duz_f * e.dN2_dz + n3.duz_f * e.dN3_dz + n4.duz_f * e.dN4_dz);
		}
	}

	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle &pcl = md.pcls[pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			double vol_de_vol_s = pcl.vol * e.de_vol_s;
			double vol_de_vol_f = pcl.vol * e.de_vol_f;
			// node 1
			Node &n1 = md.nodes[e.n1];
			n1.pcl_vol  += pcl.N1 * pcl.vol;
			n1.de_vol_s += pcl.N1 * vol_de_vol_s;
			n1.de_vol_f += pcl.N1 * vol_de_vol_f;
			// node 2
			Node &n2 = md.nodes[e.n2];
			n2.pcl_vol  += pcl.N2 * pcl.vol;
			n2.de_vol_s += pcl.N2 * vol_de_vol_s;
			n2.de_vol_f += pcl.N2 * vol_de_vol_f;
			// node 3
			Node &n3 = md.nodes[e.n3];
			n3.pcl_vol  += pcl.N3 * pcl.vol;
			n3.de_vol_s += pcl.N3 * vol_de_vol_s;
			n3.de_vol_f += pcl.N3 * vol_de_vol_f;
			// node 4
			Node &n4 = md.nodes[e.n4];
			n4.pcl_vol  += pcl.N4 * pcl.vol;
			n4.de_vol_s += pcl.N4 * vol_de_vol_s;
			n4.de_vol_f += pcl.N4 * vol_de_vol_f;
		}
	}

	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.has_mp)
		{
			n.de_vol_s /= n.pcl_vol;
			n.de_vol_f /= n.pcl_vol;
		}
	}

	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
		if (e.pcls)
		{
			Node &n1 = md.nodes[e.n1];
			Node &n2 = md.nodes[e.n2];
			Node &n3 = md.nodes[e.n3];
			Node &n4 = md.nodes[e.n4];
			// solid volumetric strain
			e.de_vol_s = (n1.de_vol_s + n2.de_vol_s + n3.de_vol_s + n4.de_vol_s) * 0.25;
			// fluid volumetric strain
			e.de_vol_f = (n1.de_vol_f + n2.de_vol_f + n3.de_vol_f + n4.de_vol_f) * 0.25;
			e.p += md.Kf * e.de_vol_f;
		}
	}

	double de12, de23, de31;
	double de_vol_s; // de_vol_f;
	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle &pcl = md.pcls[pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			Node &n1 = md.nodes[e.n1];
			Node &n2 = md.nodes[e.n2];
			Node &n3 = md.nodes[e.n3];
			Node &n4 = md.nodes[e.n4];

			// velocity
			pcl.vx_s += (n1.ax_s * pcl.N1 + n2.ax_s * pcl.N2 + n3.ax_s * pcl.N3 + n4.ax_s * pcl.N4) * self.dtime;
			pcl.vy_s += (n1.ay_s * pcl.N1 + n2.ay_s * pcl.N2 + n3.ay_s * pcl.N3 + n4.ay_s * pcl.N4) * self.dtime;
			pcl.vz_s += (n1.az_s * pcl.N1 + n2.az_s * pcl.N2 + n3.az_s * pcl.N3 + n4.az_s * pcl.N4) * self.dtime;
			pcl.vx_f += (n1.ax_f * pcl.N1 + n2.ax_f * pcl.N2 + n3.ax_f * pcl.N3 + n4.ax_f * pcl.N4) * self.dtime;
			pcl.vy_f += (n1.ay_f * pcl.N1 + n2.ay_f * pcl.N2 + n3.ay_f * pcl.N3 + n4.ay_f * pcl.N4) * self.dtime;
			pcl.vz_f += (n1.az_f * pcl.N1 + n2.az_f * pcl.N2 + n3.az_f * pcl.N3 + n4.az_f * pcl.N4) * self.dtime;

			// displacement
			pcl.ux_s += n1.dux_s * pcl.N1 + n2.dux_s * pcl.N2 + n3.dux_s * pcl.N3 + n4.dux_s * pcl.N4;
			pcl.uy_s += n1.duy_s * pcl.N1 + n2.duy_s * pcl.N2 + n3.duy_s * pcl.N3 + n4.duy_s * pcl.N4;
			pcl.uz_s += n1.duz_s * pcl.N1 + n2.duz_s * pcl.N2 + n3.duz_s * pcl.N3 + n4.duz_s * pcl.N4;
			pcl.ux_f += n1.dux_f * pcl.N1 + n2.dux_f * pcl.N2 + n3.dux_f * pcl.N3 + n4.dux_f * pcl.N4;
			pcl.uy_f += n1.duy_f * pcl.N1 + n2.duy_f * pcl.N2 + n3.duy_f * pcl.N3 + n4.duy_f * pcl.N4;
			pcl.uz_f += n1.duz_f * pcl.N1 + n2.duz_f * pcl.N2 + n3.duz_f * pcl.N3 + n4.duz_f * pcl.N4;

			// update position
			pcl.x = pcl.x_ori + pcl.ux_s;
			pcl.y = pcl.y_ori + pcl.uy_s;
			pcl.z = pcl.z_ori + pcl.uz_s;

			// strain enhancement appraoch
			//de_vol_s = n1.de_vol_s * pcl.N1 + n2.de_vol_s * pcl.N2 + n3.de_vol_s * pcl.N3 + n4.de_vol_s * pcl.N4;
			de_vol_s = e.de_vol_s;

			// strain increment
			de_vol_s /= 3.0;
			de11 = e.dde11 + de_vol_s;
			de22 = e.dde22 + de_vol_s;
			de33 = e.dde33 + de_vol_s;
			de12 = e.de12;
			de23 = e.de23;
			de31 = e.de31;
			// strain
			pcl.e11 += de11;
			pcl.e22 += de22;
			pcl.e33 += de33;
			pcl.e12 += de12;
			pcl.e23 += de23;
			pcl.e31 += de31;

			// update stress using constitutive model
			double dstrain[6] = { de11, de22, de33, de12, de23, de31 };
			pcl.mm->integrate(dstrain);
			const double *dstress = pcl.mm->get_dstress();
			pcl.s11 += dstress[0];
			pcl.s22 += dstress[1];
			pcl.s33 += dstress[2];
			pcl.s12 += dstress[3];
			pcl.s23 += dstress[4];
			pcl.s31 += dstress[5];

			// porosity
			pcl.n = (de_vol_s + pcl.n) / (1.0 + de_vol_s);
			
			// pore pressure
			//de_vol_f = n1.de_vol_f * pcl.N1 + n2.de_vol_f * pcl.N2 + n3.de_vol_f * pcl.N3;
			pcl.p = e.p;

			// fluid density
			pcl.density_f /= 1.0 - e.de_vol_f;
		}
	}
	
	return 0;
}
