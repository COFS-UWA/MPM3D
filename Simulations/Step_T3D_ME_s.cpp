#include "Simulations_pcp.h"

#include <iostream>
#include <assert.h>

#include <cmath>
#include "MaterialModel.h"

#include "Step_T3D_ME_s.h"

Step_T3D_ME_s::Step_T3D_ME_s(const char *_name) :
	Step(_name, "Step_T3D_ME_s", &solve_substep_T3D_ME_s),
	model(nullptr), damping_ratio(0.0) {}

Step_T3D_ME_s::~Step_T3D_ME_s() {}

int Step_T3D_ME_s::init_calculation()
{
	Model_T3D_ME_s &md = *model;

	if (is_first_step) {}

	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle &pcl = md.pcls[pcl_id];
		pcl.pe = md.elems;
		pcl.vol = pcl.m / pcl.density;
		pcl.x_ori = pcl.x;
		pcl.y_ori = pcl.y;
		pcl.z_ori = pcl.z;
		pcl.ux = 0.0;
		pcl.uy = 0.0;
		pcl.uz = 0.0;
	}

	return 0;
}

int Step_T3D_ME_s::finalize_calculation() { return 0; }

int solve_substep_T3D_ME_s(void *_self)
{
	typedef Model_T3D_ME_s::Particle Particle;
	typedef Model_T3D_ME_s::Element Element;
	typedef Model_T3D_ME_s::Node Node;

	Step_T3D_ME_s &self = *(Step_T3D_ME_s *)(_self);
	Model_T3D_ME_s &md = *self.model;

	// init nodes
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		n.has_mp = false;
		n.m = 0.0;
		n.vx = 0.0;
		n.vy = 0.0;
		n.vz = 0.0;
		n.fx_ext = 0.0;
		n.fy_ext = 0.0;
		n.fz_ext = 0.0;
		n.fx_int = 0.0;
		n.fy_int = 0.0;
		n.fz_int = 0.0;
		// strain enhancement
		n.pcl_vol = 0.0;
		n.de_vol = 0.0;
	}

	// init elements
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
		e.pcls = nullptr;
		e.pcl_vol = 0.0;
		e.s11 = 0.0;
		e.s22 = 0.0;
		e.s33 = 0.0;
		e.s12 = 0.0;
		e.s23 = 0.0;
		e.s31 = 0.0;
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

			// test find_in_which_element function
			//size_t id1 = md.find_in_which_element(pcl)->id;
			//size_t id2 = md.find_in_which_element_bf(pcl)->id;
			//assert(id1 == id2);
			//std::cout << id1 << ", " << id2 << "\n";

			Element &e = *pcl.pe;
			pcl.vol = pcl.m / pcl.density;
			e.pcl_vol += pcl.vol;
			e.s11 += pcl.vol * pcl.s11;
			e.s22 += pcl.vol * pcl.s22;
			e.s33 += pcl.vol * pcl.s33;
			e.s12 += pcl.vol * pcl.s12;
			e.s23 += pcl.vol * pcl.s23;
			e.s31 += pcl.vol * pcl.s31;

			double mvx = pcl.m * pcl.vx;
			double mvy = pcl.m * pcl.vy;
			double mvz = pcl.m * pcl.vz;
			// node 1
			Node &n1 = md.nodes[e.n1];
			n1.has_mp = true;
			n1.m += pcl.N1 * pcl.m;
			n1.vx += pcl.N1 * mvx;
			n1.vy += pcl.N1 * mvy;
			n1.vz += pcl.N1 * mvz;
			// node 2
			Node &n2 = md.nodes[e.n2];
			n2.has_mp = true;
			n2.m += pcl.N2 * pcl.m;
			n2.vx += pcl.N2 * mvx;
			n2.vy += pcl.N2 * mvy;
			n2.vz += pcl.N2 * mvz;
			// node 3
			Node &n3 = md.nodes[e.n3];
			n3.has_mp = true;
			n3.m += pcl.N3 * pcl.m;
			n3.vx += pcl.N3 * mvx;
			n3.vy += pcl.N3 * mvy;
			n3.vz += pcl.N3 * mvz;
			// node 4
			Node &n4 = md.nodes[e.n4];
			n4.has_mp = true;
			n4.m += pcl.N4 * pcl.m;
			n4.vx += pcl.N4 * mvx;
			n4.vy += pcl.N4 * mvy;
			n4.vz += pcl.N4 * mvz;
		}
	}

	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
		if (e.pcls)
		{
			e.s11 /= e.pcl_vol;
			e.s22 /= e.pcl_vol;
			e.s33 /= e.pcl_vol;
			e.s12 /= e.pcl_vol;
			e.s23 /= e.pcl_vol;
			e.s31 /= e.pcl_vol;
			if (e.pcl_vol > e.vol)
				e.pcl_vol = e.vol;

			Node &n1 = md.nodes[e.n1];
			Node &n2 = md.nodes[e.n2];
			Node &n3 = md.nodes[e.n3];
			Node &n4 = md.nodes[e.n4];
			// node 1
			n1.fx_int += (e.dN1_dx * e.s11 + e.dN1_dy * e.s12 + e.dN1_dz * e.s31) * e.pcl_vol;
			n1.fy_int += (e.dN1_dx * e.s12 + e.dN1_dy * e.s22 + e.dN1_dz * e.s23) * e.pcl_vol;
			n1.fz_int += (e.dN1_dx * e.s31 + e.dN1_dy * e.s23 + e.dN1_dz * e.s33) * e.pcl_vol;
			// node 2
			n2.fx_int += (e.dN2_dx * e.s11 + e.dN2_dy * e.s12 + e.dN2_dz * e.s31) * e.pcl_vol;
			n2.fy_int += (e.dN2_dx * e.s12 + e.dN2_dy * e.s22 + e.dN2_dz * e.s23) * e.pcl_vol;
			n2.fz_int += (e.dN2_dx * e.s31 + e.dN2_dy * e.s23 + e.dN2_dz * e.s33) * e.pcl_vol;
			// node 3
			n3.fx_int += (e.dN3_dx * e.s11 + e.dN3_dy * e.s12 + e.dN3_dz * e.s31) * e.pcl_vol;
			n3.fy_int += (e.dN3_dx * e.s12 + e.dN3_dy * e.s22 + e.dN3_dz * e.s23) * e.pcl_vol;
			n3.fz_int += (e.dN3_dx * e.s31 + e.dN3_dy * e.s23 + e.dN3_dz * e.s33) * e.pcl_vol;
			// node 4
			n4.fx_int += (e.dN4_dx * e.s11 + e.dN4_dy * e.s12 + e.dN4_dz * e.s31) * e.pcl_vol;
			n4.fy_int += (e.dN4_dx * e.s12 + e.dN4_dy * e.s22 + e.dN4_dz * e.s23) * e.pcl_vol;
			n4.fz_int += (e.dN4_dx * e.s31 + e.dN4_dy * e.s23 + e.dN4_dz * e.s33) * e.pcl_vol;
		}
	}

	// body force
	double bf_mag;
	for (size_t bf_id = 0; bf_id < md.bfx_num; ++bf_id)
	{
		BodyForceAtPcl &bf = md.bfxs[bf_id];
		Particle &pcl = md.pcls[bf.pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			bf_mag = pcl.m * bf.bf;
			// node 1
			Node &n1 = md.nodes[e.n1];
			n1.fx_ext += pcl.N1 * bf_mag;
			// node 2
			Node &n2 = md.nodes[e.n2];
			n2.fx_ext += pcl.N2 * bf_mag;
			// node 3
			Node &n3 = md.nodes[e.n3];
			n3.fx_ext += pcl.N3 * bf_mag;
			// node 4
			Node &n4 = md.nodes[e.n4];
			n4.fx_ext += pcl.N4 * bf_mag;
		}
	}
	for (size_t bf_id = 0; bf_id < md.bfy_num; ++bf_id)
	{
		BodyForceAtPcl &bf = md.bfys[bf_id];
		Particle &pcl = md.pcls[bf.pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			bf_mag = pcl.m * bf.bf;
			// node 1
			Node &n1 = md.nodes[e.n1];
			n1.fy_ext += pcl.N1 * bf_mag;
			// node 2
			Node &n2 = md.nodes[e.n2];
			n2.fy_ext += pcl.N2 * bf_mag;
			// node 3
			Node &n3 = md.nodes[e.n3];
			n3.fy_ext += pcl.N3 * bf_mag;
			// node 4
			Node &n4 = md.nodes[e.n4];
			n4.fy_ext += pcl.N4 * bf_mag;
		}
	}
	for (size_t bf_id = 0; bf_id < md.bfz_num; ++bf_id)
	{
		BodyForceAtPcl &bf = md.bfzs[bf_id];
		Particle &pcl = md.pcls[bf.pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			bf_mag = pcl.m * bf.bf;
			// node 1
			Node &n1 = md.nodes[e.n1];
			n1.fz_ext += pcl.N1 * bf_mag;
			// node 2
			Node &n2 = md.nodes[e.n2];
			n2.fz_ext += pcl.N2 * bf_mag;
			// node 3
			Node &n3 = md.nodes[e.n3];
			n3.fz_ext += pcl.N3 * bf_mag;
			// node 4
			Node &n4 = md.nodes[e.n4];
			n4.fz_ext += pcl.N4 * bf_mag;
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
			n1.fx_ext += pcl.N1 * tf.t;
			// node 2
			Node &n2 = md.nodes[e.n2];
			n2.fx_ext += pcl.N2 * tf.t;
			// node 3
			Node &n3 = md.nodes[e.n3];
			n3.fx_ext += pcl.N3 * tf.t;
			// node 4
			Node &n4 = md.nodes[e.n4];
			n4.fx_ext += pcl.N4 * tf.t;
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
			n1.fy_ext += pcl.N1 * tf.t;
			// node 2
			Node &n2 = md.nodes[e.n2];
			n2.fy_ext += pcl.N2 * tf.t;
			// node 3
			Node &n3 = md.nodes[e.n3];
			n3.fy_ext += pcl.N3 * tf.t;
			// node 4
			Node &n4 = md.nodes[e.n4];
			n4.fy_ext += pcl.N4 * tf.t;
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
			n1.fz_ext += pcl.N1 * tf.t;
			// node 2
			Node &n2 = md.nodes[e.n2];
			n2.fz_ext += pcl.N2 * tf.t;
			// node 3
			Node &n3 = md.nodes[e.n3];
			n3.fz_ext += pcl.N3 * tf.t;
			// node 4
			Node &n4 = md.nodes[e.n4];
			n4.fz_ext += pcl.N4 * tf.t;
		}
	}

	// update nodal acceleration of fluid pahse
	double nf, v_sign;
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.has_mp) // or n.m_f != 0.0
		{
			// fx
			if (n.vx > 0.0)
				v_sign = 1.0;
			else if (n.vx < 0.0)
				v_sign = -1.0;
			else
				v_sign = 0.0;
			nf = n.fx_ext - n.fx_int;
			n.ax = (nf - self.damping_ratio * abs(nf) * v_sign) / n.m;
			// fy
			if (n.vy > 0.0)
				v_sign = 1.0;
			else if (n.vy < 0.0)
				v_sign = -1.0;
			else
				v_sign = 0.0;
			nf = n.fy_ext - n.fy_int;
			n.ay = (nf - self.damping_ratio * abs(nf) * v_sign) / n.m;
			// fz
			if (n.vz > 0.0)
				v_sign = 1.0;
			else if (n.vz < 0.0)
				v_sign = -1.0;
			else
				v_sign = 0.0;
			nf = n.fz_ext - n.fz_int;
			n.az = (nf - self.damping_ratio * abs(nf) * v_sign) / n.m;
		}
	}

	// apply acceleration bc
	for (size_t a_id = 0; a_id < md.ax_num; ++a_id)
	{
		Node &n = md.nodes[md.axs[a_id].node_id];
		n.ax = md.axs[a_id].a;
	}
	for (size_t a_id = 0; a_id < md.ay_num; ++a_id)
	{
		Node &n = md.nodes[md.ays[a_id].node_id];
		n.ay = md.ays[a_id].a;
	}
	for (size_t a_id = 0; a_id < md.az_num; ++a_id)
	{
		Node &n = md.nodes[md.azs[a_id].node_id];
		n.az = md.azs[a_id].a;
	}

	// update nodal momentum
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.has_mp)
		{
			n.vx /= n.m;
			n.vx += n.ax * self.dtime;
			n.vy /= n.m;
			n.vy += n.ay * self.dtime;
			n.vz /= n.m;
			n.vz += n.az * self.dtime;
		}
	}

	// apply velocity bc
	for (size_t v_id = 0; v_id < md.vx_num; ++v_id)
	{
		Node &n = md.nodes[md.vxs[v_id].node_id];
		n.vx = md.vxs[v_id].v;
		n.ax = 0.0;
	}
	for (size_t v_id = 0; v_id < md.vy_num; ++v_id)
	{
		Node &n = md.nodes[md.vys[v_id].node_id];
		n.vy = md.vys[v_id].v;
		n.ay = 0.0;
	}
	for (size_t v_id = 0; v_id < md.vz_num; ++v_id)
	{
		Node &n = md.nodes[md.vzs[v_id].node_id];
		n.vz = md.vzs[v_id].v;
		n.az = 0.0;
	}

	// update displacement increment of both phases
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.has_mp)
		{
			n.dux = n.vx * self.dtime;
			n.duy = n.vy * self.dtime;
			n.duz = n.vz * self.dtime;
		}
	}

	// map variables back to particles and update their variables
	double de11, de22, de33, de_vol_over_3;
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
			de11 = n1.dux * e.dN1_dx + n2.dux * e.dN2_dx + n3.dux * e.dN3_dx + n4.dux * e.dN4_dx;
			de22 = n1.duy * e.dN1_dy + n2.duy * e.dN2_dy + n3.duy * e.dN3_dy + n4.duy * e.dN4_dy;
			de33 = n1.duz * e.dN1_dz + n2.duz * e.dN2_dz + n3.duz * e.dN3_dz + n4.duz * e.dN4_dz;
			e.de12 = (n1.dux * e.dN1_dy + n2.dux * e.dN2_dy + n3.dux * e.dN3_dy + n4.dux * e.dN4_dy
					+ n1.duy * e.dN1_dx + n2.duy * e.dN2_dx + n3.duy * e.dN3_dx + n4.duy * e.dN4_dx) * 0.5;
			e.de23 = (n1.duy * e.dN1_dz + n2.duy * e.dN2_dz + n3.duy * e.dN3_dz + n4.duy * e.dN4_dz
					+ n1.duz * e.dN1_dy + n2.duz * e.dN2_dy + n3.duz * e.dN3_dy + n4.duz * e.dN4_dy) * 0.5;
			e.de31 = (n1.duz * e.dN1_dx + n2.duz * e.dN2_dx + n3.duz * e.dN3_dx + n4.duz * e.dN4_dx
					+ n1.dux * e.dN1_dz + n2.dux * e.dN2_dz + n3.dux * e.dN3_dz + n4.dux * e.dN4_dz) * 0.5;
			e.de_vol = de11 + de22 + de33;
			de_vol_over_3 = e.de_vol / 3.0;
			e.dde11 = de11 - de_vol_over_3;
			e.dde22 = de22 - de_vol_over_3;
			e.dde33 = de33 - de_vol_over_3;
		}
	}

	// strain enhancement
	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle &pcl = md.pcls[pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			double vol_de_vol = pcl.vol * e.de_vol;
			// node 1
			Node &n1 = md.nodes[e.n1];
			n1.pcl_vol += pcl.N1 * pcl.vol;
			n1.de_vol += pcl.N1 * vol_de_vol;
			// node 2
			Node &n2 = md.nodes[e.n2];
			n2.pcl_vol += pcl.N2 * pcl.vol;
			n2.de_vol += pcl.N2 * vol_de_vol;
			// node 3
			Node &n3 = md.nodes[e.n3];
			n3.pcl_vol += pcl.N3 * pcl.vol;
			n3.de_vol += pcl.N3 * vol_de_vol;
			// node 4
			Node &n4 = md.nodes[e.n4];
			n4.pcl_vol += pcl.N4 * pcl.vol;
			n4.de_vol += pcl.N4 * vol_de_vol;
		}
	}

	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.has_mp)
		{
			n.de_vol /= n.pcl_vol;
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
			e.de_vol = (n1.de_vol + n2.de_vol + n3.de_vol + n4.de_vol) * 0.25;
		}
	}

	double ds11, ds22, ds33, ds12, ds23, ds31;
	double de12, de23, de31;
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
			pcl.vx += (n1.ax * pcl.N1 + n2.ax * pcl.N2 + n3.ax * pcl.N3 + n4.ax * pcl.N4) * self.dtime;
			pcl.vy += (n1.ay * pcl.N1 + n2.ay * pcl.N2 + n3.ay * pcl.N3 + n4.ay * pcl.N4) * self.dtime;
			pcl.vz += (n1.az * pcl.N1 + n2.az * pcl.N2 + n3.az * pcl.N3 + n4.az * pcl.N4) * self.dtime;

			// displacement
			pcl.ux += n1.dux * pcl.N1 + n2.dux * pcl.N2 + n3.dux * pcl.N3 + n4.dux * pcl.N4;
			pcl.uy += n1.duy * pcl.N1 + n2.duy * pcl.N2 + n3.duy * pcl.N3 + n4.duy * pcl.N4;
			pcl.uz += n1.duz * pcl.N1 + n2.duz * pcl.N2 + n3.duz * pcl.N3 + n4.duz * pcl.N4;

			// position
			pcl.x = pcl.x_ori + pcl.ux;
			pcl.y = pcl.y_ori + pcl.uy;
			pcl.z = pcl.z_ori + pcl.uz;

			// strain
			de_vol_over_3 = e.de_vol / 3.0;
			de11 = e.dde11 + de_vol_over_3;
			de22 = e.dde22 + de_vol_over_3;
			de33 = e.dde33 + de_vol_over_3;
			de12 = e.de12;
			de23 = e.de23;
			de31 = e.de31;
			pcl.e11 += de11;
			pcl.e22 += de22;
			pcl.e33 += de33;
			pcl.e12 += de12;
			pcl.e23 += de23;
			pcl.e31 += de31;

			// stress
			double dstrain[6] = { de11, de22, de33, de12, de23, de31 };
			pcl.mm->integrate(dstrain);
			const double *dstress = pcl.mm->get_dstress();
			pcl.s11 += dstress[0];
			pcl.s22 += dstress[1];
			pcl.s33 += dstress[2];
			pcl.s12 += dstress[3];
			pcl.s23 += dstress[4];
			pcl.s31 += dstress[5];

			// density
			pcl.density /= 1.0 + e.de_vol;
		}
	}

	return 0;
}