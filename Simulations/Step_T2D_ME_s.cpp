#include "Simulations_pcp.h"

#include <iostream>
#include <assert.h>

#include <cmath>
#include "MaterialModel.h"

#include "Step_T2D_ME_s.h"

Step_T2D_ME_s::Step_T2D_ME_s(const char* _name) :
	Step(_name, "Step_T2D_ME_s", &solve_substep_T2D_ME_s),
	model(nullptr), damping_ratio(0.0) {}

Step_T2D_ME_s::~Step_T2D_ME_s() {}

int Step_T2D_ME_s::init_calculation()
{
	Model_T2D_ME_s &md = *model;

	if (is_first_step) {}

	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle &pcl = md.pcls[pcl_id];
		pcl.pe = md.elems;
		pcl.vol = pcl.m / pcl.density;
		pcl.x_ori = pcl.x;
		pcl.y_ori = pcl.y;
		pcl.ux = 0.0;
		pcl.uy = 0.0;
	}

	return 0;
}

int Step_T2D_ME_s::finalize_calculation() { return 0; }

int solve_substep_T2D_ME_s(void *_self)
{
	typedef Model_T2D_ME_s::Particle Particle_mpm;
	typedef Model_T2D_ME_s::Element Element_mpm;
	typedef Model_T2D_ME_s::Node Node_mpm;

	Step_T2D_ME_s &self = *(Step_T2D_ME_s *)(_self);
	Model_T2D_ME_s &md = *self.model;

	// init nodes
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node_mpm &n = md.nodes[n_id];
		n.has_mp = false;
		n.m = 0.0;
		n.vx = 0.0;
		n.vy = 0.0;
		n.fx_ext = 0.0;
		n.fy_ext = 0.0;
		n.fx_int = 0.0;
		n.fy_int = 0.0;
		// strain enhancement
		n.pcl_vol = 0.0;
		n.de_vol = 0.0;
	}

	// init elements
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element_mpm &e = md.elems[e_id];
		e.pcls = nullptr;
		e.pcl_vol = 0.0;
		e.s11 = 0.0;
		e.s22 = 0.0;
		e.s12 = 0.0;
	}

	// init particles
	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle_mpm &pcl = md.pcls[pcl_id];
		if (pcl.pe)
		{
			pcl.pe = md.find_in_which_element(pcl);
			if (!pcl.pe)
				continue;
			pcl.pe->add_pcl(pcl);

			Element_mpm &e = *pcl.pe;
			pcl.vol = pcl.m / pcl.density;
			e.pcl_vol += pcl.vol;
			e.s11 += pcl.vol * pcl.s11;
			e.s22 += pcl.vol * pcl.s22;
			e.s12 += pcl.vol * pcl.s12;

			double mvx = pcl.m * pcl.vx;
			double mvy = pcl.m * pcl.vy;
			// node 1
			Node_mpm &n1 = md.nodes[e.n1];
			n1.has_mp = true;
			n1.m  += pcl.N1 * pcl.m;
			n1.vx += pcl.N1 * mvx;
			n1.vy += pcl.N1 * mvy;
			// node 2
			Node_mpm &n2 = md.nodes[e.n2];
			n2.has_mp = true;
			n2.m  += pcl.N2 * pcl.m;
			n2.vx += pcl.N2 * mvx;
			n2.vy += pcl.N2 * mvy;
			// node 3
			Node_mpm &n3 = md.nodes[e.n3];
			n3.has_mp = true;
			n3.m  += pcl.N3 * pcl.m;
			n3.vx += pcl.N3 * mvx;
			n3.vy += pcl.N3 * mvy;
		}
	}

	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element_mpm &e = md.elems[e_id];
		if (e.pcls)
		{
			e.s11 /= e.pcl_vol;
			e.s22 /= e.pcl_vol;
			e.s12 /= e.pcl_vol;
			if (e.pcl_vol > e.area)
				e.pcl_vol = e.area;

			Node_mpm &n1 = md.nodes[e.n1];
			Node_mpm &n2 = md.nodes[e.n2];
			Node_mpm &n3 = md.nodes[e.n3];
			// node 1
			n1.fx_int += (e.dN1_dx * e.s11 + e.dN1_dy * e.s12) * e.pcl_vol;
			n1.fy_int += (e.dN1_dx * e.s12 + e.dN1_dy * e.s22) * e.pcl_vol;
			// node 2
			n2.fx_int += (e.dN2_dx * e.s11 + e.dN2_dy * e.s12) * e.pcl_vol;
			n2.fy_int += (e.dN2_dx * e.s12 + e.dN2_dy * e.s22) * e.pcl_vol;
			// node 3
			n3.fx_int += (e.dN3_dx * e.s11 + e.dN3_dy * e.s12) * e.pcl_vol;
			n3.fy_int += (e.dN3_dx * e.s12 + e.dN3_dy * e.s22) * e.pcl_vol;
		}
	}

	// body force
	double bf_mag;
	for (size_t bf_id = 0; bf_id < md.bfx_num; ++bf_id)
	{
		BodyForceAtPcl &bf = md.bfxs[bf_id];
		Particle_mpm &pcl = md.pcls[bf.pcl_id];
		if (pcl.pe)
		{
			Element_mpm &e = *pcl.pe;
			bf_mag = pcl.m * bf.bf;
			// node 1
			Node_mpm &n1 = md.nodes[e.n1];
			n1.fx_ext += pcl.N1 * bf_mag;
			// node 2
			Node_mpm &n2 = md.nodes[e.n2];
			n2.fx_ext += pcl.N2 * bf_mag;
			// node 3
			Node_mpm &n3 = md.nodes[e.n3];
			n3.fx_ext += pcl.N3 * bf_mag;
		}
	}
	for (size_t bf_id = 0; bf_id < md.bfy_num; ++bf_id)
	{
		BodyForceAtPcl &bf = md.bfys[bf_id];
		Particle_mpm &pcl = md.pcls[bf.pcl_id];
		if (pcl.pe)
		{
			Element_mpm &e = *pcl.pe;
			bf_mag = pcl.m * bf.bf;
			// node 1
			Node_mpm &n1 = md.nodes[e.n1];
			n1.fy_ext += pcl.N1 * bf_mag;
			// node 2
			Node_mpm &n2 = md.nodes[e.n2];
			n2.fy_ext += pcl.N2 * bf_mag;
			// node 3
			Node_mpm &n3 = md.nodes[e.n3];
			n3.fy_ext += pcl.N3 * bf_mag;
		}
	}

	// surface force
	for (size_t tf_id = 0; tf_id < md.tx_num; ++tf_id)
	{
		TractionBCAtPcl &tf = md.txs[tf_id];
		Particle_mpm &pcl = md.pcls[tf.pcl_id];
		if (pcl.pe)
		{
			Element_mpm &e = *pcl.pe;
			// node 1
			Node_mpm &n1 = md.nodes[e.n1];
			n1.fx_ext += pcl.N1 * tf.t;
			// node 2
			Node_mpm &n2 = md.nodes[e.n2];
			n2.fx_ext += pcl.N2 * tf.t;
			// node 3
			Node_mpm &n3 = md.nodes[e.n3];
			n3.fx_ext += pcl.N3 * tf.t;
		}
	}
	for (size_t tf_id = 0; tf_id < md.ty_num; ++tf_id)
	{
		TractionBCAtPcl &tf = md.tys[tf_id];
		Particle_mpm &pcl = md.pcls[tf.pcl_id];
		if (pcl.pe)
		{
			Element_mpm &e = *pcl.pe;
			// node 1
			Node_mpm &n1 = md.nodes[e.n1];
			n1.fy_ext += pcl.N1 * tf.t;
			// node 2
			Node_mpm &n2 = md.nodes[e.n2];
			n2.fy_ext += pcl.N2 * tf.t;
			// node 3
			Node_mpm &n3 = md.nodes[e.n3];
			n3.fy_ext += pcl.N3 * tf.t;
		}
	}

	// update nodal acceleration of fluid pahse
	double nf, v_sign;
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node_mpm &n = md.nodes[n_id];
		if (n.has_mp) // or n.m_f != 0.0
		{
			// fx
			nf = n.fx_ext - n.fx_int;
			n.ax = (nf - self.damping_ratio * abs(nf) * get_sign(n.vx)) / n.m;
			// fy
			nf = n.fy_ext - n.fy_int;
			n.ay = (nf - self.damping_ratio * abs(nf) * get_sign(n.vy)) / n.m;
		}
	}

	// apply acceleration bc
	for (size_t a_id = 0; a_id < md.ax_num; ++a_id)
	{
		Node_mpm &n = md.nodes[md.axs[a_id].node_id];
		n.ax = md.axs[a_id].a;
	}
	for (size_t a_id = 0; a_id < md.ay_num; ++a_id)
	{
		Node_mpm &n = md.nodes[md.ays[a_id].node_id];
		n.ay = md.ays[a_id].a;
	}

	// update nodal momentum
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node_mpm &n = md.nodes[n_id];
		if (n.has_mp)
		{
			n.vx /= n.m;
			n.vx += n.ax * self.dtime;
			n.vy /= n.m;
			n.vy += n.ay * self.dtime;
		}
	}

	// contact detection
	//md.apply_contact_force_to_bg_mesh(self.dtime);

	// apply velocity bc
	for (size_t v_id = 0; v_id < md.vx_num; ++v_id)
	{
		Node_mpm &n = md.nodes[md.vxs[v_id].node_id];
		n.vx = md.vxs[v_id].v;
		n.ax = 0.0;
	}
	for (size_t v_id = 0; v_id < md.vy_num; ++v_id)
	{
		Node_mpm &n = md.nodes[md.vys[v_id].node_id];
		n.vy = md.vys[v_id].v;
		n.ay = 0.0;
	}
	
	// update displacement increment of both phases
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node_mpm &n = md.nodes[n_id];
		if (n.has_mp)
		{
			n.dux = n.vx * self.dtime;
			n.duy = n.vy * self.dtime;
		}
	}

	// map variables back to particles and update their variables
	double de11, de22, de12, de_vol;
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element_mpm &e = md.elems[e_id];
		if (e.pcls)
		{
			Node_mpm &n1 = md.nodes[e.n1];
			Node_mpm &n2 = md.nodes[e.n2];
			Node_mpm &n3 = md.nodes[e.n3];

			// strain increment
			de11 = n1.dux * e.dN1_dx + n2.dux * e.dN2_dx + n3.dux * e.dN3_dx;
			de22 = n1.duy * e.dN1_dy + n2.duy * e.dN2_dy + n3.duy * e.dN3_dy;
			de12 = (n1.dux * e.dN1_dy + n2.dux * e.dN2_dy + n3.dux * e.dN3_dy
				  + n1.duy * e.dN1_dx + n2.duy * e.dN2_dx + n3.duy * e.dN3_dx) * 0.5;
			de_vol = (de11 + de22) / 3.0;
			e.dde11 = de11 - de_vol;
			e.dde22 = de22 - de_vol;
			e.de12 = de12;
			e.de_vol = de_vol;
		}
	}

	// strain enhancement
	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle_mpm &pcl = md.pcls[pcl_id];
		if (pcl.pe)
		{
			Element_mpm &e = *pcl.pe;
			double vol_de_vol = pcl.vol * e.de_vol;
			// node 1
			Node_mpm &n1 = md.nodes[e.n1];
			n1.pcl_vol += pcl.N1 * pcl.vol;
			n1.de_vol += pcl.N1 * vol_de_vol;
			// node 2
			Node_mpm &n2 = md.nodes[e.n2];
			n2.pcl_vol += pcl.N2 * pcl.vol;
			n2.de_vol += pcl.N2 * vol_de_vol;
			// node 3
			Node_mpm &n3 = md.nodes[e.n3];
			n3.pcl_vol += pcl.N3 * pcl.vol;
			n3.de_vol += pcl.N3 * vol_de_vol;
		}
	}

	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node_mpm &n = md.nodes[n_id];
		if (n.has_mp)
			n.de_vol /= n.pcl_vol;
	}

	double ds11, ds22, ds12;
	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle_mpm &pcl = md.pcls[pcl_id];
		Element_mpm &e = *pcl.pe;
		if (pcl.pe)
		{
			Node_mpm &n1 = md.nodes[e.n1];
			Node_mpm &n2 = md.nodes[e.n2];
			Node_mpm &n3 = md.nodes[e.n3];

			// velocity
			pcl.vx += (n1.ax * pcl.N1 + n2.ax * pcl.N2 + n3.ax * pcl.N3) * self.dtime;
			pcl.vy += (n1.ay * pcl.N1 + n2.ay * pcl.N2 + n3.ay * pcl.N3) * self.dtime;

			// displacement
			pcl.ux += n1.dux * pcl.N1 + n2.dux * pcl.N2 + n3.dux * pcl.N3;
			pcl.uy += n1.duy * pcl.N1 + n2.duy * pcl.N2 + n3.duy * pcl.N3;

			// update position
			pcl.x = pcl.x_ori + pcl.ux;
			pcl.y = pcl.y_ori + pcl.uy;

			// strain
			//de_vol = n1.de_vol * pcl.N1 + n2.de_vol * pcl.N2 + n3.de_vol * pcl.N3;
			de_vol = (n1.de_vol + n2.de_vol + n3.de_vol) / 3.0;
			de11 = e.dde11 + de_vol;
			de22 = e.dde22 + de_vol;
			de12 = e.de12;
			pcl.e11 += de11;
			pcl.e22 += de22;
			pcl.e12 += de12;

			// stress
			// update stress using constitutive model
			double dstrain[6] = { de11, de22, 0.0, de12, 0.0, 0.0 };
			pcl.mm->integrate(dstrain);
			const double *dstress = pcl.mm->get_dstress();
			pcl.s11 += dstress[0];
			pcl.s22 += dstress[1];
			pcl.s12 += dstress[3];

			// density
			pcl.density /= 1.0 + de_vol;
		}
	}
	
	return 0;
}
