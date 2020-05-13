#include "SimulationCore_pcp.h"

#include <cmath>
#include "Step_T2D_fluid.h"

int solve_substep_T2D_fluid_MI(void *_self)
{
	typedef Model_T2D_fluid::Node Node;
	typedef Model_T2D_fluid::Element Element;
	typedef Model_T2D_fluid::Particle Particle;
	Step_T2D_fluid &self = *(Step_T2D_fluid *)(_self);
	Model_T2D_fluid &md = *self.model;

	// init nodes
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		n.has_mp = false;
		n.m = 0.0;
		n.vx = 0.0;
		n.vy = 0.0;
		n.fx_ext = 0.0;
		n.fy_ext = 0.0;
		n.fx_int = 0.0;
		n.fy_int = 0.0;
		// strain enhancement
		n.vol = 0.0;
		n.de_vol = 0.0;
		n.de_vol_dt = 0.0;
	}

	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
		e.pcls = nullptr;
		e.vol = 0.0;
		e.p = 0.0;
	}

	// init particles
	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle &pcl = md.pcls[pcl_id];
		if (pcl.pe)
		{
			pcl.pe = md.find_in_which_element(pcl);
			if (pcl.pe == nullptr) continue;
			
			Element &e = *pcl.pe;
			e.add_pcl(pcl);
			pcl.vol = pcl.m / pcl.density;
			e.vol += pcl.vol;
			e.p += pcl.p * pcl.vol;

			double mvx = pcl.m * pcl.vx;
			double mvy = pcl.m * pcl.vy;
			// node 1
			Node &n1 = md.nodes[e.n1];
			n1.has_mp = true;
			n1.m += pcl.N1 * pcl.m;
			n1.vx += pcl.N1 * mvx;
			n1.vy += pcl.N1 * mvy;
			// node 2
			Node &n2 = md.nodes[e.n2];
			n2.has_mp = true;
			n2.m += pcl.N2 * pcl.m;
			n2.vx += pcl.N2 * mvx;
			n2.vy += pcl.N2 * mvy;
			// node 3
			Node &n3 = md.nodes[e.n3];
			n3.has_mp = true;
			n3.m += pcl.N3 * pcl.m;
			n3.vx += pcl.N3 * mvx;
			n3.vy += pcl.N3 * mvy;
		}
	}

	// update nodal acceleration of fluid pahse
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.has_mp)
		{
			n.vx /= n.m;
			n.vy /= n.m;
		}
	}

	// apply velocity bc
	for (size_t v_id = 0; v_id < md.vx_num; ++v_id)
	{
		Node &n = md.nodes[md.vxs[v_id].node_id];
		n.vx = md.vxs[v_id].v;
	}
	for (size_t v_id = 0; v_id < md.vy_num; ++v_id)
	{
		Node &n = md.nodes[md.vys[v_id].node_id];
		n.vy = md.vys[v_id].v;
	}

	// internal force
	double de11_dt, de22_dt, de12_dt, de_vol_dt;
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
		if (e.pcls)
		{
			Node &n1 = md.nodes[e.n1];
			Node &n2 = md.nodes[e.n2];
			Node &n3 = md.nodes[e.n3];

			e.p /= e.vol;
			if (e.vol > e.area)
				e.vol = e.area;

			de11_dt = (e.dN1_dx * n1.vx + e.dN2_dx * n2.vx + e.dN3_dx * n3.vx) * 2.0;
			de22_dt = (e.dN1_dy * n1.vy + e.dN2_dy * n2.vy + e.dN3_dy * n3.vy) * 2.0;
			de12_dt = e.dN1_dx * n1.vy + e.dN2_dx * n2.vy + e.dN3_dx * n3.vy
					+ e.dN1_dy * n1.vx + e.dN2_dy * n2.vx + e.dN3_dy * n3.vx;
			de_vol_dt = (de11_dt + de22_dt) * 0.5;

			e.t11 = md.miu * de11_dt - md.lambda * de_vol_dt;
			e.t22 = md.miu * de22_dt - md.lambda * de_vol_dt;
			e.t12 = md.miu * de12_dt;

			// node 1
			n1.fx_int += (e.dN1_dx * (e.t11 - e.p) + e.dN1_dy * e.t12) * e.vol;
			n1.fy_int += (e.dN1_dx * e.t12 + e.dN1_dy * (e.t22 - e.p)) * e.vol;
			// node 2
			n2.fx_int += (e.dN2_dx * (e.t11 - e.p) + e.dN2_dy * e.t12) * e.vol;
			n2.fy_int += (e.dN2_dx * e.t12 + e.dN2_dy * (e.t22 - e.p)) * e.vol;
			// node 3
			n3.fx_int += (e.dN3_dx * (e.t11 - e.p) + e.dN3_dy * e.t12) * e.vol;
			n3.fy_int += (e.dN3_dx * e.t12 + e.dN3_dy * (e.t22 - e.p)) * e.vol;
		}
	}

	// external force (body force)
	double bf_mag;
	for (size_t bf_id = 0; bf_id < md.bfx_num; ++bf_id)
	{
		BodyForce &bf = md.bfxs[bf_id];
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
		}
	}
	for (size_t bf_id = 0; bf_id < md.bfy_num; ++bf_id)
	{
		BodyForce &bf = md.bfys[bf_id];
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
		}
	}

	// update nodal acceleration of fluid pahse
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.has_mp)
		{
			n.ax = (n.fx_ext - n.fx_int) / n.m;
			n.ay = (n.fy_ext - n.fy_int) / n.m;
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

	// update nodal velocity
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.has_mp)
		{
			n.vx += n.ax * self.dtime;
			n.vy += n.ay * self.dtime;
		}
	}

	// reapply velocity bc
	for (size_t v_id = 0; v_id < md.vx_num; ++v_id)
	{
		Node &n = md.nodes[md.vxs[v_id].node_id];
		n.ax = 0.0;
		n.vx = md.vxs[v_id].v;
	}
	for (size_t v_id = 0; v_id < md.vy_num; ++v_id)
	{
		Node &n = md.nodes[md.vys[v_id].node_id];
		n.ay = 0.0;
		n.vy = md.vys[v_id].v;
	}

	// update displacement increment of both phases
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.has_mp)
		{
			n.dux = n.vx * self.dtime;
			n.duy = n.vy * self.dtime;
		}
	}

	// map variables back to particles and update their variables
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
		if (e.pcls)
		{
			Node &n1 = md.nodes[e.n1];
			Node &n2 = md.nodes[e.n2];
			Node &n3 = md.nodes[e.n3];
			e.de_vol = e.dN1_dx * n1.dux
					 + e.dN2_dx * n2.dux
					 + e.dN3_dx * n3.dux
					 + e.dN1_dy * n1.duy
					 + e.dN2_dy * n2.duy
					 + e.dN3_dy * n3.duy;
		}
	}

	// strain enhancement
	double vol_N;
	for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	{
		Particle &pcl = md.pcls[pcl_id];
		if (pcl.pe)
		{
			Element &e = *pcl.pe;
			Node &n1 = md.nodes[e.n1];
			vol_N = pcl.vol * pcl.N1;
			n1.vol += vol_N;
			n1.de_vol += vol_N * e.de_vol;
			Node &n2 = md.nodes[e.n2];
			vol_N = pcl.vol * pcl.N2;
			n2.vol += vol_N;
			n2.de_vol += vol_N * e.de_vol;
			Node &n3 = md.nodes[e.n3];
			vol_N = pcl.vol * pcl.N3;
			n3.vol += vol_N;
			n3.de_vol += vol_N * e.de_vol;
		}
	}

	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		if (n.has_mp)
			n.de_vol /= n.vol;
	}

	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
		if (e.pcls)
		{
			Node &n1 = md.nodes[e.n1];
			Node &n2 = md.nodes[e.n2];
			Node &n3 = md.nodes[e.n3];
			// de_vol takes tension as positive
			// p takes compression as positive
			e.de_vol = (n1.de_vol + n2.de_vol + n3.de_vol) / 3.0;
			e.p -= md.Kf * e.de_vol;
		}
	}

	double de_vol;
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
			pcl.vx += (n1.ax * pcl.N1 + n2.ax * pcl.N2 + n3.ax * pcl.N3) * self.dtime;
			pcl.vy += (n1.ay * pcl.N1 + n2.ay * pcl.N2 + n3.ay * pcl.N3) * self.dtime;

			// displacement
			pcl.ux += n1.dux * pcl.N1 + n2.dux * pcl.N2 + n3.dux * pcl.N3;
			pcl.uy += n1.duy * pcl.N1 + n2.duy * pcl.N2 + n3.duy * pcl.N3;

			// update position
			pcl.x = pcl.x_ori + pcl.ux;
			pcl.y = pcl.y_ori + pcl.uy;

			// pore pressure
			de_vol = -(e.p - pcl.p) / md.Kf;
			pcl.p = e.p;
			
			// density
			pcl.density /= 1.0 + de_vol;
		}
	}

	return 0;
}

