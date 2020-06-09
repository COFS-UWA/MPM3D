#include "Simulations_pcp.h"

#include <iostream>
#include <assert.h>

#include <cmath>
#include "MaterialModel.h"

#include "Step_FEM_T3D_ME_s.h"

Step_FEM_T3D_ME_s::Step_FEM_T3D_ME_s(const char *_name) :
	Step(_name, "Step_FEM_T3D_ME_s", &solve_substep_FEM_T3D_ME_s),
	model(nullptr), damping_ratio(0.0) {}

Step_FEM_T3D_ME_s::~Step_FEM_T3D_ME_s() {}

int Step_FEM_T3D_ME_s::init_calculation()
{
	Model_FEM_T3D_ME_s &md = *model;

	if (is_first_step) {}

	// init nodes
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node& n = md.nodes[n_id];
		n.ux = 0.0;
		n.uy = 0.0;
		n.uz = 0.0;
	}

	return 0;
}

int Step_FEM_T3D_ME_s::finalize_calculation() { return 0; }

int solve_substep_FEM_T3D_ME_s(void *_self)
{
	typedef Model_FEM_T3D_ME_s::Particle Particle;
	typedef Model_FEM_T3D_ME_s::Element Element;
	typedef Model_FEM_T3D_ME_s::Node Node;

	Step_FEM_T3D_ME_s &self = *(Step_FEM_T3D_ME_s *)(_self);
	Model_FEM_T3D_ME_s &md = *self.model;

	// init nodes
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		n.m = 0.0;
		n.fx_ext = 0.0;
		n.fy_ext = 0.0;
		n.fz_ext = 0.0;
		n.fx_int = 0.0;
		n.fy_int = 0.0;
		n.fz_int = 0.0;
		// strain enhancement
		n.pcl_w = 0.0;
		n.de_vol = 0.0;
	}

	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element& e = md.elems[e_id];
		Node& n1 = md.nodes[e.n1];
		Node& n2 = md.nodes[e.n2];
		Node& n3 = md.nodes[e.n3];
		Node& n4 = md.nodes[e.n4];
		Particle& p1 = *e.p1;
		Particle& p2 = *e.p2;
		Particle& p3 = *e.p3;
		Particle& p4 = *e.p4;

		// cal nodal mass
		n1.m += e.density * (p1.N1*p1.w + p2.N1*p2.w + p3.N1*p3.w + p4.N1*p4.w) * e.vol;
		n2.m += e.density * (p1.N2*p1.w + p2.N2*p2.w + p3.N2*p3.w + p4.N2*p4.w) * e.vol;
		n3.m += e.density * (p1.N3 * p1.w + p2.N3 * p2.w + p3.N3 * p3.w + p4.N3 * p4.w) * e.vol;
		n4.m += e.density * (p1.N4 * p1.w + p2.N4 * p2.w + p3.N4 * p3.w + p4.N4 * p4.w) * e.vol;

		// mixed integration
		e.s11 = p1.w * p1.s11 + p2.w * p2.s11 + p3.w * p3.s11 + p4.w * p4.s11;
		e.s22 = p1.w * p1.s22 + p2.w * p2.s22 + p3.w * p3.s22 + p4.w * p4.s22;
		e.s33 = p1.w * p1.s33 + p2.w * p2.s33 + p3.w * p3.s33 + p4.w * p4.s33;
		e.s12 = p1.w * p1.s12 + p2.w * p2.s12 + p3.w * p3.s12 + p4.w * p4.s12;
		e.s23 = p1.w * p1.s23 + p2.w * p2.s23 + p3.w * p3.s23 + p4.w * p4.s23;
		e.s31 = p1.w * p1.s31 + p2.w * p2.s31 + p3.w * p3.s31 + p4.w * p4.s31;
		// internal force
		// node 1
		n1.fx_int += (e.dN1_dx * e.s11 + e.dN1_dy * e.s12 + e.dN1_dz * e.s31) * e.vol;
		n1.fy_int += (e.dN1_dx * e.s12 + e.dN1_dy * e.s22 + e.dN1_dz * e.s23) * e.vol;
		n1.fz_int += (e.dN1_dx * e.s31 + e.dN1_dy * e.s23 + e.dN1_dz * e.s33) * e.vol;
		// node 2
		n2.fx_int += (e.dN2_dx * e.s11 + e.dN2_dy * e.s12 + e.dN2_dz * e.s31) * e.vol;
		n2.fy_int += (e.dN2_dx * e.s12 + e.dN2_dy * e.s22 + e.dN2_dz * e.s23) * e.vol;
		n2.fz_int += (e.dN2_dx * e.s31 + e.dN2_dy * e.s23 + e.dN2_dz * e.s33) * e.vol;
		// node 3
		n3.fx_int += (e.dN3_dx * e.s11 + e.dN3_dy * e.s12 + e.dN3_dz * e.s31) * e.vol;
		n3.fy_int += (e.dN3_dx * e.s12 + e.dN3_dy * e.s22 + e.dN3_dz * e.s23) * e.vol;
		n3.fz_int += (e.dN3_dx * e.s31 + e.dN3_dy * e.s23 + e.dN3_dz * e.s33) * e.vol;
		// node 4
		n4.fx_int += (e.dN4_dx * e.s11 + e.dN4_dy * e.s12 + e.dN4_dz * e.s31) * e.vol;
		n4.fy_int += (e.dN4_dx * e.s12 + e.dN4_dy * e.s22 + e.dN4_dz * e.s23) * e.vol;
		n4.fz_int += (e.dN4_dx * e.s31 + e.dN4_dy * e.s23 + e.dN4_dz * e.s33) * e.vol;
	}

	// body force
	// x direction
	for (size_t bf_id = 0; bf_id < md.bfx_num; ++bf_id)
	{
		BodyForceAtElem &bf = md.bfxs[bf_id];
		Element &e = md.elems[bf.elem_id];
		Particle &p1 = *e.p1;
		Particle &p2 = *e.p2;
		Particle &p3 = *e.p3;
		Particle &p4 = *e.p4;
		// node 1
		Node &n1 = md.nodes[e.n1];
		n1.fx_ext += bf.bf * (p1.N1*p1.w + p2.N1*p2.w + p3.N1*p3.w + p4.N1*p4.w) * e.vol;
		// node 2
		Node &n2 = md.nodes[e.n2];
		n2.fx_ext += bf.bf * (p1.N2*p1.w + p2.N2*p2.w + p3.N2*p3.w + p4.N2*p4.w) * e.vol;
		// node 3
		Node &n3 = md.nodes[e.n3];
		n3.fx_ext += bf.bf * (p1.N3*p1.w + p2.N3*p2.w + p3.N3*p3.w + p4.N3*p4.w) * e.vol;
		// node 4
		Node &n4 = md.nodes[e.n4];
		n4.fx_ext += bf.bf * (p1.N4*p1.w + p2.N4*p2.w + p3.N4*p3.w + p4.N4*p4.w) * e.vol;
	}
	// y direction
	for (size_t bf_id = 0; bf_id < md.bfy_num; ++bf_id)
	{
		BodyForceAtElem& bf = md.bfys[bf_id];
		Element &e = md.elems[bf.elem_id];
		Particle &p1 = *e.p1;
		Particle &p2 = *e.p2;
		Particle &p3 = *e.p3;
		Particle &p4 = *e.p4;
		// node 1
		Node& n1 = md.nodes[e.n1];
		n1.fy_ext += bf.bf * (p1.N1 * p1.w + p2.N1 * p2.w + p3.N1 * p3.w + p4.N1 * p4.w) * e.vol;
		// node 2
		Node& n2 = md.nodes[e.n2];
		n2.fy_ext += bf.bf * (p1.N2 * p1.w + p2.N2 * p2.w + p3.N2 * p3.w + p4.N2 * p4.w) * e.vol;
		// node 3
		Node& n3 = md.nodes[e.n3];
		n3.fy_ext += bf.bf * (p1.N3 * p1.w + p2.N3 * p2.w + p3.N3 * p3.w + p4.N3 * p4.w) * e.vol;
		// node 4
		Node& n4 = md.nodes[e.n4];
		n4.fy_ext += bf.bf * (p1.N4 * p1.w + p2.N4 * p2.w + p3.N4 * p3.w + p4.N4 * p4.w) * e.vol;
	}
	// z direction
	for (size_t bf_id = 0; bf_id < md.bfz_num; ++bf_id)
	{
		BodyForceAtElem& bf = md.bfzs[bf_id];
		Element& e = md.elems[bf.elem_id];
		Particle& p1 = *e.p1;
		Particle& p2 = *e.p2;
		Particle& p3 = *e.p3;
		Particle& p4 = *e.p4;
		// node 1
		Node& n1 = md.nodes[e.n1];
		n1.fz_ext += bf.bf * (p1.N1*p1.w + p2.N1*p2.w + p3.N1*p3.w + p4.N1*p4.w) * e.vol;
		// node 2
		Node& n2 = md.nodes[e.n2];
		n2.fz_ext += bf.bf * (p1.N2*p1.w + p2.N2*p2.w + p3.N2*p3.w + p4.N2*p4.w) * e.vol;
		// node 3
		Node& n3 = md.nodes[e.n3];
		n3.fz_ext += bf.bf * (p1.N3*p1.w + p2.N3*p2.w + p3.N3*p3.w + p4.N3*p4.w) * e.vol;
		// node 4
		Node& n4 = md.nodes[e.n4];
		n4.fz_ext += bf.bf * (p1.N4*p1.w + p2.N4*p2.w + p3.N4*p3.w + p4.N4*p4.w) * e.vol;
	}

	// surface force
	double force;
	for (size_t tf_id = 0; tf_id < md.tx_num; ++tf_id)
	{
		TractionBCAtFace &tf = md.txs[tf_id];
		Element& e = md.elems[tf.elem_id];
		Node& n1 = md.nodes[e.n1];
		Node& n2 = md.nodes[e.n2];
		Node& n3 = md.nodes[e.n3];
		Node& n4 = md.nodes[e.n4];
		switch (tf.face_id)
		{
		case 0:
			force = cal_triangle_area(n1, n2, n3) * tf.t / 3.0;
			n1.fx_ext += force;
			n2.fx_ext += force;
			n3.fx_ext += force;
			break;
		case 1:
			force = cal_triangle_area(n1, n4, n2) * tf.t / 3.0;
			n1.fx_ext += force;
			n4.fx_ext += force;
			n2.fx_ext += force;
			break;
		case 2:
			force = cal_triangle_area(n1, n3, n4) * tf.t / 3.0;
			n1.fx_ext += force;
			n3.fx_ext += force;
			n4.fx_ext += force;
			break;
		case 3:
			force = cal_triangle_area(n2, n4, n3) * tf.t / 3.0;
			n2.fx_ext += force;
			n4.fx_ext += force;
			n3.fx_ext += force;
			break;
		default:
			break;
		}
	}
	for (size_t tf_id = 0; tf_id < md.ty_num; ++tf_id)
	{
		TractionBCAtFace& tf = md.tys[tf_id];
		Element& e = md.elems[tf.elem_id];
		Node& n1 = md.nodes[e.n1];
		Node& n2 = md.nodes[e.n2];
		Node& n3 = md.nodes[e.n3];
		Node& n4 = md.nodes[e.n4];
		switch (tf.face_id)
		{
		case 0:
			force = cal_triangle_area(n1, n2, n3) * tf.t / 3.0;
			n1.fy_ext += force;
			n2.fy_ext += force;
			n3.fy_ext += force;
			break;
		case 1:
			force = cal_triangle_area(n1, n4, n2) * tf.t / 3.0;
			n1.fy_ext += force;
			n4.fy_ext += force;
			n2.fy_ext += force;
			break;
		case 2:
			force = cal_triangle_area(n1, n3, n4) * tf.t / 3.0;
			n1.fy_ext += force;
			n3.fy_ext += force;
			n4.fy_ext += force;
			break;
		case 3:
			force = cal_triangle_area(n2, n4, n3) * tf.t / 3.0;
			n2.fy_ext += force;
			n4.fy_ext += force;
			n3.fy_ext += force;
			break;
		default:
			break;
		}
	}
	for (size_t tf_id = 0; tf_id < md.tz_num; ++tf_id)
	{
		TractionBCAtFace& tf = md.tzs[tf_id];
		Element& e = md.elems[tf.elem_id];
		Node& n1 = md.nodes[e.n1];
		Node& n2 = md.nodes[e.n2];
		Node& n3 = md.nodes[e.n3];
		Node& n4 = md.nodes[e.n4];
		switch (tf.face_id)
		{
		case 0:
			force = cal_triangle_area(n1, n2, n3) * tf.t / 3.0;
			n1.fz_ext += force;
			n2.fz_ext += force;
			n3.fz_ext += force;
			break;
		case 1:
			force = cal_triangle_area(n1, n4, n2) * tf.t / 3.0;
			n1.fz_ext += force;
			n4.fz_ext += force;
			n2.fz_ext += force;
			break;
		case 2:
			force = cal_triangle_area(n1, n3, n4) * tf.t / 3.0;
			n1.fz_ext += force;
			n3.fz_ext += force;
			n4.fz_ext += force;
			break;
		case 3:
			force = cal_triangle_area(n2, n4, n3) * tf.t / 3.0;
			n2.fz_ext += force;
			n4.fz_ext += force;
			n3.fz_ext += force;
			break;
		default:
			break;
		}
	}

	double nf;
	for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	{
		Node &n = md.nodes[n_id];
		// fx
		nf = n.fx_ext - n.fx_int;
		n.ax = (nf - self.damping_ratio * abs(nf) * get_sign(n.vx)) / n.m;
		// fy
		nf = n.fy_ext - n.fy_int;
		n.ay = (nf - self.damping_ratio * abs(nf) * get_sign(n.vy)) / n.m;
		// fz
		nf = n.fz_ext - n.fz_int;
		n.az = (nf - self.damping_ratio * abs(nf) * get_sign(n.vz)) / n.m;
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
		n.vx += n.ax * self.dtime;
		n.vy += n.ay * self.dtime;
		n.vz += n.az * self.dtime;
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
		n.dux = n.vx * self.dtime;
		n.ux += n.dux;
		n.duy = n.vy * self.dtime;
		n.uy += n.duy;
		n.duz = n.vz * self.dtime;
		n.uz += n.duz;
	}

	union
	{
		double dstrain[6];
		struct { double de11, de22, de33, de12, de23, de31; };
	};
	double de_vol_over_3;
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element &e = md.elems[e_id];
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

	//// strain enhancement
	//for (size_t pcl_id = 0; pcl_id < md.pcl_num; ++pcl_id)
	//{
	//	Particle &pcl = md.pcls[pcl_id];
	//	Element &e = *pcl.pe;
	//	double vol_de_vol = pcl.w * e.de_vol;
	//	// node 1
	//	Node &n1 = md.nodes[e.n1];
	//	n1.pcl_w += pcl.N1 * pcl.w;
	//	n1.de_vol += pcl.N1 * vol_de_vol;
	//	// node 2
	//	Node &n2 = md.nodes[e.n2];
	//	n2.pcl_w += pcl.N2 * pcl.w;
	//	n2.de_vol += pcl.N2 * vol_de_vol;
	//	// node 3
	//	Node &n3 = md.nodes[e.n3];
	//	n3.pcl_w += pcl.N3 * pcl.w;
	//	n3.de_vol += pcl.N3 * vol_de_vol;
	//	// node 4
	//	Node &n4 = md.nodes[e.n4];
	//	n4.pcl_w += pcl.N4 * pcl.w;
	//	n4.de_vol += pcl.N4 * vol_de_vol;
	//}

	//for (size_t n_id = 0; n_id < md.node_num; ++n_id)
	//{
	//	Node &n = md.nodes[n_id];
	//	n.de_vol /= n.pcl_w;
	//}

	//for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	//{
	//	Element &e = md.elems[e_id];
	//	Node &n1 = md.nodes[e.n1];
	//	Node &n2 = md.nodes[e.n2];
	//	Node &n3 = md.nodes[e.n3];
	//	Node &n4 = md.nodes[e.n4];
	//	e.de_vol = (n1.de_vol + n2.de_vol + n3.de_vol + n4.de_vol) * 0.25;
	//}

	const double* dstress;
	for (size_t e_id = 0; e_id < md.elem_num; ++e_id)
	{
		Element& e = md.elems[e_id];
		Node& n1 = md.nodes[e.n1];
		Node& n2 = md.nodes[e.n2];
		Node& n3 = md.nodes[e.n3];
		Node& n4 = md.nodes[e.n4];
		Particle& p1 = *e.p1;
		Particle& p2 = *e.p2;
		Particle& p3 = *e.p3;
		Particle& p4 = *e.p4;

		de_vol_over_3 = e.de_vol / 3.0;
		de11 = e.dde11 + de_vol_over_3;
		de22 = e.dde22 + de_vol_over_3;
		de33 = e.dde33 + de_vol_over_3;
		de12 = e.de12;
		de23 = e.de23;
		de31 = e.de31;

		// strain
		p1.e11 += de11;
		p1.e22 += de22;
		p1.e33 += de33;
		p1.e12 += de12;
		p1.e23 += de23;
		p1.e31 += de31;
		// stress
		p1.mm->integrate(dstrain);
		dstress = p1.mm->get_dstress();
		p1.s11 += dstress[0];
		p1.s22 += dstress[1];
		p1.s33 += dstress[2];
		p1.s12 += dstress[3];
		p1.s23 += dstress[4];
		p1.s31 += dstress[5];

		// strain
		p2.e11 += de11;
		p2.e22 += de22;
		p2.e33 += de33;
		p2.e12 += de12;
		p2.e23 += de23;
		p2.e31 += de31;
		// stress
		p2.mm->integrate(dstrain);
		dstress = p2.mm->get_dstress();
		p2.s11 += dstress[0];
		p2.s22 += dstress[1];
		p2.s33 += dstress[2];
		p2.s12 += dstress[3];
		p2.s23 += dstress[4];
		p2.s31 += dstress[5];

		// strain
		p3.e11 += de11;
		p3.e22 += de22;
		p3.e33 += de33;
		p3.e12 += de12;
		p3.e23 += de23;
		p3.e31 += de31;
		// stress
		p3.mm->integrate(dstrain);
		dstress = p3.mm->get_dstress();
		p3.s11 += dstress[0];
		p3.s22 += dstress[1];
		p3.s33 += dstress[2];
		p3.s12 += dstress[3];
		p3.s23 += dstress[4];
		p3.s31 += dstress[5];

		// strain
		p4.e11 += de11;
		p4.e22 += de22;
		p4.e33 += de33;
		p4.e12 += de12;
		p4.e23 += de23;
		p4.e31 += de31;
		// stress
		p4.mm->integrate(dstrain);
		dstress = p4.mm->get_dstress();
		p4.s11 += dstress[0];
		p4.s22 += dstress[1];
		p4.s33 += dstress[2];
		p4.s12 += dstress[3];
		p4.s23 += dstress[4];
		p4.s31 += dstress[5];
	}

	return 0;
}