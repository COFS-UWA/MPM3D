#include "SimulationsOMP_pcp.h"

#include "RoughContact2D.h"

RoughContact2D::RoughContact2D() :
	Kn_cont(0.0), Kt_cont(0.0) {}

RoughContact2D::~RoughContact2D() {}

void RoughContact2D::cal_contact_force(
	size_t substp_id,
	size_t ori_pcl_id,
	double dist,
	const Vector2D& norm,
	const Point2D& cont_pos,
	double pcl_len,
	ParticleVariablesGetter& pv_getter,
	size_t& cont_substp_id,
	Point2D& prev_cont_pos,
	Vector2D& prev_cont_tan_force,
	Vector2D& cont_force)
{
	// normal force
	// allow small overlapping to reduce oscillation
	
	const double f_cont = Kn_cont * pcl_len * dist;
	cont_force.x = f_cont * norm.x;
	cont_force.y = f_cont * norm.y;
	// tangential force
	if (cont_substp_id != substp_id)
	{
		// not in contacct
		// record contact position
		prev_cont_pos.x = cont_pos.x;
		prev_cont_pos.y = cont_pos.y;
	}
	else
	{
		// previously in contactt
		// add tangential force to contact force
		const double Kt_len = Kt_cont * pcl_len;
		//
		//cont_force.x = -(cont_pos.x - prev_cont_pos.x) * Kt_len;
		//cont_force.y = -(cont_pos.y - prev_cont_pos.y) * Kt_len;
		//
		double rx = prev_cont_pos.x - cont_pos.x;
		double ry = prev_cont_pos.y - cont_pos.y;
		const double norm_len = rx * norm.x + ry * norm.y;
		rx -= norm_len * norm.x;
		ry -= norm_len * norm.y;
		const double ctfx = rx * Kt_len;
		const double ctfy = ry * Kt_len;
		// local damping
		const double dctfx_sign = sign(ctfx - prev_cont_tan_force.x);
		const double dctfy_sign = sign(ctfy - prev_cont_tan_force.y);
		prev_cont_tan_force.x = ctfx;
		prev_cont_tan_force.y = ctfy;
		cont_force.x += ctfx + dctfx_sign * abs(ctfx) * 0.05;
		cont_force.y += ctfy + dctfy_sign * abs(ctfy) * 0.05;
	}
	cont_substp_id = substp_id + 1;
}
