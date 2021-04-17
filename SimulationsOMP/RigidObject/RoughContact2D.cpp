#include "SimulationsOMP_pcp.h"

#include "RoughContact2D.h"

RoughContact2D::RoughContact2D() :
	Kn_cont(0.0), Kt_cont(0.0) {}

RoughContact2D::~RoughContact2D() {}

void RoughContact2D::cal_contact_force(
	size_t substp_id,
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
	// allow overlapping of 0.1 pcl_len
	// reduce oscillation
	dist -= 0.1 * pcl_len;
	const double f_cont = dist > 0.0 ? Kn_cont * pcl_len * dist : 0.0;
	cont_force.x = f_cont * norm.x;
	cont_force.y = f_cont * norm.y;
	// tangential force
	if (cont_substp_id != substp_id)
	{
		if (dist > 0.0)
		{
			// not in contacct
			// record contact position
			prev_cont_pos.x = cont_pos.x;
			prev_cont_pos.y = cont_pos.y;
			cont_force.x = 0.0;
			cont_force.y = 0.0;
			cont_substp_id = substp_id + 1;
		}
	}
	else
	{
		// previously in contactt
		// add tangential force to contact force
		const double Kt_len = Kt_cont * pcl_len;
		cont_force.x = (cont_pos.x - prev_cont_pos.x) * Kt_len;
		cont_force.y = (cont_pos.y - prev_cont_pos.y) * Kt_len;
		cont_substp_id = substp_id + 1;
	}
}
