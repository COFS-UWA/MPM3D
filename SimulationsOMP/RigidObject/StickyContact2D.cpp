#include "SimulationsOMP_pcp.h"

#include "StickyContact2D.h"

StickyContact2D::StickyContact2D() :
	Kn_cont(0.0), Kt_cont(0.0), shear_strength(0.0) {}

StickyContact2D::~StickyContact2D() {}

void StickyContact2D::cal_contact_force(
	// in
	size_t substp_id,
	double dist,
	const Vector2D& norm,
	const Point2D& cont_pos,
	double pcl_len,
	ParticleVariablesGetter& pv_getter,
	// out
	size_t& cont_substp_id,
	Point2D& prev_cont_pos,
	Vector2D& prev_cont_tan_force,
	Vector2D& cont_force)
{
	// normal force
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
		// rest contact force
		prev_cont_tan_force.x = 0.0;
		prev_cont_tan_force.y = 0.0;
	}
	else
	{
		// previously in contactt
		// cal relative movement
		double rx = cont_pos.x - prev_cont_pos.x;
		double ry = cont_pos.y - prev_cont_pos.y;
		prev_cont_pos.x = cont_pos.x;
		prev_cont_pos.y = cont_pos.y;
		const double norm_len = rx * norm.x + ry * norm.y;
		rx -= norm_len * norm.x;
		ry -= norm_len * norm.y;
		const double Kt_area = Kt_cont * pcl_len;
		prev_cont_tan_force.x -= rx * Kt_area;
		prev_cont_tan_force.y -= ry * Kt_area;
		// contact consitutive model
		double tan_force = sqrt(prev_cont_tan_force.x * prev_cont_tan_force.x
							  + prev_cont_tan_force.y * prev_cont_tan_force.y);
		if (tan_force > shear_strength)
		{
			const double ratio = shear_strength / tan_force;
			prev_cont_tan_force.x *= ratio;
			prev_cont_tan_force.y *= ratio;
		}
		// add tangential force to contact force
		cont_force.x += prev_cont_tan_force.x;
		cont_force.y += prev_cont_tan_force.y;
	}
	cont_substp_id = substp_id + 1;
}
