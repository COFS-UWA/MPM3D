#include "SimulationsOMP_pcp.h"

#include "RoughContact3D.h"

RoughContact3D::RoughContact3D() :
	pv_getter(nullptr), Kn_cont(0.0), Kt_cont(0.0) {}

RoughContact3D::~RoughContact3D() {}

void RoughContact3D::cal_contact_force(
	size_t substp_id,
	double dist,
	const Vector3D& norm,
	const Point3D& cont_pos,
	ParticleVariablesGetter& pv_getter,
	size_t& cont_substp_id,
	Point3D& prev_cont_pos,
	Vector3D& prev_cont_tan_force,
	Vector3D& cont_force
	)
{
	// normal force
	// allow overlapping of 0.1 pcl_r
	// reduce oscillation
	dist -= 0.1 * pv_getter.get_square_r();
	double f_cont = dist > 0.0 ? Kn_cont * dist : 0.0;
	cont_force.x = f_cont * norm.x;
	cont_force.y = f_cont * norm.y;
	cont_force.z = f_cont * norm.z;
	// tangential force
	if (cont_substp_id != substp_id)
	{
		// not in contacct
		// record contact position
		prev_cont_pos.x = cont_pos.x;
		prev_cont_pos.y = cont_pos.y;
		prev_cont_pos.z = cont_pos.z;
	}
	else
	{
		// previously in contactt
		// add tangential force to contact force
		cont_force.x += (cont_pos.x - prev_cont_pos.x) * Kt_cont;
		cont_force.y += (cont_pos.y - prev_cont_pos.y) * Kt_cont;
		cont_force.z += (cont_pos.z - prev_cont_pos.z) * Kt_cont;
	}
	cont_substp_id = substp_id + 1;
}
