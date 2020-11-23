#include "SimulationsOMP_pcp.h"

#include "SmoothContact3D.h"

SmoothContact3D::SmoothContact3D() : Kn_cont(0.0) {}

SmoothContact3D::~SmoothContact3D() {}

void SmoothContact3D::cal_contact_force(
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
	double f_cont = Kn_cont * dist;
	cont_force.x = f_cont * norm.x;
	cont_force.y = f_cont * norm.y;
	cont_force.z = f_cont * norm.z;
}
