#include "SimulationsOMP_pcp.h"

#include "NoContact3D.h"

NoContact3D::NoContact3D() {}

NoContact3D::~NoContact3D() {}

void NoContact3D::cal_contact_force(
	size_t substp_id,
	double dist,
	const Vector3D& norm,
	const Point3D& cont_pos,
	double pcl_len,
	ParticleVariablesGetter& pv_getter,
	size_t& cont_substp_id,
	Point3D& prev_cont_pos,
	Vector3D& prev_cont_tan_force,
	Vector3D& cont_force
	)
{
	cont_force.x = 0.0;
	cont_force.y = 0.0;
	cont_force.z = 0.0;
}
