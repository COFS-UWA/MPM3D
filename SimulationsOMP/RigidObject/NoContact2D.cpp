#include "SimulationsOMP_pcp.h"

#include "NoContact2D.h"

NoContact2D::NoContact2D() {}

NoContact2D::~NoContact2D() {}

void NoContact2D::cal_contact_force(
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
	cont_force.x = 0.0;
	cont_force.y = 0.0;
}
