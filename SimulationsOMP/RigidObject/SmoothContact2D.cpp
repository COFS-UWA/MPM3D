#include "SimulationsOMP_pcp.h"

#include "SmoothContact2D.h"

SmoothContact2D::SmoothContact2D() :
	Kn_cont(0.0), prev_contact_dist(nullptr) {}

SmoothContact2D::~SmoothContact2D() {}

void SmoothContact2D::cal_contact_force(
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
	const double f_cont = Kn_cont * pcl_len * dist;
	cont_force.x = f_cont * norm.x;
	cont_force.y = f_cont * norm.y;
	double& pcd = prev_contact_dist[ori_pcl_id];
	if (cont_substp_id != substp_id) // not previously in contact
	{
		pcd = 0.0;
	}
	else // previously in contact
	{
		// local damping
		const double ddist_sign = sign(dist - pcd);
		pcd = dist;
		cont_force.x -= ddist_sign * abs(cont_force.x) * 0.02;
		cont_force.y -= ddist_sign * abs(cont_force.y) * 0.02;
	}
	cont_substp_id = substp_id + 1;
}
