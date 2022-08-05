#include "SimulationsOMP_pcp.h"

#include "SmoothContact3D.h"

SmoothContact3D::SmoothContact3D() : Kn_cont(0.0) {}

SmoothContact3D::~SmoothContact3D() {}

void SmoothContact3D::cal_contact_force(
	size_t substp_id,
	size_t ori_pcl_id,
	double dist,
	const Vector3D& norm,
	const Point3D& cont_pos,
	double pcl_len,
	ParticleVariablesGetter& pv_getter,
	Vector3D& cont_force)
{
	constexpr double K_damp_ratio = 0.2;
	// normal force
	//double f_cont = Kn_cont * pcl_len * pcl_len * dist;
	double f_cont = Kn_cont * pcl_len * pcl_len * dist * dist;
	cont_force.x = f_cont * norm.x;
	cont_force.y = f_cont * norm.y;
	cont_force.z = f_cont * norm.z;
	size_t& cont_substp_id = contact_substep_ids[ori_pcl_id];
	double& pcd = prev_contact_dists[ori_pcl_id];
	if (cont_substp_id != substp_id) // not previously in contact
	{
		pcd = 0.0;
	}
	else // previously in contact
	{
		// local damping
		const double ddist_sign = sign(dist - pcd);
		pcd = dist;
		cont_force.x -= ddist_sign * abs(cont_force.x) * K_damp_ratio;
		cont_force.y -= ddist_sign * abs(cont_force.y) * K_damp_ratio;
		cont_force.z -= ddist_sign * abs(cont_force.z) * K_damp_ratio;
	}
	cont_substp_id = substp_id + 1;
}
