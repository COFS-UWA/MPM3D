#include "SimulationsOMP_pcp.h"

#include "NonSeperateContact3D.h"

NonSeperateContact3D::NonSeperateContact3D() : Kn_cont(0.0), dist_off(0.0) {}

NonSeperateContact3D::~NonSeperateContact3D() {}

void NonSeperateContact3D::cal_contact_force(
	size_t substp_id,
	size_t ori_pcl_id,
	double dist,
	const Vector3D& norm,
	const Point3D& cont_pos,
	double pcl_len,
	ParticleVariablesGetter& pv_getter,
	Vector3D& cont_force)
{
	// normal contact force
	const double f_cont_dist = dist - dist_off;
	const double f_cont = Kn_cont * pcl_len * pcl_len * f_cont_dist;
	cont_force.x = f_cont * norm.x;
	cont_force.y = f_cont * norm.y;
	cont_force.z = f_cont * norm.z;
	// add local damping to contact force
	size_t &cont_substp_id = contact_substep_ids[ori_pcl_id];
	double& pcd = prev_contact_dists[ori_pcl_id];
	if (cont_substp_id != substp_id) // not previously in contact
	{
		pcd = -dist_off; // not damping
	}
	else // previously in contact
	{
		constexpr double Kn_damp_ratio = 0.02;
		const double ddist_sign = sign(f_cont_dist - pcd);
		pcd = f_cont_dist;
		cont_force.x -= ddist_sign * abs(cont_force.x) * Kn_damp_ratio;
		cont_force.y -= ddist_sign * abs(cont_force.y) * Kn_damp_ratio;
		cont_force.z -= ddist_sign * abs(cont_force.z) * Kn_damp_ratio;
	}

	cont_substp_id = substp_id + 1;
}
