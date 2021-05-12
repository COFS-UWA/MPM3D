#include "SimulationsOMP_pcp.h"

#include "RoughContact3D.h"

RoughContact3D::RoughContact3D() :
	Kn_cont(0.0), Kt_cont(0.0) {}

RoughContact3D::~RoughContact3D() {}

void RoughContact3D::cal_contact_force(
	size_t substp_id,
	size_t ori_pcl_id,
	double dist,
	const Vector3D& norm,
	const Point3D& cont_pos,
	double pcl_len,
	ParticleVariablesGetter& pv_getter,
	Vector3D& cont_force)
{
	constexpr double Kn_damp_ratio = 0.02;
	constexpr double Kt_damp_ratio = 0.02;
	const double pcl_area = pcl_len * pcl_len;
	// normal force
	const double f_cont = Kn_cont * pcl_area * dist;
	cont_force.x = f_cont * norm.x;
	cont_force.y = f_cont * norm.y;
	cont_force.z = f_cont * norm.z;
	size_t &cont_substp_id = contact_substep_ids[ori_pcl_id];
	Position &prev_cont_pos = prev_contact_poses[ori_pcl_id];
	Force & prev_cont_tan_force = prev_contact_tan_forces[ori_pcl_id];
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
		cont_force.x -= ddist_sign * abs(cont_force.x) * Kn_damp_ratio;
		cont_force.y -= ddist_sign * abs(cont_force.y) * Kn_damp_ratio;
		cont_force.z -= ddist_sign * abs(cont_force.z) * Kn_damp_ratio;
	}

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
		const double Kt_area = Kt_cont * pcl_area;
		double rx = prev_cont_pos.x - cont_pos.x;
		double ry = prev_cont_pos.y - cont_pos.y;
		double rz = prev_cont_pos.z - cont_pos.z;
		const double norm_len = rx * norm.x + ry * norm.y + rz * norm.z;
		rx -= norm_len * norm.x;
		ry -= norm_len * norm.y;
		rz -= norm_len * norm.z;
		const double ctfx = rx * Kt_area;
		const double ctfy = ry * Kt_area;
		const double ctfz = rz * Kt_area;
		// local damping
		const double dctfx_sign = sign(ctfx - prev_cont_tan_force.x);
		const double dctfy_sign = sign(ctfy - prev_cont_tan_force.y);
		const double dctfz_sign = sign(ctfz - prev_cont_tan_force.z);
		prev_cont_tan_force.x = ctfx;
		prev_cont_tan_force.y = ctfy;
		prev_cont_tan_force.z = ctfz;
		// add local damping to force
		cont_force.x += ctfx + dctfx_sign * abs(ctfx) * Kt_damp_ratio;
		cont_force.y += ctfy + dctfy_sign * abs(ctfy) * Kt_damp_ratio;
		cont_force.z += ctfz + dctfz_sign * abs(ctfz) * Kt_damp_ratio;
	}
	cont_substp_id = substp_id + 1;
}
