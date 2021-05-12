#include "SimulationsOMP_pcp.h"

#include "StickyContact3D.h"

StickyContact3D::StickyContact3D() : Kn_cont(0.0) {}

StickyContact3D::~StickyContact3D() {}

void StickyContact3D::cal_contact_force(
	size_t substp_id,
	size_t ori_pcl_id,
	double dist,
	const Vector3D& norm,
	const Point3D& cont_pos,
	double pcl_len,
	ParticleVariablesGetter& pv_getter,
	Vector3D& cont_force)
{
	// normal force
	const double pcl_area = pcl_len * pcl_len;
	double f_cont = Kn_cont * pcl_area * dist;
	cont_force.x = f_cont * norm.x;
	cont_force.y = f_cont * norm.y;
	cont_force.z = f_cont * norm.z;
	size_t& cont_substp_id = contact_substep_ids[ori_pcl_id];
	Position& prev_cont_pos = prev_contact_poses[ori_pcl_id];
	Force& prev_cont_tan_force = prev_contact_tan_forces[ori_pcl_id];
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
		cont_force.x -= ddist_sign * abs(cont_force.x) * 0.02;
		cont_force.y -= ddist_sign * abs(cont_force.y) * 0.02;
		cont_force.z -= ddist_sign * abs(cont_force.z) * 0.02;
	}

	if (cont_substp_id != substp_id) // not previously in contact
	{
		// not in contacct
		// record contact position
		prev_cont_pos.x = cont_pos.x;
		prev_cont_pos.y = cont_pos.y;
		prev_cont_pos.z = cont_pos.z;
		// rest contact force
		prev_cont_tan_force.x = 0.0;
		prev_cont_tan_force.y = 0.0;
		prev_cont_tan_force.z = 0.0;
	}
	else
	{
		// previously in contactt
		// cal relative movement
		double rx = cont_pos.x - prev_cont_pos.x;
		double ry = cont_pos.y - prev_cont_pos.y;
		double rz = cont_pos.z - prev_cont_pos.z;
		prev_cont_pos.x = cont_pos.x;
		prev_cont_pos.y = cont_pos.y;
		prev_cont_pos.z = cont_pos.z;
		const double norm_len = rx * norm.x + ry * norm.y + rz * norm.z;
		rx -= norm_len * norm.x;
		ry -= norm_len * norm.y;
		rz -= norm_len * norm.z;
		const double Kt_area = Kt_cont * pcl_area;
		prev_cont_tan_force.x -= rx * Kt_area;
		prev_cont_tan_force.y -= ry * Kt_area;
		prev_cont_tan_force.z -= rz * Kt_area;
		// contact consitutive model
		double tan_force = sqrt(prev_cont_tan_force.x * prev_cont_tan_force.x
							  + prev_cont_tan_force.y * prev_cont_tan_force.y
							  + prev_cont_tan_force.z * prev_cont_tan_force.z);
		const double max_shear_force = pcl_area * shear_strength;
		if (tan_force > max_shear_force)
		{
			const double ratio = max_shear_force / tan_force;
			prev_cont_tan_force.x *= ratio;
			prev_cont_tan_force.y *= ratio;
			prev_cont_tan_force.z *= ratio;
		}
		// add tangential force to contact force
		cont_force.x += prev_cont_tan_force.x;
		cont_force.y += prev_cont_tan_force.y;
		cont_force.z += prev_cont_tan_force.z;
	}
	cont_substp_id = substp_id + 1;
}
