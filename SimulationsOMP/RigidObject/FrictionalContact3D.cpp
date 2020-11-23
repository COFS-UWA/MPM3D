#include "SimulationsOMP_pcp.h"

#include "FrictionalContact3D.h"

FrictionalContact3D::FrictionalContact3D() :
	Kn_cont(0.0), Kt_cont(0.0), friction_ratio(0.0) {}

FrictionalContact3D::~FrictionalContact3D() {}

void FrictionalContact3D::cal_contact_force(
	// in
	size_t substp_id,
	double dist,
	const Vector3D& norm,
	const Point3D& cont_pos,
	ParticleVariablesGetter& pv_getter,
	// out
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
	// tangential force
	if (cont_substp_id != substp_id)
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
		double norm_len = rx * norm.x + ry * norm.y + rz * norm.z;
		rx -= norm_len * norm.x;
		ry -= norm_len * norm.y;
		rz -= norm_len * norm.z;
		prev_cont_tan_force.x += rx * Kt_cont;
		prev_cont_tan_force.y += ry * Kt_cont;
		prev_cont_tan_force.z += rz * Kt_cont;
		// contact consitutive model
		double tan_force = sqrt(prev_cont_tan_force.x * prev_cont_tan_force.x
							  + prev_cont_tan_force.y * prev_cont_tan_force.y
							  + prev_cont_tan_force.z * prev_cont_tan_force.z);
		double max_tan_force = sqrt(cont_force.x * cont_force.x
								  + cont_force.y * cont_force.y
								  + cont_force.z * cont_force.z) * friction_ratio;
		if (tan_force > max_tan_force)
		{
			double ratio = max_tan_force / tan_force;
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
