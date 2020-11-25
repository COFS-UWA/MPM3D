#ifndef __Contact_Model_3D_h__
#define __Contact_Model_3D_h__

#include "Geometry3D.h"
#include "ParticleVariablesGetter.h"

class ContactModel3D
{
public:
	explicit ContactModel3D();
	virtual ~ContactModel3D();

	virtual void cal_contact_force(
		// in
		size_t substp_id,
		double dist,
		const Vector3D& norm,
		const Point3D& cont_pos,
		double pcl_len,
		ParticleVariablesGetter &pv_getter,
		// out
		size_t& cont_substp_id,
		Point3D &prev_cont_pos,
		Vector3D &prev_cont_tan_force,
		Vector3D &cont_force
		) = 0;
};

#endif