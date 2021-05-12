#ifndef __Contact_Model_3D_h__
#define __Contact_Model_3D_h__

#include "Geometry3D.h"
#include "ParticleVariablesGetter.h"

class ContactModel3D
{
public:
	struct Position { double x, y, z; };
	struct Force { double x, y, z; };

	explicit ContactModel3D();
	virtual ~ContactModel3D();

	virtual void cal_contact_force(
		// in
		size_t substp_id,
		size_t ori_pcl_id,
		double dist,
		const Vector3D& norm,
		const Point3D& cont_pos,
		double pcl_len,
		ParticleVariablesGetter &pv_getter,
		// out
		Vector3D &cont_force
		) = 0;

	inline static double sign(double num) noexcept { return num < 0.0 ? -1.0 : 1.0; }
};

#endif