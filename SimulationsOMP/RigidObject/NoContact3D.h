#ifndef __No_Contact_3D_h__
#define __No_Contact_3D_h__

#include "ContactModel3D.h"

class NoContact3D : public ContactModel3D
{
public:
	NoContact3D();
	~NoContact3D();

	void cal_contact_force(
		// in
		size_t substp_id,
		double dist,
		const Vector3D &norm,
		const Point3D &cont_pos,
		double pcl_len,
		ParticleVariablesGetter &pv_getter,
		// out
		size_t& cont_substp_id,
		Point3D& prev_cont_pos,
		Vector3D& prev_cont_tan_force,
		Vector3D& cont_force
		) override;
};

#endif