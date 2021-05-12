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
		size_t ori_pcl_id,
		double dist,
		const Vector3D &norm,
		const Point3D &cont_pos,
		double pcl_len,
		ParticleVariablesGetter &pv_getter,
		// out
		Vector3D& cont_force
		) override;
};

#endif