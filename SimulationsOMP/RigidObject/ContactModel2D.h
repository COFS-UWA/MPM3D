#ifndef __Contact_Model_2D_h__
#define __Contact_Model_2D_h__

#include "Geometry2D.h"
#include "ParticleVariablesGetter.h"

class ContactModel2D
{
public:
	explicit ContactModel2D();
	virtual ~ContactModel2D();

	virtual void cal_contact_force(
		// in
		size_t substp_id,
		double dist,
		const Vector2D& norm,
		const Point2D& cont_pos,
		double pcl_len,
		ParticleVariablesGetter &pv_getter,
		// out
		size_t& cont_substp_id,
		Point2D &prev_cont_pos,
		Vector2D &prev_cont_tan_force,
		Vector2D &cont_force
		) = 0;
};

#endif