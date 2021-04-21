#ifndef __No_Contact_2D_h__
#define __No_Contact_2D_h__

#include "ContactModel2D.h"

class NoContact2D : public ContactModel2D
{
public:
	NoContact2D();
	~NoContact2D();

	void cal_contact_force(
		// in
		size_t substp_id,
		size_t ori_pcl_id,
		double dist,
		const Vector2D &norm,
		const Point2D &cont_pos,
		double pcl_len,
		ParticleVariablesGetter &pv_getter,
		// out
		size_t& cont_substp_id,
		Point2D& prev_cont_pos,
		Vector2D& prev_cont_tan_force,
		Vector2D& cont_force
		) override;
};

#endif