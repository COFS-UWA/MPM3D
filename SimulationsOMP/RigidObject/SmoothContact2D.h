#ifndef __Smooth_Contact_2D_h__
#define __Smooth_Contact_2D_h__

#include "ContactModel2D.h"

class SmoothContact2D : public ContactModel2D
{
protected:
	double Kn_cont;

public:
	SmoothContact2D();
	~SmoothContact2D();

	inline void set_Kn_cont(double _Kn) noexcept { Kn_cont = _Kn; }
	inline double get_Kn_cont() const noexcept { return Kn_cont; }

	void cal_contact_force(
		// in
		size_t substp_id,
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