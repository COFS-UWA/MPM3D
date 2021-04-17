#ifndef __Sticky_Contact_2D_h__
#define __Sticky_Contact_2D_h__

#include "ContactModel2D.h"

class StickyContact2D : public ContactModel2D
{
protected:
	double Kn_cont, Kt_cont;
	double shear_strength;

public:
	StickyContact2D();
	~StickyContact2D();

	inline void set_K_cont(double _Kn, double _Kt) noexcept
	{ Kn_cont = _Kn; Kt_cont = _Kt; }
	inline void set_shear_strength(double su) noexcept
	{ shear_strength = su; }
	inline double get_Kn_cont() const noexcept { return Kn_cont; }
	inline double get_Kt_cont() const noexcept { return Kt_cont; }
	inline double get_shear_strength() const noexcept { return shear_strength; }

	void cal_contact_force(
		// in
		size_t substp_id,
		double dist,
		const Vector2D& norm,
		const Point2D& cont_pos,
		double pcl_len,
		ParticleVariablesGetter& pv_getter,
		// out
		size_t& cont_substp_id,
		Point2D& prev_cont_pos,
		Vector2D& prev_cont_tan_force,
		Vector2D& cont_force
		) override;
};

#endif