#ifndef __Smooth_Contact_3D_h__
#define __Smooth_Contact_3D_h__

#include "ContactModel3D.h"

class SmoothContact3D : public ContactModel3D
{
protected:
	double Kn_cont;

public:
	SmoothContact3D();
	~SmoothContact3D();

	inline void set_Kn_cont(double _Kn) noexcept { Kn_cont = _Kn; }
	inline double get_Kn_cont() const noexcept { return Kn_cont; }

	void cal_contact_force(
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
		) override;
};

#endif