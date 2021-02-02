#ifndef __Rough_Contact_3D_h__
#define __Rough_Contact_3D_h__

#include "ContactModel3D.h"

class RoughContact3D : public ContactModel3D
{
protected:
	double Kn_cont, Kt_cont;

public:
	RoughContact3D();
	~RoughContact3D();

	inline void set_K_cont(double _Kn, double _Kt) noexcept
	{ Kn_cont = _Kn; Kt_cont = _Kt; }

	inline double get_Kn_cont() const noexcept { return Kn_cont; }
	inline double get_Kt_cont() const noexcept { return Kt_cont; }

	void cal_contact_force(
		// in
		size_t substp_id,
		double dist,
		const Vector3D& norm,
		const Point3D& cont_pos,
		double pcl_len,
		ParticleVariablesGetter& pv_getter,
		// out
		size_t& cont_substp_id,
		Point3D& prev_cont_pos,
		Vector3D& prev_cont_tan_force,
		Vector3D& cont_force
		) override;
};

#endif