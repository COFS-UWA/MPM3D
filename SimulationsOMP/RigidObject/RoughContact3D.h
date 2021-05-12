#ifndef __Rough_Contact_3D_h__
#define __Rough_Contact_3D_h__

#include "ContactModel3D.h"

class RoughContact3D : public ContactModel3D
{
protected:
	double Kn_cont, Kt_cont;

public:
	size_t* contact_substep_ids; // ori_pcl_num
	Position* prev_contact_poses; // ori_pcl_num
	Force* prev_contact_tan_forces; // ori_pcl_num
	double* prev_contact_dists; // ori_pcl_num

	RoughContact3D();
	~RoughContact3D();

	inline void set_K_cont(double _Kn, double _Kt) noexcept
	{ Kn_cont = _Kn; Kt_cont = _Kt; }

	inline double get_Kn_cont() const noexcept { return Kn_cont; }
	inline double get_Kt_cont() const noexcept { return Kt_cont; }

	void cal_contact_force(
		// in
		size_t substp_id,
		size_t ori_pcl_id,
		double dist,
		const Vector3D& norm,
		const Point3D& cont_pos,
		double pcl_len,
		ParticleVariablesGetter& pv_getter,
		// out
		Vector3D& cont_force
		) override;
};

#endif