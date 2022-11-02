#ifndef __Smooth_Contact_3D_Quad_h__
#define __Smooth_Contact_3D_Quad_h__

#include "ContactModel3D.h"

class SmoothContact3D_Quad : public ContactModel3D
{
protected:
	double Kn_cont;

public:
	size_t* contact_substep_ids; // ori_pcl_num
	double* prev_contact_dists; // ori_pcl_num

	SmoothContact3D_Quad();
	~SmoothContact3D_Quad();

	inline void set_Kn_cont(double _Kn) noexcept { Kn_cont = _Kn; }
	inline double get_Kn_cont() const noexcept { return Kn_cont; }

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