#ifndef __Non_Seperate_Contact_3D_h__
#define __Non_Seperate_Contact_3D_h__

#include "ContactModel3D.h"

// non seperate in normal direction
class NonSeperateContact3D : public ContactModel3D
{
protected:
	double Kn_cont;
	double dist_off;

public:
	size_t* contact_substep_ids; // ori_pcl_num
	double* prev_contact_dists; // ori_pcl_num

	NonSeperateContact3D();
	~NonSeperateContact3D();

	inline void set_K_cont(double _Kn) noexcept { Kn_cont = _Kn; }
	inline void set_overlap_dist_off(double _off) noexcept { dist_off = _off; }
	inline double get_Kn_cont() const noexcept { return Kn_cont; }
	inline double get_overlap_dist_off() noexcept { return dist_off; }

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