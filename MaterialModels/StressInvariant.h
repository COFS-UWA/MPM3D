#ifndef __Stress_Invariant_h__
#define __Stress_Invariant_h__

#include <math.h>
#include "MaterialModelPrecision.h"
#include "MatModelConstants.h"

// stress invariant I1, I2 and I3
inline void cal_I1_I2_I3(const __Float_Type__ stress[6], __Float_Type__ invariant[3])
{
	// I1
	invariant[0] = stress[0] + stress[1] + stress[2];
	// I2
	invariant[1] = stress[0] * stress[1]
				 + stress[1] * stress[2]
				 + stress[2] * stress[0]
				 - stress[3] * stress[3]
				 - stress[4] * stress[4]
				 - stress[5] * stress[5];
	// I3
	invariant[2] = stress[0] * stress[1] * stress[2]
				 + stress[3] * stress[4] * stress[5] * ffmat(2.0)
				 - stress[0] * stress[4] * stress[4] // sxx * syz * syz
				 - stress[1] * stress[5] * stress[5] // syy * szx * szx
				 - stress[2] * stress[3] * stress[3]; // szz * sxy * sxy
}

// stress invariant I1 and deviatoric stress invariant J2 and J3
inline void cal_I1_J2_J3(const __Float_Type__ stress[6], __Float_Type__ invariant[3])
{
	// I1
	invariant[0] = stress[0] + stress[1] + stress[2];
	// J2
	invariant[1] = ((stress[1] - stress[0]) * (stress[1] - stress[0])
		+ (stress[2] - stress[1]) * (stress[2] - stress[1])
		+ (stress[0] - stress[2]) * (stress[0] - stress[2])) / ffmat(6.0)
		+ stress[3] * stress[3]
		+ stress[4] * stress[4]
		+ stress[5] * stress[5];
	// J3
	const __Float_Type__ p = invariant[0] / ffmat(3.0);
	const __Float_Type__ sxx = stress[0] - p;
	const __Float_Type__ syy = stress[1] - p;
	const __Float_Type__ szz = stress[2] - p;
	invariant[2] = sxx * syy * szz
		+ stress[3] * stress[4] * stress[5] * ffmat(2.0)
		- sxx * stress[4] * stress[4]
		- syy * stress[5] * stress[5]
		- szz * stress[3] * stress[3];
}

// calculate invariant p and q
inline void cal_p_q(const __Float_Type__ stress[6], __Float_Type__ invariant[2])
{
	// p
	invariant[0] = (stress[0] + stress[1] + stress[2]) / ffmat(3.0);
	// q
	invariant[1] = (__Float_Type__)sqrt(
		 ((stress[1] - stress[0]) * (stress[1] - stress[0])
		+ (stress[2] - stress[1]) * (stress[2] - stress[1])
		+ (stress[0] - stress[2]) * (stress[0] - stress[2])) * ffmat(0.5)
		+ (stress[3] * stress[3]
		 + stress[4] * stress[4]
		 + stress[5] * stress[5]) * ffmat(3.0));
}

// calculate invariant p, q and lode_angle
// lode_angle in radian
inline void cal_p_q_lode_angle(const __Float_Type__ stress[6], __Float_Type__ invariant[3])
{
	// p
	const __Float_Type__ p = (stress[0] + stress[1] + stress[2]) / ffmat(3.0);
	invariant[0] = p;
	// q
	const __Float_Type__ q = (__Float_Type__)sqrt(
		((stress[1] - stress[0]) * (stress[1] - stress[0])
		+ (stress[2] - stress[1]) * (stress[2] - stress[1])
		+ (stress[0] - stress[2]) * (stress[0] - stress[2]))* ffmat(0.5)
		+ (stress[3] * stress[3]
		+ stress[4] * stress[4]
		+ stress[5] * stress[5]) * ffmat(3.0));
	invariant[1] = q;
	// J3
	const __Float_Type__ sxx = stress[0] - p;
	const __Float_Type__ syy = stress[1] - p;
	const __Float_Type__ szz = stress[2] - p;
	const __Float_Type__ J3 = sxx * syy * szz
		+ stress[3] * stress[4] * stress[5] * ffmat(2.0)
		- sxx * stress[4] * stress[4] // yz
		- syy * stress[5] * stress[5] // zx
		- szz * stress[3] * stress[3]; // xy
	// hydrostatic state taken as compression
	__Float_Type__ cos_3_theta;
	if (q < ffmat(1.0e-8))
		cos_3_theta = J3 > ffmat(0.0) ? 1.0 : -1.0;
	else
	{
		cos_3_theta = ffmat(27.0) / ffmat(2.0) * J3 / (q * q * q);
		if (cos_3_theta > ffmat(1.0))
			cos_3_theta = ffmat(1.0);
		if (cos_3_theta < ffmat(-1.0))
			cos_3_theta = ffmat(-1.0);
	}
	invariant[2] = (__Float_Type__)acos(cos_3_theta) / ffmat(3.0);
}

// principle stress s1, s2 and s3
// not so accurate as cal_sym_mat_eigen()
inline void cal_s1_s2_s3(const __Float_Type__ stress[6], __Float_Type__ ps[3])
{
	union
	{
		__Float_Type__ invars[3];
		struct { __Float_Type__ p, q, theta ; };
	};
	cal_p_q_lode_angle(stress, invars);
	ps[0] = p + ffmat(2.0) / ffmat(3.0) * q * cos(theta);
	ps[1] = p + ffmat(2.0) / ffmat(3.0) * q * cos(theta - ffmat(2.0) / ffmat(3.0) * PI);
	ps[2] = p + ffmat(2.0) / ffmat(3.0) * q * cos(theta + ffmat(2.0) / ffmat(3.0) * PI);
}

#endif