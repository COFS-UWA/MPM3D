#ifndef __Mohr_Coulomb_h__
#define __Mohr_Coulomb_h__

#include "stdint.h"
#include "MaterialModelPrecision.h"

struct MohrCoulombGlobal
{
	// params
	__Float_Type__ phi; // fricion angle
	__Float_Type__ psi; // dilation angle
	__Float_Type__ cohesion;
	__Float_Type__ E, niu;

	// cal vars
	__Float_Type__ G, lambda, lambda_2G;
	__Float_Type__ inv_E, niu_div_E, inv_G;

	__Float_Type__ half_one_minus_sin_phi;
	__Float_Type__ half_one_plus_sin_phi;
	__Float_Type__ c_cos_phi;
	__Float_Type__ c_div_tan_phi;
	__Float_Type__ one_plus_sin_phi_psi;
	__Float_Type__ one_minus_sin_psi;
	__Float_Type__ one_plus_sin_psi;

	// each row is the axises
	union
	{
		__Float_Type__ i1[2][3];
		struct { __Float_Type__ v1x, v1y, v1z; };
	};
	union
	{
		__Float_Type__ i2[2][3];
		struct { __Float_Type__ v2x, v2y, v2z; };
	};
	__Float_Type__ v1x2_v1y2_v1z2, v2x2_v2y2_v2z2;

	void set_param(__Float_Type__ _phi,
		__Float_Type__ _psi, __Float_Type__ _cohesion,
		__Float_Type__ _E, __Float_Type__ _niu);
};

struct MohrCoulomb
{
	union
	{
		__Float_Type__ stress[6];
		struct { __Float_Type__ s11, s22, s33, s12, s23, s31; };
	};
};

int32_t integrate_mohr_coulomb(
	const MohrCoulombGlobal& glb_dat,
	MohrCoulomb& mat_dat,
	const __Float_Type__ dstrain[6],
	__Float_Type__ destrain[6],
	__Float_Type__ dpstrain[6]);

#endif