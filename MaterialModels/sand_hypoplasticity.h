#ifndef __Sand_Hypoplasticity_h__
#define __Sand_Hypoplasticity_h__

#include "stdint.h"
#include "MaterialModelPrecision.h"

struct SandHypoplasticityGlobal;
struct SandHypoplasticity;
int32_t integrate_sand_hypoplasticity(
	const SandHypoplasticityGlobal &glb_dat,
	SandHypoplasticity &mat_dat,
	const __Float_Type__ dstrain[6],
	__Float_Type__ substp_size_ratio);

// ten_E << material stiffness
struct SandHypoplasticityGlobal
{
	// params
	__Float_Type__ phi, hs, n, alpha, beta;
	__Float_Type__ ei0, ec0, ed0;
	__Float_Type__ ten_E, ten_niu;

	// cal vars
	__Float_Type__ a, fbfe_coef;
	__Float_Type__ ten_G, ten_lambda, ten_lambda_2G;

	void set_param(__Float_Type__ _phi,
		__Float_Type__ _hs, __Float_Type__ _n,
		__Float_Type__ _alpha, __Float_Type__ _beta,
		__Float_Type__ _ed0, __Float_Type__ _ec0, __Float_Type__ _ei0, 
		__Float_Type__ _tE, __Float_Type__ _tniu);

protected:
	friend int32_t integrate_sand_hypoplasticity(
		const SandHypoplasticityGlobal& glb_dat,
		SandHypoplasticity& mat_dat,
		const __Float_Type__ dstrain[6],
		__Float_Type__ substp_size_ratio);

	int32_t hypoplasticity_substp(
		const __Float_Type__ stress[6],
		__Float_Type__ e,
		const __Float_Type__ dstrain[6],
		__Float_Type__ dstress[6],
		__Float_Type__& de) const;
	void elasticity(const __Float_Type__ dstrain[6],
		__Float_Type__ stress[6]) const;
};

struct SandHypoplasticity
{
	union
	{
		__Float_Type__ stress[6];
		struct { __Float_Type__ s11, s22, s33, s12, s23, s31; };
	};
	__Float_Type__ e;
	__Float_Type__ substp_size;
};

#endif