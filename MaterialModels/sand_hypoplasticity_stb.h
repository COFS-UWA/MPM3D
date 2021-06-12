#ifndef __Sand_Hypoplasticity_Stb_h__
#define __Sand_Hypoplasticity_Stb_h__

#include "stdint.h"
#include "MaterialModelPrecision.h"

struct SandHypoplasticityStbGlobal;
struct SandHypoplasticityStb;
int32_t integrate_sand_hypoplasticity_stb(
	const SandHypoplasticityStbGlobal& glb_dat,
	SandHypoplasticityStb& mat_dat,
	const __Float_Type__ dstrain[6],
	__Float_Type__ substp_size_ratio);

struct SandHypoplasticityStbGlobal
{
	// params
	// hypoplasticity
	__Float_Type__ phi, hs, n, alpha, beta;
	__Float_Type__ ei0, ec0, ed0;
	// norsand
	__Float_Type__ N, chi, H;
	// elasticity
	__Float_Type__ Ig, niu;
	// tensile
	// ten_E << material stiffness
	__Float_Type__ ten_E, ten_niu;

	// cal vars
	__Float_Type__ a, fbfe_coef;
	__Float_Type__ Mtc, N_chi_div_Mtc;
	__Float_Type__ ten_G, ten_lambda, ten_lambda_2G;

	int32_t hypoplasticity_substp(
		const __Float_Type__ stress[6],
		__Float_Type__ e,
		const __Float_Type__ dstrain[6],
		__Float_Type__ dstress[6],
		__Float_Type__& de) const;
	void ten_elasticity(const __Float_Type__ dstrain[6],
		__Float_Type__ stress[6]) const;
};

void SandHypoplasticityStbGlobal_set_param(
	SandHypoplasticityStbGlobal& dat,
	__Float_Type__ phi,
	__Float_Type__ hs, __Float_Type__ n,
	__Float_Type__ alpha, __Float_Type__ beta,
	__Float_Type__ ed0, __Float_Type__ ec0, __Float_Type__ ei0,
	__Float_Type__ N, __Float_Type__ chi, __Float_Type__ H,
	__Float_Type__ Ig, __Float_Type__ niu,
	__Float_Type__ tE, __Float_Type__ tniu);

struct SandHypoplasticityStb
{
	union
	{
		__Float_Type__ stress[6];
		struct { __Float_Type__ s11, s22, s33, s12, s23, s31; };
	};
	__Float_Type__ e;
	__Float_Type__ substp_size;
	// loading surface
	__Float_Type__ pl;
	// yield surface
	__Float_Type__ Mi, pi;
};

void SandHypoplasticityStb_set_NC_param(
	SandHypoplasticityStb& dat,
	const SandHypoplasticityStbGlobal& glb_dat,
	const __Float_Type__ stress[6],
	__Float_Type__ e,
	__Float_Type__ substp_size);

void SandHypoplasticityStb_set_OC_param(
	SandHypoplasticityStb& dat,
	const SandHypoplasticityStbGlobal& glb_dat,
	const __Float_Type__ stress[6],
	__Float_Type__ e,
	__Float_Type__ OCR,
	__Float_Type__ substp_size);

#endif