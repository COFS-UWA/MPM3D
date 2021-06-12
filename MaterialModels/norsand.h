#ifndef __Norsand_h__
#define __Norsand_h__

#include "stdint.h"
#include "MaterialModelPrecision.h"

struct NorsandGlobal
{
	// friction angle
	__Float_Type__ phi;
	// critical state
	__Float_Type__ gamma, lambda;
	// 
	__Float_Type__ N, chi;
	// hardening
	__Float_Type__ H;
	// non-linear elasticity
	__Float_Type__ Ig, niu;

	// cal var
	__Float_Type__ Mtc;
	__Float_Type__ N_chi_div_Mtc;
};

void NorsandGlobal_set_param(
	NorsandGlobal &dat,
	__Float_Type__ phi,
	__Float_Type__ gamma,
	__Float_Type__ lambda,
	__Float_Type__ N,
	__Float_Type__ chi,
	__Float_Type__ H,
	__Float_Type__ Ig,
	__Float_Type__ niu);

struct Norsand
{
	union
	{
		__Float_Type__ stress[6];
		struct { __Float_Type__ s11, s22, s33, s12, s23, s31; };
	};
	__Float_Type__ e;
	__Float_Type__ pi;
};

void Norsand_set_NC_param(
	Norsand &dat,
	const NorsandGlobal &glb_dat,
	const __Float_Type__ stress[6],
	__Float_Type__ e);

void Norsand_set_OC_param(
	Norsand &dat,
	const NorsandGlobal &glb_dat,
	const __Float_Type__ stress[6],
	__Float_Type__ e,
	__Float_Type__ OCR);

int32_t integrate_norsand(
	const NorsandGlobal &glb_dat,
	Norsand &mat_dat,
	const __Float_Type__ dstrain[6],
	__Float_Type__ destrain[6],
	__Float_Type__ dpstrain[6]);

#endif