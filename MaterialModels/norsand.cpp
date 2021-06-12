#include "MaterialModels_pcp.h"

#include <math.h>

#include "MatModelConstants.h"
#include "StressInvariant.h"
#include "SymMatEigen.h"
#include "norsand.h"

void NorsandGlobal_set_param(
	NorsandGlobal &dat,
	__Float_Type__ phi,
	__Float_Type__ gamma,
	__Float_Type__ lambda,
	__Float_Type__ N,
	__Float_Type__ chi,
	__Float_Type__ H,
	__Float_Type__ Ig,
	__Float_Type__ niu)
{
	dat.phi = phi;
	dat.gamma = gamma;
	dat.lambda = lambda;
	dat.N = N;
	dat.chi = chi;
	dat.H = H;
	dat.Ig = Ig;
	dat.niu = niu;
	//
	const __Float_Type__ sin_phi = (__Float_Type__)sin(ToRadian(phi));
	dat.Mtc = ffmat(6.0) * sin_phi / (ffmat(3.0) - sin_phi);
	dat.N_chi_div_Mtc = N * chi / dat.Mtc;
}

void Norsand_set_NC_param(
	Norsand &dat,
	const NorsandGlobal &glb_dat,
	const __Float_Type__ stress[6],
	__Float_Type__ e)
{
	dat.s11 = stress[0];
	dat.s22 = stress[1];
	dat.s33 = stress[2];
	dat.s12 = stress[3];
	dat.s23 = stress[4];
	dat.s31 = stress[5];
	dat.e = e;
	union
	{
		__Float_Type__  invars[3];
		struct { __Float_Type__  p, q, lode_angle; };
	};
	cal_p_q_lode_angle(dat.stress, invars);
	__Float_Type__ state_param = dat.e - glb_dat.gamma + glb_dat.lambda * (__Float_Type__)log(-p);
	const __Float_Type__ M_coef = (ffmat(1.0) - glb_dat.N_chi_div_Mtc * (__Float_Type__)fabs(state_param));
	const __Float_Type__ M = glb_dat.Mtc
		* (ffmat(1.0) - glb_dat.Mtc / (ffmat(3.0) + glb_dat.Mtc) * (__Float_Type__)cos(ffmat(1.5) * lode_angle));
	const __Float_Type__ Mi = M * M_coef;
	dat.pi = -p / (__Float_Type__)exp(ffmat(1.0) + q / (Mi * p));
}

void Norsand_set_OC_param(
	Norsand &dat,
	const NorsandGlobal &glb_dat,
	const __Float_Type__ stress[6],
	__Float_Type__ e,
	__Float_Type__ OCR)
{

}

int32_t integrate_norsand(
	const NorsandGlobal& glb_dat,
	Norsand& mat_dat,
	const __Float_Type__ dstrain[6],
	__Float_Type__ destrain[6],
	__Float_Type__ dpstrain[6])
{
	// save previous state
	const __Float_Type__ ori_stress[6] = {
		mat_dat.s11, mat_dat.s22, mat_dat.s33,
		mat_dat.s12, mat_dat.s23, mat_dat.s31
	};
	const __Float_Type__ ori_e = mat_dat.e;

	// non-linear elasticity
	const double I1 = -(mat_dat.s11 + mat_dat.s22 + mat_dat.s33) / 3.0;
	const double G = ffmat(2.0) * glb_dat.Ig * I1; // 2G
	const double lbd = G * glb_dat.niu / (ffmat(1.0) - glb_dat.niu - glb_dat.niu);
	const double lbd_2G = lbd + G;
	mat_dat.s11 += lbd_2G * dstrain[0] + lbd * dstrain[1] + lbd * dstrain[2];
	mat_dat.s22 += lbd * dstrain[0] + lbd_2G * dstrain[1] + lbd * dstrain[2];
	mat_dat.s33 += lbd * dstrain[0] + lbd * dstrain[1] + lbd_2G * dstrain[2];
	mat_dat.s12 += G * dstrain[3];
	mat_dat.s23 += G * dstrain[4];
	mat_dat.s31 += G * dstrain[5];
	mat_dat.e += (ffmat(1.0) + mat_dat.e) * (dstrain[0] + dstrain[1] + dstrain[2]);
	destrain[0] = dstrain[0];
	destrain[1] = dstrain[1];
	destrain[2] = dstrain[2];
	destrain[3] = dstrain[3];
	destrain[4] = dstrain[4];
	destrain[5] = dstrain[5];
	dpstrain[0] = 0.0;
	dpstrain[1] = 0.0;
	dpstrain[2] = 0.0;
	dpstrain[3] = 0.0;
	dpstrain[4] = 0.0;
	dpstrain[5] = 0.0;

	union
	{
		__Float_Type__  invars[3];
		struct { __Float_Type__  p, q, lode_angle; };
	};
	// cal M and Mi
	cal_p_q_lode_angle(ori_stress, invars);
	__Float_Type__ state_param = ori_e - glb_dat.gamma + glb_dat.lambda * (__Float_Type__)log(-p);
	const __Float_Type__ M_coef = (ffmat(1.0) - glb_dat.N_chi_div_Mtc * (__Float_Type__)fabs(state_param));
	const __Float_Type__ M = glb_dat.Mtc
		* (ffmat(1.0) - glb_dat.Mtc / (ffmat(3.0) + glb_dat.Mtc) * (__Float_Type__)cos(ffmat(1.5) * lode_angle));
	const __Float_Type__ Mi = M * M_coef;

	cal_p_q(mat_dat.stress, invars);
	__Float_Type__ f = q + Mi * p * (ffmat(1.0) - (__Float_Type__)log(-p / mat_dat.pi));
	if (f <= 0.0) // elastic
		return 0;

	const __Float_Type__ Mi_tc = glb_dat.Mtc * M_coef;
	__Float_Type__ E_dg_ds[3], dep_cor[6];
	for (int32_t iter_id = 0; iter_id < 20; ++iter_id)
	{
		const __Float_Type__ df_dp = -Mi * (__Float_Type__)log(-p/mat_dat.pi) / ffmat(3.0) * q;
		const __Float_Type__ df_ds11 = df_dp + ffmat(0.5) * (mat_dat.s11 + mat_dat.s11 - mat_dat.s22 - mat_dat.s33);
		const __Float_Type__ df_ds22 = df_dp + ffmat(0.5) * (mat_dat.s22 + mat_dat.s22 - mat_dat.s33 - mat_dat.s11);
		const __Float_Type__ df_ds33 = df_dp + ffmat(0.5) * (mat_dat.s33 + mat_dat.s33 - mat_dat.s11 - mat_dat.s22);
		const __Float_Type__ df_ds12 = ffmat(3.0) * mat_dat.s12;
		const __Float_Type__ df_ds23 = ffmat(3.0) * mat_dat.s23;
		const __Float_Type__ df_ds31 = ffmat(3.0) * mat_dat.s31;
		const __Float_Type__ dg_dp = -(Mi * (-p) - q) / ffmat(3.0) * q;
		__Float_Type__ dg_ds11 = dg_dp + (-p) * ffmat(0.5) * (mat_dat.s11 + mat_dat.s11 - mat_dat.s22 - mat_dat.s33);
		__Float_Type__ dg_ds22 = dg_dp + (-p) * ffmat(0.5) * (mat_dat.s22 + mat_dat.s22 - mat_dat.s33 - mat_dat.s11);
		__Float_Type__ dg_ds33 = dg_dp + (-p) * ffmat(0.5) * (mat_dat.s33 + mat_dat.s33 - mat_dat.s11 - mat_dat.s22);
		__Float_Type__ dg_ds12 = (-p) * ffmat(3.0) * mat_dat.s12;
		__Float_Type__ dg_ds23 = (-p) * ffmat(3.0) * mat_dat.s23;
		__Float_Type__ dg_ds31 = (-p) * ffmat(3.0) * mat_dat.s31;
		const __Float_Type__ dg_ds11_ds22 = dg_ds11 - dg_ds22;
		const __Float_Type__ dg_ds22_ds33 = dg_ds22 - dg_ds33;
		const __Float_Type__ dg_ds33_ds11 = dg_ds33 - dg_ds11;

		//__Float_Type__ dg_len = sqrt(dg_ds11*dg_ds11 + dg_ds22*dg_ds22 + dg_ds33*dg_ds33
		//							+ dg_ds12*dg_ds12 + dg_ds23*dg_ds23 + dg_ds31*dg_ds31);
		//if (dg_len != 0.0)
		//{
		//	dg_ds11 /= dg_len;
		//	dg_ds22 /= dg_len;
		//	dg_ds33 /= dg_len;
		//	dg_ds12 /= dg_len;
		//	dg_ds23 /= dg_len;
		//	dg_ds31 /= dg_len;
		//}

		// A
		state_param = mat_dat.e - glb_dat.gamma + glb_dat.lambda * (__Float_Type__)log(-p);
		const __Float_Type__ pi_max = -p * (__Float_Type__)exp(-glb_dat.chi * state_param / Mi_tc);
		const __Float_Type__ A = glb_dat.H * Mi * (-p) / (Mi_tc * mat_dat.pi) * (pi_max - mat_dat.pi)
			* (__Float_Type__)sqrt((dg_ds11_ds22 * dg_ds11_ds22 + dg_ds22_ds33 * dg_ds22_ds33 + dg_ds33_ds11 * dg_ds33_ds11) * ffmat(0.5)
								 + (dg_ds12 * dg_ds12 + dg_ds23 * dg_ds23 + dg_ds31 * dg_ds31) * ffmat(3.0)) * ffmat(2.0) / ffmat(3.0);
		E_dg_ds[0] = lbd_2G * dg_ds11 + lbd * dg_ds22 + lbd * dg_ds33;
		E_dg_ds[1] = lbd * dg_ds11 + lbd_2G * dg_ds22 + lbd * dg_ds33;
		E_dg_ds[2] = lbd * dg_ds11 + lbd * dg_ds22 + lbd_2G * dg_ds33;
		const __Float_Type__ divider = df_ds11 * E_dg_ds[0] + df_ds22 * E_dg_ds[1] + df_ds33 * E_dg_ds[2]
			+ df_ds12 * G * dg_ds12 + df_ds23 * G * dg_ds23 + df_ds31 * G * dg_ds31
			- (Mi * p / mat_dat.pi) * q * A;
		const __Float_Type__ dl = f * q / divider;
		// strain correction
		dep_cor[0] = dl * dg_ds11;
		dep_cor[1] = dl * dg_ds22;
		dep_cor[2] = dl * dg_ds33;
		dep_cor[3] = dl * dg_ds12;
		dep_cor[4] = dl * dg_ds23;
		dep_cor[5] = dl * dg_ds31;
		// elastic strain increment
		destrain[0] -= dep_cor[0];
		destrain[1] -= dep_cor[1];
		destrain[2] -= dep_cor[2];
		destrain[3] -= dep_cor[3];
		destrain[4] -= dep_cor[4];
		destrain[5] -= dep_cor[5];
		// plastic strain increment
		dpstrain[0] += dep_cor[0];
		dpstrain[1] += dep_cor[1];
		dpstrain[2] += dep_cor[2];
		dpstrain[3] += dep_cor[3];
		dpstrain[4] += dep_cor[4];
		dpstrain[5] += dep_cor[5];
		// correct stress
		mat_dat.s11 -= lbd_2G * dep_cor[0] + lbd * dep_cor[1] + lbd * dep_cor[2];
		mat_dat.s22 -= lbd * dep_cor[0] + lbd_2G * dep_cor[1] + lbd * dep_cor[2];
		mat_dat.s33 -= lbd * dep_cor[0] + lbd * dep_cor[1] + lbd_2G * dep_cor[2];
		mat_dat.s12 -= G * dep_cor[3];
		mat_dat.s23 -= G * dep_cor[4];
		mat_dat.s31 -= G * dep_cor[5];
		mat_dat.pi += A * dl;

		cal_p_q(mat_dat.stress, invars);
		f = q + Mi * p * (ffmat(1.0) - (__Float_Type__)log(-p / mat_dat.pi));
		if (f < (Mi * (__Float_Type__)fabs(p)) * ffmat(1.0e-3)) // converge
			return iter_id + 1;
	}
	
	// integration failure
	mat_dat.s11 = ori_stress[0];
	mat_dat.s22 = ori_stress[1];
	mat_dat.s33 = ori_stress[2];
	mat_dat.s12 = ori_stress[3];
	mat_dat.s23 = ori_stress[4];
	mat_dat.s31 = ori_stress[5];
	mat_dat.e = ori_e;
	return -1;
}
