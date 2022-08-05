#include "MaterialModels_pcp.h"
#include <math.h>

#include "MatModelConstants.h"
#include "StressInvariant.h"
#include "SymMatEigen.h"
#include "norsand.h"

void NorsandGlobal_set_param(
	NorsandGlobal& dat,
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
	dat.Mtc_div_3_plus_Mtc = dat.Mtc / (ffmat(3.0) + dat.Mtc);
}

void Norsand_set_NC_param(
	Norsand& dat,
	const NorsandGlobal& glb_dat,
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
	const __Float_Type__ M_lode_coef = ffmat(1.0) - glb_dat.Mtc_div_3_plus_Mtc * (__Float_Type__)cos(ffmat(1.5) * lode_angle);
	__Float_Type__ state_param = dat.e - glb_dat.gamma + glb_dat.lambda * (__Float_Type__)log(-p);
	const __Float_Type__ M_i_coef = ffmat(1.0) - glb_dat.N * (__Float_Type__)fabs(state_param);
	const __Float_Type__ Mi = glb_dat.Mtc * M_lode_coef * M_i_coef;
	dat.pi = -p / (__Float_Type__)exp(ffmat(1.0) + q / (Mi * p));
}

void Norsand_set_OC_param(
	Norsand& dat,
	const NorsandGlobal& glb_dat,
	const __Float_Type__ stress[6],
	__Float_Type__ e,
	__Float_Type__ OCR)
{}

inline double sign_(double val)
{
	if (val == 0.0)
		return 0.0;
	else if (val > 0.0)
		return 1.0;
	else
		return -1.0;
}

int32_t integrate_norsand(
	const NorsandGlobal& glb_dat,
	Norsand& mat_dat,
	const __Float_Type__ dstrain[6],
	__Float_Type__ destrain[6],
	__Float_Type__ dpstrain[6])
{
	// store previous state
	const __Float_Type__ ori_stress[6] = {
		mat_dat.s11, mat_dat.s22, mat_dat.s33,
		mat_dat.s12, mat_dat.s23, mat_dat.s31
	};
	const __Float_Type__ ori_e = mat_dat.e;

	// non-linear elasticity
	const double I1 = -(mat_dat.s11 + mat_dat.s22 + mat_dat.s33) / ffmat(3.0);
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
	const __Float_Type__ M_lode_coef = ffmat(1.0) - glb_dat.Mtc_div_3_plus_Mtc * (__Float_Type__)cos(ffmat(1.5) * lode_angle);
	__Float_Type__ state_param = ori_e - glb_dat.gamma + glb_dat.lambda * (__Float_Type__)log(-p);
	const __Float_Type__ M_i_coef = ffmat(1.0) - glb_dat.N * (__Float_Type__)fabs(state_param);
	const __Float_Type__ Mi = glb_dat.Mtc * M_lode_coef * M_i_coef;

	cal_p_q(mat_dat.stress, invars);
	__Float_Type__ f = q + Mi * p * (ffmat(1.0) - (__Float_Type__)log(-p / mat_dat.pi));
	if (f <= ffmat(0.0)) // elastic
		return 0;

	const __Float_Type__ Mi_tc = glb_dat.Mtc * M_i_coef;
	for (int32_t iter_id = 0; iter_id < 20; ++iter_id)
	{
		// df_ds
		const __Float_Type__ df_dp = -Mi * (__Float_Type__)log(-p / mat_dat.pi) / ffmat(3.0);
		const __Float_Type__ dg_dp = -(Mi * (-p) - q) / ffmat(3.0);
		const __Float_Type__ dg_dq = (-p);
		__Float_Type__ dq_ds11, dq_ds22, dq_ds33, dq_ds12, dq_ds23, dq_ds31;
		if (q > ffmat(1.0e-5))
		{
			const __Float_Type__ inv_q = ffmat(1.0) / q;
			dq_ds11 = ffmat(0.5) * (mat_dat.s11 + mat_dat.s11 - mat_dat.s22 - mat_dat.s33) * inv_q;
			dq_ds22 = ffmat(0.5) * (mat_dat.s22 + mat_dat.s22 - mat_dat.s33 - mat_dat.s11) * inv_q;
			dq_ds33 = ffmat(0.5) * (mat_dat.s33 + mat_dat.s33 - mat_dat.s11 - mat_dat.s22) * inv_q;
			dq_ds12 = ffmat(3.0) * mat_dat.s12 * inv_q;
			dq_ds23 = ffmat(3.0) * mat_dat.s23 * inv_q;
			dq_ds31 = ffmat(3.0) * mat_dat.s31 * inv_q;
		}
		else
		{
			dq_ds11 = sign_(mat_dat.s11 + mat_dat.s11 - mat_dat.s22 - mat_dat.s33);
			dq_ds22 = sign_(mat_dat.s22 + mat_dat.s22 - mat_dat.s33 - mat_dat.s11);
			dq_ds33 = sign_(mat_dat.s33 + mat_dat.s33 - mat_dat.s11 - mat_dat.s22);
			dq_ds12 = 1.732050808 * sign_(mat_dat.s12);
			dq_ds23 = 1.732050808 * sign_(mat_dat.s23);
			dq_ds31 = 1.732050808 * sign_(mat_dat.s31);
		}
		// df_ds
		const __Float_Type__ df_ds11 = df_dp + dq_ds11;
		const __Float_Type__ df_ds22 = df_dp + dq_ds22;
		const __Float_Type__ df_ds33 = df_dp + dq_ds33;
		const __Float_Type__ df_ds12 = dq_ds12;
		const __Float_Type__ df_ds23 = dq_ds23;
		const __Float_Type__ df_ds31 = dq_ds31;
		// dg_ds
		const __Float_Type__ dg_ds11 = dg_dp + dg_dq * dq_ds11;
		const __Float_Type__ dg_ds22 = dg_dp + dg_dq * dq_ds22;
		const __Float_Type__ dg_ds33 = dg_dp + dg_dq * dq_ds33;
		const __Float_Type__ dg_ds12 = dg_dq * dq_ds12;
		const __Float_Type__ dg_ds23 = dg_dq * dq_ds23;
		const __Float_Type__ dg_ds31 = dg_dq * dq_ds31;
		// A
		state_param = mat_dat.e - glb_dat.gamma + glb_dat.lambda * (__Float_Type__)log(-p);
		const __Float_Type__ pi_max = -p * (__Float_Type__)exp(-glb_dat.chi * state_param);
		const __Float_Type__ A = glb_dat.H * M_lode_coef * (pi_max - mat_dat.pi) * dg_dq;
		const __Float_Type__ E_dg_ds[3] = {
			lbd_2G * dg_ds11 + lbd * dg_ds22 + lbd * dg_ds33,
			lbd * dg_ds11 + lbd_2G * dg_ds22 + lbd * dg_ds33,
			lbd * dg_ds11 + lbd * dg_ds22 + lbd_2G * dg_ds33
		};
		// dl
		const __Float_Type__ divider = df_ds11 * E_dg_ds[0] + df_ds22 * E_dg_ds[1] + df_ds33 * E_dg_ds[2]
			+ df_ds12 * G * dg_ds12 + df_ds23 * G * dg_ds23 + df_ds31 * G * dg_ds31
			- (Mi * p / mat_dat.pi) * A;
		const __Float_Type__ dl = f / divider;

		// strain correction
		const __Float_Type__ dep_cor[6] = {
			dl * dg_ds11, dl * dg_ds22, dl * dg_ds33,
			dl * dg_ds12, dl * dg_ds23, dl * dg_ds31
		};
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
	return -1;
}
