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
	dat.Mtc_div_3_plus_Mtc = dat.Mtc / (ffmat(3.0) + dat.Mtc);
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
	const __Float_Type__ M_lode_coef = ffmat(1.0) - glb_dat.Mtc_div_3_plus_Mtc * (__Float_Type__)cos(ffmat(1.5) * lode_angle);
	__Float_Type__ state_param = dat.e - glb_dat.gamma + glb_dat.lambda * (__Float_Type__)log(-p);
	const __Float_Type__ M_i_coef = ffmat(1.0) - glb_dat.N * (__Float_Type__)fabs(state_param);
	const __Float_Type__ Mi = glb_dat.Mtc * M_lode_coef * M_i_coef;
	dat.pi = -p / (__Float_Type__)exp(ffmat(1.0) + q / (Mi * p));
}

void Norsand_set_OC_param(
	Norsand &dat,
	const NorsandGlobal &glb_dat,
	const __Float_Type__ stress[6],
	__Float_Type__ e,
	__Float_Type__ OCR)
{}

inline __Float_Type__ sign_(__Float_Type__ var)
{
	if (var > ffmat(0.0))
		return ffmat(1.0);
	else if (var < ffmat(0.0))
		return ffmat(-1.0);
	return ffmat(0.0);
}

inline void cal_dq_ds(
	__Float_Type__ q,
	const __Float_Type__ s[6],
	__Float_Type__ dq_ds[6])
{
	if (q < ffmat(1.0e-5))
	{
		dq_ds[0] = ffmat(0.0);
		dq_ds[1] = ffmat(0.0);
		dq_ds[2] = ffmat(0.0);
		dq_ds[3] = ffmat(0.0);
		dq_ds[4] = ffmat(0.0);
		dq_ds[5] = ffmat(0.0);
		return;
	}
	const __Float_Type__ inv_q = ffmat(1.0) / q;
	dq_ds[0] = ffmat(0.5) * (s[0] + s[0] - s[1] - s[2]) * inv_q;
	dq_ds[1] = ffmat(0.5) * (s[1] + s[1] - s[2] - s[0]) * inv_q;
	dq_ds[2] = ffmat(0.5) * (s[2] + s[2] - s[0] - s[1]) * inv_q;
	dq_ds[3] = ffmat(3.0) * s[3] * inv_q;
	dq_ds[4] = ffmat(3.0) * s[4] * inv_q;
	dq_ds[5] = ffmat(3.0) * s[5] * inv_q;
}

#define p_min ffmat(10.0)
#define pi_min ffmat(500.0)

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
	
	union
	{
		__Float_Type__  invars[3];
		struct { __Float_Type__  p, q, lode_angle; };
	};
	union
	{
		struct { __Float_Type__ dq_ds11, dq_ds22, dq_ds33, dq_ds12, dq_ds23, dq_ds31; };
		__Float_Type__ dq_ds[6];
	};
	__Float_Type__ df_dp, df_dq, dg_dp, dg_dq, pi_max, A, dl;
	__Float_Type__ df_ds11, df_ds22, df_ds33, df_ds12, df_ds23, df_ds31;
	__Float_Type__ dg_ds11, dg_ds22, dg_ds33, dg_ds12, dg_ds23, dg_ds31;
	__Float_Type__ dep_cor[6], stress_tmp[6], p_tmp, q_ratio, f;

	// cal M and Mi
	cal_p_q_lode_angle(ori_stress, invars);
	const __Float_Type__ M_lode_coef = ffmat(1.0) - glb_dat.Mtc_div_3_plus_Mtc * (__Float_Type__)cos(ffmat(1.5) * lode_angle);
	__Float_Type__ state_param = ori_e - glb_dat.gamma + glb_dat.lambda * (__Float_Type__)log(-p);
	const __Float_Type__ M_i_coef = ffmat(1.0) - glb_dat.N * (__Float_Type__)fabs(state_param);
	const __Float_Type__ Mi = glb_dat.Mtc * M_lode_coef * M_i_coef;
	
	// non-linear elasticity
	p_tmp = (mat_dat.s11 + mat_dat.s22 + mat_dat.s33) / ffmat(3.0);
	if (p_tmp > -p_min)
		p_tmp = -p_min;
	const double I1 = -p_tmp;
	const double G = ffmat(2.0) * glb_dat.Ig * I1; // 2G
	const double lbd = G * glb_dat.niu / (ffmat(1.0) - glb_dat.niu - glb_dat.niu);
	const double lbd_2G = lbd + G;
	const double K = (ffmat(1.0) + glb_dat.niu) * G / (ffmat(3.0) * (ffmat(1.0) - glb_dat.niu - glb_dat.niu)); // bulk modulus
	const double inv_E = ffmat(1.0) / ((ffmat(1.0) + glb_dat.niu) * G);
	const double niu_inv_E = glb_dat.niu * inv_E;
	const double inv_2G = ffmat(1.0) / G;
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
	dpstrain[0] = ffmat(0.0);
	dpstrain[1] = ffmat(0.0);
	dpstrain[2] = ffmat(0.0);
	dpstrain[3] = ffmat(0.0);
	dpstrain[4] = ffmat(0.0);
	dpstrain[5] = ffmat(0.0);
	
	cal_p_q(mat_dat.stress, invars);
	if (p > -p_min)
	{
		q_ratio = (Mi * p_min * (ffmat(1.0) - (__Float_Type__)log(p_min / mat_dat.pi))) / q;
		if (q_ratio > ffmat(1.0))
			q_ratio = ffmat(1.0);
		mat_dat.s11 = (mat_dat.s11 - p) * q_ratio + (-p_min);
		mat_dat.s22 = (mat_dat.s22 - p) * q_ratio + (-p_min);
		mat_dat.s33 = (mat_dat.s33 - p) * q_ratio + (-p_min);
		mat_dat.s12 = mat_dat.s12 * q_ratio;
		mat_dat.s23 = mat_dat.s23 * q_ratio;
		mat_dat.s31 = mat_dat.s31 * q_ratio;
		stress_tmp[0] = mat_dat.s11 - ori_stress[0];
		stress_tmp[1] = mat_dat.s22 - ori_stress[1];
		stress_tmp[2] = mat_dat.s33 - ori_stress[2];
		stress_tmp[3] = mat_dat.s12 - ori_stress[3];
		stress_tmp[4] = mat_dat.s23 - ori_stress[4];
		stress_tmp[5] = mat_dat.s31 - ori_stress[5];
		destrain[0] = inv_E * stress_tmp[0] - niu_inv_E * stress_tmp[1] - niu_inv_E * stress_tmp[2];
		destrain[1] = -niu_inv_E * stress_tmp[0] + inv_E * stress_tmp[1] - niu_inv_E * stress_tmp[2];
		destrain[2] = -niu_inv_E * stress_tmp[0] - niu_inv_E * stress_tmp[1] + inv_E * stress_tmp[2];
		destrain[3] = inv_2G * stress_tmp[3];
		destrain[4] = inv_2G * stress_tmp[4];
		destrain[5] = inv_2G * stress_tmp[5];
		dpstrain[0] = dstrain[0] - destrain[0];
		dpstrain[1] = dstrain[1] - destrain[1];
		dpstrain[2] = dstrain[2] - destrain[2];
		dpstrain[3] = dstrain[3] - destrain[3];
		dpstrain[4] = dstrain[4] - destrain[4];
		dpstrain[5] = dstrain[5] - destrain[5];
		return 0;
	}
	
	f = q + Mi * p * (ffmat(1.0) - (__Float_Type__)log(-p / mat_dat.pi));
	if (f <= ffmat(0.0)) // elastic or p < p0
		return 0;

	for (int32_t iter_id = 0; iter_id < 40; ++iter_id)
	{
		if ((-p) > mat_dat.pi * ffmat(2.71828182845904) &&
			(q * Mi) < ((-p) - mat_dat.pi * ffmat(2.71828182845904)))
		{
			// df_ds
			df_dp = p + mat_dat.pi * ffmat(2.71828182845904);
			df_dq = q;
			const __Float_Type__ df_ds_len = (__Float_Type__)sqrt(df_dp * df_dp + df_dq * df_dq);
			// f
			f = q + Mi * p * (ffmat(1.0) - (__Float_Type__)log(-p / mat_dat.pi));
			// df_ds
			df_dp *= f / df_ds_len;
			df_dq *= f / df_ds_len;
			cal_dq_ds(q, mat_dat.stress, dq_ds);
			df_ds11 = df_dp / ffmat(3.0) + df_dq * dq_ds11;
			df_ds22 = df_dp / ffmat(3.0) + df_dq * dq_ds22;
			df_ds33 = df_dp / ffmat(3.0) + df_dq * dq_ds33;
			df_ds12 = df_dq * dq_ds12;
			df_ds23 = df_dq * dq_ds23;
			df_ds31 = df_dq * dq_ds31;
		}
		else
		{
			// df_ds
			df_dp = -Mi * (__Float_Type__)log(-p / mat_dat.pi);
			cal_dq_ds(q, mat_dat.stress, dq_ds);
			df_ds11 = df_dp / ffmat(3.0) + dq_ds11;
			df_ds22 = df_dp / ffmat(3.0) + dq_ds22;
			df_ds33 = df_dp / ffmat(3.0) + dq_ds33;
			df_ds12 = dq_ds12;
			df_ds23 = dq_ds23;
			df_ds31 = dq_ds31;
		}

		// dg_ds
		dg_dp = -(Mi * (-p) - q);
		dg_dq = (-p);
		dg_ds11 = dg_dp / ffmat(3.0) + dg_dq * dq_ds11;
		dg_ds22 = dg_dp / ffmat(3.0) + dg_dq * dq_ds22;
		dg_ds33 = dg_dp / ffmat(3.0) + dg_dq * dq_ds33;
		dg_ds12 = dg_dq * dq_ds12;
		dg_ds23 = dg_dq * dq_ds23;
		dg_ds31 = dg_dq * dq_ds31;
		
		// A
		state_param = mat_dat.e - glb_dat.gamma + glb_dat.lambda * (__Float_Type__)log(-p);
		pi_max = -p * (__Float_Type__)exp(-glb_dat.chi * state_param);
		A = glb_dat.H * M_lode_coef * (pi_max - mat_dat.pi) * dg_dq;

        // dl
		const __Float_Type__ E_dg_ds[3] = {
			lbd_2G * dg_ds11 + lbd * dg_ds22 + lbd * dg_ds33,
			lbd * dg_ds11 + lbd_2G * dg_ds22 + lbd * dg_ds33,
			lbd * dg_ds11 + lbd * dg_ds22 + lbd_2G * dg_ds33
		};
		const __Float_Type__ divider = df_ds11 * E_dg_ds[0] + df_ds22 * E_dg_ds[1] + df_ds33 * E_dg_ds[2]
			+ df_ds12 * G * dg_ds12 + df_ds23 * G * dg_ds23 + df_ds31 * G * dg_ds31
			- (Mi * p / mat_dat.pi) * A;
		dl = f / divider;

		// strain correction
		dep_cor[0] = dl * dg_ds11;
		dep_cor[1] = dl * dg_ds22;
		dep_cor[2] = dl * dg_ds33,
		dep_cor[3] = dl * dg_ds12;
		dep_cor[4] = dl * dg_ds23;
		dep_cor[5] = dl * dg_ds31;
		// correct stress
		mat_dat.s11 -= lbd_2G * dep_cor[0] + lbd * dep_cor[1] + lbd * dep_cor[2];
		mat_dat.s22 -= lbd * dep_cor[0] + lbd_2G * dep_cor[1] + lbd * dep_cor[2];
		mat_dat.s33 -= lbd * dep_cor[0] + lbd * dep_cor[1] + lbd_2G * dep_cor[2];
		mat_dat.s12 -= G * dep_cor[3];
		mat_dat.s23 -= G * dep_cor[4];
		mat_dat.s31 -= G * dep_cor[5];
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
		// image stress
		mat_dat.pi += A * dl;
		if (mat_dat.pi < pi_min)
			mat_dat.pi = pi_min;

		cal_p_q(mat_dat.stress, invars);
		if (p > -p_min)
		{
			q_ratio = (Mi * p_min * (ffmat(1.0) - (__Float_Type__)log(p_min / mat_dat.pi))) / q;
			if (q_ratio > ffmat(1.0))
				q_ratio = ffmat(1.0);
			mat_dat.s11 = (mat_dat.s11 - p) * q_ratio + (-p_min);
			mat_dat.s22 = (mat_dat.s22 - p) * q_ratio + (-p_min);
			mat_dat.s33 = (mat_dat.s33 - p) * q_ratio + (-p_min);
			mat_dat.s12 = mat_dat.s12 * q_ratio;
			mat_dat.s23 = mat_dat.s23 * q_ratio;
			mat_dat.s31 = mat_dat.s31 * q_ratio;
			stress_tmp[0] = mat_dat.s11 - ori_stress[0];
			stress_tmp[1] = mat_dat.s22 - ori_stress[1];
			stress_tmp[2] = mat_dat.s33 - ori_stress[2];
			stress_tmp[3] = mat_dat.s12 - ori_stress[3];
			stress_tmp[4] = mat_dat.s23 - ori_stress[4];
			stress_tmp[5] = mat_dat.s31 - ori_stress[5];
			destrain[0] = inv_E * stress_tmp[0] - niu_inv_E * stress_tmp[1] - niu_inv_E * stress_tmp[2];
			destrain[1] = -niu_inv_E * stress_tmp[0] + inv_E * stress_tmp[1] - niu_inv_E * stress_tmp[2];
			destrain[2] = -niu_inv_E * stress_tmp[0] - niu_inv_E * stress_tmp[1] + inv_E * stress_tmp[2];
			destrain[3] = inv_2G * stress_tmp[3];
			destrain[4] = inv_2G * stress_tmp[4];
			destrain[5] = inv_2G * stress_tmp[5];
			dpstrain[0] = dstrain[0] - destrain[0];
			dpstrain[1] = dstrain[1] - destrain[1];
			dpstrain[2] = dstrain[2] - destrain[2];
			dpstrain[3] = dstrain[3] - destrain[3];
			dpstrain[4] = dstrain[4] - destrain[4];
			dpstrain[5] = dstrain[5] - destrain[5];
			return iter_id + 1;
		}
		
		f = q + Mi * p * (ffmat(1.0) - (__Float_Type__)log(-p / mat_dat.pi));
		if (f < mat_dat.pi * ffmat(1.0e-2)) // converge
			return iter_id + 1;
	}
	
	// integration failure
	destrain[0] = ffmat(0.0);
	destrain[1] = ffmat(0.0);
	destrain[2] = ffmat(0.0);
	destrain[3] = ffmat(0.0);
	destrain[4] = ffmat(0.0);
	destrain[5] = ffmat(0.0);
	dpstrain[0] = dstrain[0];
	dpstrain[1] = dstrain[1];
	dpstrain[2] = dstrain[2];
	dpstrain[3] = dstrain[3];
	dpstrain[4] = dstrain[4];
	dpstrain[5] = dstrain[5];
	mat_dat.s11 = ori_stress[0];
	mat_dat.s22 = ori_stress[1];
	mat_dat.s33 = ori_stress[2];
	mat_dat.s12 = ori_stress[3];
	mat_dat.s23 = ori_stress[4];
	mat_dat.s31 = ori_stress[5];
	mat_dat.e = ori_e;
	return -1;
}
