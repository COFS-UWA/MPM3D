#include "MaterialModels_pcp.h"

#include <math.h>

#include "sand_hypoplasticity.h"
#include "MatModelConstants.h"
#include "SymMatEigen.h"

// minimum compressive stress for
// hypoplasticity to be meaning
#define min_com_stress ffmat(1.0)

// error limit for each substep
#define error_tol ffmat(1.0e-2)

#define max_substp_num (100)
#define min_substp_size ffmat(2.0e-2)
#define tensile_min_substp_size ffmat(5.0e-2)
// < min_substp_size and tensile_min_substp_size
#define substp_tol ffmat(1.0e-3)

#define min(a, b) ((a) < (b) ? (a) : (b))
#define max(a, b) ((a) > (b) ? (a) : (b))

void SandHypoplasticityGlobal::set_param(
	__Float_Type__ _phi,
	__Float_Type__ _hs, __Float_Type__ _n,
	__Float_Type__ _alpha, __Float_Type__ _beta,
	__Float_Type__ _ed0, __Float_Type__ _ec0,
	__Float_Type__ _ei0,
	__Float_Type__ _tE, __Float_Type__ _tniu)
{
	phi = _phi;
	hs = _hs;
	n = _n;
	alpha = _alpha;
	beta = _beta;
	__Float_Type__ sin_phi = (__Float_Type__)sin(phi / ffmat(180.0) * PI);
	a = sqrt3 * (ffmat(3.0) - sin_phi) / (ffmat(2.0) * sqrt2 * sin_phi);
	
	ed0 = _ed0;
	ec0 = _ec0;
	ei0 = _ei0;
	fbfe_coef = hs * (__Float_Type__)pow(ei0/ec0, beta) / 
		(n * (ffmat(3.0) + a * a - a * sqrt3 * (__Float_Type__)pow((ei0 - ed0) / (ec0 - ed0), alpha)));
	
	ten_E = _tE;
	ten_niu = _tniu;
	ten_G = ten_E / (ffmat(1.0) + ten_niu);
	ten_lambda = ten_E * ten_niu / ((ffmat(1.0) + ten_niu) * (ffmat(1.0) - ten_niu - ten_niu));
	ten_lambda_2G = ten_lambda + ten_G;
}

int32_t SandHypoplasticityGlobal::hypoplasticity_substp(
	const __Float_Type__ stress[6],
	__Float_Type__ e,
	const __Float_Type__ dstrain[6],
	__Float_Type__ dstress[6],
	__Float_Type__& de) const
{
	const __Float_Type__ I1 = stress[0] + stress[1] + stress[2];

	const __Float_Type__ s_cap[6] = {
		stress[0] / I1,
		stress[1] / I1,
		stress[2] / I1,
		stress[3] / I1,
		stress[4] / I1,
		stress[5] / I1
	};

	const __Float_Type__ s_cap2 =
		  s_cap[0] * s_cap[0]
		+ s_cap[1] * s_cap[1]
		+ s_cap[2] * s_cap[2]
		+ s_cap[3] * s_cap[3] * ffmat(2.0)
		+ s_cap[4] * s_cap[4] * ffmat(2.0)
		+ s_cap[5] * s_cap[5] * ffmat(2.0);

	const __Float_Type__ s_star[3] = {
		s_cap[0] - ffmat(1.0) / ffmat(3.0),
		s_cap[1] - ffmat(1.0) / ffmat(3.0),
		s_cap[2] - ffmat(1.0) / ffmat(3.0)
	};

	const __Float_Type__ e_ratio = (__Float_Type__)exp(-pow(-I1/hs, n));
	const __Float_Type__ ei = ei0 * e_ratio;
	const __Float_Type__ ec = ec0 * e_ratio;
	const __Float_Type__ ed = ed0 * e_ratio;
	const __Float_Type__ fd = e > ed ? (__Float_Type__)pow((e - ed) / (ec - ed), alpha) : ffmat(0.0);
	const __Float_Type__ fbfe = fbfe_coef * (ffmat(1.0) + ei) / ei * (__Float_Type__)pow(-I1 / hs, ffmat(1.0) - n)
		* (__Float_Type__)pow(ec/e, beta);
	__Float_Type__ fb = fbfe_coef * (ffmat(1.0) + ei) / ei * (__Float_Type__)pow(-I1 / hs, ffmat(1.0) - n);
	__Float_Type__ efe = hs * (__Float_Type__)pow(ei0 / ec0, beta) /
		(n * (ffmat(3.0) + a * a - a * sqrt3 * (__Float_Type__)pow((ei0 - ed0) / (ec0 - ed0), alpha)));
	__Float_Type__ efe1 = (__Float_Type__)pow(ei0 / ec0, beta);
	__Float_Type__ efe2 = pow((ei0 - ed0) / (ec0 - ed0), alpha);
	__Float_Type__ temp = ffmat(3.0) + a * a - a * sqrt3 * (__Float_Type__)pow((ei0 - ed0) / (ec0 - ed0), alpha);
	__Float_Type__ fb1 = hs / n / (ffmat(3.0) + a * a - a * sqrt3 * (__Float_Type__)pow((ei0 - ed0) / (ec0 - ed0), alpha))
		* (ffmat(1.0) + ei) / ei * (__Float_Type__)pow(ei0 / ec0, beta) * (__Float_Type__)pow(-I1 / hs, ffmat(1.0) - n);


	const __Float_Type__ s_star_norm = (__Float_Type__)sqrt(
		  s_star[0] * s_star[0]
		+ s_star[1] * s_star[1]
		+ s_star[2] * s_star[2]
		+ s_cap[3] * s_cap[3] * ffmat(2.0)
		+ s_cap[4] * s_cap[4] * ffmat(2.0)
		+ s_cap[5] * s_cap[5] * ffmat(2.0));
	
	const __Float_Type__ TTT = (
		  s_star[0] * s_star[0] * s_star[0] // s11 * s11 * s11
		+ s_star[1] * s_star[1] * s_star[1] // s22 * s22 * s22
		+ s_star[2] * s_star[2] * s_star[2] // s33 * s33 * s33
		+ s_star[0] * s_cap[3] * s_cap[3] * ffmat(3.0) // s11 * s12 * s12
		+ s_star[0] * s_cap[5] * s_cap[5] * ffmat(3.0) // s11 * s31 * s31
		+ s_star[1] * s_cap[3] * s_cap[3] * ffmat(3.0) // s22 * s12 * s12
		+ s_star[1] * s_cap[4] * s_cap[4] * ffmat(3.0) // s22 * s23 * s23
		+ s_star[2] * s_cap[4] * s_cap[4] * ffmat(3.0) // s33 * s23 * s23
		+ s_star[2] * s_cap[5] * s_cap[5] * ffmat(3.0) // s33 * s31 * s31
		+ s_cap[3] * s_cap[4] * s_cap[5] * ffmat(6.0)); // s12 * s23 * s31
	__Float_Type__ cos_3theta;
	if ((__Float_Type__)fabs(TTT) > ffmat(1.0e-10))
	{
		cos_3theta = -sqrt6 * TTT / (s_star_norm * s_star_norm * s_star_norm);
		if (cos_3theta > ffmat(1.0))
			cos_3theta = ffmat(1.0);
		if (cos_3theta < ffmat(-1.0))
			cos_3theta = ffmat(-1.0);
	}
	else
	{
		if (TTT > ffmat(0.0))
			cos_3theta = ffmat(-1.0);
		else
			cos_3theta = ffmat(1.0);
	}

	const __Float_Type__ tan_psai = sqrt3 * s_star_norm;
	const __Float_Type__ F = (__Float_Type__)sqrt(tan_psai * tan_psai / ffmat(8.0) + (ffmat(2.0) - tan_psai * tan_psai) / (ffmat(2.0) + sqrt2 * tan_psai * cos_3theta))
			- tan_psai / (ffmat(2.0) * sqrt2);

	const __Float_Type__ fbfe_div_s_cap2 = fbfe / s_cap2;

	const __Float_Type__ dstrain_norm = (__Float_Type__)sqrt(
		  dstrain[0] * dstrain[0]
		+ dstrain[1] * dstrain[1]
		+ dstrain[2] * dstrain[2]
		+ dstrain[3] * dstrain[3] * ffmat(2.0)
		+ dstrain[4] * dstrain[4] * ffmat(2.0)
		+ dstrain[5] * dstrain[5] * ffmat(2.0));
	
	const __Float_Type__ L1_tmp = fbfe_div_s_cap2 * F * F;
	const __Float_Type__ L2_tmp = fbfe_div_s_cap2 * a * a * (
		  s_cap[0] * dstrain[0]
		+ s_cap[1] * dstrain[1]
		+ s_cap[2] * dstrain[2]
		+ s_cap[3] * dstrain[3] * ffmat(2.0)
		+ s_cap[4] * dstrain[4] * ffmat(2.0)
		+ s_cap[5] * dstrain[5] * ffmat(2.0));
	const __Float_Type__ N_tmp = fd * fbfe_div_s_cap2 * F * a * dstrain_norm;
	double N[6] = {
		N_tmp* (s_cap[0] + s_star[0]),
		N_tmp* (s_cap[1] + s_star[1]),
		N_tmp* (s_cap[2] + s_star[2]),
		N_tmp* (s_cap[3] + s_cap[3]),
		N_tmp* (s_cap[4] + s_cap[4]),
		N_tmp* (s_cap[5] + s_cap[5])
	};
	dstress[0] = L1_tmp * dstrain[0] + L2_tmp * s_cap[0] + N_tmp * (s_cap[0] + s_star[0]);
	dstress[1] = L1_tmp * dstrain[1] + L2_tmp * s_cap[1] + N_tmp * (s_cap[1] + s_star[1]);
	dstress[2] = L1_tmp * dstrain[2] + L2_tmp * s_cap[2] + N_tmp * (s_cap[2] + s_star[2]);
	dstress[3] = L1_tmp * dstrain[3] + L2_tmp * s_cap[3] + N_tmp * (s_cap[3] + s_cap[3]);
	dstress[4] = L1_tmp * dstrain[4] + L2_tmp * s_cap[4] + N_tmp * (s_cap[4] + s_cap[4]);
	dstress[5] = L1_tmp * dstrain[5] + L2_tmp * s_cap[5] + N_tmp * (s_cap[5] + s_cap[5]);

	de = (ffmat(1.0) + e) * (dstrain[0] + dstrain[1] + dstrain[2]);
	return 1;
}

void SandHypoplasticityGlobal::elasticity(
	const __Float_Type__ dstrain[6],
	__Float_Type__ stress[6]) const
{
	stress[0] += ten_lambda_2G * dstrain[0]
		+ ten_lambda * dstrain[1]
		+ ten_lambda * dstrain[2];
	stress[1] += ten_lambda * dstrain[0]
		+ ten_lambda_2G * dstrain[1]
		+ ten_lambda * dstrain[2];
	stress[2] += ten_lambda * dstrain[0]
		+ ten_lambda * dstrain[1]
		+ ten_lambda_2G * dstrain[2];
	stress[3] += ten_G * dstrain[3];
	stress[4] += ten_G * dstrain[4];
	stress[5] += ten_G * dstrain[5];
}

bool in_tensile_state(const __Float_Type__ stress[6])
{
	const __Float_Type__ I1 = stress[0] + stress[1] + stress[2];
	if (-I1 < min_com_stress * ffmat(3.0))
		return true;
	__Float_Type__ prin_s[3];
	cal_sym_mat_eigen(stress, prin_s);
	if (-prin_s[2] < min_com_stress)
		return true;
	return false;
}

int32_t integrate_sand_hypoplasticity(
	const SandHypoplasticityGlobal& glb_dat,
	SandHypoplasticity& mat_dat,
	const __Float_Type__ dstrain[6],
	__Float_Type__ substp_size_ratio)
{
	if (in_tensile_state(mat_dat.stress))
	{
		glb_dat.elasticity(dstrain, mat_dat.stress);
		return 0;
	}

	// RKF23 hypoplasticity integration
	__Float_Type__ ddstrain[6];
	__Float_Type__ dstress1[6], de1, stress1[6], e1;
	__Float_Type__ dstress2[6], de2, stress2[6], e2;
	__Float_Type__ dstress3[6], de3, stress3[6], e3;
	uint32_t substp_id = 0;
	__Float_Type__ substp_size = min(mat_dat.substp_size * substp_size_ratio, ffmat(1.0));
	__Float_Type__ total_size = ffmat(0.0);
	__Float_Type__ act_substp_size;
	while (total_size < (ffmat(1.0) - substp_tol))
	{
		if (++substp_id > max_substp_num)
			return -3;

		act_substp_size = min(substp_size, ffmat(1.0) - total_size);

		ddstrain[0] = act_substp_size * dstrain[0];
		ddstrain[1] = act_substp_size * dstrain[1];
		ddstrain[2] = act_substp_size * dstrain[2];
		ddstrain[3] = act_substp_size * dstrain[3];
		ddstrain[4] = act_substp_size * dstrain[4];
		ddstrain[5] = act_substp_size * dstrain[5];

		glb_dat.hypoplasticity_substp(mat_dat.stress,
			mat_dat.e, ddstrain, dstress1, de1);

		stress1[0] = mat_dat.stress[0] + dstress1[0] * ffmat(0.5);
		stress1[1] = mat_dat.stress[1] + dstress1[1] * ffmat(0.5);
		stress1[2] = mat_dat.stress[2] + dstress1[2] * ffmat(0.5);
		stress1[3] = mat_dat.stress[3] + dstress1[3] * ffmat(0.5);
		stress1[4] = mat_dat.stress[4] + dstress1[4] * ffmat(0.5);
		stress1[5] = mat_dat.stress[5] + dstress1[5] * ffmat(0.5);
		if (in_tensile_state(stress1))
		{
			substp_size *= 0.25;
			if (substp_size < tensile_min_substp_size)
				return 0;
			continue;
		}
		e1 = mat_dat.e + de1 * ffmat(0.5);

		glb_dat.hypoplasticity_substp(
			stress1, e1, ddstrain, dstress2, de2);

		stress2[0] = mat_dat.stress[0] - dstress1[0] + dstress2[0] + dstress2[0];
		stress2[1] = mat_dat.stress[1] - dstress1[1] + dstress2[1] + dstress2[1];
		stress2[2] = mat_dat.stress[2] - dstress1[2] + dstress2[2] + dstress2[2];
		stress2[3] = mat_dat.stress[3] - dstress1[3] + dstress2[3] + dstress2[3];
		stress2[4] = mat_dat.stress[4] - dstress1[4] + dstress2[4] + dstress2[4];
		stress2[5] = mat_dat.stress[5] - dstress1[5] + dstress2[5] + dstress2[5];
		if (in_tensile_state(stress2))
		{
			substp_size *= 0.25;
			if (substp_size < tensile_min_substp_size)
				return 0;
			continue;
		}
		e2 = mat_dat.e - de1 + de2 + de2;

		glb_dat.hypoplasticity_substp(stress2, e2, ddstrain, dstress3, de3);
			
		dstress3[0] = ffmat(1.0) / ffmat(6.0) * dstress1[0]
					+ ffmat(2.0) / ffmat(3.0) * dstress2[0]
					+ ffmat(1.0) / ffmat(6.0) * dstress3[0];
		dstress3[1] = ffmat(1.0) / ffmat(6.0) * dstress1[1]
					+ ffmat(2.0) / ffmat(3.0) * dstress2[1]
					+ ffmat(1.0) / ffmat(6.0) * dstress3[1];
		dstress3[2] = ffmat(1.0) / ffmat(6.0) * dstress1[2]
					+ ffmat(2.0) / ffmat(3.0) * dstress2[2]
					+ ffmat(1.0) / ffmat(6.0) * dstress3[2];
		dstress3[3] = ffmat(1.0) / ffmat(6.0) * dstress1[3]
					+ ffmat(2.0) / ffmat(3.0) * dstress2[3]
					+ ffmat(1.0) / ffmat(6.0) * dstress3[3];
		dstress3[4] = ffmat(1.0) / ffmat(6.0) * dstress1[4]
					+ ffmat(2.0) / ffmat(3.0) * dstress2[4]
					+ ffmat(1.0) / ffmat(6.0) * dstress3[4];
		dstress3[5] = ffmat(1.0) / ffmat(6.0) * dstress1[5]
					+ ffmat(2.0) / ffmat(3.0) * dstress2[5]
					+ ffmat(1.0) / ffmat(6.0) * dstress3[5];

		stress3[0] = mat_dat.stress[0] + dstress3[0];
		stress3[1] = mat_dat.stress[1] + dstress3[1];
		stress3[2] = mat_dat.stress[2] + dstress3[2];
		stress3[3] = mat_dat.stress[3] + dstress3[3];
		stress3[4] = mat_dat.stress[4] + dstress3[4];
		stress3[5] = mat_dat.stress[5] + dstress3[5];
		if (in_tensile_state(stress3))
		{
			substp_size *= 0.25;
			if (substp_size < tensile_min_substp_size)
				return 0;
			continue;
		}
		de3 = ffmat(1.0) / ffmat(6.0) * de1
			+ ffmat(2.0) / ffmat(3.0) * de2
			+ ffmat(1.0) / ffmat(6.0) * de3;
		e3 = mat_dat.e + de3;

		const __Float_Type__ error2 = ((de3 - de2) * (de3 - de2) 
				+ (dstress3[0] - dstress2[0]) * (dstress3[0] - dstress2[0])
				+ (dstress3[1] - dstress2[1]) * (dstress3[1] - dstress2[1])
				+ (dstress3[2] - dstress2[2]) * (dstress3[2] - dstress2[2])
				+ (dstress3[3] - dstress2[3]) * (dstress3[3] - dstress2[3])
				+ (dstress3[4] - dstress2[4]) * (dstress3[4] - dstress2[4])
				+ (dstress3[5] - dstress2[5]) * (dstress3[5] - dstress2[5]))
			/ (e3 * e3 + stress3[0] * stress3[0] + stress3[1] * stress3[1]
				+ stress3[2] * stress3[2] + stress3[3] * stress3[3]
				+ stress3[4] * stress3[4] + stress3[5] * stress3[5]);

		__Float_Type__ substp_size_adjust_ratio = ffmat(4.0);
		if (error2 != ffmat(0.0))
			substp_size_adjust_ratio = ffmat(0.9) * (__Float_Type__)pow(error_tol * error_tol / error2, ffmat(1.0) / ffmat(8.0));
			
		if (error2 < error_tol * error_tol) // accept substep
		{
			total_size += act_substp_size;
			mat_dat.stress[0] = stress3[0];
			mat_dat.stress[1] = stress3[1];
			mat_dat.stress[2] = stress3[2];
			mat_dat.stress[3] = stress3[3];
			mat_dat.stress[4] = stress3[4];
			mat_dat.stress[5] = stress3[5];
			mat_dat.e = e3;
			substp_size *= min(substp_size_adjust_ratio, ffmat(2.0));
			substp_size = min(substp_size, ffmat(1.0));
			continue;
		}
		
		// reject substep
		substp_size *= max(substp_size_adjust_ratio, ffmat(0.25));
		if (substp_size < min_substp_size)
			return -1;
	}
	mat_dat.substp_size = substp_size;
	return substp_id;
}