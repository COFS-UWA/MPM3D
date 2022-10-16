#include "MaterialModels_pcp.h"

#include <math.h>

#include "mohr_coulomb.h"
#include "MatModelConstants.h"
#include "SymMatEigen.h"

void MohrCoulombGlobal::set_param(
	__Float_Type__ _phi, __Float_Type__ _psi,
	__Float_Type__ _cohesion,
	__Float_Type__ _E, __Float_Type__ _niu)
{
	E = _E;
	niu = _niu;
	G = E / (ffmat(1.0) + niu);
	lambda = E * niu / ((ffmat(1.0) + niu) * (ffmat(1.0) - niu - niu));
	lambda_2G = lambda + G;
	inv_E = ffmat(1.0) / E;
	niu_div_E = niu / E;
	inv_G = ffmat(1.0) / G;

	phi = _phi;
	psi = _psi;
	cohesion = _cohesion;
	const __Float_Type__ sin_phi = (__Float_Type__)sin(ToRadian(phi));
	half_one_minus_sin_phi = ffmat(0.5) * (ffmat(1.0) - sin_phi);
	half_one_plus_sin_phi = ffmat(0.5) * (ffmat(1.0) + sin_phi);
	const __Float_Type__ cos_phi = (__Float_Type__)cos(ToRadian(phi));
	c_cos_phi = cohesion * cos_phi;
	c_div_tan_phi = cohesion / (__Float_Type__)tan(ToRadian(phi));
	const __Float_Type__ sin_psi = (__Float_Type__)sin(ToRadian(psi));
	one_minus_sin_psi = ffmat(1.0) - sin_psi;
	one_plus_sin_psi = ffmat(1.0) + sin_psi;
	one_plus_sin_phi_psi = ffmat(1.0) + sin_phi * sin_psi;

	__Float_Type__ len_tmp;
	i1[0][0] = one_plus_sin_phi_psi * (ffmat(1.0) + sin_psi);
	i1[0][1] = (ffmat(1.0) - sin_phi) * (ffmat(1.0) + sin_psi * sin_psi);
	i1[0][2] = one_plus_sin_phi_psi * (ffmat(1.0) - sin_psi);
	len_tmp = (__Float_Type__)sqrt(i1[0][0] * i1[0][0] + i1[0][1] * i1[0][1] + i1[0][2] * i1[0][2]);
	i1[0][0] /= len_tmp;
	i1[0][1] /= len_tmp;
	i1[0][2] /= len_tmp;
	i1[1][0] = (ffmat(1.0) - sin_phi) * (ffmat(1.0) + sin_psi);
	i1[1][1] = -ffmat(2.0) * (ffmat(1.0) + sin_phi * sin_psi);
	i1[1][2] = (ffmat(1.0) - sin_phi) * (ffmat(1.0) - sin_psi);
	len_tmp = (__Float_Type__)sqrt(i1[1][0] * i1[1][0] + i1[1][1] * i1[1][1] + i1[1][2] * i1[1][2]);
	i1[1][0] /= len_tmp;
	i1[1][1] /= len_tmp;
	i1[1][2] /= len_tmp;

	i2[0][0] = one_plus_sin_phi_psi * (ffmat(1.0) + sin_psi);
	i2[0][1] = (ffmat(1.0) + sin_phi) * (ffmat(1.0) + sin_psi * sin_psi);
	i2[0][2] = one_plus_sin_phi_psi * (ffmat(1.0) - sin_psi);
	len_tmp = (__Float_Type__)sqrt(i2[0][0] * i2[0][0] + i2[0][1] * i2[0][1] + i2[0][2] * i2[0][2]);
	i2[0][0] /= len_tmp;
	i2[0][1] /= len_tmp;
	i2[0][2] /= len_tmp;
	i2[1][0] = (ffmat(1.0) + sin_phi) * (ffmat(1.0) + sin_psi);
	i2[1][1] = -ffmat(2.0) * (ffmat(1.0) + sin_phi * sin_psi);
	i2[1][2] = (ffmat(1.0) + sin_phi) * (ffmat(1.0) - sin_psi);
	len_tmp = (__Float_Type__)sqrt(i2[1][0] * i2[1][0] + i2[1][1] * i2[1][1] + i2[1][2] * i2[1][2]);
	i2[1][0] /= len_tmp;
	i2[1][1] /= len_tmp;
	i2[1][2] /= len_tmp;

	v1x = ffmat(1.0) + sin_phi;
	v1y = ffmat(1.0) - sin_phi;
	v1z = ffmat(1.0) - sin_phi;
	v2x = ffmat(1.0) + sin_phi;
	v2y = ffmat(1.0) + sin_phi;
	v2z = ffmat(1.0) - sin_phi;
	u1x = ffmat(1.0) + sin_psi;
	u1y = ffmat(1.0) - sin_psi;
	u1z = ffmat(1.0) - sin_psi;
	u2x = ffmat(1.0) + sin_psi;
	u2y = ffmat(1.0) + sin_psi;
	u2z = ffmat(1.0) - sin_psi;
	uv1x_uv1y_uv1z = u1x * v1x + u1y * v1y + u1z * v1z;
	uv2x_uv2y_uv2z = u2x * v2x + u2y * v2y + u2z * v2z;
}

int32_t integrate_mohr_coulomb(
	const MohrCoulombGlobal& glb_dat,
	MohrCoulomb& mat_dat,
	const __Float_Type__ dstrain[6],
	__Float_Type__ destrain[6],
	__Float_Type__ dpstrain[6])
{
	mat_dat.stress[0] += glb_dat.lambda_2G * dstrain[0]
		+ glb_dat.lambda * dstrain[1]
		+ glb_dat.lambda * dstrain[2];
	mat_dat.stress[1] += glb_dat.lambda * dstrain[0]
		+ glb_dat.lambda_2G * dstrain[1]
		+ glb_dat.lambda * dstrain[2];
	mat_dat.stress[2] += glb_dat.lambda * dstrain[0]
		+ glb_dat.lambda * dstrain[1]
		+ glb_dat.lambda_2G * dstrain[2];
	mat_dat.stress[3] += glb_dat.G * dstrain[3];
	mat_dat.stress[4] += glb_dat.G * dstrain[4];
	mat_dat.stress[5] += glb_dat.G * dstrain[5];

	__Float_Type__ eigen_vals[3], eigen_vecs[3][3];
	cal_sym_mat_eigen(mat_dat.stress, eigen_vals, eigen_vecs);
	__Float_Type__ s1 = -eigen_vals[2];
	__Float_Type__ s2 = -eigen_vals[1];
	__Float_Type__ s3 = -eigen_vals[0];

	const __Float_Type__ F = glb_dat.half_one_minus_sin_phi * s1
		- s3 * glb_dat.half_one_plus_sin_phi - glb_dat.c_cos_phi;
	if (F <= ffmat(0.0)) // elasticity
	{
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
		return 0;
	}

	const __Float_Type__ s1_ctan = s1 + glb_dat.c_div_tan_phi;
	const __Float_Type__ s2_ctan = s2 + glb_dat.c_div_tan_phi;
	const __Float_Type__ s3_ctan = s3 + glb_dat.c_div_tan_phi;

	const __Float_Type__ sx_1 = glb_dat.i1[0][0] * s1_ctan + glb_dat.i1[0][1] * s2_ctan + glb_dat.i1[0][2] * s3_ctan;
	const __Float_Type__ sy_1 = glb_dat.i1[1][0] * s1_ctan + glb_dat.i1[1][1] * s2_ctan + glb_dat.i1[1][2] * s3_ctan;
	const __Float_Type__ sx_2 = glb_dat.i2[0][0] * s1_ctan + glb_dat.i2[0][1] * s2_ctan + glb_dat.i2[0][2] * s3_ctan;
	const __Float_Type__ sy_2 = glb_dat.i2[1][0] * s1_ctan + glb_dat.i2[1][1] * s2_ctan + glb_dat.i2[1][2] * s3_ctan;

	__Float_Type__ t;
	if (sy_1 <= ffmat(0.0) && sy_2 >= ffmat(0.0))
	{
		// area 0
		t = F / glb_dat.one_plus_sin_phi_psi;
		s1 -= t * glb_dat.one_minus_sin_psi;
		s3 += t * glb_dat.one_plus_sin_psi;
	}
	else if (sy_1 > ffmat(0.0) && sx_1 > ffmat(0.0))
	{
		// area 1
		t = (glb_dat.u1x * (s1 + glb_dat.c_div_tan_phi)
			+ glb_dat.u1y * (s2 + glb_dat.c_div_tan_phi)
			+ glb_dat.u1z * (s3 + glb_dat.c_div_tan_phi))
			/ glb_dat.uv1x_uv1y_uv1z;
		s1 = -glb_dat.c_div_tan_phi + glb_dat.v1x * t;
		s2 = -glb_dat.c_div_tan_phi + glb_dat.v1y * t;
		s3 = -glb_dat.c_div_tan_phi + glb_dat.v1z * t;
	}
	else if (sy_2 < ffmat(0.0) && sx_2 > ffmat(0.0))
	{
		// area 2
		t = (glb_dat.u2x * (s1 + glb_dat.c_div_tan_phi)
			+ glb_dat.u2y * (s2 + glb_dat.c_div_tan_phi)
			+ glb_dat.u2z * (s3 + glb_dat.c_div_tan_phi))
			/ glb_dat.uv2x_uv2y_uv2z;
		s1 = -glb_dat.c_div_tan_phi + glb_dat.v2x * t;
		s2 = -glb_dat.c_div_tan_phi + glb_dat.v2y * t;
		s3 = -glb_dat.c_div_tan_phi + glb_dat.v2z * t;
	}
	else
	{
		// area 3
		s1 = -glb_dat.c_div_tan_phi;
		s2 = -glb_dat.c_div_tan_phi;
		s3 = -glb_dat.c_div_tan_phi;
	}

	eigen_vals[0] = -s3;
	eigen_vals[1] = -s2;
	eigen_vals[2] = -s1;

	// rotate stress back
	__Float_Type__ s_corrected[6];
	rotate_eigen_mat_to_sym_mat(eigen_vals, eigen_vecs, s_corrected);
	//const __Float_Type__ s_tmp[3][3] = {
	//	eigen_vals[0] * eigen_vecs[0][0], eigen_vals[0] * eigen_vecs[1][0], eigen_vals[0] * eigen_vecs[2][0],
	//	eigen_vals[1] * eigen_vecs[0][1], eigen_vals[1] * eigen_vecs[1][1], eigen_vals[1] * eigen_vecs[2][1],
	//	eigen_vals[2] * eigen_vecs[0][2], eigen_vals[2] * eigen_vecs[1][2], eigen_vals[2] * eigen_vecs[2][2]
	//};
	//const __Float_Type__ s_corrected[6] = {
	//	eigen_vecs[0][0] * s_tmp[0][0] + eigen_vecs[0][1] * s_tmp[1][0] + eigen_vecs[0][2] * s_tmp[2][0],
	//	eigen_vecs[1][0] * s_tmp[0][1] + eigen_vecs[1][1] * s_tmp[1][1] + eigen_vecs[1][2] * s_tmp[2][1],
	//	eigen_vecs[2][0] * s_tmp[0][2] + eigen_vecs[2][1] * s_tmp[1][2] + eigen_vecs[2][2] * s_tmp[2][2],
	//	eigen_vecs[0][0] * s_tmp[0][1] + eigen_vecs[0][1] * s_tmp[1][1] + eigen_vecs[0][2] * s_tmp[2][1],
	//	eigen_vecs[0][0] * s_tmp[0][2] + eigen_vecs[0][1] * s_tmp[1][2] + eigen_vecs[0][2] * s_tmp[2][2],
	//	eigen_vecs[1][0] * s_tmp[0][2] + eigen_vecs[1][1] * s_tmp[1][2] + eigen_vecs[1][2] * s_tmp[2][2]
	//};

	const __Float_Type__ ds_corrected[3] = {
		mat_dat.stress[0] - s_corrected[0],
		mat_dat.stress[1] - s_corrected[1],
		mat_dat.stress[2] - s_corrected[2]
	};

	dpstrain[0] = glb_dat.inv_E * ds_corrected[0] - glb_dat.niu_div_E * ds_corrected[1] - glb_dat.niu_div_E * ds_corrected[2];
	dpstrain[1] = -glb_dat.niu_div_E * ds_corrected[0] + glb_dat.inv_E * ds_corrected[1] - glb_dat.niu_div_E * ds_corrected[2];
	dpstrain[2] = -glb_dat.niu_div_E * ds_corrected[0] - glb_dat.niu_div_E * ds_corrected[1] + glb_dat.inv_E * ds_corrected[2];
	dpstrain[3] = (mat_dat.stress[3] - s_corrected[3]) * glb_dat.inv_G;
	dpstrain[4] = (mat_dat.stress[4] - s_corrected[4]) * glb_dat.inv_G;
	dpstrain[5] = (mat_dat.stress[5] - s_corrected[5]) * glb_dat.inv_G;
	destrain[0] = dstrain[0] - dpstrain[0];
	destrain[1] = dstrain[1] - dpstrain[1];
	destrain[2] = dstrain[2] - dpstrain[2];
	destrain[3] = dstrain[3] - dpstrain[3];
	destrain[4] = dstrain[4] - dpstrain[4];
	destrain[5] = dstrain[5] - dpstrain[5];
	mat_dat.stress[0] = s_corrected[0];
	mat_dat.stress[1] = s_corrected[1];
	mat_dat.stress[2] = s_corrected[2];
	mat_dat.stress[3] = s_corrected[3];
	mat_dat.stress[4] = s_corrected[4];
	mat_dat.stress[5] = s_corrected[5];
	return 1;
}
