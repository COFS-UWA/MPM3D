#include "MaterialModels_pcp.h"

#include <stdint.h>
#include <cmath>

#include "SymMatEigen.h"
#include "MatModelConstants.h"

#define DTol (ffmat(1.0e-10))

#define cross_prod(a, b, c) \
	c[0] = a[1] * b[2] - a[2] * b[1]; \
	c[1] = a[2] * b[0] - a[0] * b[2]; \
	c[2] = a[0] * b[1] - a[1] * b[0]

#define swap(a, b, tmp) \
	(tmp) = (a), (a) = (b), (b) = (tmp)

static void cal_evec(const __Float_Type__ mat_va[6], __Float_Type__ ev, __Float_Type__ evec[3])
{
	const __Float_Type__ r0[3] = { mat_va[0] - ev, mat_va[3], mat_va[5] };
	const __Float_Type__ r1[3] = { mat_va[3], mat_va[1] - ev, mat_va[4] };
	const __Float_Type__ r2[3] = { mat_va[5], mat_va[4], mat_va[2] - ev };
	__Float_Type__ r0xr1[3], r1xr2[3], r2xr0[3];
	cross_prod(r0, r1, r0xr1);
	const __Float_Type__ r0xr1_len2 = r0xr1[0] * r0xr1[0] + r0xr1[1] * r0xr1[1] + r0xr1[2] * r0xr1[2];
	cross_prod(r1, r2, r1xr2);
	const __Float_Type__ r1xr2_len2 = r1xr2[0] * r1xr2[0] + r1xr2[1] * r1xr2[1] + r1xr2[2] * r1xr2[2];
	cross_prod(r2, r0, r2xr0);
	const __Float_Type__ r2xr0_len2 = r2xr0[0] * r2xr0[0] + r2xr0[1] * r2xr0[1] + r2xr0[2] * r2xr0[2];
	__Float_Type__ max_len2 = r0xr1_len2;
	const __Float_Type__* max_vec = r0xr1;
	if (r1xr2_len2 > max_len2)
	{
		max_len2 = r1xr2_len2;
		max_vec = r1xr2;
	}
	if (r2xr0_len2 > max_len2)
	{
		max_len2 = r2xr0_len2;
		max_vec = r2xr0;
	}
	max_len2 = (__Float_Type__)sqrt(max_len2);
	evec[0] = max_vec[0] / max_len2;
	evec[1] = max_vec[1] / max_len2;
	evec[2] = max_vec[2] / max_len2;
}

static void cal_ortho_vec(const __Float_Type__ vec0[3], __Float_Type__ vec1[3])
{
	__Float_Type__ inv_len;
	if (abs(vec0[0]) > abs(vec0[1]))
	{
		inv_len = ffmat(1.0) / (__Float_Type__)sqrt(vec0[0] * vec0[0] + vec0[2] * vec0[2]);
		vec1[0] = -vec0[2] * inv_len;
		vec1[1] = ffmat(0.0);
		vec1[2] = vec0[0] * inv_len;
	}
	else
	{
		inv_len = ffmat(1.0) / (__Float_Type__)sqrt(vec0[1] * vec0[1] + vec0[2] * vec0[2]);
		vec1[0] = ffmat(0.0);
		vec1[1] = vec0[2] * inv_len;
		vec1[2] = -vec0[1] * inv_len;
	}
}

static const __Float_Type__ identity_mat[3][3] = {
	ffmat(1.0), ffmat(0.0), ffmat(0.0),
	ffmat(0.0), ffmat(1.0), ffmat(0.0),
	ffmat(0.0), ffmat(0.0), ffmat(1.0)
};

// mat is sysmetric
// mat_va[0] = mat[0][0]
// mat_va[1] = mat[1][1]
// mat_va[2] = mat[2][2]
// mat_va[3] = mat[0][1]
// mat_va[4] = mat[1][2]
// mat_va[5] = mat[2][0]
// return: ev[0] > ev[1] > ev[2]
void cal_sym_mat_eigen(const __Float_Type__ mat_va[6], __Float_Type__ ev[3])
{
	const __Float_Type__ tmp = mat_va[3] * mat_va[3]
		+ mat_va[4] * mat_va[4] + mat_va[5] * mat_va[5];
	if (tmp < DTol * DTol)
	{
		// eigen values
		ev[0] = mat_va[0];
		ev[1] = mat_va[1];
		ev[2] = mat_va[2];
	}
	else
	{
		// cal eigenvalues
		const __Float_Type__ p = (mat_va[0] + mat_va[1] + mat_va[2]) / ffmat(3.0);
		__Float_Type__ matb_va[6];
		matb_va[0] = mat_va[0] - p;
		matb_va[1] = mat_va[1] - p;
		matb_va[2] = mat_va[2] - p;
		__Float_Type__ q = (__Float_Type__)sqrt(
			 (matb_va[0] * matb_va[0]
			+ matb_va[1] * matb_va[1]
			+ matb_va[2] * matb_va[2]
			+ tmp * ffmat(2.0)) / ffmat(6.0));
		matb_va[0] /= q;
		matb_va[1] /= q;
		matb_va[2] /= q;
		matb_va[3] = mat_va[3] / q;
		matb_va[4] = mat_va[4] / q;
		matb_va[5] = mat_va[5] / q;
		const __Float_Type__ detb = (
			  matb_va[0] * matb_va[1] * matb_va[2]
			+ matb_va[3] * matb_va[4] * matb_va[5] * ffmat(2.0)
			- matb_va[0] * matb_va[4] * matb_va[4]
			- matb_va[1] * matb_va[5] * matb_va[5]
			- matb_va[2] * matb_va[3] * matb_va[3]) * ffmat(0.5);
		__Float_Type__ phi;
		if (detb <= ffmat(-1.0))
			phi = PI / ffmat(3.0);
		else if (detb >= ffmat(1.0))
			phi = ffmat(0.0);
		else
			phi = (__Float_Type__)acos(detb) / ffmat(3.0);
		ev[0] = p + ffmat(2.0) * q * (__Float_Type__)cos(phi);
		ev[1] = p + ffmat(2.0) * q * (__Float_Type__)cos(phi + ffmat(2.0) / ffmat(3.0) * PI );
		ev[2] = p + ffmat(2.0) * q * (__Float_Type__)cos(phi + ffmat(4.0) / ffmat(3.0) * PI);
	}

	// sort eigenvalues
	// ev[0] > ev[1] > ev[2]
	__Float_Type__ sw_tmp;
	if (ev[0] < ev[1])
		swap(ev[0], ev[1], sw_tmp);
	if (ev[0] < ev[2])
		swap(ev[0], ev[2], sw_tmp);
	if (ev[1] < ev[2])
		swap(ev[1], ev[2], sw_tmp);
}

// mat is sysmetric
// mat_va[0] = mat[0][0]
// mat_va[1] = mat[1][1]
// mat_va[2] = mat[2][2]
// mat_va[3] = mat[0][1]
// mat_va[4] = mat[1][2]
// mat_va[5] = mat[2][0]
// return: ev[0] > ev[1] > ev[2]
void cal_sym_mat_eigen(
	const __Float_Type__ mat_va[6],
	__Float_Type__ ev[3],
	__Float_Type__ evecs[3][3])
{
	const __Float_Type__ tmp = mat_va[3] * mat_va[3]
		+ mat_va[4] * mat_va[4]	+ mat_va[5] * mat_va[5];
	if (tmp < DTol * DTol)
	{
		// eigen values
		ev[0] = mat_va[0];
		ev[1] = mat_va[1];
		ev[2] = mat_va[2];
		// sort eigenvalues
		// ev[0] > ev[1] > ev[2]
		uint8_t ev0_id = 0;
		uint8_t ev1_id = 1;
		uint8_t ev2_id = 2;
		__Float_Type__ sw_tmp_d;
		uint8_t sw_tmp_uint8;
		if (ev[0] < ev[1])
		{
			swap(ev[0], ev[1], sw_tmp_d);
			swap(ev0_id, ev1_id, sw_tmp_uint8);
		}
		if (ev[0] < ev[2])
		{
			swap(ev[0], ev[2], sw_tmp_d);
			swap(ev0_id, ev2_id, sw_tmp_uint8);
		}
		if (ev[1] < ev[2])
		{
			swap(ev[1], ev[2], sw_tmp_d);
			swap(ev1_id, ev2_id, sw_tmp_uint8);
		}
		// eigen vectors are identity matrix
		evecs[0][0] = identity_mat[0][ev0_id];
		evecs[1][0] = identity_mat[1][ev0_id];
		evecs[2][0] = identity_mat[2][ev0_id];
		evecs[0][1] = identity_mat[0][ev1_id];
		evecs[1][1] = identity_mat[1][ev1_id];
		evecs[2][1] = identity_mat[2][ev1_id];
		evecs[0][2] = identity_mat[0][ev2_id];
		evecs[1][2] = identity_mat[1][ev2_id];
		evecs[2][2] = identity_mat[2][ev2_id];
	}
	else
	{
		// cal eigenvalues
		const __Float_Type__ p = (mat_va[0] + mat_va[1] + mat_va[2]) / ffmat(3.0);
		__Float_Type__ matb_va[6];
		matb_va[0] = mat_va[0] - p;
		matb_va[1] = mat_va[1] - p;
		matb_va[2] = mat_va[2] - p;
		__Float_Type__ q = (__Float_Type__)sqrt(
			 (matb_va[0] * matb_va[0]
			+ matb_va[1] * matb_va[1]
			+ matb_va[2] * matb_va[2]
			+ tmp * ffmat(2.0)) / ffmat(6.0));
		matb_va[0] /= q;
		matb_va[1] /= q;
		matb_va[2] /= q;
		matb_va[3] = mat_va[3] / q;
		matb_va[4] = mat_va[4] / q;
		matb_va[5] = mat_va[5] / q;
		const __Float_Type__ detb = (
			  matb_va[0] * matb_va[1] * matb_va[2]
			+ matb_va[3] * matb_va[4] * matb_va[5] * ffmat(2.0)
			- matb_va[0] * matb_va[4] * matb_va[4]
			- matb_va[1] * matb_va[5] * matb_va[5]
			- matb_va[2] * matb_va[3] * matb_va[3]) * ffmat(0.5);
		__Float_Type__ phi;
		if (detb <= ffmat(-1.0))
			phi = PI / ffmat(3.0);
		else if (detb >= ffmat(1.0))
			phi = ffmat(0.0);
		else
			phi = acos(detb) / ffmat(3.0);
		ev[0] = p + ffmat(2.0) * q * cos(phi);
		ev[1] = p + ffmat(2.0) * q * cos(phi + ffmat(2.0) * PI / ffmat(3.0));
		ev[2] = p + ffmat(2.0) * q * cos(phi + ffmat(4.0) * PI / ffmat(3.0));

		// sort eigenvalues
		// ev[0] > ev[1] > ev[2]
		double sw_tmp_d;
		if (ev[0] < ev[1])
			swap(ev[0], ev[1], sw_tmp_d);
		if (ev[0] < ev[2])
			swap(ev[0], ev[2], sw_tmp_d);
		if (ev[1] < ev[2])
			swap(ev[1], ev[2], sw_tmp_d);

		// cal eigen vector
		const __Float_Type__ d_tol0 = (ev[0] - ev[2]) * ffmat(0.5) * DTol;
		const __Float_Type__ d_tol1 = (ev[0] + ev[1] + ev[2]) * ffmat(0.333333) * DTol;
		__Float_Type__ d_tol = DTol;
		if (d_tol < d_tol0)
			d_tol = d_tol0;
		if (d_tol < d_tol1)
			d_tol = d_tol1;

		__Float_Type__ evec_tmp0[3], evec_tmp1[3], evec_tmp2[3];
		if ((__Float_Type__)abs(ev[0] - ev[1]) > d_tol) // ev[0] != ev[1]
		{
			const __Float_Type__ diff12 = abs(ev[1] - ev[2]);
			if ((__Float_Type__)abs(ev[1] - ev[2]) > d_tol) // ev[0] != ev[1] != ev[2]
			{
				cal_evec(mat_va, ev[0], evec_tmp0);
				evecs[0][0] = evec_tmp0[0];
				evecs[1][0] = evec_tmp0[1];
				evecs[2][0] = evec_tmp0[2];
				cal_evec(mat_va, ev[1], evec_tmp1);
				evecs[0][1] = evec_tmp1[0];
				evecs[1][1] = evec_tmp1[1];
				evecs[2][1] = evec_tmp1[2];
				cross_prod(evec_tmp0, evec_tmp1, evec_tmp2);
				evecs[0][2] = evec_tmp2[0];
				evecs[1][2] = evec_tmp2[1];
				evecs[2][2] = evec_tmp2[2];
			}
			else // ev[0] != ev[1] == ev[2]
			{
				cal_evec(mat_va, ev[0], evec_tmp0);
				evecs[0][0] = evec_tmp0[0];
				evecs[1][0] = evec_tmp0[1];
				evecs[2][0] = evec_tmp0[2];
				cal_ortho_vec(evec_tmp0, evec_tmp1);
				evecs[0][1] = evec_tmp1[0];
				evecs[1][1] = evec_tmp1[1];
				evecs[2][1] = evec_tmp1[2];
				cross_prod(evec_tmp0, evec_tmp1, evec_tmp2);
				evecs[0][2] = evec_tmp2[0];
				evecs[1][2] = evec_tmp2[1];
				evecs[2][2] = evec_tmp2[2];
			}
		}
		else // ev[0] == ev[1]
		{
			if (abs(ev[1] - ev[2]) > d_tol) // ev[0] == ev[1] != ev[2]
			{
				cal_evec(mat_va, ev[2], evec_tmp0);
				evecs[0][2] = evec_tmp0[0];
				evecs[1][2] = evec_tmp0[1];
				evecs[2][2] = evec_tmp0[2];
				cal_ortho_vec(evec_tmp0, evec_tmp1);
				evecs[0][0] = evec_tmp1[0];
				evecs[1][0] = evec_tmp1[1];
				evecs[2][0] = evec_tmp1[2];
				cross_prod(evec_tmp0, evec_tmp1, evec_tmp2);
				evecs[0][1] = evec_tmp2[0];
				evecs[1][1] = evec_tmp2[1];
				evecs[2][1] = evec_tmp2[2];
			}
			else // ev[0] == ev[1] == ev[2]
			{
				// eigen values
				ev[0] = mat_va[0];
				ev[1] = mat_va[1];
				ev[2] = mat_va[2];
				// eigen vectors are identity matrix
				evecs[0][0] = ffmat(1.0);
				evecs[1][0] = ffmat(0.0);
				evecs[2][0] = ffmat(0.0);
				evecs[0][1] = ffmat(0.0);
				evecs[1][1] = ffmat(1.0);
				evecs[2][1] = ffmat(0.0);
				evecs[0][2] = ffmat(0.0);
				evecs[1][2] = ffmat(0.0);
				evecs[2][2] = ffmat(1.0);
			}
		}
	}
}

void rotate_sym_mat_to_eigen_mat(
	const __Float_Type__ mat[6],
	const __Float_Type__ evecs[3][3],
	__Float_Type__ evs[3])
{
	const double tmp[3][3] = {
		mat[0] * evecs[0][0] + mat[3] * evecs[1][0] + mat[5] * evecs[2][0], // 00
		mat[0] * evecs[0][1] + mat[3] * evecs[1][1] + mat[5] * evecs[2][1], // 01
		mat[0] * evecs[0][2] + mat[3] * evecs[1][2] + mat[5] * evecs[2][2], // 02
		mat[3] * evecs[0][0] + mat[1] * evecs[1][0] + mat[4] * evecs[2][0], // 10
		mat[3] * evecs[0][1] + mat[1] * evecs[1][1] + mat[4] * evecs[2][1], // 11
		mat[3] * evecs[0][2] + mat[1] * evecs[1][2] + mat[4] * evecs[2][2], // 12
		mat[5] * evecs[0][0] + mat[4] * evecs[1][0] + mat[2] * evecs[2][0], // 20
		mat[5] * evecs[0][1] + mat[4] * evecs[1][1] + mat[2] * evecs[2][1], // 21
		mat[5] * evecs[0][2] + mat[4] * evecs[1][2] + mat[2] * evecs[2][2]  // 22
	};
	evs[0] = evecs[0][0] * tmp[0][0] + evecs[1][0] * tmp[1][0] + evecs[2][0] * tmp[2][0];
	evs[1] = evecs[0][1] * tmp[0][1] + evecs[1][1] * tmp[1][1] + evecs[2][1] * tmp[2][1];
	evs[2] = evecs[0][2] * tmp[0][2] + evecs[1][2] * tmp[1][2] + evecs[2][2] * tmp[2][2];
}

void rotate_eigen_mat_to_sym_mat(
	const __Float_Type__ evs[3],
	const __Float_Type__ evecs[3][3],
	__Float_Type__ mat[6])
{
	const double tmp[3][3] = {
		evs[0] * evecs[0][0], evs[0] * evecs[1][0], evs[0] * evecs[2][0],
		evs[1] * evecs[0][1], evs[1] * evecs[1][1], evs[1] * evecs[2][1],
		evs[2] * evecs[0][2], evs[2] * evecs[1][2], evs[2] * evecs[2][2]
	};
	mat[0] = evecs[0][0] * tmp[0][0] + evecs[0][1] * tmp[1][0] + evecs[0][2] * tmp[2][0];
	mat[1] = evecs[1][0] * tmp[0][1] + evecs[1][1] * tmp[1][1] + evecs[1][2] * tmp[2][1];
	mat[2] = evecs[2][0] * tmp[0][2] + evecs[2][1] * tmp[1][2] + evecs[2][2] * tmp[2][2];
	mat[3] = evecs[0][0] * tmp[0][1] + evecs[0][1] * tmp[1][1] + evecs[0][2] * tmp[2][1];
	mat[4] = evecs[1][0] * tmp[0][2] + evecs[1][1] * tmp[1][2] + evecs[1][2] * tmp[2][2];
	mat[5] = evecs[0][0] * tmp[0][2] + evecs[0][1] * tmp[1][2] + evecs[0][2] * tmp[2][2];
}
