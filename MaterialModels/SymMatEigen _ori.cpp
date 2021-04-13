#include "MaterialModels_pcp.h"

#include <cmath>

#include "SymMatEigen.h"

#define DTol (1.0e-10)
#define PI (3.14159265359)

#define cross_prod(a, b, c) \
	c[0] = a[1] * b[2] - a[2] * b[1]; \
	c[1] = a[2] * b[0] - a[0] * b[2]; \
	c[2] = a[0] * b[1] - a[1] * b[0]

#define swap_ll(a, b)    \
		(a) = (a) ^ (b); \
		(b) = (a) ^ (b); \
		(a) = (a) ^ (b)

static void cal_evec(const double mat_va[6], double ev, double evec[3])
{
	double r0[3] = { mat_va[0] - ev, mat_va[3], mat_va[5] };
	double r1[3] = { mat_va[3], mat_va[1] - ev, mat_va[4] };
	double r2[3] = { mat_va[5], mat_va[4], mat_va[2] - ev };
	double r0xr1[3], r1xr2[3], r2xr0[3];
	cross_prod(r0, r1, r0xr1);
	double r0xr1_len2 = r0xr1[0] * r0xr1[0] + r0xr1[1] * r0xr1[1] + r0xr1[2] * r0xr1[2];
	cross_prod(r1, r2, r1xr2);
	double r1xr2_len2 = r1xr2[0] * r1xr2[0] + r1xr2[1] * r1xr2[1] + r1xr2[2] * r1xr2[2];
	cross_prod(r2, r0, r2xr0);
	double r2xr0_len2 = r2xr0[0] * r2xr0[0] + r2xr0[1] * r2xr0[1] + r2xr0[2] * r2xr0[2];
	double max_len2 = r0xr1_len2;
	double* max_vec = r0xr1;
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
	max_len2 = sqrt(max_len2);
	evec[0] = max_vec[0] / max_len2;
	evec[1] = max_vec[1] / max_len2;
	evec[2] = max_vec[2] / max_len2;
}

static void cal_ortho_vec(const double vec0[3], double vec1[3])
{
	double inv_len;
	if (abs(vec0[0]) > abs(vec0[1]))
	{
		inv_len = 1.0 / sqrt(vec0[0] * vec0[0] + vec0[2] * vec0[2]);
		vec1[0] = -vec0[2] * inv_len;
		vec1[1] = 0.0;
		vec1[2] = vec0[0] * inv_len;
	}
	else
	{
		inv_len = 1.0 / sqrt(vec0[1] * vec0[1] + vec0[2] * vec0[2]);
		vec1[0] = 0.0;
		vec1[1] = vec0[2] * inv_len;
		vec1[2] = -vec0[1] * inv_len;
	}
}

static const double identity_mat[3][3] = {
	1.0, 0.0, 0.0,
	0.0, 1.0, 0.0,
	0.0, 0.0, 1.0
};

// mat is sysmetric
// mat_va[0] = mat[0][0]
// mat_va[1] = mat[1][1]
// mat_va[2] = mat[2][2]
// mat_va[3] = mat[0][1]
// mat_va[4] = mat[1][2]
// mat_va[5] = mat[2][0]
void cal_sym_mat_eigen(const double mat_va[6], double ev[3], double evecs[3][3])
{
#define ev_ll(id) (((long long *)(ev))[(id)])
	double tmp = mat_va[3] * mat_va[3]
		+ mat_va[4] * mat_va[4]
		+ mat_va[5] * mat_va[5];
	if (tmp < DTol * DTol)
	{
		// eigen values
		ev[0] = mat_va[0];
		ev[1] = mat_va[1];
		ev[2] = mat_va[2];
		// sort eigenvalues
		// ev[0] > ev[1] > ev[2]
		size_t ev0_id = 0;
		size_t ev1_id = 1;
		size_t ev2_id = 2;
		if (ev[0] < ev[1])
		{
			swap_ll(ev_ll(0), ev_ll(1));
			swap_ll(ev0_id, ev1_id);
		}
		if (ev[0] < ev[2])
		{
			swap_ll(ev_ll(0), ev_ll(2));
			swap_ll(ev0_id, ev2_id);
		}
		if (ev[1] < ev[2])
		{
			swap_ll(ev_ll(1), ev_ll(2));
			swap_ll(ev1_id, ev2_id);
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
		double p = (mat_va[0] + mat_va[1] + mat_va[2]) / 3.0;
		double matb_va[6];
		matb_va[0] = mat_va[0] - p;
		matb_va[1] = mat_va[1] - p;
		matb_va[2] = mat_va[2] - p;
		double q = sqrt((matb_va[0] * matb_va[0] + matb_va[1] * matb_va[1]
					   + matb_va[2] * matb_va[2] + 2.0 * tmp) / 6.0);
		matb_va[0] /= q;
		matb_va[1] /= q;
		matb_va[2] /= q;
		matb_va[3] = mat_va[3] / q;
		matb_va[4] = mat_va[4] / q;
		matb_va[5] = mat_va[5] / q;
		double detb = (matb_va[0] * matb_va[1] * matb_va[2]
			+ 2.0 * matb_va[3] * matb_va[4] * matb_va[5]
			- matb_va[0] * matb_va[4] * matb_va[4]
			- matb_va[1] * matb_va[5] * matb_va[5]
			- matb_va[2] * matb_va[3] * matb_va[3]) * 0.5;
		double phi;
		if (detb <= -1.0)
			phi = PI / 3.0;
		else if (detb >= 1.0)
			phi = 0.0;
		else
			phi = acos(detb) / 3.0;
		ev[0] = p + 2.0 * q * cos(phi);
		ev[1] = p + 2.0 * q * cos(phi + 2.0 * PI / 3.0);
		ev[2] = p + 2.0 * q * cos(phi + 4.0 * PI / 3.0);

		// sort eigenvalues
		// ev[0] > ev[1] > ev[2]
		if (ev[0] < ev[1]) { swap_ll(ev_ll(0), ev_ll(1)); }
		if (ev[0] < ev[2]) { swap_ll(ev_ll(0), ev_ll(2)); }
		if (ev[1] < ev[2]) { swap_ll(ev_ll(1), ev_ll(2)); }

		// cal eigen vector
		double d_tol0 = (ev[0] - ev[2]) * 0.5 * 1.0e-10;
		double d_tol1 = (ev[0] + ev[1] + ev[2]) * 0.333333 * 1.0e-10;
		double d_tol = 1.0e-10;
		if (d_tol < d_tol0)
			d_tol = d_tol0;
		if (d_tol < d_tol1)
			d_tol = d_tol1;

		double evec_tmp0[3], evec_tmp1[3], evec_tmp2[3];
		if (abs(ev[0] - ev[1]) > d_tol) // ev[0] != ev[1]
		{
			double diff12 = abs(ev[1] - ev[2]);
			if (abs(ev[1] - ev[2]) > d_tol) // ev[0] != ev[1] != ev[2]
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
				evecs[0][0] = 1.0;
				evecs[1][0] = 0.0;
				evecs[2][0] = 0.0;
				evecs[0][1] = 0.0;
				evecs[1][1] = 1.0;
				evecs[2][1] = 0.0;
				evecs[0][2] = 0.0;
				evecs[1][2] = 0.0;
				evecs[2][2] = 1.0;
			}
		}
	}
#undef ev_ll
}

#undef DTol
#undef PI
#undef cross_prod
#undef swap_ll
