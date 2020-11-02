#include "MaterialModels_pcp.h"

#include "SymMatEigen.h"
#include "Tresca.h"

namespace MatModel
{
	using namespace MatModel_Internal;

#define swap_size_t(a, b) \
	(a) = (a) ^ (b); \
	(b) = (a) ^ (b); \
	(a) = (a) ^ (b)
	
	// return value:
	//     = 0 - elstic
	//     > 0 - plastic
	//     < 0 - convergence failure
	int tresca_integration_function(MaterialModel* _self, double dstrain[6])
	{
		Tresca& self = static_cast<Tresca&>(*_self);
		double* dstress = self.dstress;
		double* stress = self.stress;
		double* dstrain_e = self.dstrain_e;
		double* dstrain_p = self.dstrain_p;
		double(*De_mat)[6] = self.De_mat;
		double(*Dep_mat)[6] = self.Dep_mat;
		double e_stress[6];
		union
		{
			double stress_ev[6];
			struct { double s1, s2, s3; };
		};
		double stress_evec[3][3];

		matrix6x6_prod_vector6(De_mat, dstrain, dstress);
		vector6_add(stress, dstress, e_stress);
		vector6_copy(dstrain, dstrain_e);
		self.copy_De_to_Dep();

		// cal principle stress
		cal_sym_mat_eigen(e_stress, stress_ev, stress_evec);

		// criteria0
		if ((s1 - s3 - self.two_c) <= 0.0)
		{
			// remain elastic
			vector6_copy(e_stress, stress);
			vector6_zeros(dstrain_p);
			return 0;
		}

		double old_stress[6];
		vector6_copy(stress, old_stress);

		// correct stress back to yield surface
		double s_mean, lambda;
		double criteria1 = s1 + s3 - s2 - s2;
		double Dep_tmp = De_mat[0][0] - De_mat[0][1];
		if (criteria1 < -self.two_c) // return to line1
		{
			// stress
			s_mean = (s1 + s2 + s3) / 3.0;
			s1 = s_mean + 2.0 / 3.0 * self.cohesion;
			s2 = s1;
			s3 = s_mean - 4.0 / 3.0 * self.cohesion;
			// Dep
			Dep_tmp *= 1.0 / 6.0;
			Dep_mat[0][0] -= Dep_tmp;
			Dep_mat[0][1] -= Dep_tmp;
			Dep_mat[0][2] -= -2.0 * Dep_tmp;
			Dep_mat[1][0] -= Dep_tmp;
			Dep_mat[1][1] -= Dep_tmp;
			Dep_mat[1][2] -= -2.0 * Dep_tmp;
			Dep_mat[2][0] -= -2.0 * Dep_tmp;
			Dep_mat[2][1] -= -2.0 * Dep_tmp;
			Dep_mat[2][2] -= 4.0 * Dep_tmp;
		}
		else if (criteria1 > self.two_c) // return to line2
		{
			// stress
			s_mean = (s1 + s2 + s3) / 3.0;
			s1 = s_mean + 4.0 / 3.0 * self.cohesion;
			s2 = s_mean - 2.0 / 3.0 * self.cohesion;
			s3 = s2;
			// Dep
			Dep_tmp *= 1.0 / 6.0;
			Dep_mat[0][0] -= 4.0 * Dep_tmp;
			Dep_mat[0][1] -= -2.0 * Dep_tmp;
			Dep_mat[0][2] -= -2.0 * Dep_tmp;
			Dep_mat[1][0] -= -2.0 * Dep_tmp;
			Dep_mat[1][1] -= Dep_tmp;
			Dep_mat[1][2] -= Dep_tmp;
			Dep_mat[2][0] -= -2.0 * Dep_tmp;
			Dep_mat[2][1] -= Dep_tmp;
			Dep_mat[2][2] -= Dep_tmp;
		}
		else // return to plane
		{
			// stress
			lambda = (s1 - s3 - self.two_c) / (Dep_tmp + Dep_tmp);
			s1 -= (De_mat[0][0] - De_mat[0][2]) * lambda;
			s2 -= (De_mat[1][0] - De_mat[1][2]) * lambda;
			s3 -= (De_mat[2][0] - De_mat[2][2]) * lambda;
			// Dep
			Dep_tmp *= 0.5;
			Dep_mat[0][0] -= Dep_tmp;
			Dep_mat[0][2] += Dep_tmp;
			Dep_mat[2][0] += Dep_tmp;
			Dep_mat[2][2] -= Dep_tmp;
		}

		// transform matrix from principle space to current space
		double tmat6x6[6][6];
		tmat6x6[0][0] = stress_evec[0][0] * stress_evec[0][0];
		tmat6x6[0][1] = stress_evec[0][1] * stress_evec[0][1];
		tmat6x6[0][2] = stress_evec[0][2] * stress_evec[0][2];
		tmat6x6[0][3] = 2.0 * stress_evec[0][0] * stress_evec[0][1];
		tmat6x6[0][4] = 2.0 * stress_evec[0][1] * stress_evec[0][2];
		tmat6x6[0][5] = 2.0 * stress_evec[0][2] * stress_evec[0][0];
		tmat6x6[1][0] = stress_evec[1][0] * stress_evec[1][0];
		tmat6x6[1][1] = stress_evec[1][1] * stress_evec[1][1];
		tmat6x6[1][2] = stress_evec[1][2] * stress_evec[1][2];
		tmat6x6[1][3] = 2.0 * stress_evec[1][0] * stress_evec[1][1];
		tmat6x6[1][4] = 2.0 * stress_evec[1][1] * stress_evec[1][2];
		tmat6x6[1][5] = 2.0 * stress_evec[1][2] * stress_evec[1][0];
		tmat6x6[2][0] = stress_evec[2][0] * stress_evec[2][0];
		tmat6x6[2][1] = stress_evec[2][1] * stress_evec[2][1];
		tmat6x6[2][2] = stress_evec[2][2] * stress_evec[2][2];
		tmat6x6[2][3] = 2.0 * stress_evec[2][0] * stress_evec[2][1];
		tmat6x6[2][4] = 2.0 * stress_evec[2][1] * stress_evec[2][2];
		tmat6x6[2][5] = 2.0 * stress_evec[2][2] * stress_evec[2][0];
		tmat6x6[3][0] = stress_evec[0][0] * stress_evec[1][0];
		tmat6x6[3][1] = stress_evec[0][1] * stress_evec[1][1];
		tmat6x6[3][2] = stress_evec[0][2] * stress_evec[1][2];
		tmat6x6[3][3] = stress_evec[0][0] * stress_evec[1][1] + stress_evec[0][1] * stress_evec[1][0];
		tmat6x6[3][4] = stress_evec[0][1] * stress_evec[1][2] + stress_evec[0][2] * stress_evec[1][1];
		tmat6x6[3][5] = stress_evec[0][0] * stress_evec[1][2] + stress_evec[0][2] * stress_evec[1][0];
		tmat6x6[4][0] = stress_evec[1][0] * stress_evec[2][0];
		tmat6x6[4][1] = stress_evec[1][1] * stress_evec[2][1];
		tmat6x6[4][2] = stress_evec[1][2] * stress_evec[2][2];
		tmat6x6[4][3] = stress_evec[1][0] * stress_evec[2][1] + stress_evec[1][1] * stress_evec[2][0];
		tmat6x6[4][4] = stress_evec[1][1] * stress_evec[2][2] + stress_evec[1][2] * stress_evec[2][1];
		tmat6x6[4][5] = stress_evec[1][0] * stress_evec[2][2] + stress_evec[1][2] * stress_evec[2][0];
		tmat6x6[5][0] = stress_evec[2][0] * stress_evec[0][0];
		tmat6x6[5][1] = stress_evec[2][1] * stress_evec[0][1];
		tmat6x6[5][2] = stress_evec[2][2] * stress_evec[0][2];
		tmat6x6[5][3] = stress_evec[2][0] * stress_evec[0][1] + stress_evec[2][1] * stress_evec[0][0];
		tmat6x6[5][4] = stress_evec[2][1] * stress_evec[0][2] + stress_evec[2][2] * stress_evec[0][1];
		tmat6x6[5][5] = stress_evec[2][0] * stress_evec[0][2] + stress_evec[2][2] * stress_evec[0][0];
		// stress
		stress_ev[3] = 0.0;
		stress_ev[4] = 0.0;
		stress_ev[5] = 0.0;
		matrix6x6_prod_vector6(tmat6x6, stress_ev, stress);

		// stress increment
		vector6_subtract(stress, old_stress, dstress);

		// plastic strain
		double stress_diff[6];
		vector6_subtract(e_stress, stress, stress_diff);
		double coef1 = 1.0 / self.E;
		double coef2 = -self.niu / self.E;
		double coef3 = (1.0 + self.niu) / self.E;
		dstrain_p[0] = coef1 * stress_diff[0] + coef2 * stress_diff[1] + coef2 * stress_diff[2];
		dstrain_p[1] = coef2 * stress_diff[0] + coef1 * stress_diff[1] + coef2 * stress_diff[2];
		dstrain_p[2] = coef2 * stress_diff[0] + coef2 * stress_diff[1] + coef1 * stress_diff[2];
		dstrain_p[3] = coef3 * stress_diff[3];
		dstrain_p[4] = coef3 * stress_diff[4];
		dstrain_p[5] = coef3 * stress_diff[5];
		vector6_subtract(dstrain_p, dstrain_e);

		// currently no Dep to save time
		self.copy_De_to_Dep();
		//// transpose of transformation matrix
		//double tmat6x6T[6][6]; // transpose of tmat6x6
		//tmat6x6T[0][0] = tmat3x3(0, 0) * tmat3x3(0, 0);
		//tmat6x6T[0][1] = tmat3x3(1, 0) * tmat3x3(1, 0);
		//tmat6x6T[0][2] = tmat3x3(2, 0) * tmat3x3(2, 0);
		//tmat6x6T[0][3] = 2.0 * tmat3x3(0, 0) * tmat3x3(1, 0);
		//tmat6x6T[0][4] = 2.0 * tmat3x3(1, 0) * tmat3x3(2, 0);
		//tmat6x6T[0][5] = 2.0 * tmat3x3(2, 0) * tmat3x3(0, 0);
		//tmat6x6T[1][0] = tmat3x3(0, 1) * tmat3x3(0, 1);
		//tmat6x6T[1][1] = tmat3x3(1, 1) * tmat3x3(1, 1);
		//tmat6x6T[1][2] = tmat3x3(2, 1) * tmat3x3(2, 1);
		//tmat6x6T[1][3] = 2.0 * tmat3x3(0, 1) * tmat3x3(1, 1);
		//tmat6x6T[1][4] = 2.0 * tmat3x3(1, 1) * tmat3x3(2, 1);
		//tmat6x6T[1][5] = 2.0 * tmat3x3(2, 1) * tmat3x3(0, 1);
		//tmat6x6T[2][0] = tmat3x3(0, 2) * tmat3x3(0, 2);
		//tmat6x6T[2][1] = tmat3x3(1, 2) * tmat3x3(1, 2);
		//tmat6x6T[2][2] = tmat3x3(2, 2) * tmat3x3(2, 2);
		//tmat6x6T[2][3] = 2.0 * tmat3x3(0, 2) * tmat3x3(1, 2);
		//tmat6x6T[2][4] = 2.0 * tmat3x3(1, 2) * tmat3x3(2, 2);
		//tmat6x6T[2][5] = 2.0 * tmat3x3(2, 2) * tmat3x3(0, 2);
		//tmat6x6T[3][0] = tmat3x3(0, 0) * tmat3x3(0, 1);
		//tmat6x6T[3][1] = tmat3x3(1, 0) * tmat3x3(1, 1);
		//tmat6x6T[3][2] = tmat3x3(2, 0) * tmat3x3(2, 1);
		//tmat6x6T[3][3] = tmat3x3(0, 0) * tmat3x3(1, 1) + tmat3x3(1, 0) * tmat3x3(0, 1);
		//tmat6x6T[3][4] = tmat3x3(1, 0) * tmat3x3(2, 1) + tmat3x3(2, 0) * tmat3x3(1, 1);
		//tmat6x6T[3][5] = tmat3x3(0, 0) * tmat3x3(2, 1) + tmat3x3(2, 0) * tmat3x3(0, 1);
		//tmat6x6T[4][0] = tmat3x3(0, 1) * tmat3x3(0, 2);
		//tmat6x6T[4][1] = tmat3x3(1, 1) * tmat3x3(1, 2);
		//tmat6x6T[4][2] = tmat3x3(2, 1) * tmat3x3(2, 2);
		//tmat6x6T[4][3] = tmat3x3(0, 1) * tmat3x3(1, 2) + tmat3x3(1, 1) * tmat3x3(0, 2);
		//tmat6x6T[4][4] = tmat3x3(1, 1) * tmat3x3(2, 2) + tmat3x3(2, 1) * tmat3x3(1, 2);
		//tmat6x6T[4][5] = tmat3x3(0, 1) * tmat3x3(2, 2) + tmat3x3(2, 1) * tmat3x3(0, 2);
		//tmat6x6T[5][0] = tmat3x3(0, 2) * tmat3x3(0, 0);
		//tmat6x6T[5][1] = tmat3x3(1, 2) * tmat3x3(1, 0);
		//tmat6x6T[5][2] = tmat3x3(2, 2) * tmat3x3(2, 0);
		//tmat6x6T[5][3] = tmat3x3(0, 2) * tmat3x3(1, 0) + tmat3x3(1, 2) * tmat3x3(0, 0);
		//tmat6x6T[5][4] = tmat3x3(1, 2) * tmat3x3(2, 0) + tmat3x3(2, 2) * tmat3x3(1, 0);
		//tmat6x6T[5][5] = tmat3x3(0, 2) * tmat3x3(2, 0) + tmat3x3(2, 2) * tmat3x3(0, 0);
		//// trasform Dep from principle space to current space
		//double tmp_mat[6][6];
		//matrix6x6_prod_matrix6x6(Dep_mat, tmat6x6T, tmp_mat);
		//matrix6x6_prod_matrix6x6(tmat6x6, tmp_mat, Dep_mat);

		return 0;
	}
}
