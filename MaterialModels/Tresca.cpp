#include "MaterialModels_pcp.h"

#include <Eigen/Dense>
#include "Tresca.h"

namespace MatModel
{
	using namespace MatModel_Internal;

	template <typename Item>
	inline void swap(Item& i1, Item& i2)
	{
		Item tmp = i1;
		i1 = i2;
		i2 = tmp;
	}

	// return value:
	//     = 0 - elstic
	//     > 0 - plastic
	//     < 0 - convergence failure
	int tresca_integration_function(MaterialModel* _self, double dstrain[6])
	{
		Tresca& self = *static_cast<Tresca *>(_self);

		double ori_stress[6];
		double *dstress = self.dstress;
		double *stress = self.stress;
		double *dstrain_e = self.dstrain_e;
		double *dstrain_p = self.dstrain_p;
		double (*De_mat)[6] = self.De_mat;
		double (*Dep_mat)[6] = self.Dep_mat;
		vector6_copy(stress, ori_stress);
		matrix6x6_prod_vector6(De_mat, dstrain, dstress);
		vector6_add(dstress, stress);
		vector6_copy(dstrain, dstrain_e);
		vector6_zeros(dstrain_p);
		matrix6x6_copy(De_mat, Dep_mat);

		// cal principle stress
		Eigen::Matrix3d stress_mat;
		stress_mat << stress[0], stress[3], stress[5],
					  stress[3], stress[1], stress[4],
					  stress[5], stress[4], stress[2];
		Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(stress_mat);
		//if (eigen_solver.info() != Eigen::Success)
		//	assert(0);
		const Eigen::Vector3d &prin_stress = eigen_solver.eigenvalues();

		// sorts principle stress
		union
		{
			struct { double s1, s2, s3; };
			struct { size_t s1_ui, s2_ui, s3_ui; };
		};
		s1 = prin_stress[0];
		s2 = prin_stress[1];
		s3 = prin_stress[2];
		size_t s1_id, s2_id, s3_id;
		s1_id = 0;
		s2_id = 1;
		s3_id = 2;
		if (s1 < s2)
		{
			swap(s1, s2);
			swap(s1_id, s2_id);
		}
		if (s1 < s3)
		{
			swap(s1, s3);
			swap(s1_id, s3_id);
		}
		if (s2 < s3)
		{
			swap(s2, s3);
			swap(s2_id, s3_id);
		}

		// criteria0
		if ((s1 - s3 - self.two_c) < 0.0)
			return 0;

		double s_mean, lambda;
		double criteria1 = s1 + s3 - s2 - s2;
		double Dep_tmp = De_mat[0][0] - De_mat[0][1];
		if (criteria1 < -self.two_c) // return to line1
		{
			s_mean = (s1 + s2 + s3) / 3.0;
			s1 = s_mean + 2.0 / 3.0 * self.cohesion;
			s2 = s1;
			s3 = s_mean - 4.0 / 3.0 * self.cohesion;
			// Dep
			Dep_tmp *= 1.0 / 6.0;
			Dep_mat[s1_id][s1_id] -= Dep_tmp;
			Dep_mat[s1_id][s2_id] -= Dep_tmp;
			Dep_mat[s1_id][s3_id] -= -2.0 * Dep_tmp;
			Dep_mat[s2_id][s1_id] -= Dep_tmp;
			Dep_mat[s2_id][s2_id] -= Dep_tmp;
			Dep_mat[s2_id][s3_id] -= -2.0 * Dep_tmp;
			Dep_mat[s3_id][s1_id] -= -2.0 * Dep_tmp;
			Dep_mat[s3_id][s2_id] -= -2.0 * Dep_tmp;
			Dep_mat[s3_id][s3_id] -= 4.0 * Dep_tmp;
		}
		else if (criteria1 > self.two_c) // return to line2
		{
			s_mean = (s1 + s2 + s3) / 3.0;
			s1 = s_mean + 4.0 / 3.0 * self.cohesion;
			s2 = s_mean - 2.0 / 3.0 * self.cohesion;
			s3 = s2;
			// Dep
			Dep_tmp *= 1.0 / 6.0;
			Dep_mat[s1_id][s1_id] -= 4.0 * Dep_tmp;
			Dep_mat[s1_id][s2_id] -= -2.0 * Dep_tmp;
			Dep_mat[s1_id][s3_id] -= -2.0 * Dep_tmp;
			Dep_mat[s2_id][s1_id] -= -2.0 * Dep_tmp;
			Dep_mat[s2_id][s2_id] -= Dep_tmp;
			Dep_mat[s2_id][s3_id] -= Dep_tmp;
			Dep_mat[s3_id][s1_id] -= -2.0 * Dep_tmp;
			Dep_mat[s3_id][s2_id] -= Dep_tmp;
			Dep_mat[s3_id][s3_id] -= Dep_tmp;
		}
		else // return to plane
		{
			lambda = (self.two_c - (s1 - s3)) / (2.0 * (De_mat[0][1] - De_mat[0][0]));
			s1 -= (De_mat[0][0] - De_mat[0][2]) * lambda;
			s2 -= (De_mat[1][0] - De_mat[1][2]) * lambda;
			s3 -= (De_mat[2][0] - De_mat[2][2]) * lambda;
			// update Dep
			Dep_tmp *= 0.5;
			Dep_mat[s1_id][s1_id] -= Dep_tmp;
			Dep_mat[s1_id][s3_id] += Dep_tmp;
			Dep_mat[s3_id][s1_id] += Dep_tmp;
			Dep_mat[s3_id][s3_id] -= Dep_tmp;
		}

		// transfer back to global coordinates
		//Eigen::Matrix3d cor_pstress_mat;
		//cor_pstress_mat.setZero();
		//cor_pstress_mat(s1_id, s1_id) = s1;
		//cor_pstress_mat(s2_id, s2_id) = s2;
		//cor_pstress_mat(s3_id, s3_id) = s3;
		//Eigen::Matrix3d cor_stress = eigen_vector
		//	* cor_pstress_mat * eigen_vector.transpose();
		//std::cout << cor_stress << "\n";

		const Eigen::Matrix3d& tmat3x3 = eigen_solver.eigenvectors();
		double tmat6x6[6][6];
		tmat6x6[0][0] = tmat3x3(0, 0) * tmat3x3(0, 0);
		tmat6x6[0][1] = tmat3x3(0, 1) * tmat3x3(0, 1);
		tmat6x6[0][2] = tmat3x3(0, 2) * tmat3x3(0, 2);
		tmat6x6[0][3] = 2.0 * tmat3x3(0, 0) * tmat3x3(0, 1);
		tmat6x6[0][4] = 2.0 * tmat3x3(0, 1) * tmat3x3(0, 2);
		tmat6x6[0][5] = 2.0 * tmat3x3(0, 2) * tmat3x3(0, 0);
		tmat6x6[1][0] = tmat3x3(1, 0) * tmat3x3(1, 0);
		tmat6x6[1][1] = tmat3x3(1, 1) * tmat3x3(1, 1);
		tmat6x6[1][2] = tmat3x3(1, 2) * tmat3x3(1, 2);
		tmat6x6[1][3] = 2.0 * tmat3x3(1, 0) * tmat3x3(1, 1);
		tmat6x6[1][4] = 2.0 * tmat3x3(1, 1) * tmat3x3(1, 2);
		tmat6x6[1][5] = 2.0 * tmat3x3(1, 2) * tmat3x3(1, 0);
		tmat6x6[2][0] = tmat3x3(2, 0) * tmat3x3(2, 0);
		tmat6x6[2][1] = tmat3x3(2, 1) * tmat3x3(2, 1);
		tmat6x6[2][2] = tmat3x3(2, 2) * tmat3x3(2, 2);
		tmat6x6[2][3] = 2.0 * tmat3x3(2, 0) * tmat3x3(2, 1);
		tmat6x6[2][4] = 2.0 * tmat3x3(2, 1) * tmat3x3(2, 2);
		tmat6x6[2][5] = 2.0 * tmat3x3(2, 2) * tmat3x3(2, 0);
		tmat6x6[3][0] = tmat3x3(0, 0) * tmat3x3(1, 0);
		tmat6x6[3][1] = tmat3x3(0, 1) * tmat3x3(1, 1);
		tmat6x6[3][2] = tmat3x3(0, 2) * tmat3x3(1, 2);
		tmat6x6[3][3] = tmat3x3(0, 0) * tmat3x3(1, 1) + tmat3x3(0, 1) * tmat3x3(1, 0);
		tmat6x6[3][4] = tmat3x3(0, 1) * tmat3x3(1, 2) + tmat3x3(0, 2) * tmat3x3(1, 1);
		tmat6x6[3][5] = tmat3x3(0, 0) * tmat3x3(1, 2) + tmat3x3(0, 2) * tmat3x3(1, 0);
		tmat6x6[4][0] = tmat3x3(1, 0) * tmat3x3(2, 0);
		tmat6x6[4][1] = tmat3x3(1, 1) * tmat3x3(2, 1);
		tmat6x6[4][2] = tmat3x3(1, 2) * tmat3x3(2, 2);
		tmat6x6[4][3] = tmat3x3(1, 0) * tmat3x3(2, 1) + tmat3x3(1, 1) * tmat3x3(2, 0);
		tmat6x6[4][4] = tmat3x3(1, 1) * tmat3x3(2, 2) + tmat3x3(1, 2) * tmat3x3(2, 1);
		tmat6x6[4][5] = tmat3x3(1, 0) * tmat3x3(2, 2) + tmat3x3(1, 2) * tmat3x3(2, 0);
		tmat6x6[5][0] = tmat3x3(2, 0) * tmat3x3(0, 0);
		tmat6x6[5][1] = tmat3x3(2, 1) * tmat3x3(0, 1);
		tmat6x6[5][2] = tmat3x3(2, 2) * tmat3x3(0, 2);
		tmat6x6[5][3] = tmat3x3(2, 0) * tmat3x3(0, 1) + tmat3x3(2, 1) * tmat3x3(0, 0);
		tmat6x6[5][4] = tmat3x3(2, 1) * tmat3x3(0, 2) + tmat3x3(2, 2) * tmat3x3(0, 1);
		tmat6x6[5][5] = tmat3x3(2, 0) * tmat3x3(0, 2) + tmat3x3(2, 2) * tmat3x3(0, 0);
		double pstress_cor[6];
		pstress_cor[s1_id] = s1;
		pstress_cor[s2_id] = s2;
		pstress_cor[s3_id] = s3;
		pstress_cor[3] = 0.0;
		pstress_cor[4] = 0.0;
		pstress_cor[5] = 0.0;
		double stress_cor[6];
		matrix6x6_prod_vector6(tmat6x6, pstress_cor, stress_cor);
		//std::cout << stress_cor[0] << ", " << stress_cor[1] << ", "
		//		  << stress_cor[2] << ", " << stress_cor[3] << ", "
		//		  << stress_cor[4] << ", " << stress_cor[5] << "\n";

		double tmat6x6T[6][6]; // transpose of tmat6x6
		tmat6x6T[0][0] = tmat3x3(0, 0) * tmat3x3(0, 0);
		tmat6x6T[0][1] = tmat3x3(1, 0) * tmat3x3(1, 0);
		tmat6x6T[0][2] = tmat3x3(2, 0) * tmat3x3(2, 0);
		tmat6x6T[0][3] = 2.0 * tmat3x3(0, 0) * tmat3x3(1, 0);
		tmat6x6T[0][4] = 2.0 * tmat3x3(1, 0) * tmat3x3(2, 0);
		tmat6x6T[0][5] = 2.0 * tmat3x3(2, 0) * tmat3x3(0, 0);
		tmat6x6T[1][0] = tmat3x3(0, 1) * tmat3x3(0, 1);
		tmat6x6T[1][1] = tmat3x3(1, 1) * tmat3x3(1, 1);
		tmat6x6T[1][2] = tmat3x3(2, 1) * tmat3x3(2, 1);
		tmat6x6T[1][3] = 2.0 * tmat3x3(0, 1) * tmat3x3(1, 1);
		tmat6x6T[1][4] = 2.0 * tmat3x3(1, 1) * tmat3x3(2, 1);
		tmat6x6T[1][5] = 2.0 * tmat3x3(2, 1) * tmat3x3(0, 1);
		tmat6x6T[2][0] = tmat3x3(0, 2) * tmat3x3(0, 2);
		tmat6x6T[2][1] = tmat3x3(1, 2) * tmat3x3(1, 2);
		tmat6x6T[2][2] = tmat3x3(2, 2) * tmat3x3(2, 2);
		tmat6x6T[2][3] = 2.0 * tmat3x3(0, 2) * tmat3x3(1, 2);
		tmat6x6T[2][4] = 2.0 * tmat3x3(1, 2) * tmat3x3(2, 2);
		tmat6x6T[2][5] = 2.0 * tmat3x3(2, 2) * tmat3x3(0, 2);
		tmat6x6T[3][0] = tmat3x3(0, 0) * tmat3x3(0, 1);
		tmat6x6T[3][1] = tmat3x3(1, 0) * tmat3x3(1, 1);
		tmat6x6T[3][2] = tmat3x3(2, 0) * tmat3x3(2, 1);
		tmat6x6T[3][3] = tmat3x3(0, 0) * tmat3x3(1, 1) + tmat3x3(1, 0) * tmat3x3(0, 1);
		tmat6x6T[3][4] = tmat3x3(1, 0) * tmat3x3(2, 1) + tmat3x3(2, 0) * tmat3x3(1, 1);
		tmat6x6T[3][5] = tmat3x3(0, 0) * tmat3x3(2, 1) + tmat3x3(2, 0) * tmat3x3(0, 1);
		tmat6x6T[4][0] = tmat3x3(0, 1) * tmat3x3(0, 2);
		tmat6x6T[4][1] = tmat3x3(1, 1) * tmat3x3(1, 2);
		tmat6x6T[4][2] = tmat3x3(2, 1) * tmat3x3(2, 2);
		tmat6x6T[4][3] = tmat3x3(0, 1) * tmat3x3(1, 2) + tmat3x3(1, 1) * tmat3x3(0, 2);
		tmat6x6T[4][4] = tmat3x3(1, 1) * tmat3x3(2, 2) + tmat3x3(2, 1) * tmat3x3(1, 2);
		tmat6x6T[4][5] = tmat3x3(0, 1) * tmat3x3(2, 2) + tmat3x3(2, 1) * tmat3x3(0, 2);
		tmat6x6T[5][0] = tmat3x3(0, 2) * tmat3x3(0, 0);
		tmat6x6T[5][1] = tmat3x3(1, 2) * tmat3x3(1, 0);
		tmat6x6T[5][2] = tmat3x3(2, 2) * tmat3x3(2, 0);
		tmat6x6T[5][3] = tmat3x3(0, 2) * tmat3x3(1, 0) + tmat3x3(1, 2) * tmat3x3(0, 0);
		tmat6x6T[5][4] = tmat3x3(1, 2) * tmat3x3(2, 0) + tmat3x3(2, 2) * tmat3x3(1, 0);
		tmat6x6T[5][5] = tmat3x3(0, 2) * tmat3x3(2, 0) + tmat3x3(2, 2) * tmat3x3(0, 0);
		//double tmat6x6_prod[6][6];
		//for (size_t i = 0; i < 6; i++)
		//{
		//	for (size_t j = 0; j < 6; j++)
		//	{
		//		tmat6x6_prod[i][j] = 0.0;
		//		for (size_t k = 0; k < 6; k++)
		//			tmat6x6_prod[i][j] += tmat6x6[i][k] * tmat6x6T[k][j];
		//		std::cout << tmat6x6_prod[i][j] << " ";
		//	}
		//	std::cout << "\n";
		//}

		// dep = D^-1 * (s_p - s_c)
		//Eigen::Matrix<double, 6, 6> De;
		//De << De_mat[0][0], De_mat[0][1], De_mat[0][2], De_mat[0][3], De_mat[0][4], De_mat[0][5],
		//	  De_mat[1][0], De_mat[1][1], De_mat[1][2], De_mat[1][3], De_mat[1][4], De_mat[1][5],
		//	  De_mat[2][0], De_mat[2][1], De_mat[2][2], De_mat[2][3], De_mat[2][4], De_mat[2][5], 
		//	  De_mat[3][0], De_mat[3][1], De_mat[3][2], De_mat[3][3], De_mat[3][4], De_mat[3][5], 
		//	  De_mat[4][0], De_mat[4][1], De_mat[4][2], De_mat[4][3], De_mat[4][4], De_mat[4][5], 
		//	  De_mat[5][0], De_mat[5][1], De_mat[5][2], De_mat[5][3], De_mat[5][4], De_mat[5][5];
		//Eigen::Matrix<double, 6, 1> s_diff;
		//s_diff[0] = stress[0] - cor_stress(0, 0);
		//s_diff[1] = stress[1] - cor_stress(1, 1);
		//s_diff[2] = stress[2] - cor_stress(2, 2);
		//s_diff[3] = stress[3] - cor_stress(0, 1);
		//s_diff[4] = stress[4] - cor_stress(1, 2);
		//s_diff[5] = stress[5] - cor_stress(0, 2);
		//Eigen::Matrix<double, 6, 1> dstrain_p_vec = De.partialPivLu().solve(s_diff);
		double stress_diff[6];
		stress_diff[0] = stress[0] - stress_cor[0];
		stress_diff[1] = stress[1] - stress_cor[1];
		stress_diff[2] = stress[2] - stress_cor[2];
		stress_diff[3] = stress[3] - stress_cor[3];
		stress_diff[4] = stress[4] - stress_cor[4];
		stress_diff[5] = stress[5] - stress_cor[5];
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

		stress[0] = stress_cor[0];
		stress[1] = stress_cor[1];
		stress[2] = stress_cor[2];
		stress[3] = stress_cor[3];
		stress[4] = stress_cor[4];
		stress[5] = stress_cor[5];
		vector6_subtract(stress, ori_stress, dstress);

		// update Dep
		double tmp_mat[6][6];
		matrix6x6_prod_matrix6x6(Dep_mat, tmat6x6T, tmp_mat);
		matrix6x6_prod_matrix6x6(tmat6x6, tmp_mat, Dep_mat);
		//for (size_t i = 0; i < 6; i++)
		//{
		//	for (size_t j = 0; j < 6; j++)
		//		std::cout << Dep_mat[i][j] << ", ";
		//	std::cout << "\n";
		//}

		return 0;
	}
}
