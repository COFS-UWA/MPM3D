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
		Tresca& self = static_cast<Tresca&>(*_self);

		double ori_stress[6];
		double *dstress = self.dstress;
		double *stress = self.stress;
		double *dstrain_e = self.dstrain_e;
		double *dstrain_p = self.dstrain_p;
		double (*De_mat)[6] = self.De_mat;
		vector6_copy(stress, ori_stress);
		matrix6x6_prod_vector6(De_mat, dstrain, dstress);
		vector6_add(dstress, stress);
		vector6_copy(dstrain, dstrain_e);
		vector6_zeros(dstrain_p);

		Eigen::Matrix3d stress_mat;
		stress_mat << stress[0], stress[3], stress[5],
					  stress[3], stress[1], stress[4],
					  stress[5], stress[4], stress[2];

		// ??
		Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(stress_mat);
		if (eigen_solver.info() != Eigen::Success)
			assert(0);
		
		const Eigen::Vector3d &prin_stress = eigen_solver.eigenvalues();
		const Eigen::Matrix3d &eigen_vector = eigen_solver.eigenvectors();

		size_t s1_id, s2_id, s3_id;
		double s1, s2, s3;
		s1_id = 0;
		s2_id = 1;
		s3_id = 2;
		s1 = prin_stress[0];
		s2 = prin_stress[1];
		s3 = prin_stress[2];

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
		if (criteria1 < -self.two_c) // return to line1
		{
			s_mean = (s1 + s2 + s3) / 3.0;
			s1 = s_mean + 2.0 / 3.0 * self.cohesion;
			s2 = s1;
			s3 = s_mean - 4.0 / 3.0 * self.cohesion;
		}
		else if (criteria1 > self.two_c) // return to line2
		{
			s_mean = (s1 + s2 + s3) / 3.0;
			s1 = s_mean + 4.0 / 3.0 * self.cohesion;
			s2 = s_mean - 2.0 / 3.0 * self.cohesion;
			s3 = s2;
		}
		else // return to plane
		{
			lambda = (self.two_c - (s1 - s3)) / (2.0 * (De_mat[0][1] - De_mat[0][0]));
			s1 -= (De_mat[0][0] - De_mat[0][2]) * lambda;
			s2 -= (De_mat[1][0] - De_mat[1][2]) * lambda;
			s3 -= (De_mat[2][0] - De_mat[2][2]) * lambda;
		}

		// transfer back to global coordinates
		Eigen::Matrix3d cor_pstress_mat;
		cor_pstress_mat.setZero();
		cor_pstress_mat(s1_id, s1_id) = s1;
		cor_pstress_mat(s2_id, s2_id) = s2;
		cor_pstress_mat(s3_id, s3_id) = s3;
		Eigen::Matrix3d cor_stress = eigen_vector.transpose() * cor_pstress_mat * eigen_vector;

		// dep = D^-1 * (s_p - s_c)
		Eigen::Matrix<double, 6, 6> De;
		De << De_mat[0][0], De_mat[0][1], De_mat[0][2], De_mat[0][3], De_mat[0][4], De_mat[0][5],
			  De_mat[1][0], De_mat[1][1], De_mat[1][2], De_mat[1][3], De_mat[1][4], De_mat[1][5],
			  De_mat[2][0], De_mat[2][1], De_mat[2][2], De_mat[2][3], De_mat[2][4], De_mat[2][5], 
			  De_mat[3][0], De_mat[3][1], De_mat[3][2], De_mat[3][3], De_mat[3][4], De_mat[3][5], 
			  De_mat[4][0], De_mat[4][1], De_mat[4][2], De_mat[4][3], De_mat[4][4], De_mat[4][5], 
			  De_mat[5][0], De_mat[5][1], De_mat[5][2], De_mat[5][3], De_mat[5][4], De_mat[5][5];
		Eigen::Matrix<double, 6, 1> s_diff;
		s_diff[0] = stress[0] - cor_stress(0, 0);
		s_diff[1] = stress[1] - cor_stress(1, 1);
		s_diff[2] = stress[2] - cor_stress(2, 2);
		s_diff[3] = stress[3] - cor_stress(0, 1);
		s_diff[4] = stress[4] - cor_stress(1, 2);
		s_diff[5] = stress[5] - cor_stress(0, 2);
		Eigen::Matrix<double, 6, 1> dstrain_p_vec = De.partialPivLu().solve(s_diff);
		dstrain_p[0] = dstrain_p_vec[0];
		dstrain_p[1] = dstrain_p_vec[1];
		dstrain_p[2] = dstrain_p_vec[2];
		dstrain_p[3] = dstrain_p_vec[3];
		dstrain_p[4] = dstrain_p_vec[4];
		dstrain_p[5] = dstrain_p_vec[5];
		vector6_subtract(dstrain_p, dstrain_e);
		stress[0] = cor_stress(0, 0);
		stress[1] = cor_stress(1, 1);
		stress[2] = cor_stress(2, 2);
		stress[3] = cor_stress(0, 1);
		stress[4] = cor_stress(1, 2);
		stress[5] = cor_stress(0, 2);
		vector6_subtract(stress, ori_stress, dstress);

		// update De here...

		return 0;
	}
}
