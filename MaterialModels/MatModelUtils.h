#ifndef __Mat_Model_Utils_h__
#define __Mat_Model_Utils_h__

#include "MaterialModel.h"

namespace MatModel
{
	namespace MatModel_Internal
	{
		inline void form_linear_elastic_stiffness_matrix(
			const double E,
			const double niu,
			double De_mat[6][6]
			)
		{
			double coef;
			coef = (1.0 - niu) / ((1.0 + niu) * (1.0 - 2.0 * niu)) * E;
			De_mat[0][0] = coef;
			De_mat[1][1] = coef;
			De_mat[2][2] = coef;

			coef = niu / ((1.0 + niu) * (1.0 - 2.0 * niu)) * E;
			De_mat[0][1] = coef;
			De_mat[0][2] = coef;
			De_mat[1][0] = coef;
			De_mat[1][2] = coef;
			De_mat[2][0] = coef;
			De_mat[2][1] = coef;

			De_mat[0][3] = 0.0;
			De_mat[0][4] = 0.0;
			De_mat[0][5] = 0.0;
			De_mat[1][3] = 0.0;
			De_mat[1][4] = 0.0;
			De_mat[1][5] = 0.0;
			De_mat[2][3] = 0.0;
			De_mat[2][4] = 0.0;
			De_mat[2][5] = 0.0;

			De_mat[3][0] = 0.0;
			De_mat[3][1] = 0.0;
			De_mat[3][2] = 0.0;
			De_mat[4][0] = 0.0;
			De_mat[4][1] = 0.0;
			De_mat[4][2] = 0.0;
			De_mat[5][0] = 0.0;
			De_mat[5][1] = 0.0;
			De_mat[5][2] = 0.0;

			coef = E / (1.0 + niu);
			De_mat[3][3] = coef;
			De_mat[4][4] = coef;
			De_mat[5][5] = coef;

			De_mat[3][4] = 0.0;
			De_mat[3][5] = 0.0;
			De_mat[4][3] = 0.0;
			De_mat[4][5] = 0.0;
			De_mat[5][3] = 0.0;
			De_mat[5][4] = 0.0;
		}

		inline void vector6_zeros(double res[6])
		{
			res[0] = 0.0;
			res[1] = 0.0;
			res[2] = 0.0;
			res[3] = 0.0;
			res[4] = 0.0;
			res[5] = 0.0;
		}

		inline void vector6_copy(
			const double vec[6],
			double res[6]
			)
		{
			res[0] = vec[0];
			res[1] = vec[1];
			res[2] = vec[2];
			res[3] = vec[3];
			res[4] = vec[4];
			res[5] = vec[5];
		}
		
		inline void vector6_add(
			const double vec[6],
			double res[6]
			)
		{
			res[0] += vec[0];
			res[1] += vec[1];
			res[2] += vec[2];
			res[3] += vec[3];
			res[4] += vec[4];
			res[5] += vec[5];
		}

		inline void vector6_add(
			const double vec1[6],
			const double vec2[6],
			double res[6]
		)
		{
			res[0] = vec1[0] + vec2[0];
			res[1] = vec1[1] + vec2[1];
			res[2] = vec1[2] + vec2[2];
			res[3] = vec1[3] + vec2[3];
			res[4] = vec1[4] + vec2[4];
			res[5] = vec1[5] + vec2[5];
		}
		
		inline void vector6_subtract(
			const double vec[6],
			double res[6]
			)
		{
			res[0] -= vec[0];
			res[1] -= vec[1];
			res[2] -= vec[2];
			res[3] -= vec[3];
			res[4] -= vec[4];
			res[5] -= vec[5];
		}

		inline void vector6_subtract(
			const double vec1[6],
			const double vec2[6],
			double res[6]
			)
		{
			res[0] = vec1[0] - vec2[0];
			res[1] = vec1[1] - vec2[1];
			res[2] = vec1[2] - vec2[2];
			res[3] = vec1[3] - vec2[3];
			res[4] = vec1[4] - vec2[4];
			res[5] = vec1[5] - vec2[5];
		}

		inline void vector6_scale(
			const double vec[6],
			const double factor,
			double res[6]
			)
		{
			res[0] = factor * vec[0];
			res[1] = factor * vec[1];
			res[2] = factor * vec[2];
			res[3] = factor * vec[3];
			res[4] = factor * vec[4];
			res[5] = factor * vec[5];
		}

		inline void matrix6x6_subtract(
			const double mat1[6][6],
			const double mat2[6][6],
			double res[6][6]
		)
		{
			res[0][0] = mat1[0][0] - mat2[0][0];
			res[0][1] = mat1[0][1] - mat2[0][1];
			res[0][2] = mat1[0][2] - mat2[0][2];
			res[0][3] = mat1[0][3] - mat2[0][3];
			res[0][4] = mat1[0][4] - mat2[0][4];
			res[0][5] = mat1[0][5] - mat2[0][5];
			res[1][0] = mat1[1][0] - mat2[1][0];
			res[1][1] = mat1[1][1] - mat2[1][1];
			res[1][2] = mat1[1][2] - mat2[1][2];
			res[1][3] = mat1[1][3] - mat2[1][3];
			res[1][4] = mat1[1][4] - mat2[1][4];
			res[1][5] = mat1[1][5] - mat2[1][5];
			res[2][0] = mat1[2][0] - mat2[2][0];
			res[2][1] = mat1[2][1] - mat2[2][1];
			res[2][2] = mat1[2][2] - mat2[2][2];
			res[2][3] = mat1[2][3] - mat2[2][3];
			res[2][4] = mat1[2][4] - mat2[2][4];
			res[2][5] = mat1[2][5] - mat2[2][5];
			res[3][0] = mat1[3][0] - mat2[3][0];
			res[3][1] = mat1[3][1] - mat2[3][1];
			res[3][2] = mat1[3][2] - mat2[3][2];
			res[3][3] = mat1[3][3] - mat2[3][3];
			res[3][4] = mat1[3][4] - mat2[3][4];
			res[3][5] = mat1[3][5] - mat2[3][5];
			res[4][0] = mat1[4][0] - mat2[4][0];
			res[4][1] = mat1[4][1] - mat2[4][1];
			res[4][2] = mat1[4][2] - mat2[4][2];
			res[4][3] = mat1[4][3] - mat2[4][3];
			res[4][4] = mat1[4][4] - mat2[4][4];
			res[4][5] = mat1[4][5] - mat2[4][5];
			res[5][0] = mat1[5][0] - mat2[5][0];
			res[5][1] = mat1[5][1] - mat2[5][1];
			res[5][2] = mat1[5][2] - mat2[5][2];
			res[5][3] = mat1[5][3] - mat2[5][3];
			res[5][4] = mat1[5][4] - mat2[5][4];
			res[5][5] = mat1[5][5] - mat2[5][5];
		}

		inline void matrix6x6_prod_vector6(
			const double mat[6][6],
			const double vec[6],
			double res[6]
			)
		{
			res[0] = mat[0][0] * vec[0] + mat[0][1] * vec[1]
					+ mat[0][2] * vec[2] + mat[0][3] * vec[3]
					+ mat[0][4] * vec[4] + mat[0][5] * vec[5];
			res[1] = mat[1][0] * vec[0] + mat[1][1] * vec[1]
					+ mat[1][2] * vec[2] + mat[1][3] * vec[3]
					+ mat[1][4] * vec[4] + mat[1][5] * vec[5];
			res[2] = mat[2][0] * vec[0] + mat[2][1] * vec[1]
					+ mat[2][2] * vec[2] + mat[2][3] * vec[3]
					+ mat[2][4] * vec[4] + mat[2][5] * vec[5];
			res[3] = mat[3][0] * vec[0] + mat[3][1] * vec[1]
					+ mat[3][2] * vec[2] + mat[3][3] * vec[3]
					+ mat[3][4] * vec[4] + mat[3][5] * vec[5];
			res[4] = mat[4][0] * vec[0] + mat[4][1] * vec[1]
					+ mat[4][2] * vec[2] + mat[4][3] * vec[3]
					+ mat[4][4] * vec[4] + mat[4][5] * vec[5];
			res[5] = mat[5][0] * vec[0] + mat[5][1] * vec[1]
					+ mat[5][2] * vec[2] + mat[5][3] * vec[3]
					+ mat[5][4] * vec[4] + mat[5][5] * vec[5];
		}

		inline double vector6_prod_matrix6x6_prod_vector6(
			const double vec1[6],
			const double mat[6][6],
			const double vec2[6]
			)
		{
			return vec1[0] * mat[0][0] * vec2[0] + vec1[0] * mat[0][1] * vec2[1]
				+ vec1[0] * mat[0][2] * vec2[2] + vec1[0] * mat[0][3] * vec2[3]
				+ vec1[0] * mat[0][4] * vec2[4] + vec1[0] * mat[0][5] * vec2[5]
				+ vec1[1] * mat[1][0] * vec2[0] + vec1[1] * mat[1][1] * vec2[1]
				+ vec1[1] * mat[1][2] * vec2[2] + vec1[1] * mat[1][3] * vec2[3]
				+ vec1[1] * mat[1][4] * vec2[4] + vec1[1] * mat[1][5] * vec2[5]
				+ vec1[2] * mat[2][0] * vec2[0] + vec1[2] * mat[2][1] * vec2[1]
				+ vec1[2] * mat[2][2] * vec2[2] + vec1[2] * mat[2][3] * vec2[3]
				+ vec1[2] * mat[2][4] * vec2[4] + vec1[2] * mat[2][5] * vec2[5]
				+ vec1[3] * mat[3][0] * vec2[0] + vec1[3] * mat[3][1] * vec2[1]
				+ vec1[3] * mat[3][2] * vec2[2] + vec1[3] * mat[3][3] * vec2[3]
				+ vec1[3] * mat[3][4] * vec2[4] + vec1[3] * mat[3][5] * vec2[5]
				+ vec1[4] * mat[4][0] * vec2[0] + vec1[4] * mat[4][1] * vec2[1]
				+ vec1[4] * mat[4][2] * vec2[2] + vec1[4] * mat[4][3] * vec2[3]
				+ vec1[4] * mat[4][4] * vec2[4] + vec1[4] * mat[4][5] * vec2[5]
				+ vec1[5] * mat[5][0] * vec2[0] + vec1[5] * mat[5][1] * vec2[1]
				+ vec1[5] * mat[5][2] * vec2[2] + vec1[5] * mat[5][3] * vec2[3]
				+ vec1[5] * mat[5][4] * vec2[4] + vec1[5] * mat[5][5] * vec2[5];
		}

		inline void matrix6x6_copy(
			const double mat[6][6],
			double res[6][6]
			)
		{
			res[0][0] = mat[0][0];
			res[0][1] = mat[0][1];
			res[0][2] = mat[0][2];
			res[0][3] = mat[0][3];
			res[0][4] = mat[0][4];
			res[0][5] = mat[0][5];
			res[1][0] = mat[1][0];
			res[1][1] = mat[1][1];
			res[1][2] = mat[1][2];
			res[1][3] = mat[1][3];
			res[1][4] = mat[1][4];
			res[1][5] = mat[1][5];
			res[2][0] = mat[2][0];
			res[2][1] = mat[2][1];
			res[2][2] = mat[2][2];
			res[2][3] = mat[2][3];
			res[2][4] = mat[2][4];
			res[2][5] = mat[2][5];
			res[3][0] = mat[3][0];
			res[3][1] = mat[3][1];
			res[3][2] = mat[3][2];
			res[3][3] = mat[3][3];
			res[3][4] = mat[3][4];
			res[3][5] = mat[3][5];
			res[4][0] = mat[4][0];
			res[4][1] = mat[4][1];
			res[4][2] = mat[4][2];
			res[4][3] = mat[4][3];
			res[4][4] = mat[4][4];
			res[4][5] = mat[4][5];
			res[5][0] = mat[5][0];
			res[5][1] = mat[5][1];
			res[5][2] = mat[5][2];
			res[5][3] = mat[5][3];
			res[5][4] = mat[5][4];
			res[5][5] = mat[5][5];
		}

		inline void matrix6x6_prod_matrix6x6(
			const double mat1[6][6],
			const double mat2[6][6],
			double res[6][6]
			)
		{
#define mat_prod(i, j) res[i][j] = \
		  mat1[i][0] * mat2[0][j] + mat1[i][1] * mat2[1][j] \
		+ mat1[i][2] * mat2[2][j] + mat1[i][3] * mat2[3][j] \
		+ mat1[i][4] * mat2[4][j] + mat1[i][5] * mat2[5][j]
			mat_prod(0, 0);
			mat_prod(0, 1);
			mat_prod(0, 2);
			mat_prod(0, 3);
			mat_prod(0, 4);
			mat_prod(0, 5);
			mat_prod(1, 0);
			mat_prod(1, 1);
			mat_prod(1, 2);
			mat_prod(1, 3);
			mat_prod(1, 4);
			mat_prod(1, 5);
			mat_prod(2, 0);
			mat_prod(2, 1);
			mat_prod(2, 2);
			mat_prod(2, 3);
			mat_prod(2, 4);
			mat_prod(2, 5);
			mat_prod(3, 0);
			mat_prod(3, 1);
			mat_prod(3, 2);
			mat_prod(3, 3);
			mat_prod(3, 4);
			mat_prod(3, 5);
			mat_prod(4, 0);
			mat_prod(4, 1);
			mat_prod(4, 2);
			mat_prod(4, 3);
			mat_prod(4, 4);
			mat_prod(4, 5);
			mat_prod(5, 0);
			mat_prod(5, 1);
			mat_prod(5, 2);
			mat_prod(5, 3);
			mat_prod(5, 4);
			mat_prod(5, 5);
#undef mat_prod
		}
	}

	bool integrate_drained_triaxial_test(
		MaterialModel& model,
		double de,
		double dstrain[6],
		double tol = 1.0e-3,
		size_t max_iter_num = 100);
}

#endif