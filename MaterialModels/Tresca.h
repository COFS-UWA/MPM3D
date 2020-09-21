#ifndef __Tresca_h__
#define __Tresca_h__

#include <string>
#include <cmath>
#include "MatModelUtils.h"
#include "MaterialModel.h"

namespace MatModel
{
	int tresca_integration_function(MaterialModel* _self, double dstrain[6]);

	class Tresca : public MaterialModel
	{
		friend int tresca_integration_function(MaterialModel* _self, double dstrain[6]);

	protected:
		double E, niu;
		double cohesion;
		double two_c;

	public:
		Tresca() : MaterialModel(tresca_integration_function, "Tresca"),
			E(0.0), niu(0.0), cohesion(1.0), two_c(cohesion+cohesion)
		{
			for (size_t i = 0; i < 6; ++i)
			{
				stress[i] = 0.0;
				dstress[i] = 0.0;
				dstrain_e[i] = 0.0;
				dstrain_p[i] = 0.0;
			}
		}
		~Tresca() {}

		inline void set_param(double _E, double _niu, double _cohesion)
		{
			E = _E;
			niu = _niu;
			cohesion = _cohesion;
			two_c = cohesion + cohesion;
			form_De_mat();
		}

	protected:
		inline void form_De_mat()
		{
			using MatModel_Internal::form_linear_elastic_stiffness_matrix;
			form_linear_elastic_stiffness_matrix(E, niu, De_mat);
		}
		inline void form_Dep_mat() noexcept
		{
			memcpy(Dep_mat, De_mat, 6 * 6 * sizeof(double));
		}
		inline void form_Dep_mat(double dg_ds[6], double divider) noexcept
		{
			double De_dg_ds[6];
			for (size_t i = 0; i < 6; ++i)
				De_dg_ds[i] = De_mat[i][0] * dg_ds[0] + De_mat[i][1] * dg_ds[1]
				+ De_mat[i][2] * dg_ds[2] + De_mat[i][3] * dg_ds[3]
				+ De_mat[i][4] * dg_ds[4] + De_mat[i][5] * dg_ds[5];
			for (size_t i = 0; i < 6; ++i)
				for (size_t j = 0; j < 6; ++j)
					Dep_mat[i][j] = De_mat[i][j] - De_dg_ds[i] * De_dg_ds[j] / divider;
		}
	};
}

#endif