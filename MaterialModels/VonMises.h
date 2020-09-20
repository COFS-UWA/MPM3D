#ifndef __Von_Mises_h__
#define __Von_Mises_h__

#include <string>
#include <cmath>
#include "MatModelUtils.h"
#include "MaterialModel.h"

namespace Model_hdf5_utilities
{
	struct VonMisesStateData;
}

namespace MatModel
{
	int von_mises_integration_function(MaterialModel* _self, double dstrain[6]);

	class VonMises : public MaterialModel
	{
		friend int von_mises_integration_function(MaterialModel* _self, double dstrain[6]);
		friend struct Model_hdf5_utilities::VonMisesStateData;
	protected:
		double E, niu;
		double cohesion;

	public:
		VonMises() :
			MaterialModel(von_mises_integration_function, "VonMises"),
			E(0.0), niu(0.0), cohesion(1.0)
		{
			for (size_t i = 0; i < 6; ++i)
			{
				stress[i] = 0.0;
				dstress[i] = 0.0;
				dstrain_e[i] = 0.0;
				dstrain_p[i] = 0.0;
			}
		}
		~VonMises() {}

		inline void set_param(
			double _E,
			double _niu,
			double _cohesion,
			const double ini_stress[6] = nullptr
			)
		{
			E = _E;
			niu = _niu;
			cohesion = _cohesion;
			if (ini_stress)
			{
				stress[0] = ini_stress[0];
				stress[1] = ini_stress[1];
				stress[2] = ini_stress[2];
				stress[3] = ini_stress[3];
				stress[4] = ini_stress[4];
				stress[5] = ini_stress[5];
			}
			else
			{
				stress[0] = 0.0;
				stress[1] = 0.0;
				stress[2] = 0.0;
				stress[3] = 0.0;
				stress[4] = 0.0;
				stress[5] = 0.0;
			}
			form_De_mat();
		}

		inline double get_p() const noexcept { return cal_p(); }
		inline double get_q() const noexcept { return cal_q(); }
		inline double get_f() const noexcept { return cal_q() - cohesion; }
		inline double get_norm_f() const noexcept { return cal_q() / cohesion - 1.0; }

	protected:
		// form elastic stiffness matrix
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

		inline double cal_p() const noexcept
		{ return (s11 + s22 + s33) / 3.0; }
		inline double cal_q() const noexcept
		{
			double s11_s22_diff = s11 - s22;
			double s22_s33_diff = s22 - s33;
			double s33_s11_diff = s33 - s11;
			double q2 = (s11_s22_diff * s11_s22_diff
					   + s22_s33_diff * s22_s33_diff
					   + s33_s11_diff * s33_s11_diff) * 0.5
					   + (s12 * s12 + s23 * s23 + s31 * s31) * 3.0;
			return sqrt(q2);
		}
		inline double cal_norm_f(double f) const noexcept { return f / cohesion; }

		void cal_dg_stress(double dg_ds[6]) const noexcept
		{
			//dg_ds[0] = (s11 + s11 - s22 - s33) * 0.5 / q;
			//dg_ds[1] = (s22 + s22 - s33 - s11) * 0.5 / q;
			//dg_ds[2] = (s33 + s33 - s11 - s22) * 0.5 / q;
			//dg_ds[3] = 3.0 * s12 / q;
			//dg_ds[4] = 3.0 * s23 / q;
			//dg_ds[5] = 3.0 * s31 / q;
			dg_ds[0] = 0.5 * (s11 + s11 - s22 - s33);
			dg_ds[1] = 0.5 * (s22 + s22 - s33 - s11);
			dg_ds[2] = 0.5 * (s33 + s33 - s11 - s22);
			dg_ds[3] = 3.0 * s12;
			dg_ds[4] = 3.0 * s23;
			dg_ds[5] = 3.0 * s31;
		}
		inline double cal_divider(double dg_ds[6]) const noexcept
		{
			double divider = 0.0;
			for (size_t i = 0; i < 6; i++)
				for (size_t j = 0; j < 6; j++)
					divider += dg_ds[i] * De_mat[i][j] * dg_ds[j];
			return divider;
		}
	};
}

#endif