#include "MaterialModels_pcp.h"

#include "MatModelUtils.h"

namespace MatModel
{
	bool integrate_drained_triaxial_test(
		MaterialModel& model,
		double de,
		double dstrain[6],
		double tol,
		size_t max_iter_num)
	{
		if (de == 0.0)
			return true;

		const double* stress = model.get_stress();
		const double init_s11 = stress[0];
		const double init_s22 = stress[1];
		const double (*D)[6] = (const double(*)[6])model.get_Dep_mat();
		double ddstrain[6] = {
			0.0, 0.0,
			// 0.7 is an empirical factor
			//0.7 * (D[1][2] * D[0][1] - D[1][1] * D[0][2]) * de / (D[0][0] * D[1][1] - D[1][0] * D[0][1]), // e11
			//0.7 * (D[0][2] * D[1][0] - D[0][0] * D[1][2]) * de / (D[0][0] * D[1][1] - D[1][0] * D[0][1]), // e22
			de,	0.0, 0.0, 0.0
		};
		int res = model.integrate(ddstrain);
		dstrain[0] = ddstrain[0];
		dstrain[1] = ddstrain[1];
		dstrain[2] = ddstrain[2];
		dstrain[3] = ddstrain[3];
		dstrain[4] = ddstrain[4];
		dstrain[5] = ddstrain[5];

		ddstrain[2] = 0.0;
		const double two_de_tol2 = 2.0 * de * de * tol * tol;
		const double two_stress_tol2 = 2.0 * tol * tol;
		double ds11 = stress[0] - init_s11;
		double ds22 = stress[1] - init_s22;
		for (size_t iter_id = 0; iter_id < max_iter_num; ++iter_id)
		{
			ddstrain[0] = (D[0][1] * ds22 - D[1][1] * ds11) / (D[0][0] * D[1][1] - D[1][0] * D[0][1]);
			ddstrain[1] = (D[1][0] * ds11 - D[0][0] * ds22) / (D[0][0] * D[1][1] - D[1][0] * D[0][1]);
			res = model.integrate(ddstrain);
			dstrain[0] += ddstrain[0];
			dstrain[1] += ddstrain[1];
			ds11 = stress[0] - init_s11;
			ds22 = stress[1] - init_s22;
			if ((ddstrain[0] * ddstrain[0] + ddstrain[1] * ddstrain[1]) < two_de_tol2 &&
				(ds11 * ds11 + ds22 * ds22) < two_stress_tol2)
				return true;
		}
		return false;
	}
}
