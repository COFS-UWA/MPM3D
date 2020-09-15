#include "MaterialModels_pcp.h"

#include <Eigen/Eigen>

#include "Tresca.h"

namespace MatModel
{
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
		double f;
		for (size_t i = 0; i < 6; ++i)
		{
			// stress
			dstress[i] = De_mat[i][0] * dstrain[0] + De_mat[i][1] * dstrain[1]
				+ De_mat[i][2] * dstrain[2] + De_mat[i][3] * dstrain[3]
				+ De_mat[i][4] * dstrain[4] + De_mat[i][5] * dstrain[5];
			stress[i] += dstress[i];
			// strain
			dstrain_e[i] = dstrain[i];
			dstrain_p[i] = 0.0;
		}

		f = self.get_f();
		if (f <= 0.0)
		{
			self.form_Dep_mat();
			return 0;
		}

		double dg_ds[6], divider, dl, dep_cor[6], ds_cor[6];
		for (size_t iter_id = 0; iter_id < 20; ++iter_id)
		{
			self.cal_dg_stress(dg_ds);
			divider = self.cal_divider(dg_ds);
			dl = f / divider;
			// strain correction
			dep_cor[0] = dl * dg_ds[0];
			dep_cor[1] = dl * dg_ds[1];
			dep_cor[2] = dl * dg_ds[2];
			dep_cor[3] = dl * dg_ds[3];
			dep_cor[4] = dl * dg_ds[4];
			dep_cor[5] = dl * dg_ds[5];
			// elastic strain
			dstrain_e[0] -= dep_cor[0];
			dstrain_e[1] -= dep_cor[1];
			dstrain_e[2] -= dep_cor[2];
			dstrain_e[3] -= dep_cor[3];
			dstrain_e[4] -= dep_cor[4];
			dstrain_e[5] -= dep_cor[5];
			// plastic strain
			dstrain_p[0] += dep_cor[0];
			dstrain_p[1] += dep_cor[1];
			dstrain_p[2] += dep_cor[2];
			dstrain_p[3] += dep_cor[3];
			dstrain_p[4] += dep_cor[4];
			dstrain_p[5] += dep_cor[5];
			// correct stress
			for (size_t i = 0; i < 6; ++i)
				ds_cor[i] = De_mat[i][0] * dep_cor[0] + De_mat[i][1] * dep_cor[1]
				+ De_mat[i][2] * dep_cor[2] + De_mat[i][3] * dep_cor[3]
				+ De_mat[i][4] * dep_cor[4] + De_mat[i][5] * dep_cor[5];
			dstress[0] -= ds_cor[0];
			dstress[1] -= ds_cor[1];
			dstress[2] -= ds_cor[2];
			dstress[3] -= ds_cor[3];
			dstress[4] -= ds_cor[4];
			dstress[5] -= ds_cor[5];
			stress[0] -= ds_cor[0];
			stress[1] -= ds_cor[1];
			stress[2] -= ds_cor[2];
			stress[3] -= ds_cor[3];
			stress[4] -= ds_cor[4];
			stress[5] -= ds_cor[5];

			if (self.get_norm_f() < 1.0e-8)
			{
				// form new stiffness matrix
				self.cal_dg_stress(dg_ds);
				divider = self.cal_divider(dg_ds);
				self.form_Dep_mat(dg_ds, divider);
				return iter_id + 1;
			}
		}

		// form new stiffness matrix
		self.cal_dg_stress(dg_ds);
		divider = self.cal_divider(dg_ds);
		self.form_Dep_mat(dg_ds, divider);
		return -1;
	}
}
