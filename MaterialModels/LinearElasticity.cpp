#include "MaterialModels_pcp.h"

#include "LinearElasticity.h"

namespace MatModel
{
	int linear_elasticity_integration_function(MaterialModel *_self, double dstrain[6])
	{
		LinearElasticity &self = static_cast<LinearElasticity &>(*_self);

		double *dstrain_e = self.dstrain_e;
		double *dstrain_p = self.dstrain_p;
		double *dstress = self.dstress;
		double(*De_mat)[6] = self.De_mat;

		dstrain_e[0] = dstrain[0];
		dstrain_e[1] = dstrain[1];
		dstrain_e[2] = dstrain[2];
		dstrain_e[3] = dstrain[3];
		dstrain_e[4] = dstrain[4];
		dstrain_e[5] = dstrain[5];

		dstrain_p[0] = 0.0;
		dstrain_p[1] = 0.0;
		dstrain_p[2] = 0.0;
		dstrain_p[3] = 0.0;
		dstrain_p[4] = 0.0;
		dstrain_p[5] = 0.0;

		dstress[0] = De_mat[0][0] * dstrain[0] + De_mat[0][1] * dstrain[1]
			+ De_mat[0][2] * dstrain[2];
		dstress[1] = De_mat[1][0] * dstrain[0] + De_mat[1][1] * dstrain[1]
			+ De_mat[1][2] * dstrain[2];
		dstress[2] = De_mat[2][0] * dstrain[0] + De_mat[2][1] * dstrain[1]
			+ De_mat[2][2] * dstrain[2];
		dstress[3] = De_mat[3][3] * dstrain[3];
		dstress[4] = De_mat[4][4] * dstrain[4];
		dstress[5] = De_mat[5][5] * dstrain[5];

		return 0;
	}

	void LinearElasticity::set_param(double _E, double _niu)
	{
		E = _E;
		niu = _niu;

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

		memcpy(Dep_mat, De_mat, sizeof(double) * 36);
	}

};