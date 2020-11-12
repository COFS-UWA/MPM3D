#ifndef __Linear_Elasticity_h__
#define __Linear_Elasticity_h__

#include "MaterialModel.h"

namespace MatModel
{
	int linear_elasticity_integration_function(MaterialModel *_self, double dstrain[6]);

	class LinearElasticity : public MaterialModel
	{
		friend int linear_elasticity_integration_function(MaterialModel *_self, double dstrain[6]);
	public:
		double E, niu;

		LinearElasticity() :
			MaterialModel(linear_elasticity_integration_function, "LinearElasticity")
		{
			dstrain_p[0] = 0.0;
			dstrain_p[1] = 0.0;
			dstrain_p[2] = 0.0;
			dstrain_p[3] = 0.0;
			dstrain_p[4] = 0.0;
			dstrain_p[5] = 0.0;
			stress[0] = 0.0;
			stress[1] = 0.0;
			stress[2] = 0.0;
			stress[3] = 0.0;
			stress[4] = 0.0;
			stress[5] = 0.0;
		}
		~LinearElasticity() {}

		void set_param(double _E, double _niu);
		void set_stress(double s[6]);
	};

};

#endif