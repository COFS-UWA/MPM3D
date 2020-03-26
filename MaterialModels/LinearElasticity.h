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
			MaterialModel(linear_elasticity_integration_function, MaterialModelType::LinearElasticity) {}
		~LinearElasticity() {}

		void set_param(double _E, double _niu);
	};

};

#endif