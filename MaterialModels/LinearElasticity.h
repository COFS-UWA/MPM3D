#ifndef __Linear_Elasticity_h__
#define __Linear_Elasticity_h__

#include "MaterialModel.h"

namespace MatModel
{
	int linear_elasticity_integration_function(MaterialModel *_self, double dstrain[6]);
	void linear_elasticity_store_to_function(const MaterialModel* _self, char* mat_model_mem);
	void linear_elasticity_retrieve_from_function(MaterialModel* _self, const char* mat_model_mem);

	class LinearElasticity : public MaterialModel
	{
		friend int linear_elasticity_integration_function(MaterialModel *_self, double dstrain[6]);
		friend void linear_elasticity_store_to_function(const MaterialModel* _self, char* mat_model_mem);
		friend void linear_elasticity_retrieve_from_function(MaterialModel* _self, const char* mat_model_mem);

	public:
		double E, niu;

		LinearElasticity() :
			MaterialModel(&linear_elasticity_integration_function,
				"LinearElasticity",
				&linear_elasticity_store_to_function,
				&linear_elasticity_retrieve_from_function)
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