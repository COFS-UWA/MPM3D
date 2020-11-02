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
			MatModel_Internal::vector6_zeros(stress);
			MatModel_Internal::vector6_zeros(dstress);
			MatModel_Internal::vector6_zeros(dstrain_e);
			MatModel_Internal::vector6_zeros(dstrain_p);
		}
		~Tresca() {}

		inline double get_E() const noexcept { return E; }
		inline double get_niu() const noexcept { return niu; }
		inline double get_cohesion() const noexcept { return cohesion; }

		inline void set_param(
			double _E,
			double _niu,
			double _cohesion,
			double ini_stress[6] = nullptr
			)
		{
			E = _E;
			niu = _niu;
			cohesion = _cohesion;
			two_c = cohesion + cohesion;
			if (ini_stress)
				MatModel_Internal::vector6_copy(ini_stress, stress);
			else
				MatModel_Internal::vector6_zeros(stress);
			form_De_mat();
			copy_De_to_Dep();
		}

	protected:
		inline void form_De_mat() noexcept
		{
			MatModel_Internal::form_linear_elastic_stiffness_matrix(E, niu, De_mat);
		}
		inline void copy_De_to_Dep() noexcept
		{
			MatModel_Internal::matrix6x6_copy(De_mat, Dep_mat);
		}
	};
}

#endif