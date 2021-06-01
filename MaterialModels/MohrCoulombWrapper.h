#ifndef __Mohr_Coulomb_Wrapper_h__
#define __Mohr_Coulomb_Wrapper_h__

#include "MaterialModel.h"
#include "mohr_coulomb.h"

namespace MatModel
{
	int mohr_coulomb_wrapper_integration_function(MaterialModel* _self, double dstrain[6]);

	class MohrCoulombWrapper : public MaterialModel
	{
		friend int mohr_coulomb_wrapper_integration_function(MaterialModel* _self, double dstrain[6]);
	
	protected:
		MohrCoulombGlobal glb;
		MohrCoulomb mat;
		int32_t status_code;

	public:
		MohrCoulombWrapper() :
			MaterialModel(mohr_coulomb_wrapper_integration_function, "MohrCoulomb") {}
		MohrCoulombWrapper(const MohrCoulombWrapper& other);
		~MohrCoulombWrapper() {}

		void set_param(const __Float_Type__ ini_stress[6],
			__Float_Type__ _phi, __Float_Type__ _psi, __Float_Type__ _cohesion,
			__Float_Type__ _E, __Float_Type__ _niu);
		const __Float_Type__* get_stress() const noexcept { return mat.stress; }
		
		double get_phi() const noexcept { return glb.phi; }
		double get_psi() const noexcept { return glb.psi; }
		double get_cohesion() const noexcept { return glb.cohesion; }
		double get_E() const noexcept { return glb.E; }
		double get_niu() const noexcept { return glb.niu; }
		int32_t get_status_code() const noexcept { return status_code; }
	};
}

#endif