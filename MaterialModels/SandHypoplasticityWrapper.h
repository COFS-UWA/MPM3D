#ifndef __Sand_Hypoplasticity_Wrapper_h__
#define __Sand_Hypoplasticity_Wrapper_h__

#include "MaterialModel.h"
#include "sand_hypoplasticity.h"

namespace MatModel
{
	int sand_hypoplasticity_wrapper_integration_function(MaterialModel* _self, double dstrain[6]);

	class SandHypoplasticityWrapper : public MaterialModel
	{
		friend int sand_hypoplasticity_wrapper_integration_function(MaterialModel* _self, double dstrain[6]);
	
	protected:
		SandHypoplasticityGlobal glb;
		SandHypoplasticity mat;
		int32_t status_code;

	public:
		SandHypoplasticityWrapper();
		SandHypoplasticityWrapper(const SandHypoplasticityWrapper& other);
		~SandHypoplasticityWrapper() {}
		void set_substep_size(__Float_Type__ substp_size) noexcept
		{ mat.substp_size = substp_size; };
		void set_param(const __Float_Type__* ini_stress,
			__Float_Type__ e0,
			__Float_Type__ phi,
			__Float_Type__ hs, __Float_Type__ n,
			__Float_Type__ ed0, __Float_Type__ ec0, __Float_Type__ ei0,
			__Float_Type__ alpha, __Float_Type__ beta);
		__Float_Type__ get_phi() const noexcept { return glb.phi; }
		__Float_Type__ get_hs() const noexcept { return glb.hs; }
		__Float_Type__ get_n() const noexcept { return glb.n; }
		__Float_Type__ get_ei0() const noexcept { return glb.ei0; }
		__Float_Type__ get_ec0() const noexcept { return glb.ec0; }
		__Float_Type__ get_ed0() const noexcept { return glb.ed0; }
		__Float_Type__ get_alpha() const noexcept { return glb.alpha; }
		__Float_Type__ get_beta() const noexcept { return glb.beta; }
		const __Float_Type__* get_stress() const noexcept { return mat.stress; }
		__Float_Type__ get_e() const noexcept { return mat.e; }
		__Float_Type__ get_substp_size() const noexcept { return mat.substp_size; }
		int32_t get_status_code() const noexcept { return status_code; }
	};
}

#endif