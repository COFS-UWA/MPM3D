#ifndef __Sand_Hypoplasticity_Stb_Wrapper_h__
#define __Sand_Hypoplasticity_Stb_Wrapper_h__

#include "MaterialModel.h"
#include "sand_hypoplasticity_stb.h"

namespace MatModel
{
	int sand_hypoplasticity_stb_wrapper_integration_function(MaterialModel* _self, double dstrain[6]);

	class SandHypoplasticityStbWrapper : public MaterialModel
	{
		friend int sand_hypoplasticity_stb_wrapper_integration_function(MaterialModel* _self, double dstrain[6]);

	protected:
		SandHypoplasticityStbGlobal glb;
		SandHypoplasticityStb mat;
		int32_t status_code;

	public:
		SandHypoplasticityStbWrapper();
		SandHypoplasticityStbWrapper(const SandHypoplasticityStbWrapper& other);
		~SandHypoplasticityStbWrapper() {}
		void set_substep_size(__Float_Type__ substp_size) noexcept
		{
			mat.substp_size = substp_size;
		};
		void set_param(const __Float_Type__* ini_stress,
			__Float_Type__ e0,
			__Float_Type__ phi,
			__Float_Type__ hs, __Float_Type__ n,
			__Float_Type__ alpha, __Float_Type__ beta,
			__Float_Type__ ed0, __Float_Type__ ec0, __Float_Type__ ei0,
			__Float_Type__ N, __Float_Type__ chi, __Float_Type__ H,
			__Float_Type__ Ig, __Float_Type__ niu);
		void set_yield_surface(__Float_Type__ Mi, __Float_Type__ pi, __Float_Type__ pl);
		__Float_Type__ get_phi() const noexcept { return glb.phi; }
		__Float_Type__ get_hs() const noexcept { return glb.hs; }
		__Float_Type__ get_n() const noexcept { return glb.n; }
		__Float_Type__ get_alpha() const noexcept { return glb.alpha; }
		__Float_Type__ get_beta() const noexcept { return glb.beta; }
		__Float_Type__ get_ei0() const noexcept { return glb.ei0; }
		__Float_Type__ get_ec0() const noexcept { return glb.ec0; }
		__Float_Type__ get_ed0() const noexcept { return glb.ed0; }
		__Float_Type__ get_N() const noexcept { return glb.N; }
		__Float_Type__ get_chi() const noexcept { return glb.chi; }
		__Float_Type__ get_H() const noexcept { return glb.H; }
		__Float_Type__ get_Ig() const noexcept { return glb.Ig; }
		__Float_Type__ get_niu() const noexcept { return glb.niu; }
		const __Float_Type__* get_stress() const noexcept { return mat.stress; }
		__Float_Type__ get_e() const noexcept { return mat.e; }
		__Float_Type__ get_pi() const noexcept { return mat.pi; }
		__Float_Type__ get_Mi() const noexcept { return mat.Mi; }
		__Float_Type__ get_pl() const noexcept { return mat.pl; }
		__Float_Type__ get_substp_size() const noexcept { return mat.substp_size; }
		int32_t get_status_code() const noexcept { return status_code; }
	};
}

#endif