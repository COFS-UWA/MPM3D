#include "MaterialModels_pcp.h"

#include "MatModelConstants.h"
#include "SandHypoplasticityWrapper.h"

namespace MatModel
{
	SandHypoplasticityWrapper::SandHypoplasticityWrapper(
		const SandHypoplasticityWrapper& other)
	{
		const auto& o_glb = other.glb;
		glb.set_param(o_glb.phi, o_glb.hs, o_glb.n,
			o_glb.alpha, o_glb.beta,
			o_glb.ei0, o_glb.ec0, o_glb.ed0,
			o_glb.ten_E, o_glb.ten_niu);
		const auto& o_mat = other.mat;
		mat.stress[0] = o_mat.stress[0];
		mat.stress[1] = o_mat.stress[1];
		mat.stress[2] = o_mat.stress[2];
		mat.stress[3] = o_mat.stress[3];
		mat.stress[4] = o_mat.stress[4];
		mat.stress[5] = o_mat.stress[5];
		mat.e = o_mat.e;
		mat.substp_size = o_mat.substp_size;
	}

	void SandHypoplasticityWrapper::set_param(
		const __Float_Type__* ini_stress,
		__Float_Type__ e0,
		__Float_Type__ phi,
		__Float_Type__ hs, __Float_Type__ n,
		__Float_Type__ ei0, __Float_Type__ ec0, __Float_Type__ ed0,
		__Float_Type__ alpha, __Float_Type__ beta)
	{
		glb.set_param(phi, hs, n, alpha, beta, ei0, ec0, ed0, 1000.0, 0.48);
		mat.stress[0] = ini_stress[0];
		mat.stress[1] = ini_stress[1];
		mat.stress[2] = ini_stress[2];
		mat.stress[3] = ini_stress[3];
		mat.stress[4] = ini_stress[4];
		mat.stress[5] = ini_stress[5];
		mat.e = e0;
		mat.substp_size = ffmat(1.0);
	}

	int sand_hypoplasticity_wrapper_integration_function(MaterialModel* _self, double dstrain[6])
	{
		SandHypoplasticityWrapper& self = *static_cast<SandHypoplasticityWrapper*>(_self);
		self.status_code = integrate_sand_hypoplasticity(self.glb, self.mat, dstrain, ffmat(1.0));
		return self.status_code;
	}
}
