#include "MaterialModels_pcp.h"

#include "MatModelConstants.h"
#include "SandHypoplasticityStbWrapper.h"

namespace MatModel
{
	SandHypoplasticityStbWrapper::SandHypoplasticityStbWrapper() :
		MaterialModel(&sand_hypoplasticity_stb_wrapper_integration_function,
			"SandHypoplasticityStbWrapper",
			&sand_hypoplasticity_stb_wrapper_store_to_function,
			&sand_hypoplasticity_stb_wrapper_retrieve_from_function)
	{
		dee11 = 0.0;
		dee22 = 0.0;
		dee33 = 0.0;
		dee12 = 0.0;
		dee23 = 0.0;
		dee31 = 0.0;
	}
	
	SandHypoplasticityStbWrapper::SandHypoplasticityStbWrapper(
		const SandHypoplasticityStbWrapper& other) :
		MaterialModel(sand_hypoplasticity_stb_wrapper_integration_function,
			"SandHypoplasticityStbWrapper")
	{
		const auto& o_glb = other.glb;
		SandHypoplasticityStbGlobal_set_param(glb,
			o_glb.phi, o_glb.hs, o_glb.n,
			o_glb.alpha, o_glb.beta,
			o_glb.ei0, o_glb.ec0, o_glb.ed0,
			o_glb.N, glb.chi, glb.H,
			o_glb.Ig, o_glb.niu,
			o_glb.ten_E, o_glb.ten_niu, 1.0);
		const auto& o_mat = other.mat;
		mat.stress[0] = o_mat.stress[0];
		mat.stress[1] = o_mat.stress[1];
		mat.stress[2] = o_mat.stress[2];
		mat.stress[3] = o_mat.stress[3];
		mat.stress[4] = o_mat.stress[4];
		mat.stress[5] = o_mat.stress[5];
		mat.e = o_mat.e;
		mat.substp_size = o_mat.substp_size;
		mat.Mi = o_mat.Mi;
		mat.pi = o_mat.pi;
		mat.pl = o_mat.pl;
		mat.strain[0] = 0.0;
		mat.strain[1] = 0.0;
		mat.strain[2] = 0.0;
		mat.strain[3] = 0.0;
		mat.strain[4] = 0.0;
		mat.strain[5] = 0.0;
	}

	void SandHypoplasticityStbWrapper::set_param(
		const __Float_Type__* ini_stress,
		__Float_Type__ e0,
		__Float_Type__ phi,
		__Float_Type__ hs, __Float_Type__ n,
		__Float_Type__ alpha, __Float_Type__ beta,
		__Float_Type__ ed0, __Float_Type__ ec0, __Float_Type__ ei0,
		__Float_Type__ N, __Float_Type__ chi, __Float_Type__ H,
		__Float_Type__ Ig, __Float_Type__ niu,
		__Float_Type__ alpha_vol,
		const __Float_Type__* ini_strain)
	{
		SandHypoplasticityStbGlobal_set_param(
			glb, phi, hs, n,
			alpha, beta,
			ed0, ec0, ei0,
			N, chi, H,
			Ig, niu,
			1000.0, 0.0,
			alpha_vol);
		SandHypoplasticityStb_set_NC_param(
			mat, glb,
			ini_stress,
			e0, ffmat(1.0));
		if (ini_strain)
		{
			mat.e11 = ini_strain[0];
			mat.e22 = ini_strain[1];
			mat.e33 = ini_strain[2];
			mat.e12 = ini_strain[3];
			mat.e23 = ini_strain[4];
			mat.e31 = ini_strain[5];
		}
		else
		{
			mat.e11 = 0.0; mat.e22 = 0.0; mat.e33 = 0.0;
			mat.e12 = 0.0; mat.e23 = 0.0; mat.e31 = 0.0;
		}
	}

	void SandHypoplasticityStbWrapper::set_yield_surface(
		__Float_Type__ Mi, __Float_Type__ pi, __Float_Type__ pl)
	{
		mat.Mi = Mi; mat.pi = pi; mat.pl = pl;
	}

	int sand_hypoplasticity_stb_wrapper_integration_function(MaterialModel* _self, double dstrain[6])
	{
		SandHypoplasticityStbWrapper& self = *static_cast<SandHypoplasticityStbWrapper*>(_self);
		const double ori_stress[6] = {
			self.mat.s11,
			self.mat.s22,
			self.mat.s33,
			self.mat.s12,
			self.mat.s23,
			self.mat.s31
		};
		self.status_code = integrate_sand_hypoplasticity_stb(self.glb, self.mat, dstrain, ffmat(1.0));
		self.ds11 = self.mat.s11 - ori_stress[0];
		self.ds22 = self.mat.s22 - ori_stress[1];
		self.ds33 = self.mat.s33 - ori_stress[2];
		self.ds12 = self.mat.s12 - ori_stress[3];
		self.ds23 = self.mat.s23 - ori_stress[4];
		self.ds31 = self.mat.s31 - ori_stress[5];
		self.dep11 = dstrain[0];
		self.dep22 = dstrain[1];
		self.dep33 = dstrain[2];
		self.dep12 = dstrain[3];
		self.dep23 = dstrain[4];
		self.dep31 = dstrain[5];
		return self.status_code;
	}

	void sand_hypoplasticity_stb_wrapper_store_to_function(const MaterialModel* _self, char* mat_model_mem)
	{
		memcpy(mat_model_mem, _self, sizeof(SandHypoplasticityStbWrapper));
	}

	void sand_hypoplasticity_stb_wrapper_retrieve_from_function(MaterialModel* _self, const char* mat_model_mem)
	{
		memcpy(_self, mat_model_mem, sizeof(SandHypoplasticityStbWrapper));
	}
}
