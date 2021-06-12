#include "MaterialModels_pcp.h"

#include "MatModelConstants.h"
#include "NorsandWrapper.h"

namespace MatModel
{
	NorsandWrapper::NorsandWrapper(
		const NorsandWrapper& other)
	{
		const auto& o_glb = other.glb;
		NorsandGlobal_set_param(glb, o_glb.phi,
			o_glb.gamma, o_glb.lambda,
			o_glb.N, o_glb.chi, o_glb.H,
			o_glb.Ig, o_glb.niu);
		const auto& o_mat = other.mat;
		mat.stress[0] = o_mat.stress[0];
		mat.stress[1] = o_mat.stress[1];
		mat.stress[2] = o_mat.stress[2];
		mat.stress[3] = o_mat.stress[3];
		mat.stress[4] = o_mat.stress[4];
		mat.stress[5] = o_mat.stress[5];
		mat.e = o_mat.e;
		mat.pi = o_mat.pi;
	}

	void NorsandWrapper::set_param(
		const double* ini_stress,
		double e0,
		double phi,
		double gamma, double lambda,
		double N, double chi, double H,
		double Ig, double niu)
	{
		NorsandGlobal_set_param(glb, phi, gamma, lambda, N, chi, H, Ig, niu);
		Norsand_set_NC_param(mat, glb, ini_stress, e0);
	}

	int nor_sand_wrapper_integration_function(
		MaterialModel* _self,
		double dstrain[6])
	{
		NorsandWrapper &self = *static_cast<NorsandWrapper *>(_self);
		const double ori_stress[6] = {
			self.mat.s11,
			self.mat.s22,
			self.mat.s33,
			self.mat.s12,
			self.mat.s23,
			self.mat.s31
		};
		const int status_code = integrate_norsand(
			self.glb,
			self.mat,
			dstrain,
			self.dstrain_e,
			self.dstrain_p);
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
		return status_code;
	}
}
