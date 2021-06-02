#include "MaterialModels_pcp.h"

#include "MohrCoulombWrapper.h"

namespace MatModel
{
	MohrCoulombWrapper::MohrCoulombWrapper(const MohrCoulombWrapper& other) :
		MaterialModel(mohr_coulomb_wrapper_integration_function, "MohrCoulomb")
	{
		const auto& o_glb = other.glb;
		glb.set_param(o_glb.phi, o_glb.psi,
			o_glb.cohesion, o_glb.E, o_glb.niu);
		const auto& o_mat = other.mat;
		mat.stress[0] = o_mat.stress[0];
		mat.stress[1] = o_mat.stress[1];
		mat.stress[2] = o_mat.stress[2];
		mat.stress[3] = o_mat.stress[3];
		mat.stress[4] = o_mat.stress[4];
		mat.stress[5] = o_mat.stress[5];
		stress[0] = other.stress[0];
		stress[1] = other.stress[1];
		stress[2] = other.stress[2];
		stress[3] = other.stress[3];
		stress[4] = other.stress[4];
		stress[5] = other.stress[5];
	}

	void MohrCoulombWrapper::set_param(
		const __Float_Type__ ini_stress[6],
		__Float_Type__ _phi,
		__Float_Type__ _psi,
		__Float_Type__ _cohesion,
		__Float_Type__ _E, __Float_Type__ _niu)
	{
		glb.set_param(_phi, _psi, _cohesion, _E, _niu);
		mat.stress[0] = ini_stress[0];
		mat.stress[1] = ini_stress[1];
		mat.stress[2] = ini_stress[2];
		mat.stress[3] = ini_stress[3];
		mat.stress[4] = ini_stress[4];
		mat.stress[5] = ini_stress[5];
		stress[0] = ini_stress[0];
		stress[1] = ini_stress[1];
		stress[2] = ini_stress[2];
		stress[3] = ini_stress[3];
		stress[4] = ini_stress[4];
		stress[5] = ini_stress[5];
	}

	int mohr_coulomb_wrapper_integration_function(MaterialModel* _self, double dstrain[6])
	{
		MohrCoulombWrapper& self = *static_cast<MohrCoulombWrapper *>(_self);
		self.status_code = integrate_mohr_coulomb(
			self.glb,
			self.mat,
			dstrain,
			self.dstrain_e,
			self.dstrain_p);
		self.ds11 = self.mat.s11 - self.s11;
		self.ds22 = self.mat.s22 - self.s22;
		self.ds33 = self.mat.s33 - self.s33;
		self.ds12 = self.mat.s12 - self.s12;
		self.ds23 = self.mat.s23 - self.s23;
		self.ds31 = self.mat.s31 - self.s31;
		self.s11 = self.mat.s11;
		self.s22 = self.mat.s22;
		self.s33 = self.mat.s33;
		self.s12 = self.mat.s12;
		self.s23 = self.mat.s23;
		self.s31 = self.mat.s31;
		return self.status_code;
	}
}
