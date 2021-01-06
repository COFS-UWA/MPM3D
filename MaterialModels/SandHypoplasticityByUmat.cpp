#include "MaterialModels_pcp.h"

#include <exception>
#include <windows.h>

#include "SandHypoplasticityByUmat.h"

extern "C"
{
	typedef void (_cdecl *SandHypoIntegrateFunc)(
		double stress[6], double ddsdde[6][6], double statev[13],
		double dstran[6], double props[16]);
}

static HINSTANCE umat_dll = nullptr;
static SandHypoIntegrateFunc dll_inte_func = nullptr;

namespace MatModel
{

	int SandHypoplasticityByUmat::init(const char* dll_path)
	{
		if (umat_dll && dll_inte_func)
			return 0;

		umat_dll = LoadLibrary(dll_path);
		if (!umat_dll)
			throw std::exception("SandHypoplasticityByUmat::"
				"can't load library\n");

		// resolve function address here
		dll_inte_func = (SandHypoIntegrateFunc)GetProcAddress(
			umat_dll, "SAND_HYPOPLASTICITY_INTEGRATION");
		if (!dll_inte_func)
		{
			FreeLibrary(umat_dll);
			umat_dll = nullptr;
			throw std::exception("SandHypoplasticityByUmat::"
				"can't load function SAND_HYPOPLASTICITY_INTEGRATION\n");
		}
		return 0;
	}

	void SandHypoplasticityByUmat::free()
	{
		dll_inte_func = nullptr;
		if (umat_dll)
		{
			FreeLibrary(umat_dll);
			umat_dll = nullptr;
		}
	}

	SandHypoplasticityByUmat::SandHypoplasticityByUmat() :
		MaterialModel(&SandHypoplasticityByUmat_integration_func,
			"SandHypoplasticityByUmat")
	{
		// Assume all strain are plastic strain
		dstrain_e[0] = 0.0;	dstrain_e[1] = 0.0;
		dstrain_e[2] = 0.0;	dstrain_e[3] = 0.0;
		dstrain_e[4] = 0.0;	dstrain_e[5] = 0.0;

		for (size_t i = 0; i < 6; ++i)
			for (size_t j = 0; j < 6; ++j)
				De_mat[i][j] = 0.0;
		De_mat[0][0] = 1.0;	De_mat[1][1] = 1.0;
		De_mat[2][2] = 1.0;	De_mat[3][3] = 1.0;
		De_mat[4][4] = 1.0;	De_mat[5][5] = 1.0;
	}

	SandHypoplasticityByUmat::~SandHypoplasticityByUmat() {}

	void SandHypoplasticityByUmat::set_param(
		double _stress[6], double _e,
		double _fric_ang,
		double _hs, double _en,
		double _ed0, double _ec0, double _ei0,
		double _alpha, double _beta,
		double _m_R, double _m_T, double _R, double _beta_r, double _chi)
	{
		stress[0] = _stress[0]; stress[1] = _stress[1];
		stress[2] = _stress[2];	stress[3] = _stress[3];
		stress[4] = _stress[4];	stress[5] = _stress[5];

		stress_umat[0] = _stress[0];
		stress_umat[1] = _stress[1];
		stress_umat[2] = _stress[2];
		stress_umat[3] = _stress[3];
		stress_umat[4] = _stress[5]; // s31
		stress_umat[5] = _stress[4]; // s23

		// material properties
		// initial void ratio
		e = _e;
		props[15] = _e;

		//hypoplasticity parameters
		props[0] = _fric_ang;
		props[2] = _hs;
		props[3] = _en;
		props[4] = _ed0;
		props[5] = _ec0;
		props[6] = _ei0;
		props[7] = _alpha;
		props[8] = _beta;

		// intergranular parameters
		props[9] = _m_R;
		props[10] = _m_T;
		// maximum intergranular strain
		props[11] = _R;
		props[12] = _beta_r;
		props[13] = _chi;

		// pore fluid parameters
		// initial hydrostatic pore pressure
		props[1] = 0.0;
		// bulk modulus of pore fluid
		props[14] = 0.0;

		// integranular strain
		statev[0] = 0.0;
		statev[1] = 0.0;
		statev[2] = 0.0;
		statev[3] = 0.0;
		statev[4] = 0.0;
		statev[5] = 0.0;
		// void ratio
		statev[6] = 0.0;
		statev[7] = 0.0;
		statev[8] = 0.0;
		statev[9] = 0.0;
		statev[10] = 0.0;
		statev[11] = 0.0;
		// substep size for integration
		statev[12] = 1.0;

		// get initial stiffness matrix
		integrate(dstrain_e);
	}

	namespace
	{
		inline void swap(double& a, double& b)
		{
			double tmp = a;
			a = b;
			b = tmp;
		}
	}

	int SandHypoplasticityByUmat_integration_func(MaterialModel* _self, double dstrain[6])
	{
		SandHypoplasticityByUmat& self = *static_cast<SandHypoplasticityByUmat*>(_self);

		self.dstrain_p[0] = dstrain[0];
		self.dstrain_p[1] = dstrain[1];
		self.dstrain_p[2] = dstrain[2];
		self.dstrain_p[3] = dstrain[3];
		self.dstrain_p[4] = dstrain[4];
		self.dstrain_p[5] = dstrain[5];

		// convert dstrain to Abaqus standard format
		dstrain[3] *= 2.0;
		dstrain[4] *= 2.0;
		dstrain[5] *= 2.0;
		swap(dstrain[4], dstrain[5]);

		double ddsdde[6][6];
		(*dll_inte_func)(
			self.stress_umat,
			ddsdde,
			self.statev,
			dstrain,
			self.props);

		self.delta[0] = self.statev[0];
		self.delta[1] = self.statev[1];
		self.delta[2] = self.statev[2];
		self.delta[3] = self.statev[3];
		self.delta[4] = self.statev[4];
		self.delta[5] = self.statev[5];
		self.e = self.statev[6];

		self.dstress[0] = self.stress_umat[0] - self.stress[0];
		self.dstress[1] = self.stress_umat[1] - self.stress[1];
		self.dstress[2] = self.stress_umat[2] - self.stress[2];
		self.dstress[3] = self.stress_umat[3] - self.stress[3];
		self.dstress[4] = self.stress_umat[5] - self.stress[4];
		self.dstress[5] = self.stress_umat[4] - self.stress[5];

		self.stress[0] = self.stress_umat[0];
		self.stress[1] = self.stress_umat[1];
		self.stress[2] = self.stress_umat[2];
		self.stress[3] = self.stress_umat[3];
		self.stress[4] = self.stress_umat[5];
		self.stress[5] = self.stress_umat[4];

		for (size_t i = 0; i < 6; ++i)
			for (size_t j = 0; j < 6; j++)
				self.Dep_mat[i][j] = ddsdde[j][i];
		for (size_t i = 0; i < 6; ++i)
			swap(self.Dep_mat[4][i], self.Dep_mat[5][i]);
		for (size_t i = 0; i < 6; ++i)
			swap(self.Dep_mat[i][4], self.Dep_mat[i][5]);

		return 0;
	}

}