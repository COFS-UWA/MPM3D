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
		double const _stress[6], double _e,
		double _fric_ang,
		double _hs, double _en,
		double _ed0, double _ec0, double _ei0,
		double _alpha, double _beta,
		double _m_R, double _m_T, double _R, double _beta_r, double _chi,
		double const _ig_strain[6],
		double _Kw, double _pore_pressure)
	{
		stress[0] = _stress[0];
		stress[1] = _stress[1];
		stress[2] = _stress[2];
		stress[3] = _stress[3];
		stress[4] = _stress[4];
		stress[5] = _stress[5];

		// material properties: props[16]
		// --- hypoplasticity parameters ---
		friction_angle = _fric_ang; // 0
		pore_pressure_on_stress = 0.0; // 1
		hs = _hs; // 2
		en = _en; // 3
		ed0 = _ed0; // 4
		ec0 = _ec0; // 5
		ei0 = _ei0; // 6
		alpha = _alpha; // 7
		beta = _beta; // 8
		// --- intergranular parameters ---
		m_R = _m_R; // 9
		m_T = _m_T; // 10
		// maximum intergranular strain
		R = _R; // 11
		beta_r = _beta_r; // 12
		chi = _chi; // 13
		// --- bulk modulus of pore fluid ---
		Kw = _Kw; // 14
		// --- initial void ratio ---
		init_e = _e; // 15

		// state variables: statev[13]
		// --- integranular strain ---
		if (_ig_strain)
		{
			statev[0] = _ig_strain[0]; // 0
			statev[1] = _ig_strain[1]; // 1
			statev[2] = _ig_strain[2]; // 2
			statev[3] = _ig_strain[3]; // 3
			statev[4] = _ig_strain[4]; // 4
			statev[5] = _ig_strain[5]; // 5
		}
		else
		{
			statev[0] = 0.0; // 0
			statev[1] = 0.0; // 1
			statev[2] = 0.0; // 2
			statev[3] = 0.0; // 3
			statev[4] = 0.0; // 4
			statev[5] = 0.0; // 5
		}
		// void ratio
		e = _e; // 6
		// e = 0.0;
		// pore pressure
		negative_pore_pressure = -_pore_pressure; // 7

		mean_stress = 0.0; // 8
		func_eval_num = 0.0; // 9
		Matsuka_Nakai_friction_angle = 0.0; // 10
		intergranular_strain_norm_len = 0.0; // 11

		// substep size for integration
		integration_time_step = 1.0; // 12

		// get initial stiffness matrix Dep
		// dstrain_e = { 0, 0, 0, 0, 0, 0 }
		integrate(dstrain_e);
	}

	inline static void swap(double& a, double& b)
	{ double tmp = a; a = b; b = tmp; }

	int SandHypoplasticityByUmat_integration_func(
		MaterialModel* _self,
		double dstrain[6])
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

		// integrate the model
		double stress_umat[6] = {
			self.stress[0],
			self.stress[1],
			self.stress[2],
			self.stress[3],
			self.stress[5],
			self.stress[4]
		};
		double ddsdde[6][6];
		swap(self.statev[4], self.statev[5]);
		(*dll_inte_func)(
			stress_umat,
			ddsdde,
			self.statev,
			dstrain,
			self.props);

		swap(self.statev[4], self.statev[5]);

		self.dstress[0] = stress_umat[0] - self.stress[0];
		self.dstress[1] = stress_umat[1] - self.stress[1];
		self.dstress[2] = stress_umat[2] - self.stress[2];
		self.dstress[3] = stress_umat[3] - self.stress[3];
		self.dstress[4] = stress_umat[5] - self.stress[4];
		self.dstress[5] = stress_umat[4] - self.stress[5];

		self.stress[0] = stress_umat[0];
		self.stress[1] = stress_umat[1];
		self.stress[2] = stress_umat[2];
		self.stress[3] = stress_umat[3];
		self.stress[4] = stress_umat[5];
		self.stress[5] = stress_umat[4];

		// cal Dep
		for (size_t i = 0; i < 6; ++i)
			for (size_t j = 0; j < 6; j++)
				self.Dep_mat[i][j] = ddsdde[j][i];
		for (size_t i = 0; i < 6; ++i)
			swap(self.Dep_mat[4][i], self.Dep_mat[5][i]);
		for (size_t i = 0; i < 6; ++i)
			swap(self.Dep_mat[i][4], self.Dep_mat[i][5]);
		for (size_t i = 0; i < 6; ++i)
			for (size_t j = 3; j < 6; ++j)
				self.Dep_mat[i][j] *= 2.0; // ????

		return 0;
	}
}
