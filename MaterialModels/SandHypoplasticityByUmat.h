#ifndef __Sand_Hypoplasticity_By_Umat_h__
#define __Sand_Hypoplasticity_By_Umat_h__

#include "MaterialModel.h"

namespace MatModel
{
	int SandHypoplasticityByUmat_integration_func(MaterialModel* _self, double dstrain[6]);

	/* ===== model parameters =====
		fric_angle (in degree)
		hs, en
		ed0, ec0, ei0
		alpha, beta
		m_R, m_T, R,
		beta_r, chi
	   ============================ */
	class SandHypoplasticityByUmat : public MaterialModel
	{
		friend int SandHypoplasticityByUmat_integration_func(MaterialModel* _self, double dstrain[6]);
	
	protected:
		// state variables
		// intergranular strain
		double delta[6];
		// void ratio
		double e;

		// used by umat
		double stress_umat[6];
		double statev[13];
		double props[16];

	public:
		SandHypoplasticityByUmat();
		~SandHypoplasticityByUmat();

		void set_param(double _stress[6], double _e,
			double _fric_ang /* in degree */,
			double _hs, double _en,
			double _ed0, double _ec0, double _ei0,
			double _alpha, double _beta,
			double _m_R, double _m_T, double _R, double _beta_r, double _chi);

		inline double get_e() const { return e; }

	public:
		static int init(const char* dll_path);
		static void free();
	};
}

#endif