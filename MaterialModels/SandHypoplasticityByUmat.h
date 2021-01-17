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
		union
		{
			struct
			{
				// intergranular strain
				union
				{
					struct
					{
						double delta11; // 0
						double delta22; // 1
						double delta33; // 2
						double delta12; // 3
						double delta23; // 5
						double delta31; // 4
					};
					double delta[6];
				};
				// void ratio
				double e; // 6
				double negative_pore_pressure; // 7
				// 8 ~ 11 parameters needn't be specified
				// (s11 + s22 + s33) / 3.0
				double mean_stress; // 8
				// number of function evaluated (nfev)
				// not clear about its function, may relate to debugging
				double func_eval_num; // 9
				// friction angle of Matsuka ¨C Nakai yield surface
				double Matsuka_Nakai_friction_angle; // 10
				// normalized length of intergranular strain ¦Ñ = ¡¬¦Ä¡¬ / R
				double intergranular_strain_norm_len; // 11
				// time step of constitutive model integration
				double integration_time_step; // 12
			};
			// state variables
			double statev[13];
		};

		union
		{
			struct
			{
				double friction_angle; // 0
				// pore pressure posed on stress
				// used as apparent cohesion
				// compressive as positive
				double pore_pressure_on_stress; // 1
				double hs; // 2
				double en; // 3
				double ed0; // 4
				double ec0; // 5
				double ei0; // 6
				double alpha; // 7
				double beta; // 8
				double m_R; // 9
				double m_T; // 10
				// maximum intergranular strain
				double R; // 11
				double beta_r; // 12
				double chi; // 13
				// bulk modulus of water
				double Kw; // 14
				double init_e; // 15 (not used here)
			};
			double props[16];
		};

	public:
		SandHypoplasticityByUmat();
		~SandHypoplasticityByUmat();

		void set_param(const double _stress[6],
			double _e,
			double _fric_ang /* in degree */,
			double _hs, double _en,
			double _ed0, double _ec0, double _ei0,
			double _alpha, double _beta,
			double _m_R, double _m_T, double _R, double _beta_r, double _chi,
			double const _ig_strain[6] = nullptr, /* intergranular strain */
			double _Kw = 0.0, double _pore_pressure = 0.0);

		inline double get_friction_angle() const noexcept { return friction_angle; }
		inline double get_apparent_cohesion() const noexcept { return pore_pressure_on_stress; }
		inline double get_hs() const noexcept { return hs; }
		inline double get_en() const noexcept { return en; }
		inline double get_ed0() const noexcept { return ed0; }
		inline double get_ec0() const noexcept { return ec0; }
		inline double get_ei0() const noexcept { return ei0; }
		inline double get_alpha() const noexcept { return alpha; }
		inline double get_beta() const noexcept { return beta; }
		inline double get_m_R() const noexcept { return m_R; }
		inline double get_m_T() const noexcept { return m_T; }
		inline double get_R() const noexcept { return R; }
		inline double get_beta_r() const noexcept { return beta_r; }
		inline double get_chi() const noexcept { return chi; }
		inline double get_Kw() const noexcept { return Kw; }
		inline double get_e() const noexcept { return e != 0.0 ? e : init_e; }
		inline const double* get_delta() const noexcept { return delta; }
		inline double get_delta11() const noexcept { return delta11; }
		inline double get_delta22() const noexcept { return delta22; }
		inline double get_delta33() const noexcept { return delta33; }
		inline double get_delta12() const noexcept { return delta12; }
		inline double get_delta31() const noexcept { return delta31; }
		inline double get_delta23() const noexcept { return delta23; }
		inline double get_pore_pressure() const noexcept { return -negative_pore_pressure; }
		inline double get_integration_step_ratio() const noexcept { return integration_time_step; }

		inline void set_apparent_cohesion(double _c) noexcept
		{ pore_pressure_on_stress = _c; }
		inline void set_integration_step_ratio(double _ratio) noexcept
		{ 
			if (_ratio > 1.0)
				_ratio = 1.0;
			if (_ratio < 0.0)
				_ratio = 0.0;
			integration_time_step = _ratio;
		}

	public:
		static int init(const char* dll_path);
		static void free();
	};
}

#endif