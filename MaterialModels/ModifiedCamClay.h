#ifndef __Modified_Cam_Clay_h__
#define __Modified_Cam_Clay_h__

#include <cmath>

#include "MaterialModel.h"

namespace MatModel
{

	int modified_cam_clay_integration_function(MaterialModel* _self, double dstrain[6]);

	// This model uses constant stiffness instead of
	// deducing it from recompression line
	class ModifiedCamClay : public MaterialModel
	{
		friend int modified_cam_clay_integration_function(MaterialModel* _self, double dstrain[6]);
	public:
		double niu; // possion ratio
		double kappa; // logrithmic recompression modulus
		double lambda; // logrithmic compression modulus
		double N;
		double fric_angle, M, M2; // critical state line parameters
		double Gamma;
		double pc; // pre-consolidation stress
		double e; // void ratio

		ModifiedCamClay() :
			MaterialModel(modified_cam_clay_integration_function, "ModifiedCamClay"),
			e(0.0), pc(0.0), niu(0.0),
			kappa(0.0), lambda(0.0), fric_angle(0.0), M2(0.0)
		{
			for (size_t i = 0; i < 6; ++i)
			{
				stress[i] = 0.0;
				dstress[i] = 0.0;
				dstrain_e[i] = 0.0;
				dstrain_p[i] = 0.0;
			}
		}
		~ModifiedCamClay() {}

		// angle unit is degree
		void set_param_NC(
			double _niu,
			double _kappa,
			double _lambda,
			double _fric_ang,
			double _N, // void ratio at p = kPa on NCL
			double _s[6]
		)
		{
			// elastic
			niu = _niu;
			kappa = _kappa;
			// NCL
			lambda = _lambda;
			N = _N;
			// CSL
			Gamma = N - (lambda - kappa) * log(2.0);
			fric_angle = _fric_ang;
			double fric_ang_tmp = _fric_ang / 180.0 * 3.14159265359;
			M = 6.0 * sin(fric_ang_tmp) / (3.0 - sin(fric_ang_tmp));
			M2 = M * M;
			// stress
			memcpy(stress, _s, sizeof(double) * 6);
			// pre-consolidation stress
			double p = cal_p();
			double q = cal_q();
			pc = cal_pc(p, q);
			// void ratio
			e = N + (kappa - lambda) * log(pc) - kappa * log(-p);
			// stiffness mat
			form_De_mat(p);
			//form_Dep_mat();
			double dg_ds[6], dg_dpc;
			double A_pc, divider;
			cal_dg_stress(p, q, dg_ds, dg_dpc);
			A_pc = cal_A_pc(dg_ds);
			divider = cal_divider(dg_ds, dg_dpc, A_pc);
			form_Dep_mat(dg_ds, divider);
		}

		void set_param_OC(double _niu,
			double _kappa, double _lambda, double _fric_ang,
			double _N, double _s[6], double _pc)
		{
			// elastic
			niu = _niu;
			kappa = _kappa;
			// NCL
			lambda = _lambda;
			N = _N;
			// CSL
			Gamma = N - (lambda - kappa) * log(2.0);
			fric_angle = _fric_ang;
			double fric_ang_tmp = _fric_ang / 180.0 * 3.14159265359;
			M = 6.0 * sin(fric_ang_tmp) / (3.0 - sin(fric_ang_tmp));
			M2 = M * M;
			// stress
			memcpy(stress, _s, sizeof(double) * 6);
			// pre-consolidation stress
			double p = cal_p();
			double q = cal_q();
			pc = cal_pc(p, q);
			if (_pc > pc) // overconsolidated
			{
				pc = _pc;
				// void ratio
				e = N + (kappa - lambda) * log(pc) - kappa * log(-p);
				// stiffness mat
				form_De_mat(p);
				form_Dep_mat();
			}
			else // normally consolidated
			{
				// void ratio
				//e = N - lambda * log(pc);
				e = N + (kappa - lambda) * log(pc) - kappa * log(-p);
				// stiffness mat
				form_De_mat(p);
				double dg_ds[6], dg_dpc;
				double A_pc, divider;
				cal_dg_stress(p, q, dg_ds, dg_dpc);
				A_pc = cal_A_pc(dg_ds);
				divider = cal_divider(dg_ds, dg_dpc, A_pc);
				form_Dep_mat(dg_ds, divider);
			}
		}

		inline double get_p() const noexcept { return cal_p(); }
		inline double get_q() const noexcept { return cal_q(); }
		inline double get_pc() const noexcept { return pc; }
		inline double get_e_by_strain() const noexcept { return e; }
		inline double get_e_by_model() const noexcept
		{
			return N - lambda * log(pc) + kappa * log(-pc / cal_p());
		}
		inline double get_f() const noexcept { return cal_f(cal_p(), cal_q()); }
		inline double get_norm_f() const noexcept { return cal_norm_f(cal_p(), cal_q()); }

	protected:
		// form elastic stiffness matrix
		inline void form_De_mat(double p)
		{
			double K_mod = (1.0 + e) * -p / kappa;
			double E_mod = 3.0 * (1.0 - 2.0 * niu) * K_mod;

			double coef;
			coef = (1.0 - niu) / ((1.0 + niu) * (1.0 - 2.0 * niu)) * E_mod;
			De_mat[0][0] = coef;
			De_mat[1][1] = coef;
			De_mat[2][2] = coef;

			coef = niu / ((1.0 + niu) * (1.0 - 2.0 * niu)) * E_mod;
			De_mat[0][1] = coef;
			De_mat[0][2] = coef;
			De_mat[1][0] = coef;
			De_mat[1][2] = coef;
			De_mat[2][0] = coef;
			De_mat[2][1] = coef;

			De_mat[0][3] = 0.0;
			De_mat[0][4] = 0.0;
			De_mat[0][5] = 0.0;
			De_mat[1][3] = 0.0;
			De_mat[1][4] = 0.0;
			De_mat[1][5] = 0.0;
			De_mat[2][3] = 0.0;
			De_mat[2][4] = 0.0;
			De_mat[2][5] = 0.0;

			De_mat[3][0] = 0.0;
			De_mat[3][1] = 0.0;
			De_mat[3][2] = 0.0;
			De_mat[4][0] = 0.0;
			De_mat[4][1] = 0.0;
			De_mat[4][2] = 0.0;
			De_mat[5][0] = 0.0;
			De_mat[5][1] = 0.0;
			De_mat[5][2] = 0.0;

			coef = E_mod / (1.0 + niu);
			De_mat[3][3] = coef;
			De_mat[4][4] = coef;
			De_mat[5][5] = coef;

			De_mat[3][4] = 0.0;
			De_mat[3][5] = 0.0;
			De_mat[4][3] = 0.0;
			De_mat[4][5] = 0.0;
			De_mat[5][3] = 0.0;
			De_mat[5][4] = 0.0;
		}
		inline void form_Dep_mat() noexcept
		{
			memcpy(Dep_mat, De_mat, 6 * 6 * sizeof(double));
		}
		inline void form_Dep_mat(double dg_ds[6], double divider) noexcept
		{
			double De_dg_ds[6];
			for (size_t i = 0; i < 6; ++i)
				De_dg_ds[i] = De_mat[i][0] * dg_ds[0] + De_mat[i][1] * dg_ds[1]
				+ De_mat[i][2] * dg_ds[2] + De_mat[i][3] * dg_ds[3]
				+ De_mat[i][4] * dg_ds[4] + De_mat[i][5] * dg_ds[5];
			for (size_t i = 0; i < 6; ++i)
				for (size_t j = 0; j < 6; ++j)
					Dep_mat[i][j] = De_mat[i][j] - De_dg_ds[i] * De_dg_ds[j] / divider;
		}

		inline double cal_p() const noexcept
		{
			return (s11 + s22 + s33) / 3.0;
		}
		inline double cal_q() const noexcept
		{
			double s11_s22_diff = s11 - s22;
			double s22_s33_diff = s22 - s33;
			double s33_s11_diff = s33 - s11;
			double q2 = (s11_s22_diff * s11_s22_diff
				+ s22_s33_diff * s22_s33_diff
				+ s33_s11_diff * s33_s11_diff) * 0.5
				+ (s12 * s12 + s23 * s23 + s31 * s31) * 3.0;
			return sqrt(q2);
		}
		// equivalent pre-consolidation stress
		inline double cal_pc(double p, double q) const noexcept
		{
			return -(q * q / (M2 * p) + p);
		}
		inline double cal_f(double p, double q) const noexcept
		{
			return q * q + M2 * p * (p + pc);
		}
		// normalized f - covergence craterion
		inline double cal_norm_f(double f, double p) const noexcept
		{
			return -f / (M2 * p * pc);
		}
		void cal_dg_stress(double p, double q, double dg_ds[6], double& dg_dpc) const noexcept
		{
			double dg_dp, dg_dq;
			dg_dp = M2 * (2.0 * p + pc);
			dg_dq = 2.0; // *q in dq_ds
			dg_dpc = M2 * p;
			// dq_ds
			double dq_ds11, dq_ds22, dq_ds33;
			double dq_ds12, dq_ds23, dq_ds31;
			dq_ds11 = 0.5 * (s11 + s11 - s22 - s33);
			dq_ds22 = 0.5 * (s22 + s22 - s33 - s11);
			dq_ds33 = 0.5 * (s33 + s33 - s11 - s22);
			dq_ds12 = 3.0 * s12;
			dq_ds23 = 3.0 * s23;
			dq_ds31 = 3.0 * s31;
			// dg_ds
			dg_ds[0] = dg_dp / 3.0 + dg_dq * dq_ds11; // dg_ds11
			dg_ds[1] = dg_dp / 3.0 + dg_dq * dq_ds22; // dg_ds22
			dg_ds[2] = dg_dp / 3.0 + dg_dq * dq_ds33; // dg_ds33
			dg_ds[3] = dg_dq * dq_ds12;               // dg_ds12
			dg_ds[4] = dg_dq * dq_ds23;               // dg_ds23
			dg_ds[5] = dg_dq * dq_ds31;               // dg_ds31
		}
		inline double cal_A_pc(double dg_ds[6]) const noexcept
		{
			return -(1.0 + e) * pc / (lambda - kappa) * (dg_ds[0] + dg_ds[1] + dg_ds[2]);
		}
		inline double cal_divider(double dg_ds[6], double dg_dpc, double A_pc) const noexcept
		{
			double divider = 0.0;
			for (size_t i = 0; i < 6; i++)
				for (size_t j = 0; j < 6; j++)
					divider += dg_ds[i] * De_mat[i][j] * dg_ds[j];
			divider -= dg_dpc * A_pc;
			return divider;
		}
	};

}

#endif