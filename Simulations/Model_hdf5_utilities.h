#ifndef __Model_hdf5_utilities_h__
#define __Model_hdf5_utilities_h__

#include "hdf5.h"
#include "BCs.h"
#include "MatModelContainer.h"
#include "ResultFile_hdf5.h"

namespace Model_hdf5_utilities
{
	// 2D Background mesh
	struct Node2DData
	{
		unsigned long long id;
		double x, y;
	};

	inline hid_t get_nd_2d_dt_id()
	{
		hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(Node2DData));
		H5Tinsert(res, "id", HOFFSET(Node2DData, id), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "x", HOFFSET(Node2DData, x), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "y", HOFFSET(Node2DData, y), H5T_NATIVE_DOUBLE);
		return res;
	}

	struct Elem2DData
	{
		unsigned long long id;
		unsigned long long n1, n2, n3;
	};

	inline hid_t get_ed_2d_dt_id()
	{
		hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(Elem2DData));
		H5Tinsert(res, "id", HOFFSET(Elem2DData, id), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "n1", HOFFSET(Elem2DData, n1), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "n2", HOFFSET(Elem2DData, n2), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "n3", HOFFSET(Elem2DData, n3), H5T_NATIVE_ULLONG);
		return res;
	}

	// 3D Background mesh
	struct Node3DData
	{
		unsigned long long id;
		double x, y, z;
	};

	inline hid_t get_nd_3d_dt_id()
	{
		hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(Node3DData));
		H5Tinsert(res, "id", HOFFSET(Node3DData, id), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "x", HOFFSET(Node3DData, x), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "y", HOFFSET(Node3DData, y), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "z", HOFFSET(Node3DData, z), H5T_NATIVE_DOUBLE);
		return res;
	}

	struct Elem3DData
	{
		unsigned long long id;
		unsigned long long n1, n2, n3, n4;
	};

	inline hid_t get_ed_3d_dt_id()
	{
		hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(Elem3DData));
		H5Tinsert(res, "id", HOFFSET(Elem3DData, id), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "n1", HOFFSET(Elem3DData, n1), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "n2", HOFFSET(Elem3DData, n2), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "n3", HOFFSET(Elem3DData, n3), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "n4", HOFFSET(Elem3DData, n4), H5T_NATIVE_ULLONG);
		return res;
	}

	// Boundary condition
	struct BodyForceAtPclData
	{
		unsigned long long pcl_id;
		double bf;
		inline void from_bf(BodyForceAtPcl& _bf)
		{
			pcl_id = _bf.pcl_id;
			bf = _bf.bf;
		}
		inline void to_bf(BodyForceAtPcl& _bf)
		{
			_bf.pcl_id = pcl_id;
			_bf.bf = bf;
		}
	};

	inline hid_t get_bf_pcl_dt_id()
	{
		hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(BodyForceAtPclData));
		H5Tinsert(res, "pcl_id", HOFFSET(BodyForceAtPclData, pcl_id), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "bf", HOFFSET(BodyForceAtPclData, bf), H5T_NATIVE_DOUBLE);
		return res;
	}

	struct BodyForceAtElemData
	{
		unsigned long long elem_id;
		double bf;
	};

	inline hid_t get_bf_elem_dt_id()
	{
		hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(BodyForceAtElemData));
		H5Tinsert(res, "elem_id", HOFFSET(BodyForceAtElemData, elem_id), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "bf", HOFFSET(BodyForceAtElemData, bf), H5T_NATIVE_DOUBLE);
		return res;
	}

	struct TractionBCAtPclData
	{
		unsigned long long pcl_id;
		double t;
		inline void from_tbc(TractionBCAtPcl& tbc)
		{
			pcl_id = tbc.pcl_id;
			t = tbc.t;
		}
		inline void to_tbc(TractionBCAtPcl &tbc)
		{
			tbc.pcl_id = pcl_id;
			tbc.t = t;
		}
	};

	inline hid_t get_tbc_pcl_dt_id()
	{
		hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(TractionBCAtPclData));
		H5Tinsert(res, "pcl_id", HOFFSET(TractionBCAtPclData, pcl_id), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "t", HOFFSET(TractionBCAtPclData, t), H5T_NATIVE_DOUBLE);
		return res;
	}

	struct TractionBCAtFaceData
	{
		unsigned long long elem_id;
		unsigned long long face_id;
		double t;
	};

	inline hid_t get_tbc_face_dt_id()
	{
		hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(TractionBCAtFaceData));
		H5Tinsert(res, "elem_id", HOFFSET(TractionBCAtFaceData, elem_id), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "face_id", HOFFSET(TractionBCAtFaceData, face_id), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "t", HOFFSET(TractionBCAtFaceData, t), H5T_NATIVE_DOUBLE);
		return res;
	}

	struct AccelerationBCData
	{
		unsigned long long node_id;
		double a;
		inline void from_abc(AccelerationBC &abc)
		{
			node_id = abc.node_id;
			a = abc.a;
		}
		inline void to_abc(AccelerationBC &abc)
		{
			abc.node_id = node_id;
			abc.a = a;
		}
	};

	inline hid_t get_abc_dt_id()
	{
		hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(AccelerationBCData));
		H5Tinsert(res, "node_id", HOFFSET(AccelerationBCData, node_id), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "a", HOFFSET(AccelerationBCData, a), H5T_NATIVE_DOUBLE);
		return res;
	}

	struct VelocityBCData
	{
		unsigned long long node_id;
		double v;
		inline void from_vbc(VelocityBC &vbc)
		{
			node_id = vbc.node_id;
			v = vbc.v;
		}
		inline void to_vbc(VelocityBC &vbc)
		{
			vbc.node_id = node_id;
			vbc.v = v;
		}
	};

	inline hid_t get_vbc_dt_id()
	{
		hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(VelocityBCData));
		H5Tinsert(res, "node_id", HOFFSET(VelocityBCData, node_id), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "v", HOFFSET(VelocityBCData, v), H5T_NATIVE_DOUBLE);
		return res;
	}

	// Material model
	// Linear elasticity
	struct LinearElasticityStateData
	{
		unsigned long long id;
		double E; // Young's modulus
		double niu; // possion ratio
		double s11, s22, s33, s12, s23, s31;
		inline void from_mm(MatModel::LinearElasticity &mm)
		{
			id = mm.get_id();
			E = mm.E;
			niu = mm.niu;
			const double* stress = mm.get_stress();
			s11 = stress[0];
			s22 = stress[1];
			s33 = stress[2];
			s12 = stress[3];
			s23 = stress[4];
			s31 = stress[5];
		}
		inline void to_mm(MatModel::LinearElasticity &mm)
		{
			mm.set_id(id);
			mm.set_param(E, niu);
			double stress[6] = { s11, s22, s33, s12, s23, s31 };
			mm.set_stress(stress);
		}
	};

	inline hid_t get_le_hdf5_dt_id()
	{
		hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(LinearElasticityStateData));
		H5Tinsert(res, "id", HOFFSET(LinearElasticityStateData, id), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "E", HOFFSET(LinearElasticityStateData, E), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "niu", HOFFSET(LinearElasticityStateData, niu), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s11", HOFFSET(LinearElasticityStateData, s11), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s22", HOFFSET(LinearElasticityStateData, s22), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s33", HOFFSET(LinearElasticityStateData, s33), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s12", HOFFSET(LinearElasticityStateData, s12), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s23", HOFFSET(LinearElasticityStateData, s23), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s31", HOFFSET(LinearElasticityStateData, s31), H5T_NATIVE_DOUBLE);
		return res;
	}

	// ModifiedCamClay
	struct ModifiedCamClayStateData
	{
		unsigned long long id;
		double niu; // possion ratio
		double kappa; // logrithmic recompression modulus
		double lambda; // logrithmic compression modulus
		double fric_angle; // friction angle
		double e; // void ratio
		double pc; // pre-consolidation stress
		double N; // normal consolidation line
		double Gamma; // critical state line
		double M; // critical state line q = Mp
		union // stress state
		{
			double stress[6];
			struct { double s11, s22, s33, s12, s23, s31; };
		};

		inline void from_mm(MatModel::ModifiedCamClay &mm)
		{
			id = mm.get_id();
			niu = mm.niu;
			kappa = mm.kappa;
			lambda = mm.lambda;
			fric_angle = mm.fric_angle;
			e = mm.e;
			pc = mm.pc;
			N = mm.N;
			Gamma = mm.Gamma;
			M = mm.M;
			const double * mm_stress = mm.get_stress();
			s11 = mm_stress[0];
			s22 = mm_stress[1];
			s33 = mm_stress[2];
			s12 = mm_stress[3];
			s23 = mm_stress[4];
			s31 = mm_stress[5];
		}

		inline void to_mm(MatModel::ModifiedCamClay &mm)
		{
			mm.set_id(id);
			mm.set_param_OC(niu, kappa, lambda, fric_angle, N, stress, pc);
			mm.e = e;
			mm.pc = pc;
		}
	};

	inline hid_t get_mcc_hdf5_dt_id()
	{
		hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(ModifiedCamClayStateData));
		H5Tinsert(res, "id", HOFFSET(ModifiedCamClayStateData, id), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "niu", HOFFSET(ModifiedCamClayStateData, niu), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "kappa", HOFFSET(ModifiedCamClayStateData, kappa), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "lambda", HOFFSET(ModifiedCamClayStateData, lambda), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "fric_angle", HOFFSET(ModifiedCamClayStateData, fric_angle), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "e", HOFFSET(ModifiedCamClayStateData, e), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "pc", HOFFSET(ModifiedCamClayStateData, pc), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "N", HOFFSET(ModifiedCamClayStateData, N), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "Gamma", HOFFSET(ModifiedCamClayStateData, Gamma), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "M", HOFFSET(ModifiedCamClayStateData, M), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s11", HOFFSET(ModifiedCamClayStateData, s11), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s22", HOFFSET(ModifiedCamClayStateData, s22), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s33", HOFFSET(ModifiedCamClayStateData, s33), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s12", HOFFSET(ModifiedCamClayStateData, s12), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s23", HOFFSET(ModifiedCamClayStateData, s23), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s31", HOFFSET(ModifiedCamClayStateData, s31), H5T_NATIVE_DOUBLE);
		return res;
	}

	// UndrainedModifiedCamClay
	struct UndrainedModifiedCamClayStateData
	{
		unsigned long long id;
		double niu; // possion ratio
		double kappa; // logrithmic recompression modulus
		double lambda; // logrithmic compression modulus
		double fric_angle; // friction angle
		double e; // void ratio
		double pc; // pre-consolidation stress
		double N; // normal consolidation line
		double Gamma; // critical state line
		double M; // critical state line q = Mp 
		union // effective stress state
		{
			double stress[6];
			struct { double s11, s22, s33, s12, s23, s31; };
		};
		double Kw; // bulk modulus
		double pore_pressure; // pore pressure

		inline void from_mm(MatModel::UndrainedModifiedCamClay& mm)
		{
			id = mm.get_id();
			MatModel::ModifiedCamClay& mcc = mm.mcc;
			niu = mcc.niu;
			kappa = mcc.kappa;
			lambda = mcc.lambda;
			fric_angle = mcc.fric_angle;
			e = mcc.e;
			pc = mcc.pc;
			N = mcc.N;
			Gamma = mcc.Gamma;
			M = mcc.M;
			const double* mm_stress = mm.get_stress();
			s11 = mm_stress[0];
			s22 = mm_stress[1];
			s33 = mm_stress[2];
			s12 = mm_stress[3];
			s23 = mm_stress[4];
			s31 = mm_stress[5];
			Kw = mm.Kw;
			pore_pressure = mm.pore_pressure;
		}

		inline void to_mm(MatModel::UndrainedModifiedCamClay& mm)
		{
			mm.set_id(id);
			mm.set_param_OC(niu, kappa, lambda, fric_angle, N, stress, pc,
							Kw, pore_pressure);
		}
	};

	inline hid_t get_undrained_mcc_hdf5_dt_id()
	{
		hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(UndrainedModifiedCamClayStateData));
		H5Tinsert(res, "id", HOFFSET(UndrainedModifiedCamClayStateData, id), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "niu", HOFFSET(UndrainedModifiedCamClayStateData, niu), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "kappa", HOFFSET(UndrainedModifiedCamClayStateData, kappa), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "lambda", HOFFSET(UndrainedModifiedCamClayStateData, lambda), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "fric_angle", HOFFSET(UndrainedModifiedCamClayStateData, fric_angle), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "e", HOFFSET(UndrainedModifiedCamClayStateData, e), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "pc", HOFFSET(UndrainedModifiedCamClayStateData, pc), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "N", HOFFSET(UndrainedModifiedCamClayStateData, N), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "Gamma", HOFFSET(UndrainedModifiedCamClayStateData, Gamma), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "M", HOFFSET(UndrainedModifiedCamClayStateData, M), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s11", HOFFSET(UndrainedModifiedCamClayStateData, s11), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s22", HOFFSET(UndrainedModifiedCamClayStateData, s22), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s33", HOFFSET(UndrainedModifiedCamClayStateData, s33), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s12", HOFFSET(UndrainedModifiedCamClayStateData, s12), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s23", HOFFSET(UndrainedModifiedCamClayStateData, s23), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s31", HOFFSET(UndrainedModifiedCamClayStateData, s31), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "Kw", HOFFSET(UndrainedModifiedCamClayStateData, Kw), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "pore_pressure", HOFFSET(UndrainedModifiedCamClayStateData, pore_pressure), H5T_NATIVE_DOUBLE);
		return res;
	}

	// Von Mises
	struct VonMisesStateData
	{
		unsigned long long id;
		double E; // Young's modulus
		double niu; // possion ratio
		double cohesion;
		union
		{
			double stress[6];
			struct { double s11, s22, s33, s12, s23, s31; };
		};
		inline void from_mm(MatModel::VonMises &mm)
		{
			id = mm.get_id();
			E = mm.E;
			niu = mm.niu;
			cohesion = mm.cohesion;
			const double* mm_stress = mm.get_stress();
			s11 = mm_stress[0];
			s22 = mm_stress[1];
			s33 = mm_stress[2];
			s12 = mm_stress[3];
			s23 = mm_stress[4];
			s31 = mm_stress[5];
		}
		inline void to_mm(MatModel::VonMises &mm)
		{
			mm.set_id(id);
			mm.set_param(E, niu, cohesion, stress);
		}
	};

	inline hid_t get_von_mises_hdf5_dt_id()
	{
		hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(VonMisesStateData));
		H5Tinsert(res, "id", HOFFSET(VonMisesStateData, id), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "E", HOFFSET(VonMisesStateData, E), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "niu", HOFFSET(VonMisesStateData, niu), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "cohesion", HOFFSET(VonMisesStateData, cohesion), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s11", HOFFSET(VonMisesStateData, s11), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s22", HOFFSET(VonMisesStateData, s22), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s33", HOFFSET(VonMisesStateData, s33), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s12", HOFFSET(VonMisesStateData, s12), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s23", HOFFSET(VonMisesStateData, s23), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s31", HOFFSET(VonMisesStateData, s31), H5T_NATIVE_DOUBLE);
		return res;
	}

	// Tresca
	struct TrescaStateData
	{
		unsigned long long id;
		double E; // Young's modulus
		double niu; // possion ratio
		double cohesion;
		union
		{
			double stress[6];
			struct { double s11, s22, s33, s12, s23, s31; };
		};
		inline void from_mm(MatModel::Tresca& mm)
		{
			id = mm.get_id();
			E = mm.get_E();
			niu = mm.get_niu();
			cohesion = mm.get_cohesion();
			const double* mm_stress = mm.get_stress();
			s11 = mm_stress[0];
			s22 = mm_stress[1];
			s33 = mm_stress[2];
			s12 = mm_stress[3];
			s23 = mm_stress[4];
			s31 = mm_stress[5];
		}
		inline void to_mm(MatModel::Tresca& mm)
		{
			mm.set_id(id);
			mm.set_param(E, niu, cohesion, stress);
		}
	};

	inline hid_t get_tresca_hdf5_dt_id()
	{
		hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(TrescaStateData));
		H5Tinsert(res, "id", HOFFSET(TrescaStateData, id), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "E", HOFFSET(TrescaStateData, E), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "niu", HOFFSET(TrescaStateData, niu), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "cohesion", HOFFSET(TrescaStateData, cohesion), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s11", HOFFSET(TrescaStateData, s11), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s22", HOFFSET(TrescaStateData, s22), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s33", HOFFSET(TrescaStateData, s33), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s12", HOFFSET(TrescaStateData, s12), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s23", HOFFSET(TrescaStateData, s23), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s31", HOFFSET(TrescaStateData, s31), H5T_NATIVE_DOUBLE);
		return res;
	}

	// MohrCoulomb
	struct MohrCoulombStateData
	{
		unsigned long long id;
		double phi, psi, cohesion;
		double E, niu;
		union
		{
			double stress[6];
			struct { double s11, s22, s33, s12, s23, s31; };
		};
		int status_code;

		inline void from_mm(MatModel::MohrCoulombWrapper& mm)
		{
			id = mm.get_id();
			phi = mm.get_phi();
			psi = mm.get_psi();
			cohesion = mm.get_cohesion();
			E = mm.get_E();
			niu = mm.get_niu();
			const double* mm_stress = mm.get_stress();
			s11 = mm_stress[0];
			s22 = mm_stress[1];
			s33 = mm_stress[2];
			s12 = mm_stress[3];
			s23 = mm_stress[4];
			s31 = mm_stress[5];
			status_code = mm.get_status_code();
		}

		inline void to_mm(MatModel::MohrCoulombWrapper& mm)
		{
			mm.set_id(id);
			mm.set_param(stress, phi, psi, cohesion, E, niu);
			status_code = 0;
		}
	};

	inline hid_t get_mohr_coulomb_hdf5_dt_id()
	{
		hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(MohrCoulombStateData));
		H5Tinsert(res, "id", HOFFSET(MohrCoulombStateData, id), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "phi", HOFFSET(MohrCoulombStateData, phi), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "psi", HOFFSET(MohrCoulombStateData, psi), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "cohesion", HOFFSET(MohrCoulombStateData, cohesion), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "E", HOFFSET(MohrCoulombStateData, E), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "niu", HOFFSET(MohrCoulombStateData, niu), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "status_code", HOFFSET(MohrCoulombStateData, status_code), H5T_NATIVE_INT);
		return res;
	}

	// SandHypoplasticityByUmat
	struct SandHypoplasticityByUmatStateData
	{
		unsigned long long id;
		double friction_angle;
		double apparent_coehsion;
		double hs, en;
		double ed0, ec0, ei0;
		double alpha, beta;
		double m_R, m_T, R, beta_r, chi;
		union
		{
			struct
			{
				double delta11; // 0
				double delta22; // 1
				double delta33; // 2
				double delta12; // 3
				double delta31; // 4
				double delta23; // 5
			};
			double delta[6];
		};
		double e;
		union
		{
			double stress[6];
			struct { double s11, s22, s33, s12, s23, s31; };
		};
		double pore_pressure, Kw;
		double integration_time_step;
		int status_code;

		inline void from_mm(MatModel::SandHypoplasticityByUmat&mm)
		{
			id = mm.get_id();
			friction_angle = mm.get_friction_angle();
			apparent_coehsion = mm.get_apparent_cohesion();
			hs = mm.get_hs();
			en = mm.get_en();
			ed0 = mm.get_ed0();
			ec0 = mm.get_ec0();
			ei0 = mm.get_ei0();
			alpha = mm.get_alpha();
			beta = mm.get_beta();
			m_R = mm.get_m_R();
			m_T = mm.get_m_T();
			R = mm.get_R();
			beta_r = mm.get_beta_r();
			chi = mm.get_chi();
			const double* mm_delta = mm.get_delta();
			delta[0] = mm_delta[0];
			delta[1] = mm_delta[1];
			delta[2] = mm_delta[2];
			delta[3] = mm_delta[3];
			delta[4] = mm_delta[4];
			delta[5] = mm_delta[5];
			e = mm.get_e();
			const double* mm_stress = mm.get_stress();
			s11 = mm_stress[0];
			s22 = mm_stress[1];
			s33 = mm_stress[2];
			s12 = mm_stress[3];
			s23 = mm_stress[4];
			s31 = mm_stress[5];
			pore_pressure = mm.get_pore_pressure();
			Kw = mm.get_Kw();
			integration_time_step = mm.get_integration_step_ratio();
			status_code = mm.get_status_code();
		}

		inline void to_mm(MatModel::SandHypoplasticityByUmat &mm)
		{
			mm.set_id(id);
			mm.set_apparent_cohesion(apparent_coehsion);
			mm.set_integration_step_ratio(integration_time_step);
			mm.set_param(stress, e, friction_angle, hs, en,
				ed0, ec0, ei0, alpha, beta, m_R, m_T, R, beta_r, chi,
				delta, Kw, pore_pressure);
		}
	};

	inline hid_t get_sand_hypoplasticity_by_umat_hdf5_dt_id()
	{
		hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(SandHypoplasticityByUmatStateData));
		H5Tinsert(res, "id", HOFFSET(SandHypoplasticityByUmatStateData, id), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "friction_angle", HOFFSET(SandHypoplasticityByUmatStateData, friction_angle), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "apparent_coehsion", HOFFSET(SandHypoplasticityByUmatStateData, apparent_coehsion), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "hs", HOFFSET(SandHypoplasticityByUmatStateData, hs), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "en", HOFFSET(SandHypoplasticityByUmatStateData, en), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "ed0", HOFFSET(SandHypoplasticityByUmatStateData, ed0), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "ec0", HOFFSET(SandHypoplasticityByUmatStateData, ec0), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "ei0", HOFFSET(SandHypoplasticityByUmatStateData, ei0), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "alpha", HOFFSET(SandHypoplasticityByUmatStateData, alpha), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "beta", HOFFSET(SandHypoplasticityByUmatStateData, beta), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "m_R", HOFFSET(SandHypoplasticityByUmatStateData, m_R), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "m_T", HOFFSET(SandHypoplasticityByUmatStateData, m_T), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "R", HOFFSET(SandHypoplasticityByUmatStateData, R), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "beta_r", HOFFSET(SandHypoplasticityByUmatStateData, beta_r), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "chi", HOFFSET(SandHypoplasticityByUmatStateData, chi), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "delta11", HOFFSET(SandHypoplasticityByUmatStateData, delta11), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "delta22", HOFFSET(SandHypoplasticityByUmatStateData, delta22), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "delta33", HOFFSET(SandHypoplasticityByUmatStateData, delta33), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "delta12", HOFFSET(SandHypoplasticityByUmatStateData, delta12), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "delta31", HOFFSET(SandHypoplasticityByUmatStateData, delta31), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "delta23", HOFFSET(SandHypoplasticityByUmatStateData, delta23), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "e", HOFFSET(SandHypoplasticityByUmatStateData, e), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s11", HOFFSET(SandHypoplasticityByUmatStateData, s11), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s22", HOFFSET(SandHypoplasticityByUmatStateData, s22), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s33", HOFFSET(SandHypoplasticityByUmatStateData, s33), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s12", HOFFSET(SandHypoplasticityByUmatStateData, s12), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s23", HOFFSET(SandHypoplasticityByUmatStateData, s23), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s31", HOFFSET(SandHypoplasticityByUmatStateData, s31), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "pore_pressure", HOFFSET(SandHypoplasticityByUmatStateData, pore_pressure), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "Kw", HOFFSET(SandHypoplasticityByUmatStateData, Kw), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "integration_time_step", HOFFSET(SandHypoplasticityByUmatStateData, integration_time_step), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "status_code", HOFFSET(SandHypoplasticityByUmatStateData, status_code), H5T_NATIVE_INT);
		return res;
	}

	// SandHypoplasticity
	struct SandHypoplasticityStateData
	{
		unsigned long long id;
		double phi, hs, n;
		double alpha, beta;
		double ed0, ec0, ei0;
		union
		{
			double stress[6];
			struct { double s11, s22, s33, s12, s23, s31; };
		};
		double e;
		double substep_size;
		int status_code;

		inline void from_mm(MatModel::SandHypoplasticityWrapper &mm)
		{
			id = mm.get_id();
			phi = mm.get_phi();
			hs = mm.get_hs();
			n = mm.get_n();
			alpha = mm.get_alpha();
			beta = mm.get_beta();
			ed0 = mm.get_ed0();
			ec0 = mm.get_ec0();
			ei0 = mm.get_ei0();
			const double* mm_stress = mm.get_stress();
			s11 = mm_stress[0];
			s22 = mm_stress[1];
			s33 = mm_stress[2];
			s12 = mm_stress[3];
			s23 = mm_stress[4];
			s31 = mm_stress[5];
			e = mm.get_e();
			substep_size = mm.get_substp_size();
			status_code = mm.get_status_code();
		}

		inline void to_mm(MatModel::SandHypoplasticityWrapper& mm)
		{
			mm.set_id(id);
			mm.set_param(stress, e, phi, hs, n,	ed0, ec0, ei0, alpha, beta);
			mm.set_substep_size(substep_size);
			status_code = 0;
		}
	};
	
	inline hid_t get_sand_hypoplasticity_hdf5_dt_id()
	{
		hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(SandHypoplasticityStateData));
		H5Tinsert(res, "id", HOFFSET(SandHypoplasticityStateData, id), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "phi", HOFFSET(SandHypoplasticityStateData, phi), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "hs", HOFFSET(SandHypoplasticityStateData, hs), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "n", HOFFSET(SandHypoplasticityStateData, n), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "alpha", HOFFSET(SandHypoplasticityStateData, alpha), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "beta", HOFFSET(SandHypoplasticityStateData, beta), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "ed0", HOFFSET(SandHypoplasticityStateData, ed0), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "ec0", HOFFSET(SandHypoplasticityStateData, ec0), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "ei0", HOFFSET(SandHypoplasticityStateData, ei0), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s11", HOFFSET(SandHypoplasticityStateData, s11), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s22", HOFFSET(SandHypoplasticityStateData, s22), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s33", HOFFSET(SandHypoplasticityStateData, s33), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s12", HOFFSET(SandHypoplasticityStateData, s12), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s23", HOFFSET(SandHypoplasticityStateData, s23), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s31", HOFFSET(SandHypoplasticityStateData, s31), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "e", HOFFSET(SandHypoplasticityStateData, e), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "substep_size", HOFFSET(SandHypoplasticityStateData, substep_size), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "status_code", HOFFSET(SandHypoplasticityStateData, status_code), H5T_NATIVE_INT);
		return res;
	}

	// SandHypoplasticityStb
	struct SandHypoplasticityStbStateData
	{
		unsigned long long id;
		double phi, hs, n;
		double alpha, beta;
		double ed0, ec0, ei0;
		double N, chi, H, alpha_vol;
		double Ig, niu;
		union
		{
			double stress[6];
			struct { double s11, s22, s33, s12, s23, s31; };
		};
		double e;
		double substep_size;
		double Mi, pi;
		double pl;
		union
		{
			double strain[6];
			struct { double e11, e22, e33, e12, e23, e31; };
		};
		int status_code;

		inline void from_mm(MatModel::SandHypoplasticityStbWrapper& mm)
		{
			id = mm.get_id();
			phi = mm.get_phi();
			hs = mm.get_hs();
			n = mm.get_n();
			alpha = mm.get_alpha();
			beta = mm.get_beta();
			ed0 = mm.get_ed0();
			ec0 = mm.get_ec0();
			ei0 = mm.get_ei0();
			N = mm.get_N();
			chi = mm.get_chi();
			H = mm.get_H();
			alpha_vol = mm.get_alpha_vol();
			Ig = mm.get_Ig();
			niu = mm.get_niu();
			const double* mm_stress = mm.get_stress();
			s11 = mm_stress[0];
			s22 = mm_stress[1];
			s33 = mm_stress[2];
			s12 = mm_stress[3];
			s23 = mm_stress[4];
			s31 = mm_stress[5];
			e = mm.get_e();
			substep_size = mm.get_substp_size();
			status_code = mm.get_status_code();
			Mi = mm.get_Mi();
			pi = mm.get_pi();
			pl = mm.get_pl();
			const double *mm_strain = mm.get_strain();
			e11 = mm_strain[0];
			e22 = mm_strain[1];
			e33 = mm_strain[2];
			e12 = mm_strain[3];
			e23 = mm_strain[4];
			e31 = mm_strain[5];
		}

		inline void to_mm(MatModel::SandHypoplasticityStbWrapper& mm)
		{
			mm.set_id(id);
			mm.set_param(stress, e,
				phi, hs, n,
				alpha, beta,
				ed0, ec0, ei0,
				N, chi, H,
				Ig, niu, alpha_vol, strain);
			mm.set_substep_size(substep_size);
			mm.set_yield_surface(Mi, pi, pl);
			status_code = 0;
		}
	};

	inline hid_t get_sand_hypoplasticity_stb_hdf5_dt_id()
	{
		hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(SandHypoplasticityStbStateData));
		H5Tinsert(res, "id", HOFFSET(SandHypoplasticityStbStateData, id), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "phi", HOFFSET(SandHypoplasticityStbStateData, phi), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "hs", HOFFSET(SandHypoplasticityStbStateData, hs), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "n", HOFFSET(SandHypoplasticityStbStateData, n), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "alpha", HOFFSET(SandHypoplasticityStbStateData, alpha), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "beta", HOFFSET(SandHypoplasticityStbStateData, beta), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "ed0", HOFFSET(SandHypoplasticityStbStateData, ed0), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "ec0", HOFFSET(SandHypoplasticityStbStateData, ec0), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "ei0", HOFFSET(SandHypoplasticityStbStateData, ei0), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "N", HOFFSET(SandHypoplasticityStbStateData, N), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "chi", HOFFSET(SandHypoplasticityStbStateData, chi), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "H", HOFFSET(SandHypoplasticityStbStateData, H), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "alpha_vol", HOFFSET(SandHypoplasticityStbStateData, alpha_vol), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "Ig", HOFFSET(SandHypoplasticityStbStateData, Ig), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "niu", HOFFSET(SandHypoplasticityStbStateData, niu), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s11", HOFFSET(SandHypoplasticityStbStateData, s11), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s22", HOFFSET(SandHypoplasticityStbStateData, s22), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s33", HOFFSET(SandHypoplasticityStbStateData, s33), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s12", HOFFSET(SandHypoplasticityStbStateData, s12), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s23", HOFFSET(SandHypoplasticityStbStateData, s23), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s31", HOFFSET(SandHypoplasticityStbStateData, s31), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "e", HOFFSET(SandHypoplasticityStbStateData, e), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "substep_size", HOFFSET(SandHypoplasticityStbStateData, substep_size), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "Mi", HOFFSET(SandHypoplasticityStbStateData, Mi), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "pi", HOFFSET(SandHypoplasticityStbStateData, pi), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "pl", HOFFSET(SandHypoplasticityStbStateData, pl), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "e11", HOFFSET(SandHypoplasticityStbStateData, e11), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "e22", HOFFSET(SandHypoplasticityStbStateData, e22), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "e33", HOFFSET(SandHypoplasticityStbStateData, e33), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "e12", HOFFSET(SandHypoplasticityStbStateData, e12), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "e23", HOFFSET(SandHypoplasticityStbStateData, e23), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "e31", HOFFSET(SandHypoplasticityStbStateData, e31), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "status_code", HOFFSET(SandHypoplasticityStbStateData, status_code), H5T_NATIVE_INT);
		return res;
	}

	// Norsand
	struct NorsandStateData
	{
		unsigned long long id;
		double phi;
		double gamma, lambda;
		double N, chi, H;
		double Ig, niu;	
		union
		{
			double stress[6];
			struct { double s11, s22, s33, s12, s23, s31; };
		};
		double e;
		double pi;

		inline void from_mm(MatModel::NorsandWrapper& mm)
		{
			id = mm.get_id();
			phi = mm.get_phi();
			gamma = mm.get_gamma();
			lambda = mm.get_lambda();
			N = mm.get_N();
			chi = mm.get_chi();
			H = mm.get_H();
			Ig = mm.get_Ig();
			niu = mm.get_niu();
			const double* mm_stress = mm.get_stress();
			s11 = mm_stress[0];
			s22 = mm_stress[1];
			s33 = mm_stress[2];
			s12 = mm_stress[3];
			s23 = mm_stress[4];
			s31 = mm_stress[5];
			e = mm.get_e();
			pi = mm.get_pi();
		}

		inline void to_mm(MatModel::NorsandWrapper& mm)
		{
			mm.set_id(id);
			mm.set_param(
				stress, e,
				phi,
				gamma, lambda,
				N, chi, H,
				Ig, niu);
			mm.set_pi(pi);
		}
	};

	inline hid_t get_norsand_dt_id()
	{
		hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(NorsandStateData));
		H5Tinsert(res, "id", HOFFSET(NorsandStateData, id), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "phi", HOFFSET(NorsandStateData, phi), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "gamma", HOFFSET(NorsandStateData, gamma), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "lambda", HOFFSET(NorsandStateData, lambda), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "N", HOFFSET(NorsandStateData, N), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "chi", HOFFSET(NorsandStateData, chi), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "H", HOFFSET(NorsandStateData, H), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "Ig", HOFFSET(NorsandStateData, Ig), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "niu", HOFFSET(NorsandStateData, niu), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s11", HOFFSET(NorsandStateData, s11), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s22", HOFFSET(NorsandStateData, s22), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s33", HOFFSET(NorsandStateData, s33), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s12", HOFFSET(NorsandStateData, s12), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s23", HOFFSET(NorsandStateData, s23), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "s31", HOFFSET(NorsandStateData, s31), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "e", HOFFSET(NorsandStateData, e), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "pi", HOFFSET(NorsandStateData, pi), H5T_NATIVE_DOUBLE);
		return res;
	}

	// material model container
	int output_material_model_container_to_hdf5_file(
		MatModel::MatModelContainer &mc, ResultFile_hdf5 &rf, hid_t mc_grp_id);
	int load_material_model_container_from_hdf5_file(
		MatModel::MatModelContainer &mc, ResultFile_hdf5 &rf, hid_t mc_grp_id);
};

#endif