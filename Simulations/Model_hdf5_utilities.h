#ifndef __Model_hdf5_utilities_H__
#define __Model_hdf5_utilities_H__

#include "hdf5.h"
#include "MatModelContainer.h"

namespace Model_hdf5_utilities
{
// 2D Background mesh
struct Node2DData
{
	unsigned long long id;
	double x;
	double y;
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
	unsigned long long n1;
	unsigned long long n2;
	unsigned long long n3;
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
	double x;
	double y;
	double z;
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
	unsigned long long n1;
	unsigned long long n2;
	unsigned long long n3;
	unsigned long long n4;
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
};

inline hid_t get_bf_pcl_dt_id()
{
	hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(BodyForceAtPclData));
	H5Tinsert(res, "pcl_id", HOFFSET(BodyForceAtPclData, pcl_id), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "bf", HOFFSET(BodyForceAtPclData, bf), H5T_NATIVE_DOUBLE);
	return res;
}

struct TractionBCAtPclData
{
	unsigned long long pcl_id;
	double t;
};

inline hid_t get_tbc_pcl_dt_id()
{
	hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(TractionBCAtPclData));
	H5Tinsert(res, "pcl_id", HOFFSET(TractionBCAtPclData, pcl_id), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "t", HOFFSET(TractionBCAtPclData, t), H5T_NATIVE_DOUBLE);
	return res;
}

struct AccelerationBCData
{
	unsigned long long node_id;
	double a;
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
	unsigned long long pcl_id;
	double E; // Young's modulus
	double niu; // possion ratio
	inline void from_mm(MatModel::LinearElasticity &mm)
	{
		E = mm.E;
		niu = mm.niu;
	}
	inline void to_mm(MatModel::LinearElasticity &mm)
	{
		mm.set_param(E, niu);
	}
};

inline hid_t get_le_hdf5_dt_id()
{
	hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(LinearElasticityStateData));
	H5Tinsert(res, "pcl_id", HOFFSET(LinearElasticityStateData, pcl_id), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "E", HOFFSET(LinearElasticityStateData, E), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "niu", HOFFSET(LinearElasticityStateData, niu), H5T_NATIVE_DOUBLE);
	return res;
}

// ModifiedCamClay
struct ModifiedCamClayStateData
{
	unsigned long long pcl_id;
	double niu; // possion ratio
	double kappa; // logrithmic recompression modulus
	double lambda; // logrithmic compression modulus
	double fric_angle; // friction angle
	double e; // void ratio
	double pc; // pre-consolidation stress
	double N; // normal consolidation line
	double Gamma; // critical state line
	double M; // critical state line q = Mp 
	double s11, s22, s33, s12, s23, s31; // stress state
	inline void from_mm(MatModel::ModifiedCamClay &mm)
	{
		niu = mm.niu;
		kappa = mm.kappa;
		lambda = mm.lambda;
		fric_angle = mm.fric_angle;
		e = mm.e;
		pc = mm.pc;
		N = mm.N;
		Gamma = mm.Gamma;
		M = mm.M;
		const double *stress = mm.get_stress();
		s11 = stress[0];
		s22 = stress[1];
		s33 = stress[2];
		s12 = stress[3];
		s23 = stress[4];
		s31 = stress[5];
	}
	inline void to_mm(MatModel::ModifiedCamClay &mm)
	{
		double stress[6] = { s11, s22, s33, s12, s23, s31 };
		fric_angle = fric_angle / 3.14159265359 * 180.0;
		mm.set_param_OC(niu, kappa, lambda, fric_angle, e, stress, pc);
	}
};

inline hid_t get_mcc_hdf5_dt_id()
{
	hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(ModifiedCamClayStateData));
	H5Tinsert(res, "pcl_id", HOFFSET(ModifiedCamClayStateData, pcl_id), H5T_NATIVE_ULLONG);
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

};

#endif