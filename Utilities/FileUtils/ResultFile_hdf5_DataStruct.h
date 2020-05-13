#ifndef __Result_File_hdf5_DataStruct_H__
#define __Result_File_hdf5_DataStruct_H__

#include "hdf5.h"
#include "ModelContainer.h"

namespace ResultFile_hdf5_DataStruct
{

// Background mesh
struct NodeData
{
	unsigned long long id;
	double x;
	double y;
};

inline hid_t get_nd_dt_id(void)
{
	hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(NodeData));
	H5Tinsert(res, "id", HOFFSET(NodeData, id), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "x", HOFFSET(NodeData, x), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "y", HOFFSET(NodeData, y), H5T_NATIVE_DOUBLE);
	return res;
}

struct ElemData
{
	unsigned long long id;
	unsigned long long n1;
	unsigned long long n2;
	unsigned long long n3;
};

inline hid_t get_ed_dt_id(void)
{
	hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(ElemData));
	H5Tinsert(res, "id", HOFFSET(ElemData, id), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "n1", HOFFSET(ElemData, n1), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "n2", HOFFSET(ElemData, n2), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "n3", HOFFSET(ElemData, n3), H5T_NATIVE_ULLONG);
	return res;
}

// constitutive model
// Linear Elasticity
struct LinearElasticityStateData
{
	unsigned long long pcl_id;
	double E; // Young's modulus
	double niu; // possion ratio
	void from_cm(LinearElasticity &cm)
	{
		E = cm.E;
		niu = cm.niu;
	}
	void to_cm(LinearElasticity &cm)
	{
		cm.set_param(E, niu);
	}
};

inline hid_t get_le_hdf5_dt_id(void)
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
	void from_cm(ModifiedCamClay &cm)
	{
		niu = cm.niu;
		kappa = cm.kappa;
		lambda = cm.lambda;
		fric_angle = cm.fric_angle;
		e = cm.e;
		pc = cm.pc;
		N = cm.N;
		Gamma = cm.Gamma;
		M = cm.M;
		s11 = cm.s11;
		s22 = cm.s22;
		s33 = cm.s33;
		s12 = cm.s12;
		s23 = cm.s23;
		s31 = cm.s31;
	}
	void to_cm(ModifiedCamClay &cm)
	{
		double stress[6] = { s11, s22, s33, s12, s23, s31 };
		fric_angle = fric_angle / 3.14159265359 * 180.0;
		cm.set_param_OC(niu, kappa, lambda, fric_angle, e, stress, pc);
	}
};

inline hid_t get_mcc_hdf5_dt_id(void)
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

// rigid circle
struct RigidBodyParticleData
{
	double xr;
	double yr;
	double vol;
};

inline hid_t get_rc_dt_id(void)
{
	hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(RigidBodyParticleData));
	H5Tinsert(res, "xr",  HOFFSET(RigidBodyParticleData, xr),  H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "yr",  HOFFSET(RigidBodyParticleData, yr),  H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "vol", HOFFSET(RigidBodyParticleData, vol), H5T_NATIVE_DOUBLE);
	return res;
}

};

#endif