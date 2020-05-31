#ifndef __Model_T3D_CHM_s_hdf5_utilities_h__
#define __Model_T3D_CHM_s_hdf5_utilities_h__

#include "hdf5.h"
#include "ResultFile_hdf5.h"
#include "Model_T3D_CHM_s.h"

namespace Model_T3D_CHM_s_hdf5_utilities
{
struct ParticleData
{
	unsigned long long id;
	double n;
	double m_s;
	double density_s;
	double density_f;
	double x;
	double y;
	double z;
	double vol;
	double vx_s;
	double vy_s;
	double vz_s;
	double vx_f;
	double vy_f;
	double vz_f;
	double s11;
	double s22;
	double s33;
	double s12;
	double s23;
	double s31;
	double p;
	double e11;
	double e22;
	double e33;
	double e12;
	double e23;
	double e31;
	void from_pcl(Model_T3D_CHM_s::Particle &pcl)
	{
		id = pcl.id;
		n = pcl.n;
		m_s = pcl.m_s;
		density_s = pcl.density_s;
		density_f = pcl.density_f;
		x = pcl.x;
		y = pcl.y;
		z = pcl.z;
		vol = pcl.get_vol();
		vx_s = pcl.vx_s;
		vy_s = pcl.vy_s;
		vz_s = pcl.vz_s;
		vx_f = pcl.vx_f;
		vy_f = pcl.vy_f;
		vz_f = pcl.vz_f;
		s11 = pcl.s11;
		s22 = pcl.s22;
		s33 = pcl.s33;
		s12 = pcl.s12;
		s23 = pcl.s23;
		s31 = pcl.s31;
		p = pcl.p;
		e11 = pcl.e11;
		e22 = pcl.e22;
		e33 = pcl.e33;
		e12 = pcl.e12;
		e23 = pcl.e23;
		e31 = pcl.e31;
	}
	void to_pcl(Model_T3D_CHM_s::Particle &pcl)
	{
		pcl.id = id;
		pcl.n = n;
		pcl.m_s = m_s;
		pcl.density_s = density_s;
		pcl.density_f = density_f;
		pcl.x = x;
		pcl.y = y;
		pcl.z = z;
		pcl.vx_s = vx_s;
		pcl.vy_s = vy_s;
		pcl.vz_s = vz_s;
		pcl.vx_f = vx_f;
		pcl.vy_f = vy_f;
		pcl.vz_f = vz_f;
		pcl.s11 = s11;
		pcl.s22 = s22;
		pcl.s33 = s33;
		pcl.s12 = s12;
		pcl.s23 = s23;
		pcl.s31 = s31;
		pcl.p = p;
		pcl.e11 = e11;
		pcl.e22 = e22;
		pcl.e33 = e33;
		pcl.e12 = e12;
		pcl.e23 = e23;
		pcl.e31 = e31;
	}
};

inline hid_t get_pcl_dt_id(void)
{
	hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(ParticleData));
	H5Tinsert(res, "id", HOFFSET(ParticleData, id), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "n", HOFFSET(ParticleData, n), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "m_s", HOFFSET(ParticleData, m_s), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "density_s", HOFFSET(ParticleData, density_s), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "density_f", HOFFSET(ParticleData, density_f), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "x", HOFFSET(ParticleData, x), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "y", HOFFSET(ParticleData, y), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "z", HOFFSET(ParticleData, z), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "vol", HOFFSET(ParticleData, vol), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "vx_s", HOFFSET(ParticleData, vx_s), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "vy_s", HOFFSET(ParticleData, vy_s), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "vz_s", HOFFSET(ParticleData, vz_s), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "vx_f", HOFFSET(ParticleData, vx_f), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "vy_f", HOFFSET(ParticleData, vy_f), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "vz_f", HOFFSET(ParticleData, vz_f), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "s11", HOFFSET(ParticleData, s11), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "s22", HOFFSET(ParticleData, s22), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "s33", HOFFSET(ParticleData, s33), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "s12", HOFFSET(ParticleData, s12), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "s23", HOFFSET(ParticleData, s23), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "s31", HOFFSET(ParticleData, s31), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "p", HOFFSET(ParticleData, p), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "e11", HOFFSET(ParticleData, e11), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "e22", HOFFSET(ParticleData, e22), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "e33", HOFFSET(ParticleData, e33), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "e12", HOFFSET(ParticleData, e12), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "e23", HOFFSET(ParticleData, e23), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "e31", HOFFSET(ParticleData, e31), H5T_NATIVE_DOUBLE);
	return res;
}

int output_background_mesh_to_hdf5_file(Model_T3D_CHM_s& md, ResultFile_hdf5& rf, hid_t grp_id);
int load_background_mesh_from_hdf5_file(Model_T3D_CHM_s& md, ResultFile_hdf5& rf, hid_t grp_id);

int output_boundary_condition_to_hdf5_file(Model_T3D_CHM_s& md, ResultFile_hdf5& rf, hid_t grp_id);
int load_boundary_condition_from_hdf5_file(Model_T3D_CHM_s& md, ResultFile_hdf5& rf, hid_t grp_id);

int output_pcl_data_to_hdf5_file(Model_T3D_CHM_s& md, ResultFile_hdf5& rf, hid_t grp_id);
int load_pcl_data_from_hdf5_file(Model_T3D_CHM_s& md, ResultFile_hdf5& rf, hid_t grp_id);

int output_material_model_to_hdf5_file(Model_T3D_CHM_s& md, ResultFile_hdf5& rf, hid_t grp_id);
int load_material_model_from_hdf5_file(Model_T3D_CHM_s& md, ResultFile_hdf5& rf, hid_t grp_id);

// output the whole model to ModelData
int output_model_to_hdf5_file(Model_T3D_CHM_s& md, ResultFile_hdf5& rf);

// output the particle data and material models to hdf5 (used by time history)
int time_history_complete_output_to_hdf5_file(Model_T3D_CHM_s& md, ResultFile_hdf5& rf, hid_t frame_grp_id);

// load model data from hdf5 to model data
int load_model_from_hdf5_file(Model_T3D_CHM_s& md, const char* hdf5_name, const char* th_name, size_t frame_id);
};

#endif