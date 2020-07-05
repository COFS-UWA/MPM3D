#ifndef __Model_T2D_ME_s_hdf5_io_utilities_H__
#define __Model_T2D_ME_s_hdf5_io_utilities_H__

#include "hdf5.h"
#include "ResultFile_hdf5.h"
#include "MatModelContainer.h"
#include "RigidCircle.h"
#include "Model_T2D_ME_s.h"

namespace Model_T2D_ME_s_hdf5_utilities
{
struct ParticleData
{
	unsigned long long id;
	double m;
	double density;
	double x;
	double y;
	double vol;
	double vx;
	double vy;
	double s11;
	double s22;
	double s12;
	double e11;
	double e22;
	double e12;
	void from_pcl(Model_T2D_ME_s::Particle &pcl)
	{
		id = pcl.id;
		m = pcl.m;
		density = pcl.density;
		x = pcl.x;
		y = pcl.y;
		vol = pcl.get_vol();
		vx = pcl.vx;
		vy = pcl.vy;
		s11 = pcl.s11;
		s22 = pcl.s22;
		s12 = pcl.s12;
		e11 = pcl.e11;
		e22 = pcl.e22;
		e12 = pcl.e12;
	}
	void to_pcl(Model_T2D_ME_s::Particle &pcl)
	{
		pcl.id = id;
		pcl.m = m;
		pcl.density = density;
		pcl.x = x;
		pcl.y = y;
		pcl.vx = vx;
		pcl.vy = vy;
		pcl.s11 = s11;
		pcl.s22 = s22;
		pcl.s12 = s12;
		pcl.e11 = e11;
		pcl.e22 = e22;
		pcl.e12 = e12;
	}
};

inline hid_t get_pcl_dt_id()
{
	hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(ParticleData));
	H5Tinsert(res, "id", HOFFSET(ParticleData, id), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "m", HOFFSET(ParticleData, m), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "density", HOFFSET(ParticleData, density), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "x", HOFFSET(ParticleData, x), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "y", HOFFSET(ParticleData, y), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "vol", HOFFSET(ParticleData, vol), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "vx", HOFFSET(ParticleData, vx), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "vy", HOFFSET(ParticleData, vy), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "s11", HOFFSET(ParticleData, s11), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "s22", HOFFSET(ParticleData, s22), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "s12", HOFFSET(ParticleData, s12), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "e11", HOFFSET(ParticleData, e11), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "e22", HOFFSET(ParticleData, e22), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "e12", HOFFSET(ParticleData, e12), H5T_NATIVE_DOUBLE);
	return res;
}

int output_background_mesh_to_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);
int load_background_mesh_from_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);

int output_boundary_condition_to_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);
int load_boundary_condition_from_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);

int output_pcl_data_to_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);
int load_pcl_data_from_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);

int output_material_model_to_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);
int load_material_model_from_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);

int output_rigid_circle_to_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);
int load_rigid_circle_to_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);

// output the whole model to ModelData
int output_model_to_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf);

// output the particle data and material models to hdf5 (used by time history)
int time_history_complete_output_to_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t frame_grp_id);

// load model data from hdf5 to model data
int load_me_s_model_from_hdf5_file(Model_T2D_ME_s &md, const char *hdf5_name, const char *th_name,	size_t frame_id);
};

#endif