#ifndef __Model_FEM_T3D_ME_s_hdf5_utilities_h__
#define __Model_FEM_T3D_ME_s_hdf5_utilities_h__

#include "hdf5.h"
#include "ResultFile_hdf5.h"
#include "Model_FEM_T3D_ME_s.h"

namespace Model_FEM_T3D_ME_s_hdf5_utilities
{
struct NodeData
{
	unsigned long long id;
	double x, y, z;
	double ux, uy, uz;
	double vx, vy, vz;
	void from_node(Model_FEM_T3D_ME_s::Node &n)
	{
		id = n.id;
		x = n.x;
		y = n.y;
		z = n.z;
		ux = n.ux;
		uy = n.uy;
		uz = n.uz;
		vx = n.vx;
		vy = n.vy;
		vz = n.vz;
	}
	void to_node(Model_FEM_T3D_ME_s::Node &n)
	{
		n.id = id;
		n.x = x;
		n.y = y;
		n.z = z;
		n.ux = ux;
		n.uy = uy;
		n.uz = uz;
		n.vx = vx;
		n.vy = vy;
		n.vz = vz;
	}
};

inline hid_t get_node_dt_id()
{
	hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(NodeData));
	H5Tinsert(res, "id", HOFFSET(NodeData, id), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "x", HOFFSET(NodeData, x), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "y", HOFFSET(NodeData, y), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "z", HOFFSET(NodeData, z), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "ux", HOFFSET(NodeData, ux), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "uy", HOFFSET(NodeData, uy), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "uz", HOFFSET(NodeData, uz), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "vx", HOFFSET(NodeData, vx), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "vy", HOFFSET(NodeData, vy), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "vz", HOFFSET(NodeData, vz), H5T_NATIVE_DOUBLE);
	return res;
}

struct ElemData
{
	unsigned long long id;
	unsigned long long n1, n2, n3, n4;
	double vol, density;
	// gauss points
	unsigned long long p1_id, p2_id, p3_id, p4_id;
	void from_elem(Model_FEM_T3D_ME_s::Element& e)
	{
		id = e.id;
		n1 = e.n1;
		n2 = e.n2;
		n3 = e.n3;
		n4 = e.n4;
		vol = e.vol;
		density = e.density;
		p1_id = e.p1->id;
		p2_id = e.p2->id;
		p3_id = e.p3->id;
		p4_id = e.p4->id;
	}
	void to_elem(
		Model_FEM_T3D_ME_s::Element &e,
		Model_FEM_T3D_ME_s::Particle *pcls)
	{
		e.id = id;
		e.n1 = n1;
		e.n2 = n2;
		e.n3 = n3;
		e.n4 = n4;
		e.vol = e.vol;
		e.density = e.density;
		e.p1 = pcls + p1_id;
		e.p2 = pcls + p2_id;
		e.p3 = pcls + p3_id;
		e.p4 = pcls + p4_id;
	}
};

inline hid_t get_elem_dt_id()
{
	hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(ElemData));
	H5Tinsert(res, "id", HOFFSET(ElemData, id), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "n1", HOFFSET(ElemData, n1), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "n2", HOFFSET(ElemData, n2), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "n3", HOFFSET(ElemData, n3), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "n4", HOFFSET(ElemData, n4), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "vol", HOFFSET(ElemData, vol), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "density", HOFFSET(ElemData, density), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "p1_id", HOFFSET(ElemData, p1_id), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "p2_id", HOFFSET(ElemData, p2_id), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "p3_id", HOFFSET(ElemData, p3_id), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "p4_id", HOFFSET(ElemData, p4_id), H5T_NATIVE_ULLONG);
	return res;
}

struct ParticleData
{
	unsigned long long id;
	double s11;
	double s22;
	double s33;
	double s12;
	double s23;
	double s31;
	double e11;
	double e22;
	double e33;
	double e12;
	double e23;
	double e31;
	void from_pcl(Model_FEM_T3D_ME_s::Particle &pcl)
	{
		id = pcl.id;
		s11 = pcl.s11;
		s22 = pcl.s22;
		s33 = pcl.s33;
		s12 = pcl.s12;
		s23 = pcl.s23;
		s31 = pcl.s31;
		e11 = pcl.e11;
		e22 = pcl.e22;
		e33 = pcl.e33;
		e12 = pcl.e12;
		e23 = pcl.e23;
		e31 = pcl.e31;
	}
	void to_pcl(Model_FEM_T3D_ME_s::Particle &pcl)
	{
		pcl.id = id;
		pcl.s11 = s11;
		pcl.s22 = s22;
		pcl.s33 = s33;
		pcl.s12 = s12;
		pcl.s23 = s23;
		pcl.s31 = s31;
		pcl.e11 = e11;
		pcl.e22 = e22;
		pcl.e33 = e33;
		pcl.e12 = e12;
		pcl.e23 = e23;
		pcl.e31 = e31;
	}
};

inline hid_t get_pcl_dt_id()
{
	hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(ParticleData));
	H5Tinsert(res, "id", HOFFSET(ParticleData, id), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "s11", HOFFSET(ParticleData, s11), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "s22", HOFFSET(ParticleData, s22), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "s33", HOFFSET(ParticleData, s33), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "s12", HOFFSET(ParticleData, s12), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "s23", HOFFSET(ParticleData, s23), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "s31", HOFFSET(ParticleData, s31), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "e11", HOFFSET(ParticleData, e11), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "e22", HOFFSET(ParticleData, e22), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "e33", HOFFSET(ParticleData, e33), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "e12", HOFFSET(ParticleData, e12), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "e23", HOFFSET(ParticleData, e23), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "e31", HOFFSET(ParticleData, e31), H5T_NATIVE_DOUBLE);
	return res;
}

int output_mesh_to_hdf5_file(Model_FEM_T3D_ME_s &md, ResultFile_hdf5 &rf, hid_t grp_id);
int load_mesh_from_hdf5_file(Model_FEM_T3D_ME_s &md, ResultFile_hdf5 &rf, hid_t grp_id);

int output_boundary_condition_to_hdf5_file(Model_FEM_T3D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);
int load_boundary_condition_from_hdf5_file(Model_FEM_T3D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);

int output_material_model_to_hdf5_file(Model_FEM_T3D_ME_s &md, ResultFile_hdf5 &rf, hid_t grp_id);
int load_material_model_from_hdf5_file(Model_FEM_T3D_ME_s &md, ResultFile_hdf5 &rf, hid_t grp_id);

// output the whole model to ModelData
int output_model_to_hdf5_file(Model_FEM_T3D_ME_s& md, ResultFile_hdf5 &rf);

// output the particle data and material models to hdf5 (used by time history)
int time_history_complete_output_to_hdf5_file(Model_FEM_T3D_ME_s& md, ResultFile_hdf5 &rf, hid_t frame_grp_id);

// load model data from hdf5 to model data
int load_model_from_hdf5_file(Model_FEM_T3D_ME_s &md, const char *hdf5_name, const char *th_name, size_t frame_id);
};

#endif