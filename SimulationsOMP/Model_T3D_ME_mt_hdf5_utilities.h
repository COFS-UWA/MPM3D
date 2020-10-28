#ifndef __Model_T3D_ME_mt_hdf5_utilities_h__
#define __Model_T3D_ME_mt_hdf5_utilities_h__

#include "ResultFile_hdf5.h"
#include "MatModelContainer.h"
#include "Model_T3D_ME_mt.h"
#include "Step_T3D_ME_mt.h"

namespace Model_T3D_ME_mt_hdf5_utilities
{
struct ParticleData
{
	size_t id;
	double x, y, z;
	double bfx, bfy, bfz;
	double tx, ty, tz;
	double m, density, vol;
	double vx, vy, vz;
	double s11, s22, s33, s12, s23, s31;
	size_t mat_id; // material model id

	void from_pcl(
		Model_T3D_ME_mt& md,
		size_t pcl_offset,
		size_t pcl_sorted_var_id
		)
	{
		Model_T3D_ME_mt::PclSortedVarArray& psva
			= md.pcl_sorted_var_array[pcl_sorted_var_id];
		id = psva.pcl_index[pcl_offset];
		Model_T3D_ME_mt::PclPos& pcl_pos = md.pcl_pos[id];
		Model_T3D_ME_mt::PclDisp& pcl_disp = psva.pcl_disp[pcl_offset];
		x = pcl_pos.x + pcl_disp.ux;
		y = pcl_pos.y + pcl_disp.uy;
		z = pcl_pos.z + pcl_disp.uz;
		Model_T3D_ME_mt::PclBodyForce& pcl_bf = md.pcl_bf[id];
		bfx = pcl_bf.bfx;
		bfy = pcl_bf.bfy;
		bfz = pcl_bf.bfz;
		Model_T3D_ME_mt::PclTraction& pcl_t = md.pcl_t[id];
		tx = pcl_t.tx;
		ty = pcl_t.ty;
		tz = pcl_t.tz;
		m = md.pcl_m[id];
		density = psva.pcl_density[pcl_offset];
		vol = m / density;
		Model_T3D_ME_mt::PclV& pcl_v = psva.pcl_v[pcl_offset];
		vx = pcl_v.vx;
		vy = pcl_v.vy;
		vz = pcl_v.vz;
		Model_T3D_ME_mt::PclStress& pcl_stress = psva.pcl_stress[pcl_offset];
		s11 = pcl_stress.s11;
		s22 = pcl_stress.s22;
		s33 = pcl_stress.s33;
		s12 = pcl_stress.s12;
		s23 = pcl_stress.s23;
		s31 = pcl_stress.s31;
		mat_id = md.pcl_mat_model[id]->get_id();
	}
	
	void from_pcl(
		Model_T3D_ME_mt &md,
		size_t pcl_offset,
		size_t pcl_sorted_var_id,
		const size_t *new_to_ori_pcl_map
		)
	{
		Model_T3D_ME_mt::PclSortedVarArray &psva
			= md.pcl_sorted_var_array[pcl_sorted_var_id ^ 1];
		size_t ori_pcl_id = new_to_ori_pcl_map[pcl_offset];
		id = psva.pcl_index[ori_pcl_id];
		Model_T3D_ME_mt::PclPos& pcl_pos = md.pcl_pos[id];
		Model_T3D_ME_mt::PclDisp& pcl_disp = psva.pcl_disp[ori_pcl_id];
		x = pcl_pos.x + pcl_disp.ux;
		y = pcl_pos.y + pcl_disp.uy;
		z = pcl_pos.z + pcl_disp.uz;
		Model_T3D_ME_mt::PclBodyForce& pcl_bf = md.pcl_bf[id];
		bfx = pcl_bf.bfx;
		bfy = pcl_bf.bfy;
		bfz = pcl_bf.bfz;
		Model_T3D_ME_mt::PclTraction& pcl_t = md.pcl_t[id];
		tx = pcl_t.tx;
		ty = pcl_t.ty;
		tz = pcl_t.tz;
		m = md.pcl_m[id];
		density = psva.pcl_density[ori_pcl_id];
		vol = m / density;
		Model_T3D_ME_mt::PclV& pcl_v = psva.pcl_v[ori_pcl_id];
		vx = pcl_v.vx;
		vy = pcl_v.vy;
		vz = pcl_v.vz;
		Model_T3D_ME_mt::PclStress& pcl_stress = psva.pcl_stress[ori_pcl_id];
		s11 = pcl_stress.s11;
		s22 = pcl_stress.s22;
		s33 = pcl_stress.s33;
		s12 = pcl_stress.s12;
		s23 = pcl_stress.s23;
		s31 = pcl_stress.s31;
		mat_id = md.pcl_mat_model[id]->get_id();
	}

	void to_pcl(
		Model_T3D_ME_mt &md,
		size_t pcl_offset,
		size_t pcl_sorted_var_id,
		MatModel::MaterialModel &mm
		)
	{
		Model_T3D_ME_mt::PclSortedVarArray& psva
			= md.pcl_sorted_var_array[pcl_sorted_var_id];
		psva.pcl_index[pcl_offset] = id;
		Model_T3D_ME_mt::PclPos &pcl_pos = md.pcl_pos[id];
		pcl_pos.x = x;
		pcl_pos.y = y;
		pcl_pos.z = z;
		Model_T3D_ME_mt::PclBodyForce& pcl_bf = md.pcl_bf[id];
		pcl_bf.bfx = bfx;
		pcl_bf.bfy = bfy;
		pcl_bf.bfz = bfz;
		Model_T3D_ME_mt::PclTraction& pcl_t = md.pcl_t[id];
		pcl_t.tx = tx;
		pcl_t.ty = ty;
		pcl_t.tz = tz;
		md.pcl_m[id] = m;
		psva.pcl_density[pcl_offset] = density;
		Model_T3D_ME_mt::PclV &pcl_v = psva.pcl_v[pcl_offset];
		pcl_v.vx = vx;
		pcl_v.vy = vy;
		pcl_v.vz = vz;
		Model_T3D_ME_mt::PclStress& pcl_stress
			= psva.pcl_stress[pcl_offset];
		pcl_stress.s11 = s11;
		pcl_stress.s22 = s22;
		pcl_stress.s33 = s33;
		pcl_stress.s12 = s12;
		pcl_stress.s23 = s23;
		pcl_stress.s31 = s31;
		md.pcl_mat_model[id] = &mm;
	}
};

inline hid_t get_pcl_dt_id()
{
	hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(ParticleData));
	H5Tinsert(res, "id", HOFFSET(ParticleData, id), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "x", HOFFSET(ParticleData, x), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "y", HOFFSET(ParticleData, y), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "z", HOFFSET(ParticleData, z), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "bfx", HOFFSET(ParticleData, bfx), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "bfy", HOFFSET(ParticleData, bfy), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "bfz", HOFFSET(ParticleData, bfz), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "tx", HOFFSET(ParticleData, tx), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "ty", HOFFSET(ParticleData, ty), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "tz", HOFFSET(ParticleData, tz), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "m", HOFFSET(ParticleData, m), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "density", HOFFSET(ParticleData, density), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "vol", HOFFSET(ParticleData, vol), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "vx", HOFFSET(ParticleData, vx), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "vy", HOFFSET(ParticleData, vy), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "vz", HOFFSET(ParticleData, vz), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "s11", HOFFSET(ParticleData, s11), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "s22", HOFFSET(ParticleData, s22), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "s33", HOFFSET(ParticleData, s33), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "s12", HOFFSET(ParticleData, s12), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "s23", HOFFSET(ParticleData, s23), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "s31", HOFFSET(ParticleData, s31), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "mat_id", HOFFSET(ParticleData, mat_id), H5T_NATIVE_ULLONG);
	return res;
}

struct NodeData
{
	size_t id;
	double x, y, z;
};

inline hid_t get_node_dt_id()
{
	hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(NodeData));
	H5Tinsert(res, "id", HOFFSET(NodeData, id), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "x", HOFFSET(NodeData, x), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "y", HOFFSET(NodeData, y), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "z", HOFFSET(NodeData, z), H5T_NATIVE_DOUBLE);
	return res;
}

struct ElementData
{
	size_t id;
	size_t n1, n2, n3, n4;
	double vol;
	double a1, b1, c1, d1;
	double a2, b2, c2, d2;
	double a3, b3, c3, d3;
	double a4, b4, c4, d4;
};

inline hid_t get_element_dt_id()
{
	hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(ElementData));
	H5Tinsert(res, "id", HOFFSET(ElementData, id), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "n1", HOFFSET(ElementData, n1), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "n2", HOFFSET(ElementData, n2), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "n3", HOFFSET(ElementData, n3), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "n4", HOFFSET(ElementData, n4), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "vol", HOFFSET(ElementData, vol), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "a1", HOFFSET(ElementData, a1), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "b1", HOFFSET(ElementData, b1), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "c1", HOFFSET(ElementData, c1), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "d1", HOFFSET(ElementData, d1), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "a2", HOFFSET(ElementData, a2), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "b2", HOFFSET(ElementData, b2), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "c2", HOFFSET(ElementData, c2), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "d2", HOFFSET(ElementData, d2), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "a3", HOFFSET(ElementData, a3), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "b3", HOFFSET(ElementData, b3), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "c3", HOFFSET(ElementData, c3), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "d3", HOFFSET(ElementData, d3), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "a4", HOFFSET(ElementData, a4), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "b4", HOFFSET(ElementData, b4), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "c4", HOFFSET(ElementData, c4), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "d4", HOFFSET(ElementData, d4), H5T_NATIVE_DOUBLE);
	return res;
}

struct NodeVBCData { bool has_vx_bc, has_vy_bc, has_vz_bc; };

inline hid_t get_node_vbc_dt_id()
{
	hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(NodeVBCData));
	H5Tinsert(res, "has_vx_bc", HOFFSET(NodeVBCData, has_vx_bc), H5T_NATIVE_HBOOL);
	H5Tinsert(res, "has_vy_bc", HOFFSET(NodeVBCData, has_vy_bc), H5T_NATIVE_HBOOL);
	H5Tinsert(res, "has_vz_bc", HOFFSET(NodeVBCData, has_vz_bc), H5T_NATIVE_HBOOL);
	return res;
}

int output_background_mesh_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
int load_background_mesh_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);

int output_boundary_condition_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
int load_boundary_condition_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);

int output_ori_pcl_data_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
int output_pcl_data_to_hdf5_file(Model_T3D_ME_mt& md, Step_T3D_ME_mt &stp, ResultFile_hdf5& rf, hid_t grp_id);
int load_pcl_data_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);

int output_material_model_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
int load_material_model_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);

//int output_rigid_circle_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
//int load_rigid_circle_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);

int output_rigid_rect_to_hdf5_file(Model_T3D_ME_mt& md, Step_T3D_ME_mt &stp, ResultFile_hdf5& rf, hid_t grp_id);
int output_rigid_rect_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
int load_rigid_rect_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);

// output the whole model to ModelData
int output_model_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf);

// output the particle data and material models to hdf5 (used by time history)
int time_history_complete_output_to_hdf5_file(Model_T3D_ME_mt& md, Step_T3D_ME_mt& stp, ResultFile_hdf5& rf, hid_t frame_grp_id);

// load model data from hdf5 to model data
int load_me_mt_model_from_hdf5_file(Model_T3D_ME_mt &md, Step_T3D_ME_mt& step, const char *hdf5_name, const char *th_name,	size_t frame_id);
};

#endif