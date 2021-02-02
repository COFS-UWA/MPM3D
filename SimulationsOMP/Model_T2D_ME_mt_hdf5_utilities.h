#ifndef __Model_T2D_ME_mt_hdf5_utilities_h__
#define __Model_T2D_ME_mt_hdf5_utilities_h__

#include "ResultFile_hdf5.h"
#include "MatModelContainer.h"
#include "Model_T2D_ME_mt.h"
#include "Step_T2D_ME_mt.h"

namespace Model_T2D_ME_mt_hdf5_utilities
{
struct ParticleData
{
	size_t id;
	double x, y;
	double bfx, bfy;
	double tx, ty;
	double m, density, vol;
	double vx, vy;
	double s11, s22, s12;
	double e11, e22, e12;
	double ee11, ee22, ee12;
	double pe11, pe22, pe12;
	size_t mat_id; // material model id
	size_t elem_id;

	void from_pcl(
		Model_T2D_ME_mt& md,
		size_t pcl_offset,
		size_t pcl_sorted_var_id)
	{
		Model_T2D_ME_mt::SortedPclVarArrays &spva
			= md.sorted_pcl_var_arrays[pcl_sorted_var_id];
		id = spva.pcl_index[pcl_offset];
		Model_T2D_ME_mt::Position& pcl_pos = md.pcl_pos[id];
		Model_T2D_ME_mt::Displacement& pcl_disp = spva.pcl_disp[pcl_offset];
		x = pcl_pos.x + pcl_disp.ux;
		y = pcl_pos.y + pcl_disp.uy;
		Model_T2D_ME_mt::Force& pcl_bf = md.pcl_bf[id];
		bfx = pcl_bf.fx;
		bfy = pcl_bf.fy;
		Model_T2D_ME_mt::Force & pcl_t = md.pcl_t[id];
		tx = pcl_t.fx;
		ty = pcl_t.fy;
		m = md.pcl_m[id];
		density = spva.pcl_density[pcl_offset];
		vol = m / density;
		Model_T2D_ME_mt::Velocity &pcl_v = spva.pcl_v[pcl_offset];
		vx = pcl_v.vx;
		vy = pcl_v.vy;
		Model_T2D_ME_mt::Stress& pcl_stress = spva.pcl_stress[pcl_offset];
		s11 = pcl_stress.s11;
		s22 = pcl_stress.s22;
		s12 = pcl_stress.s12;
		Model_T2D_ME_mt::Strain& pcl_e = spva.pcl_strain[pcl_offset];
		e11 = pcl_e.e11;
		e22 = pcl_e.e22;
		e12 = pcl_e.e12;
		Model_T2D_ME_mt::Strain& pcl_ee = spva.pcl_estrain[pcl_offset];
		ee11 = pcl_ee.e11;
		ee22 = pcl_ee.e22;
		ee12 = pcl_ee.e12;
		Model_T2D_ME_mt::Strain& pcl_pe = spva.pcl_pstrain[pcl_offset];
		pe11 = pcl_pe.e11;
		pe22 = pcl_pe.e22;
		pe12 = pcl_pe.e12;
		mat_id = md.pcl_mat_model[id]->get_id();
		elem_id = SIZE_MAX;
	}
	
	void from_pcl(
		Step_T2D_ME_mt& stp,
		size_t sorted_var_id,
		size_t pcl_offset)
	{
		typedef Step_T2D_ME_mt::SortedPclVarArrays SortedPclVarArrays;
		sorted_var_id ^= 1;
		SortedPclVarArrays& spva = stp.sorted_pcl_var_arrays[sorted_var_id];
		id = spva.pcl_index[pcl_offset];
		Model_T2D_ME_mt::Position &pcl_pos = stp.pcl_pos[id];
		Model_T2D_ME_mt::Displacement& pcl_disp = spva.pcl_disp[pcl_offset];
		x = pcl_pos.x + pcl_disp.ux;
		y = pcl_pos.y + pcl_disp.uy;
		Model_T2D_ME_mt::Force& pcl_bf = stp.pcl_bf[id];
		bfx = pcl_bf.fx;
		bfy = pcl_bf.fy;
		Model_T2D_ME_mt::Force &pcl_t = stp.pcl_t[id];
		tx = pcl_t.fx;
		ty = pcl_t.fy;
		m = stp.pcl_m[id];
		density = spva.pcl_density[pcl_offset];
		vol = m / density;
		Model_T2D_ME_mt::Velocity &pcl_v = spva.pcl_v[pcl_offset];
		vx = pcl_v.vx;
		vy = pcl_v.vy;
		Model_T2D_ME_mt::Stress& pcl_stress = spva.pcl_stress[pcl_offset];
		s11 = pcl_stress.s11;
		s22 = pcl_stress.s22;
		s12 = pcl_stress.s12;
		Model_T2D_ME_mt::Strain& pcl_e = spva.pcl_strain[pcl_offset];
		e11 = pcl_e.e11;
		e22 = pcl_e.e22;
		e12 = pcl_e.e12;
		Model_T2D_ME_mt::Strain& pcl_ee = spva.pcl_estrain[pcl_offset];
		ee11 = pcl_ee.e11;
		ee22 = pcl_ee.e22;
		ee12 = pcl_ee.e12;
		Model_T2D_ME_mt::Strain& pcl_pe = spva.pcl_pstrain[pcl_offset];
		pe11 = pcl_pe.e11;
		pe22 = pcl_pe.e22;
		pe12 = pcl_pe.e12;
		mat_id = stp.pcl_mat_model[id]->get_id();
		elem_id = stp.get_pcl_in_elem()[pcl_offset];
	}

	void to_pcl(
		Model_T2D_ME_mt &md,
		size_t pcl_offset,
		size_t pcl_sorted_var_id,
		MatModel::MaterialModel &mm)
	{
		Model_T2D_ME_mt::SortedPclVarArrays& spva
			= md.sorted_pcl_var_arrays[pcl_sorted_var_id];
		spva.pcl_index[pcl_offset] = id;
		Model_T2D_ME_mt::Position &pcl_pos = md.pcl_pos[id];
		pcl_pos.x = x;
		pcl_pos.y = y;
		Model_T2D_ME_mt::Force &pcl_bf = md.pcl_bf[id];
		pcl_bf.fx = bfx;
		pcl_bf.fy = bfy;
		Model_T2D_ME_mt::Force &pcl_t = md.pcl_t[id];
		pcl_t.fx = tx;
		pcl_t.fy = ty;
		md.pcl_m[id] = m;
		spva.pcl_density[pcl_offset] = density;
		Model_T2D_ME_mt::Velocity &pcl_v = spva.pcl_v[pcl_offset];
		pcl_v.vx = vx;
		pcl_v.vy = vy;
		Model_T2D_ME_mt::Stress& pcl_stress = spva.pcl_stress[pcl_offset];
		pcl_stress.s11 = s11;
		pcl_stress.s22 = s22;
		pcl_stress.s12 = s12;
		Model_T2D_ME_mt::Strain& pcl_strain = spva.pcl_strain[pcl_offset];
		pcl_strain.e11 = e11;
		pcl_strain.e22 = e22;
		pcl_strain.e12 = e12;
		Model_T2D_ME_mt::Strain& pcl_estrain = spva.pcl_estrain[pcl_offset];
		pcl_estrain.e11 = ee11;
		pcl_estrain.e22 = ee22;
		pcl_estrain.e12 = ee12;
		Model_T2D_ME_mt::Strain& pcl_pstrain = spva.pcl_pstrain[pcl_offset];
		pcl_pstrain.e11 = pe11;
		pcl_pstrain.e22 = pe22;
		pcl_pstrain.e12 = pe12;
		md.pcl_mat_model[id] = &mm;
	}
};

inline hid_t get_pcl_dt_id()
{
	hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(ParticleData));
	H5Tinsert(res, "id", HOFFSET(ParticleData, id), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "x", HOFFSET(ParticleData, x), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "y", HOFFSET(ParticleData, y), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "bfx", HOFFSET(ParticleData, bfx), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "bfy", HOFFSET(ParticleData, bfy), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "tx", HOFFSET(ParticleData, tx), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "ty", HOFFSET(ParticleData, ty), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "m", HOFFSET(ParticleData, m), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "density", HOFFSET(ParticleData, density), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "vol", HOFFSET(ParticleData, vol), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "vx", HOFFSET(ParticleData, vx), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "vy", HOFFSET(ParticleData, vy), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "s11", HOFFSET(ParticleData, s11), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "s22", HOFFSET(ParticleData, s22), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "s12", HOFFSET(ParticleData, s12), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "e11", HOFFSET(ParticleData, e11), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "e22", HOFFSET(ParticleData, e22), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "e12", HOFFSET(ParticleData, e12), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "ee11", HOFFSET(ParticleData, ee11), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "ee22", HOFFSET(ParticleData, ee22), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "ee12", HOFFSET(ParticleData, ee12), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "pe11", HOFFSET(ParticleData, pe11), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "pe22", HOFFSET(ParticleData, pe22), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "pe12", HOFFSET(ParticleData, pe12), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "mat_id", HOFFSET(ParticleData, mat_id), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "elem_id", HOFFSET(ParticleData, elem_id), H5T_NATIVE_ULLONG);
	return res;
}

struct NodeData
{
	size_t id;
	double x, y;
};

inline hid_t get_node_dt_id()
{
	hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(NodeData));
	H5Tinsert(res, "id", HOFFSET(NodeData, id), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "x", HOFFSET(NodeData, x), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "y", HOFFSET(NodeData, y), H5T_NATIVE_DOUBLE);
	return res;
}

struct ElementData
{
	size_t id;
	size_t n1, n2, n3;
	double area;
	double a1, b1, c1;
	double a2, b2, c2;
	double a3, b3, c3;
};

inline hid_t get_element_dt_id()
{
	hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(ElementData));
	H5Tinsert(res, "id", HOFFSET(ElementData, id), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "n1", HOFFSET(ElementData, n1), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "n2", HOFFSET(ElementData, n2), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "n3", HOFFSET(ElementData, n3), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "area", HOFFSET(ElementData, area), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "a1", HOFFSET(ElementData, a1), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "b1", HOFFSET(ElementData, b1), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "c1", HOFFSET(ElementData, c1), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "a2", HOFFSET(ElementData, a2), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "b2", HOFFSET(ElementData, b2), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "c2", HOFFSET(ElementData, c2), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "a3", HOFFSET(ElementData, a3), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "b3", HOFFSET(ElementData, b3), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "c3", HOFFSET(ElementData, c3), H5T_NATIVE_DOUBLE);
	return res;
}

struct NodeVBCData { bool has_vx_bc, has_vy_bc; };

inline hid_t get_node_vbc_dt_id()
{
	hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(NodeVBCData));
	H5Tinsert(res, "has_vx_bc", HOFFSET(NodeVBCData, has_vx_bc), H5T_NATIVE_HBOOL);
	H5Tinsert(res, "has_vy_bc", HOFFSET(NodeVBCData, has_vy_bc), H5T_NATIVE_HBOOL);
	return res;
}

int output_background_mesh_to_hdf5_file(Model_T2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
int load_background_mesh_from_hdf5_file(Model_T2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);

int output_boundary_condition_to_hdf5_file(Model_T2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
int load_boundary_condition_from_hdf5_file(Model_T2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);

int output_ori_pcl_data_to_hdf5_file(Model_T2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
int output_pcl_data_to_hdf5_file(Model_T2D_ME_mt& md, Step_T2D_ME_mt &stp, ResultFile_hdf5& rf, hid_t grp_id);
int load_pcl_data_from_hdf5_file(Model_T2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);

int output_material_model_to_hdf5_file(Model_T2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
int load_material_model_from_hdf5_file(Model_T2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);

//int output_rigid_circle_to_hdf5_file(Model_T2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
//int load_rigid_circle_from_hdf5_file(Model_T2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);

int output_rigid_rect_to_hdf5_file(Model_T2D_ME_mt& md, Step_T2D_ME_mt &stp, ResultFile_hdf5& rf, hid_t grp_id);
int output_rigid_rect_to_hdf5_file(Model_T2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
int load_rigid_rect_from_hdf5_file(Model_T2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);

// output the whole model to ModelData
int output_model_to_hdf5_file(Model_T2D_ME_mt& md, ResultFile_hdf5& rf);

// output the particle data and material models to hdf5 (used by time history)
int time_history_complete_output_to_hdf5_file(Model_T2D_ME_mt& md, Step_T2D_ME_mt& stp, ResultFile_hdf5& rf, hid_t frame_grp_id);

// load model data from hdf5 to model data
int load_me_mt_model_from_hdf5_file(Model_T2D_ME_mt &md, Step_T2D_ME_mt& step, const char *hdf5_name, const char *th_name,	size_t frame_id);
};

#endif