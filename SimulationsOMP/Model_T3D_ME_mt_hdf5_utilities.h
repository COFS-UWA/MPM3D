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
	double e11, e22, e33, e12, e23, e31;
	double ee11, ee22, ee33, ee12, ee23, ee31;
	double pe11, pe22, pe33, pe12, pe23, pe31;
	size_t mat_id; // material model id
	size_t elem_id;

	void from_pcl(
		Model_T3D_ME_mt& md,
		size_t sorted_var_id,
		size_t pcl_offset
		)
	{
		typedef Model_T3D_ME_mt::SortedPclVarArrays SortedPclVarArrays;
		SortedPclVarArrays& spva = md.sorted_pcl_var_arrays[sorted_var_id];
		id = spva.pcl_index[pcl_offset];
		Model_T3D_ME_mt::Position &p_p = md.pcl_pos[id];
		Model_T3D_ME_mt::Displacement &p_d = spva.pcl_disp[pcl_offset];
		x = p_p.x + p_d.ux;
		y = p_p.y + p_d.uy;
		z = p_p.z + p_d.uz;
		Model_T3D_ME_mt::Force &p_bf = md.pcl_bf[id];
		bfx = p_bf.fx;
		bfy = p_bf.fy;
		bfz = p_bf.fz;
		Model_T3D_ME_mt::Force &p_t = md.pcl_t[id];
		tx = p_t.fx;
		ty = p_t.fy;
		tz = p_t.fz;
		m = md.pcl_m[id];
		density = spva.pcl_density[pcl_offset];
		vol = m / density;
		Model_T3D_ME_mt::Velocity &p_v = spva.pcl_v[pcl_offset];
		vx = p_v.vx;
		vy = p_v.vy;
		vz = p_v.vz;
		Model_T3D_ME_mt::Stress &p_s = spva.pcl_stress[pcl_offset];
		s11 = p_s.s11;
		s22 = p_s.s22;
		s33 = p_s.s33;
		s12 = p_s.s12;
		s23 = p_s.s23;
		s31 = p_s.s31;
		Model_T3D_ME_mt::Strain &p_e = spva.pcl_strain[pcl_offset];
		e11 = p_e.e11;
		e22 = p_e.e22;
		e33 = p_e.e33;
		e12 = p_e.e12;
		e23 = p_e.e23;
		e31 = p_e.e31;
		Model_T3D_ME_mt::Strain& p_ee = spva.pcl_estrain[pcl_offset];
		ee11 = p_ee.e11;
		ee22 = p_ee.e22;
		ee33 = p_ee.e33;
		ee12 = p_ee.e12;
		ee23 = p_ee.e23;
		ee31 = p_ee.e31;
		Model_T3D_ME_mt::Strain& p_pe = spva.pcl_pstrain[pcl_offset];
		pe11 = p_pe.e11;
		pe22 = p_pe.e22;
		pe33 = p_pe.e33;
		pe12 = p_pe.e12;
		pe23 = p_pe.e23;
		pe31 = p_pe.e31;
		mat_id = md.pcl_mat_model[id]->get_id();
		elem_id = SIZE_MAX;
	}
	
	void from_pcl(
		Step_T3D_ME_mt &stp,
		size_t sorted_var_id,
		size_t pcl_offset
		)
	{
		typedef Step_T3D_ME_mt::SortedPclVarArrays SortedPclVarArrays;
		sorted_var_id ^= 1;
		SortedPclVarArrays &spva = stp.sorted_pcl_var_arrays[sorted_var_id];
		id = spva.pcl_index[pcl_offset];
		Step_T3D_ME_mt::Position &p_p = stp.pcl_pos[id];
		Step_T3D_ME_mt::Displacement& p_d = spva.pcl_disp[pcl_offset];
		x = p_p.x + p_d.ux;
		y = p_p.y + p_d.uy;
		z = p_p.z + p_d.uz;
		Step_T3D_ME_mt::Force& p_bf = stp.pcl_bf[id];
		bfx = p_bf.fx;
		bfy = p_bf.fy;
		bfz = p_bf.fz;
		Step_T3D_ME_mt::Force &p_t = stp.pcl_t[id];
		tx = p_t.fx;
		ty = p_t.fy;
		tz = p_t.fz;
		m = stp.pcl_m[id];
		density = spva.pcl_density[pcl_offset];
		vol = m / density;
		Step_T3D_ME_mt::Velocity &p_v = spva.pcl_v[pcl_offset];
		vx = p_v.vx;
		vy = p_v.vy;
		vz = p_v.vz;
		Step_T3D_ME_mt::Stress &p_s = spva.pcl_stress[pcl_offset];
		s11 = p_s.s11;
		s22 = p_s.s22;
		s33 = p_s.s33;
		s12 = p_s.s12;
		s23 = p_s.s23;
		s31 = p_s.s31;
		Step_T3D_ME_mt::Strain& p_e = spva.pcl_strain[pcl_offset];
		e11 = p_e.e11;
		e22 = p_e.e22;
		e33 = p_e.e33;
		e12 = p_e.e12;
		e23 = p_e.e23;
		e31 = p_e.e31;
		Step_T3D_ME_mt::Strain& p_ee = spva.pcl_estrain[pcl_offset];
		ee11 = p_ee.e11;
		ee22 = p_ee.e22;
		ee33 = p_ee.e33;
		ee12 = p_ee.e12;
		ee23 = p_ee.e23;
		ee31 = p_ee.e31;
		Step_T3D_ME_mt::Strain& p_pe = spva.pcl_pstrain[pcl_offset];
		pe11 = p_pe.e11;
		pe22 = p_pe.e22;
		pe33 = p_pe.e33;
		pe12 = p_pe.e12;
		pe23 = p_pe.e23;
		pe31 = p_pe.e31;
		mat_id = stp.pcl_mat_model[id]->get_id();
		elem_id = stp.get_pcl_in_elem()[pcl_offset];
	}

	void to_pcl(
		Model_T3D_ME_mt &md,
		size_t sorted_pcl_var_id,
		size_t pcl_offset,
		MatModel::MaterialModel &mm
		)
	{
		typedef Model_T3D_ME_mt::SortedPclVarArrays SortedPclVarArrays;
		SortedPclVarArrays& spva = md.sorted_pcl_var_arrays[sorted_pcl_var_id];
		spva.pcl_index[pcl_offset] = id;
		Model_T3D_ME_mt::Position &p_p = md.pcl_pos[id];
		p_p.x = x;
		p_p.y = y;
		p_p.z = z;
		Model_T3D_ME_mt::Force& p_bf = md.pcl_bf[id];
		p_bf.fx = bfx;
		p_bf.fy = bfy;
		p_bf.fz = bfz;
		Model_T3D_ME_mt::Force &p_t = md.pcl_t[id];
		p_t.fx = tx;
		p_t.fy = ty;
		p_t.fz = tz;
		md.pcl_m[id] = m;
		Model_T3D_ME_mt::Displacement& p_d = spva.pcl_disp[pcl_offset];
		p_d.ux = 0.0;
		p_d.uy = 0.0;
		p_d.uz = 0.0;
		spva.pcl_density[pcl_offset] = density;
		Model_T3D_ME_mt::Velocity &p_v = spva.pcl_v[pcl_offset];
		p_v.vx = vx;
		p_v.vy = vy;
		p_v.vz = vz;
		Model_T3D_ME_mt::Stress& p_s = spva.pcl_stress[pcl_offset];
		p_s.s11 = s11;
		p_s.s22 = s22;
		p_s.s33 = s33;
		p_s.s12 = s12;
		p_s.s23 = s23;
		p_s.s31 = s31;
		Model_T3D_ME_mt::Strain& p_e = spva.pcl_strain[pcl_offset];
		p_e.e11 = e11;
		p_e.e22 = e22;
		p_e.e33 = e33;
		p_e.e12 = e12;
		p_e.e23 = e23;
		p_e.e31 = e31;
		Model_T3D_ME_mt::Strain& p_ee = spva.pcl_estrain[pcl_offset];
		p_ee.e11 = ee11;
		p_ee.e22 = ee22;
		p_ee.e33 = ee33;
		p_ee.e12 = ee12;
		p_ee.e23 = ee23;
		p_ee.e31 = ee31;
		Model_T3D_ME_mt::Strain& p_pe = spva.pcl_pstrain[pcl_offset];
		p_pe.e11 = pe11;
		p_pe.e22 = pe22;
		p_pe.e33 = pe33;
		p_pe.e12 = pe12;
		p_pe.e23 = pe23;
		p_pe.e31 = pe31;
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
	H5Tinsert(res, "e11", HOFFSET(ParticleData, e11), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "e22", HOFFSET(ParticleData, e22), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "e33", HOFFSET(ParticleData, e33), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "e12", HOFFSET(ParticleData, e12), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "e23", HOFFSET(ParticleData, e23), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "e31", HOFFSET(ParticleData, e31), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "ee11", HOFFSET(ParticleData, ee11), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "ee22", HOFFSET(ParticleData, ee22), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "ee33", HOFFSET(ParticleData, ee33), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "ee12", HOFFSET(ParticleData, ee12), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "ee23", HOFFSET(ParticleData, ee23), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "ee31", HOFFSET(ParticleData, ee31), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "pe11", HOFFSET(ParticleData, pe11), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "pe22", HOFFSET(ParticleData, pe22), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "pe33", HOFFSET(ParticleData, pe33), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "pe12", HOFFSET(ParticleData, pe12), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "pe23", HOFFSET(ParticleData, pe23), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "pe31", HOFFSET(ParticleData, pe31), H5T_NATIVE_DOUBLE);
	H5Tinsert(res, "mat_id", HOFFSET(ParticleData, mat_id), H5T_NATIVE_ULLONG);
	H5Tinsert(res, "elem_id", HOFFSET(ParticleData, elem_id), H5T_NATIVE_ULLONG);
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

int output_search_mesh_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
int load_search_mesh_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);

int output_ori_pcl_data_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
int output_pcl_data_to_hdf5_file(Model_T3D_ME_mt& md, Step_T3D_ME_mt &stp, ResultFile_hdf5& rf, hid_t grp_id);
int load_pcl_data_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);

int output_material_model_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
int load_material_model_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);

int output_rigid_cylinder_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
int load_rigid_cylinder_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);

int output_rigid_cone_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
int load_rigid_cone_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);

int output_rigid_cube_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
int load_rigid_cube_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);

int output_t3d_rigid_mesh_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
int output_t3d_rigid_mesh_state_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
int load_t3d_rigid_mesh_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);

// output the whole model to ModelData
int output_model_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf);

// output the particle data and material models to hdf5 (used by time history)
int time_history_complete_output_to_hdf5_file(Model_T3D_ME_mt& md, Step_T3D_ME_mt& stp, ResultFile_hdf5& rf, hid_t frame_grp_id);

// load model data from hdf5 to model data
int load_me_mt_model_from_hdf5_file(Model_T3D_ME_mt& md, const char* hdf5_name);
int load_me_mt_model_from_hdf5_file(Model_T3D_ME_mt &md, Step_T3D_ME_mt& step, const char *hdf5_name, const char *th_name,	size_t frame_id);

};

#endif