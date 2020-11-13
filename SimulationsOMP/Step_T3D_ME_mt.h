#ifndef __Step_T3D_ME_mt_h__
#define __Step_T3D_ME_mt_h__

#include <fstream>
#include <iomanip>

#include "CacheAlignedMem.h"
#include "Step_OMP.h"
#include "Model_T3D_ME_mt.h"

class Model_T3D_ME_mt;
class Step_T3D_ME_mt;

class ResultFile_hdf5;
namespace Model_T3D_ME_mt_hdf5_utilities
{
	struct ParticleData;
	int output_background_mesh_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_background_mesh_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_boundary_condition_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_boundary_condition_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_pcl_data_to_hdf5_file(Model_T3D_ME_mt& md, Step_T3D_ME_mt& stp, ResultFile_hdf5& rf, hid_t grp_id);
	int load_pcl_data_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_material_model_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_material_model_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_me_mt_model_from_hdf5_file(Model_T3D_ME_mt& md, Step_T3D_ME_mt& step, const char* hdf5_name, const char* th_name, size_t frame_id);
}

int substep_func_omp_T3D_ME_mt(void* _self, size_t my_th_id,
	double dt, double cur_time, size_t substp_id);

class Step_T3D_ME_mt : public Step_OMP
{
public:
	typedef Model_T3D_ME_mt::ShapeFunc ShapeFunc;
	typedef Model_T3D_ME_mt::DShapeFuncABC DShapeFuncABC;
	typedef Model_T3D_ME_mt::DShapeFuncD DShapeFuncD;
	typedef Model_T3D_ME_mt::Force Force;
	typedef Model_T3D_ME_mt::Position Position;
	typedef Model_T3D_ME_mt::Displacement Displacement;
	typedef Model_T3D_ME_mt::Velocity Velocity;
	typedef Model_T3D_ME_mt::Acceleration Acceleration;
	typedef Model_T3D_ME_mt::Stress Stress;
	typedef Model_T3D_ME_mt::Strain Strain;
	typedef Model_T3D_ME_mt::StrainInc StrainInc;
	typedef Model_T3D_ME_mt::ElemNodeIndex ElemNodeIndex;
	typedef Model_T3D_ME_mt::ElemNodeVM ElemNodeVM;
	typedef Model_T3D_ME_mt::ElemNodeForce ElemNodeForce;
	typedef Model_T3D_ME_mt::NodeHasVBC NodeHasVBC;

	struct SortedPclVarArrays
	{
		size_t* pcl_index; // ori_pcl_num
		double* pcl_density; // ori_pcl_num
		Displacement* pcl_disp; // ori_pcl_num
		Velocity* pcl_v; // ori_pcl_num
		Stress* pcl_stress; // ori_pcl_num
		Strain* pcl_strain; // ori_pcl_num
		Strain* pcl_estrain; // ori_pcl_num
		Strain* pcl_pstrain; // ori_pcl_num
		size_t* pcl_in_elem; // ori_pcl_num
		size_t *contact_substep_id; // ori_pcl_num
		Position *prev_contact_pos; // ori_pcl_num
		Force *prev_contact_force; // ori_pcl_num
	};

protected:
	double* pcl_m;
	Force *pcl_bf;
	Force *pcl_t;
	Position* pcl_pos;
	double* pcl_vol;
	ShapeFunc* pcl_N;
	MatModel::MaterialModel** pcl_mat_model;

	SortedPclVarArrays sorted_pcl_var_arrays[2];
	size_t *pcl_in_elem_tmp;

	size_t elem_num;

	ElemNodeIndex* elem_node_id;
	size_t* elem_id_array;
	size_t* node_elem_id_array;
	size_t* node_elem_list;
	DShapeFuncABC* elem_dN_abc;
	DShapeFuncD* elem_dN_d;
	double* elem_vol;

	size_t* elem_substep_id;
	double* elem_density;
	double* elem_pcl_m;
	double* elem_pcl_vol;
	StrainInc* elem_de;
	double* elem_m_de_vol;

	ElemNodeVM *elem_node_vm;
	ElemNodeForce *elem_node_force;

	size_t *node_substep_id;
	Acceleration *node_a;
	Velocity *node_v;
	NodeHasVBC* node_has_vbc;
	double *node_am;
	double *node_de_vol;

	// task division
	union ThreadData
	{
		struct
		{
			size_t sorted_pcl_var_id;
		};
		char padding[Cache_Alignment];
	};
	ThreadData* thread_datas;
	size_t *node_range;
	size_t *node_elem_range;

	// radix sort
	size_t* elem_count_bin;
	size_t* elem_sum_bin;

	RigidCylinder *prc;
	double Kn_cont, Kt_cont;

	size_t pcl_num;
	double rr_fx_cont, rr_fy_cont, rr_fz_cont;
	double rr_mx_cont, rr_my_cont, rr_mz_cont;

	CacheAlignedMem task_range_mem;
	CacheAlignedMem radix_sort_var_mem;
	CacheAlignedMem elem_bin_mem;

	RigidCylinder rigid_cylinder;

	int apply_rigid_cylinder(
		size_t p_id0, size_t p_id1,
		size_t *pcl_in_elem,
		SortedPclVarArrays &cur_spva,
		ContactForce3D &rc_cf);

public:
	int init_calculation() override;
	friend int substep_func_omp_T3D_ME_mt(void* _self,
		size_t my_th_id, double dt, double cur_time, size_t substp_id);
	int finalize_calculation() override;

public:
	Step_T3D_ME_mt(const char* _name);
	~Step_T3D_ME_mt();

	inline size_t get_pcl_num() const noexcept { return pcl_num; }
	inline size_t get_sorted_pcl_var_id() const noexcept { return thread_datas[0].sorted_pcl_var_id; }
	inline double get_rr_fx_contact() const noexcept { return rr_fx_cont; }
	inline double get_rr_fy_contact() const noexcept { return rr_fy_cont; }
	inline double get_rr_fz_contact() const noexcept { return rr_fz_cont; }
	inline double get_rr_mx_contact() const noexcept { return rr_mx_cont; }
	inline double get_rr_my_contact() const noexcept { return rr_my_cont; }
	inline double get_rr_mz_contact() const noexcept { return rr_mz_cont; }

	friend struct Model_T3D_ME_mt_hdf5_utilities::ParticleData;
	friend int Model_T3D_ME_mt_hdf5_utilities::output_background_mesh_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::load_background_mesh_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::output_boundary_condition_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::load_boundary_condition_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::output_pcl_data_to_hdf5_file(Model_T3D_ME_mt& md, Step_T3D_ME_mt& stp, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::load_pcl_data_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::output_material_model_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::load_material_model_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::load_me_mt_model_from_hdf5_file(Model_T3D_ME_mt& md, Step_T3D_ME_mt& step, const char* hdf5_name, const char* th_name, size_t frame_id);
};

#endif