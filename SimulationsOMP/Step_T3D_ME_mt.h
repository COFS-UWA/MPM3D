#ifndef __Step_T3D_ME_mt_h__
#define __Step_T3D_ME_mt_h__

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
	int load_me_mt_model_from_hdf5_file(Model_T3D_ME_mt& md, const char* hdf5_name);
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
	typedef Model_T3D_ME_mt::NodeHasVBC NodeHasVBC;
	typedef Model_T3D_ME_mt::SortedPclVarArrays SortedPclVarArrays;
	
protected:
	double* pcl_m;
	Force *pcl_bf;
	Force *pcl_t;
	Position* pcl_pos;
	double* pcl_vol;
	MatModel::MaterialModel** pcl_mat_model;
	SortedPclVarArrays sorted_pcl_var_arrays[2];

	size_t* contact_substep_id; // ori_pcl_num
	Position* prev_contact_pos; // ori_pcl_num
	Force* prev_contact_tan_force; // ori_pcl_num

	size_t elem_num;
	size_t node_num;
	
	ElemNodeIndex* elem_node_id;
	DShapeFuncABC* elem_dN_abc;
	DShapeFuncD* elem_dN_d;
	double* elem_vol;

	double* elem_density;
	double* elem_pcl_m;
	StrainInc* elem_de;
	double* elem_m_de_vol;

	ElemNodeVM *elem_node_vm;
	Force *elem_node_force;

	Acceleration *node_a;
	Velocity *node_v;
	NodeHasVBC* node_has_vbc;
	double *node_am;
	double *node_de_vol;

	// thread-wise data
	union ThreadData
	{
		struct
		{
			size_t sorted_pcl_var_id;
			size_t sorted_pcl_in_elem_id;
			double max_pcl_vol;
			PclVar_T3D_ME_mt pcl_var_getter;
#ifdef _DEBUG
			// time profiling
			std::chrono::nanoseconds total;
			std::chrono::nanoseconds d0;
			std::chrono::nanoseconds d1, d10[8], d11[8], d12[8], d13[8];
			std::chrono::nanoseconds d20;
			std::chrono::nanoseconds d3, d30[8], d31[8], d32[8], d33[8];
			std::chrono::nanoseconds d4, d40, d49;
			std::chrono::nanoseconds d5, d59;
			std::chrono::nanoseconds d6, d69;
			std::chrono::nanoseconds d7;
#endif
		};
		char padding[Cache_Alignment * 20];
		ThreadData() : pcl_var_getter() {}
		~ThreadData() {}
	};
	ThreadData* thread_datas;
	
	size_t *pcl_in_elems[2]; // pcl_num
	size_t *prev_pcl_ids[2]; // pcl_num
	size_t* valid_elem_id; // elem_num
	size_t *node_has_elems[2]; // elem_num * 4
	size_t *node_elem_pairs[2]; // elem_num * 4
	// radix sort
	size_t* elem_count_bin;
	size_t* elem_sum_bin;

	RigidCylinder* prcy;
	RigidCone* prco;
	RigidCube* prcu;
	RigidObjectByT3DMesh* prmesh;

	ContactModel3D *pcf;

	size_t prev_valid_pcl_num, valid_pcl_num;
	size_t valid_elem_num;
	Force3D cf_tmp;
	
	CacheAlignedMem thread_mem;
	CacheAlignedMem cal_mem;

	int apply_rigid_cylinder(
		size_t p_id0, size_t p_id1,
		size_t *pcl_in_elem,
		SortedPclVarArrays &cur_spva,
		Force3D &rc_cf,
		size_t substp_id,
		ThreadData &thd) noexcept;

	int apply_rigid_cube(
		size_t p_id0, size_t p_id1,
		size_t* pcl_in_elem,
		SortedPclVarArrays& cur_spva,
		Force3D& rc_cf,
		size_t substp_id,
		ThreadData& thd) noexcept;

	int apply_t3d_rigid_object(
		size_t p_id0, size_t p_id1,
		size_t* pcl_in_elem,
		SortedPclVarArrays& cur_spva,
		Force3D& rc_cf,
		size_t substp_id,
		ThreadData& thd) noexcept;

#ifdef _DEBUG
	size_t prev_valid_pcl_num_tmp;
#endif

public:
	int init_calculation() override;
	friend int substep_func_omp_T3D_ME_mt(void* _self,
		size_t my_th_id, double dt, double cur_time, size_t substp_id);
	int finalize_calculation() override;

public:
	Step_T3D_ME_mt(const char* _name);
	~Step_T3D_ME_mt();

	inline size_t get_pcl_num() const noexcept { return prev_valid_pcl_num; }
	inline size_t get_sorted_pcl_var_id() const noexcept { return thread_datas[0].sorted_pcl_var_id; }
	inline size_t *get_pcl_in_elem() const noexcept { return pcl_in_elems[thread_datas[0].sorted_pcl_in_elem_id]; }

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