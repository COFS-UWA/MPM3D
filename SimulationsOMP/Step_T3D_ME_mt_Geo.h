#ifndef __Step_T3D_ME_mt_Geo_h__
#define __Step_T3D_ME_mt_Geo_h__

#include "CacheAlignedMem.h"
#include "Step_OMP.h"
#include "Model_T3D_ME_mt.h"

class Model_T3D_ME_mt;
class Step_T3D_ME_mt_Geo;

class ResultFile_hdf5;
namespace Model_T3D_ME_mt_hdf5_utilities
{
	struct ParticleData;
	int output_background_mesh_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_background_mesh_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_boundary_condition_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_boundary_condition_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_pcl_data_to_hdf5_file(Model_T3D_ME_mt& md, Step_T3D_ME_mt_Geo& stp, ResultFile_hdf5& rf, hid_t grp_id);
	int load_pcl_data_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_material_model_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_material_model_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_me_mt_model_from_hdf5_file(Model_T3D_ME_mt& md, const char* hdf5_name);
	int load_me_mt_model_from_hdf5_file(Model_T3D_ME_mt& md, Step_T3D_ME_mt_Geo& step, const char* hdf5_name, const char* th_name, size_t frame_id);
}

int substep_func_omp_T3D_ME_mt_Geo(void *_self, size_t my_th_id,
	double dt, double cur_time, size_t substp_id);

class Step_T3D_ME_mt_Geo : public Step_OMP
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
	typedef Model_T3D_ME_mt::NodeVBCVec NodeVBCVec;
	typedef Model_T3D_ME_mt::SortedPclVarArrays SortedPclVarArrays;
	
protected:
	double* pcl_m;
	Force *pcl_bf;
	Force *pcl_t;
	Position* pcl_pos;
	double* pcl_vol;
	MatModel::MaterialModel** pcl_mat_model;
	SortedPclVarArrays sorted_pcl_var_arrays[2];

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
	NodeVBCVec* node_vbc_vec;
	double *node_am;
	double *node_de_vol;

	// thread-wise data
	union ThreadData
	{
		struct
		{
			size_t sorted_pcl_var_id;
			size_t sorted_pcl_in_elem_id;
		};
		char padding[Cache_Alignment];
	};
	ThreadData* thread_datas;

	size_t *pcl_in_elems[2];
	size_t *prev_pcl_ids[2];
	size_t* valid_elem_id;
	size_t *node_has_elems[2];
	size_t *node_elem_pairs[2];
	// radix sort
	size_t* elem_count_bin;
	size_t* elem_sum_bin;

	size_t prev_valid_pcl_num, valid_pcl_num;
	size_t valid_elem_num;
	Force3D cf_tmp;
	
	CacheAlignedMem cal_mem;
	CacheAlignedMem thread_mem;

#ifdef _DEBUG
	size_t prev_valid_pcl_num_tmp;
#endif

	// 0 continue
	// 1 reset velocity
	// 2 exit
	unsigned char cal_status;
	bool init_f_ub_is_init;
	bool e_kin_max_is_init;
	double e_kin_prev;
	double f_ub_ratio_limit;
	double e_kin_ratio_limit;

	double init_f_ub;
	double f_ub;
	double f_ub_ratio;

	double e_kin_max;
	double e_kin;
	double e_kin_ratio;

public:
	int init_calculation() override;
	friend int substep_func_omp_T3D_ME_mt_Geo(void* _self,
		size_t my_th_id, double dt, double cur_time, size_t substp_id);
	int finalize_calculation() override;

public:
	Step_T3D_ME_mt_Geo(const char* _name);
	~Step_T3D_ME_mt_Geo();

	inline size_t get_pcl_num() const noexcept { return prev_valid_pcl_num; }
	inline size_t get_sorted_pcl_var_id() const noexcept { return thread_datas[0].sorted_pcl_var_id; }
	inline size_t *get_pcl_in_elem() const noexcept { return pcl_in_elems[thread_datas[0].sorted_pcl_in_elem_id]; }
	inline double get_f_ub() const noexcept { return f_ub; }
	inline double get_e_kin() const noexcept { return e_kin; }
	inline double get_f_ub_ratio() const noexcept { return f_ub_ratio; }
	inline double get_e_kin_ratio() const noexcept { return e_kin_ratio; }
	
	inline void set_f_ub_ratio_limit(double _lim) noexcept { f_ub_ratio_limit = _lim; }
	inline void set_e_kin_ratio_limit(double _lim) noexcept { e_kin_ratio_limit = _lim; }

	friend struct Model_T3D_ME_mt_hdf5_utilities::ParticleData;
	friend int Model_T3D_ME_mt_hdf5_utilities::output_background_mesh_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::load_background_mesh_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::output_boundary_condition_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::load_boundary_condition_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::output_pcl_data_to_hdf5_file(Model_T3D_ME_mt& md, Step_T3D_ME_mt_Geo& stp, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::load_pcl_data_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::output_material_model_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::load_material_model_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::load_me_mt_model_from_hdf5_file(Model_T3D_ME_mt& md, Step_T3D_ME_mt_Geo& step, const char* hdf5_name, const char* th_name, size_t frame_id);
};

#endif