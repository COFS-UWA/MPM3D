#ifndef __Step_T2D_CHM_mt_h__
#define __Step_T2D_CHM_mt_h__

#include "CacheAlignedMem.h"
#include "Step_OMP.h"
#include "Model_T2D_CHM_mt.h"
#include "RigidObject/Force2D.h"

class Model_T2D_CHM_mt;
class Step_T2D_CHM_mt;
//namespace Model_T2D_CHM_mt_hdf5_utilities
//{
//	int load_me_mt_model_from_hdf5_file(Model_T2D_CHM_mt& md, Step_T2D_CHM_mt& step, const char* hdf5_name, const char* th_name, size_t frame_id);
//}

int substep_func_omp_T2D_CHM_mt(void* _self, size_t my_th_id,
	double dt, double cur_time, size_t substp_id);

class Step_T2D_CHM_mt : public Step_OMP
{
protected:
	typedef Model_T2D_CHM_mt::ShapeFunc ShapeFunc;
	typedef Model_T2D_CHM_mt::DShapeFuncAB DShapeFuncAB;
	typedef Model_T2D_CHM_mt::DShapeFuncC DShapeFuncC;
	typedef Model_T2D_CHM_mt::Force Force;
	typedef Model_T2D_CHM_mt::Position Position;
	typedef Model_T2D_CHM_mt::Displacement Displacement;
	typedef Model_T2D_CHM_mt::Velocity Velocity;
	typedef Model_T2D_CHM_mt::Acceleration Acceleration;
	typedef Model_T2D_CHM_mt::Stress Stress;
	typedef Model_T2D_CHM_mt::Strain Strain;
	typedef Model_T2D_CHM_mt::StrainInc StrainInc;
	typedef Model_T2D_CHM_mt::SortedPclVarArrays SortedPclVarArrays;
	typedef Model_T2D_CHM_mt::ElemNodeIndex ElemNodeIndex;
	typedef Model_T2D_CHM_mt::ElemNodeVM ElemNodeVM;
	typedef Model_T2D_CHM_mt::NodeHasVBC NodeHasVBC;

	size_t elem_num;
	size_t node_num;

	double* pcl_m_s; // ori_pcl_num
	double* pcl_density_s; // ori_pcl_num
	Force* pcl_bf; // ori_pcl_num
	Force* pcl_t; // ori_pcl_num
	Position* pcl_pos; // ori_pcl_num
	double* pcl_vol; // ori_pcl_num
	MatModel::MaterialModel** pcl_mat_model; // ori_pcl_num

	SortedPclVarArrays sorted_pcl_var_arrays[2];

	ElemNodeIndex *elem_node_id;
	DShapeFuncAB *elem_N_ab;
	DShapeFuncC *elem_N_c;
	double* elem_area;

	double* elem_density_f; // elem_num
	double* elem_pcl_n; // elem_num
	double* elem_pcl_m_s; // elem_num
	double* elem_pcl_m_f; // elem_num
	StrainInc* elem_de; // elem_num
	double* elem_p; // elem_num
	double* elem_m_de_vol_s; // elem_num
	double* elem_m_de_vol_f; // elem_num

	// element-node data
	ElemNodeVM* elem_node_vm_s; // elem_num * 3
	ElemNodeVM* elem_node_vm_f; // elem_num * 3
	Force* elem_node_force_s; // elem_num * 3
	Force* elem_node_force_f; // elem_num * 3

	Acceleration* node_a_s; // node_num
	Acceleration* node_a_f; // node_num
	Velocity* node_v_s; // node_num
	Velocity* node_v_f; // node_num
	NodeHasVBC* node_has_vbc_s;
	NodeHasVBC* node_has_vbc_f; // node_num
	double* node_am_s; // node_num
	double* node_am_f; // node_num
	double* node_de_vol_s; // node_num
	double* node_de_vol_f; // node_num

	union ThreadData
	{
		size_t sorted_pcl_var_id;
		size_t sorted_pcl_in_elem_id;
		char padding[Cache_Alignment];
	};
	ThreadData *thread_datas;

	// radix sort
	size_t* elem_count_bin;
	size_t* elem_sum_bin;
	
	size_t *prev_pcl_ids[2];
	size_t *pcl_in_elems[2];
	size_t* node_has_elems[2];
	size_t* node_elem_pairs[2];

	size_t pcl_num;
	size_t valid_pcl_num;
	size_t valid_elem_num;

	CacheAlignedMem elem_bin_mem;
	CacheAlignedMem task_datas_mem;
	CacheAlignedMem radix_sort_var_mem;

	double Kn_cont;
	Force2D cf_tmp;

	int apply_rigid_circle(
		size_t p_id0, size_t p_id1,
		size_t* pcl_in_elem,
		SortedPclVarArrays& cur_spva,
		Force2D& rc_cf,
		size_t substp_id,
		ThreadData& thd) noexcept;

public:
	int init_calculation() override;
	friend int substep_func_omp_T2D_CHM_mt(void* _self,
		size_t my_th_id, double dt, double cur_time, size_t substp_id);
	int finalize_calculation() override;

public:
	Step_T2D_CHM_mt(const char* _name);
	~Step_T2D_CHM_mt();

	inline size_t get_pcl_num() const noexcept { return pcl_num; }
};

#endif