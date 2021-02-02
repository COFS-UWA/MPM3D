#ifndef __Step_T2D_ME_TBB_h__
#define __Step_T2D_ME_TBB_h__

#include "CacheAlignedMem.h"
#include "Model_T2D_ME_mt.h"
#include "SortTask.h"
#include "Step_T2D_ME_TBB_task.h"
#include "Step_TBB.h"

class Step_T2D_ME_TBB;
namespace Model_T2D_ME_mt_hdf5_utilities
{
	int load_me_mt_model_from_hdf5_file(Model_T2D_ME_mt& md, Step_T2D_ME_TBB& step, const char* hdf5_name, const char* th_name, size_t frame_id);
}

int cal_substep_func_T2D_ME_TBB(void* _self);

class Step_T2D_ME_TBB : public Step_TBB
{
	friend class Step_T2D_ME_TBB_task::CalData_T2D_ME_TBB;
	friend class Step_T2D_ME_TBB_task::MapPclToBgMeshTask;

protected:
	typedef Model_T2D_ME_mt::SortedPclVarArrays SortedPclVarArrays;
	typedef Model_T2D_ME_mt::Force Force;
	typedef Model_T2D_ME_mt::Acceleration Acceleration;
	typedef Model_T2D_ME_mt::Velocity Velocity;
	typedef Model_T2D_ME_mt::Displacement Displacement;
	typedef Model_T2D_ME_mt::Position Position;
	typedef Model_T2D_ME_mt::ShapeFunc ShapeFunc;
	typedef Model_T2D_ME_mt::ShapeFuncAB ShapeFuncAB;
	typedef Model_T2D_ME_mt::ShapeFuncC ShapeFuncC;
	typedef Model_T2D_ME_mt::Stress Stress;
	typedef Model_T2D_ME_mt::Strain Strain;
	typedef Model_T2D_ME_mt::StrainInc StrainInc;
	typedef Model_T2D_ME_mt::ElemNodeIndex ElemNodeIndex;
	typedef Model_T2D_ME_mt::ElemNodeVM ElemNodeVM;
	typedef Model_T2D_ME_mt::NodeHasVBC NodeHasVBC;
	typedef Step_T2D_ME_TBB_task::Range Range;

	size_t sorted_pcl_var_id;
	size_t prev_valid_pcl_num, valid_pcl_num;
#ifdef _DEBUG
	size_t prev_valid_pcl_num_tmp;
#endif
	size_t valid_elem_num;
	Force2D cf_tmp;
	
	size_t pcl_num_per_map_pcl_to_bg_mesh_task;
	size_t pcl_num_per_map_bg_mesh_to_pcl_task;

	SortUtils::SortMem pcl_sort_mem;
	SortUtils::SortMem node_sort_mem;
	
	// thread wise data
	Range **pcl_ranges;
	Range **node_ranges;
	Range **elem_ranges;
	size_t** valid_elem_arrays; // thread_num

	//int apply_rigid_rect(
	//	size_t p_id0, size_t p_id1,
	//	size_t* pcl_in_elem,
	//	SortedPclVarArrays& cur_spva,
	//	Force2D& rc_cf,
	//	size_t substp_id,
	//	ThreadData& thd) noexcept;

	int init_calculation() override;
	friend int cal_substep_func_T2D_ME_TBB(void* _self);
	int finalize_calculation() override;

public:
	Step_T2D_ME_TBB(const char* _name);
	~Step_T2D_ME_TBB();

	//inline size_t get_pcl_num() const noexcept { return prev_valid_pcl_num; }
	//inline size_t get_sorted_pcl_var_id() const noexcept { return thread_datas[0].sorted_pcl_var_id; }
	//inline size_t* get_pcl_in_elem() const noexcept { return pcl_in_elems[thread_datas[0].sorted_pcl_in_elem_id]; }

	friend int Model_T2D_ME_mt_hdf5_utilities::load_me_mt_model_from_hdf5_file(Model_T2D_ME_mt& md, Step_T2D_ME_TBB& step, const char* hdf5_name, const char* th_name, size_t frame_id);
};

#endif