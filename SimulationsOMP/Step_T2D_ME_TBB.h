#ifndef __Step_T2D_ME_TBB_h__
#define __Step_T2D_ME_TBB_h__

#include "tbb/task_scheduler_init.h"

#include "Model_T2D_ME_mt.h"
#include "Step_T2D_ME_Task.h"
#include "Step_TBB.h"

class Step_T2D_ME_TBB;
namespace Model_T2D_ME_mt_hdf5_utilities
{
	struct ParticleData;
	int output_pcl_data_to_hdf5_file(
		Model_T2D_ME_mt& md, Step_T2D_ME_TBB& stp,
		ResultFile_hdf5& rf, hid_t grp_id);
	int time_history_complete_output_to_hdf5_file(
		Model_T2D_ME_mt& md, Step_T2D_ME_TBB& stp,
		ResultFile_hdf5& rf, hid_t frame_grp_id);
}

int cal_substep_func_T2D_ME_TBB(void* _self);

class Step_T2D_ME_TBB : public Step_TBB
{
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

	Step_T2D_ME_Task::CalData cal_data;
	Step_T2D_ME_Task::InitPcl init_pcl;
	Step_T2D_ME_Task::MapPclToBgMesh map_pcl_to_mesh;
	Step_T2D_ME_Task::UpdateAccelerationAndVelocity update_a_and_v;
	Step_T2D_ME_Task::CalElemDeAndMapToNode cal_elem_de;
	Step_T2D_ME_Task::CalNodeDe cal_node_de;
	Step_T2D_ME_Task::MapBgMeshToPcl map_mesh_to_pcl;
	
	tbb::task_scheduler_init sche_init;

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

	friend struct Model_T2D_ME_mt_hdf5_utilities::ParticleData;
	friend int Model_T2D_ME_mt_hdf5_utilities::output_pcl_data_to_hdf5_file(
		Model_T2D_ME_mt& md, Step_T2D_ME_TBB& stp, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_ME_mt_hdf5_utilities::time_history_complete_output_to_hdf5_file(
		Model_T2D_ME_mt& md, Step_T2D_ME_TBB& stp, ResultFile_hdf5& rf, hid_t frame_grp_id);
};

#endif