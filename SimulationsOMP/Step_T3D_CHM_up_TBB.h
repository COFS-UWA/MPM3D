#ifndef __Step_T3D_CHM_up_TBB_h__
#define __Step_T3D_CHM_up_TBB_h__

#include "PclSort.h"
#include "NodeElemSort.hpp"
#include "Model_T3D_CHM_up_mt.h"
#include "Step_T3D_CHM_up_TBB_Task.h"
#include "Step_TBB.h"

class Step_T3D_CHM_up_TBB;

namespace Model_T3D_CHM_up_mt_hdf5_utilities
{
	struct ParticleData;
	int output_pcl_data_to_hdf5_file(Model_T3D_CHM_up_mt& md, Step_T3D_CHM_up_TBB& stp, ResultFile_hdf5& rf, hid_t grp_id);
	int time_history_complete_output_to_hdf5_file(Model_T3D_CHM_up_mt& md, Step_T3D_CHM_up_TBB& stp, ResultFile_hdf5& rf, hid_t frame_grp_id);
	int load_model_from_hdf5_file(Model_T3D_CHM_up_mt& md, const char* hdf5_name);
	int load_model_from_hdf5_file(Model_T3D_CHM_up_mt& md, Step_T3D_CHM_up_TBB& step, const char* hdf5_name, const char* th_name, size_t frame_id);
}

int cal_substep_func_T3D_CHM_up_TBB(void *_self);

class Step_T3D_CHM_up_TBB : public Step_TBB
{
protected:
	typedef Model_T3D_CHM_up_mt::SortedPclVarArrays SortedPclVarArrays;
	typedef Model_T3D_CHM_up_mt::Force Force;
	typedef Model_T3D_CHM_up_mt::Acceleration Acceleration;
	typedef Model_T3D_CHM_up_mt::Velocity Velocity;
	typedef Model_T3D_CHM_up_mt::Displacement Displacement;
	typedef Model_T3D_CHM_up_mt::Position Position;
	typedef Model_T3D_CHM_up_mt::ShapeFunc ShapeFunc;
	typedef Model_T3D_CHM_up_mt::DShapeFuncABC ShapeFuncABC;
	typedef Model_T3D_CHM_up_mt::DShapeFuncD ShapeFuncD;
	typedef Model_T3D_CHM_up_mt::Stress Stress;
	typedef Model_T3D_CHM_up_mt::Strain Strain;
	typedef Model_T3D_CHM_up_mt::StrainInc StrainInc;
	typedef Model_T3D_CHM_up_mt::ElemNodeIndex ElemNodeIndex;
	typedef Model_T3D_CHM_up_mt::ElemNodeVM ElemNodeVM;
	typedef Model_T3D_CHM_up_mt::NodeHasVBC NodeHasVBC;
	// ParaUtil
	typedef Step_T3D_CHM_up_TBB_Task::InitPcl InitPcl;
	typedef Step_T3D_CHM_up_TBB_Task::MapPclToBgMesh MapPclToBgMesh;
	typedef Step_T3D_CHM_up_TBB_Task::FindSoilSurface FindSoilSurface;
	typedef Step_T3D_CHM_up_TBB_Task::UpdateAccelerationAndVelocity UpdateAccelerationAndVelocity;
	typedef Step_T3D_CHM_up_TBB_Task::CalElemDeAndMapToNode CalElemDeAndMapToNode;
	typedef Step_T3D_CHM_up_TBB_Task::CalNodeDe CalNodeDe;
	typedef Step_T3D_CHM_up_TBB_Task::MapBgMeshToPcl MapBgMeshToPcl;
	typedef Step_T3D_CHM_up_TBB_Task::ContactRigidBody ContactRigidBody;
	// reduce res
	typedef Step_T3D_CHM_up_TBB_Task::InitPclRes InitPclRes;
	typedef Step_T3D_CHM_up_TBB_Task::MapPclToBgMeshRes MapPclToBgMeshRes;
	typedef Step_T3D_CHM_up_TBB_Task::MapBgMeshToPclRes MapBgMeshToPclRes;
	//
	friend class InitPcl;
	friend class MapPclToBgMesh;
	friend class FindSoilSurface;
	friend class ContactRigidBody;
	friend class UpdateAccelerationAndVelocity;
	friend class CalElemDeAndMapToNode;
	friend class CalNodeDe;
	friend class MapBgMeshToPcl;

	tbb::task_scheduler_init sche_init;
	
	Model_T3D_CHM_up_mt* pmodel;

	Model_T3D_CHM_up_mt::SortedPclVarArrays spvas[2];

	ParaUtil::PclSort pcl_sort;
	ParaUtil::NodeElemSort<4> ne_sort;

	// for sorting pcls
	size_t* in_pcl_in_elems, * in_prev_pcl_ids;
	const size_t* pcl_in_elems, * prev_pcl_ids;
	// for sorting nodes and elements
	const size_t *elem_ids, *node_ids, *node_elem_offs;

	InitPcl init_pcl;
	MapPclToBgMesh map_pcl_to_mesh;
	FindSoilSurface find_soil_surface;
	UpdateAccelerationAndVelocity update_a_and_v;
	CalElemDeAndMapToNode cal_elem_de;
	CalNodeDe cal_node_de;
	MapBgMeshToPcl map_mesh_to_pcl;
	ContactRigidBody cont_rigid_body;

	// data changed during computation
	InitPclRes init_pcl_res;
	union
	{
		Force3D react_force;
		MapPclToBgMeshRes map_pcl_to_mesh_res;
	};
	union
	{
		size_t valid_pcl_num;
		MapBgMeshToPclRes map_mesh_to_pcl_res;
	};
	size_t prev_valid_pcl_num;
	size_t valid_elem_num;
#ifdef _DEBUG
	size_t prev_valid_pcl_num_tmp;
	size_t ori_pcl_num, elem_num, node_num;
#endif
	
	int init_calculation() override;
	friend int cal_substep_func_T3D_CHM_up_TBB(void* _self);
	int finalize_calculation() override;

public:
	Step_T3D_CHM_up_TBB(const char* _name);
	~Step_T3D_CHM_up_TBB();

	inline size_t get_pcl_num() const noexcept { return prev_valid_pcl_num; }
	inline size_t next_spva_id() const noexcept { return (substep_index + 1) & 1; }
	inline size_t prev_spva_id() const noexcept { return substep_index & 1; }

	friend struct Model_T3D_CHM_up_mt_hdf5_utilities::ParticleData;
	friend int Model_T3D_CHM_up_mt_hdf5_utilities::output_pcl_data_to_hdf5_file(
		Model_T3D_CHM_up_mt& md, Step_T3D_CHM_up_TBB& stp, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_CHM_up_mt_hdf5_utilities::time_history_complete_output_to_hdf5_file(
		Model_T3D_CHM_up_mt& md, Step_T3D_CHM_up_TBB& stp, ResultFile_hdf5& rf, hid_t frame_grp_id);
	friend int Model_T3D_CHM_up_mt_hdf5_utilities::load_model_from_hdf5_file(
		Model_T3D_CHM_up_mt& md, const char* hdf5_name);
	friend int Model_T3D_CHM_up_mt_hdf5_utilities::load_model_from_hdf5_file(
		Model_T3D_CHM_up_mt& md, Step_T3D_CHM_up_TBB& step,
		const char* hdf5_name, const char* th_name, size_t frame_id);
};

#endif