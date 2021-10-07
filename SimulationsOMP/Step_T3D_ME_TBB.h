#ifndef __Step_T3D_ME_TBB_h__
#define __Step_T3D_ME_TBB_h__

#include "PclSort.h"
#include "NodeElemSort.hpp"
#include "Model_T3D_ME_mt.h"
#include "Step_T3D_ME_TBB_Task.h"
#include "Step_TBB.h"

namespace Model_T3D_ME_mt_hdf5_utilities
{
	struct ParticleData;
	int output_pcl_data_to_hdf5_file(
		Model_T3D_ME_mt& md, Step_T3D_ME_TBB& stp,
		ResultFile_hdf5& rf, hid_t grp_id);
	int time_history_complete_output_to_hdf5_file(
		Model_T3D_ME_mt& md, Step_T3D_ME_TBB& stp,
		ResultFile_hdf5& rf, hid_t frame_grp_id);
}

int cal_substep_func_T3D_ME_TBB(void *_self);

class Step_T3D_ME_TBB : public Step_TBB
{
protected:
	typedef Model_T3D_ME_mt::SortedPclVarArrays SortedPclVarArrays;
	typedef Model_T3D_ME_mt::Force Force;
	typedef Model_T3D_ME_mt::Acceleration Acceleration;
	typedef Model_T3D_ME_mt::Velocity Velocity;
	typedef Model_T3D_ME_mt::Displacement Displacement;
	typedef Model_T3D_ME_mt::Position Position;
	typedef Model_T3D_ME_mt::ShapeFunc ShapeFunc;
	typedef Model_T3D_ME_mt::DShapeFuncABC ShapeFuncABC;
	typedef Model_T3D_ME_mt::DShapeFuncD ShapeFuncD;
	typedef Model_T3D_ME_mt::Stress Stress;
	typedef Model_T3D_ME_mt::Strain Strain;
	typedef Model_T3D_ME_mt::StrainInc StrainInc;
	typedef Model_T3D_ME_mt::ElemNodeIndex ElemNodeIndex;
	typedef Model_T3D_ME_mt::ElemNodeVM ElemNodeVM;
	typedef Model_T3D_ME_mt::NodeHasVBC NodeHasVBC;
	typedef Step_T3D_ME_TBB_Task::PclRange PclRange;
	typedef Step_T3D_ME_TBB_Task::NodeElemRange NodeElemRange;
	typedef Step_T3D_ME_TBB_Task::InitPclRes InitPclRes;
	typedef Step_T3D_ME_TBB_Task::InitPcl InitPcl;
	typedef Step_T3D_ME_TBB_Task::ContactForceRes ContactForceRes;
	typedef Step_T3D_ME_TBB_Task::MapPclToBgMesh MapPclToBgMesh;
	typedef Step_T3D_ME_TBB_Task::ContactRigidBody ContactRigidBody;
	typedef Step_T3D_ME_TBB_Task::UpdateAccelerationAndVelocity UpdateAccelerationAndVelocity;
	typedef Step_T3D_ME_TBB_Task::CalElemDeAndMapToNode CalElemDeAndMapToNode;
	typedef Step_T3D_ME_TBB_Task::CalNodeDe CalNodeDe;
	typedef Step_T3D_ME_TBB_Task::NewValidPclNum NewValidPclNum;
	typedef Step_T3D_ME_TBB_Task::MapBgMeshToPcl MapBgMeshToPcl;
	friend class InitPclRes;
	friend class InitPcl;
	friend class ContactForceRes;
	friend class MapPclToBgMesh;
	friend class ContactRigidBody;
	friend class UpdateAccelerationAndVelocity;
	friend class CalElemDeAndMapToNode;
	friend class CalNodeDe;
	friend class NewValidPclNum;
	friend class MapBgMeshToPcl;

	Model_T3D_ME_mt* pmodel;

	// pcl data
	const double* pcl_m;
	const Model_T3D_ME_mt::Force* pcl_bf;
	const Model_T3D_ME_mt::Force* pcl_t;
	const Model_T3D_ME_mt::Position* pcl_pos;
	double* pcl_vol;
	MatModel::MaterialModel** pcl_mat_model;

	inline size_t prev_spva_id() const noexcept { return substep_index & 1; }
	inline size_t next_spva_id() const noexcept { return (substep_index + 1) & 1; }
	Model_T3D_ME_mt::SortedPclVarArrays spvas[2];

	// elem data
	const Model_T3D_ME_mt::ElemNodeIndex* elem_node_id;
	const Model_T3D_ME_mt::DShapeFuncABC* elem_dN_abc;
	const Model_T3D_ME_mt::DShapeFuncD* elem_dN_d;
	const double* elem_vol;

	double* elem_pcl_m;
	double* elem_density;
	Model_T3D_ME_mt::StrainInc* elem_de;
	double* elem_m_de_vol;

	// elem node data
	Model_T3D_ME_mt::ElemNodeVM* elem_node_vm;
	Model_T3D_ME_mt::Force* elem_node_force;

	// node data
	Model_T3D_ME_mt::Acceleration* node_a;
	Model_T3D_ME_mt::Velocity* node_v;
	Model_T3D_ME_mt::NodeHasVBC* node_has_vbc;
	double* node_am;
	double* node_de_vol;

#ifdef _DEBUG
	size_t ori_pcl_num, elem_num, node_num;
#endif

	tbb::task_scheduler_init sche_init;

	size_t* in_pcl_in_elems, * in_prev_pcl_ids;
	const size_t* pcl_in_elems, * prev_pcl_ids;
	const size_t* elem_ids, * node_ids, * node_elem_offs;

	ParaUtil::PclSort pcl_sort;
	ParaUtil::NodeElemSort<4> ne_sort;

	Util::DataMem range_mem;
	PclRange* pcl_ranges;
	NodeElemRange* node_elem_ranges;
	
	InitPclRes init_pcl_res;
	InitPcl init_pcl;
	ContactForceRes cont_force_res;
	MapPclToBgMesh map_pcl_to_mesh;
	ContactRigidBody cont_rigid_body;
	UpdateAccelerationAndVelocity update_a_and_v;
	CalElemDeAndMapToNode cal_elem_de;
	CalNodeDe cal_node_de;
	NewValidPclNum new_valid_pcl_num;
	MapBgMeshToPcl map_mesh_to_pcl;

	// data changed during computation
	size_t prev_valid_pcl_num;
	size_t valid_pcl_num;
	size_t valid_elem_num;
	
	int init_calculation() override;
	friend int cal_substep_func_T3D_ME_TBB(void* _self);
	int finalize_calculation() override;

public:
	Step_T3D_ME_TBB(const char* _name);
	~Step_T3D_ME_TBB();

	friend struct Model_T3D_ME_mt_hdf5_utilities::ParticleData;
	friend int Model_T3D_ME_mt_hdf5_utilities::output_pcl_data_to_hdf5_file(
		Model_T3D_ME_mt& md, Step_T3D_ME_TBB& stp, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::time_history_complete_output_to_hdf5_file(
		Model_T3D_ME_mt& md, Step_T3D_ME_TBB& stp, ResultFile_hdf5& rf, hid_t frame_grp_id);
};

#endif