#ifndef __Step_T3D_CHM_TBB_h__
#define __Step_T3D_CHM_TBB_h__

#include "PclSort.h"
#include "NodeElemSort.hpp"
#include "Model_T3D_CHM_mt.h"
#include "Step_T3D_CHM_TBB_Task.h"
#include "Step_TBB.h"

class Step_T3D_CHM_TBB;

namespace Model_T3D_CHM_mt_hdf5_utilities
{
	struct ParticleData;
	int output_pcl_data_to_hdf5_file(
		Model_T3D_CHM_mt& md, Step_T3D_CHM_TBB& stp,
		ResultFile_hdf5& rf, hid_t grp_id);
	int time_history_complete_output_to_hdf5_file(
		Model_T3D_CHM_mt& md, Step_T3D_CHM_TBB& stp,
		ResultFile_hdf5& rf, hid_t frame_grp_id);
	int load_model_from_hdf5_file(Model_T3D_CHM_mt& md, const char* hdf5_name);
	int load_model_from_hdf5_file(Model_T3D_CHM_mt& md, Step_T3D_CHM_TBB& step,
		const char* hdf5_name, const char* th_name, size_t frame_id);
}

int substep_func_T3D_CHM_TBB(void *_self);

class Step_T3D_CHM_TBB : public Step_TBB
{
protected:
	typedef Model_T3D_CHM_mt::SortedPclVarArrays SortedPclVarArrays;
	typedef Model_T3D_CHM_mt::Force Force;
	typedef Model_T3D_CHM_mt::Acceleration Acceleration;
	typedef Model_T3D_CHM_mt::Velocity Velocity;
	typedef Model_T3D_CHM_mt::Displacement Displacement;
	typedef Model_T3D_CHM_mt::Position Position;
	typedef Model_T3D_CHM_mt::ShapeFunc ShapeFunc;
	typedef Model_T3D_CHM_mt::DShapeFuncABC DShapeFuncABC;
	typedef Model_T3D_CHM_mt::DShapeFuncD DShapeFuncD;
	typedef Model_T3D_CHM_mt::Stress Stress;
	typedef Model_T3D_CHM_mt::Strain Strain;
	typedef Model_T3D_CHM_mt::StrainInc StrainInc;
	typedef Model_T3D_CHM_mt::ElemNodeIndex ElemNodeIndex;
	typedef Model_T3D_CHM_mt::ElemNodeVM ElemNodeVM;
	typedef Model_T3D_CHM_mt::NodeHasVBC NodeHasVBC;
	typedef Model_T3D_CHM_mt::NodeVBCVec NodeVBCVec;
	typedef Step_T3D_CHM_TBB_Task::InitPcl InitPcl;
	typedef Step_T3D_CHM_TBB_Task::MapPclToBgMesh MapPclToBgMesh;
	typedef Step_T3D_CHM_TBB_Task::ContactRigidBody ContactRigidBody;
	typedef Step_T3D_CHM_TBB_Task::UpdateAccelerationAndVelocity UpdateAccelerationAndVelocity;
	typedef Step_T3D_CHM_TBB_Task::CalElemDeAndMapToNode CalElemDeAndMapToNode;
	typedef Step_T3D_CHM_TBB_Task::CalNodeDe CalNodeDe;
	typedef Step_T3D_CHM_TBB_Task::MapBgMeshToPcl MapBgMeshToPcl;
	typedef Step_T3D_CHM_TBB_Task::InitPclRes InitPclRes;
	typedef Step_T3D_CHM_TBB_Task::MapPclToBgMeshRes MapPclToBgMeshRes;
	typedef Step_T3D_CHM_TBB_Task::MapBgMeshToPclRes MapBgMeshToPclRes;
	typedef Step_T3D_CHM_TBB_Task::InitPclTbb InitPclTbb;
	typedef Step_T3D_CHM_TBB_Task::MapPclToBgMeshTbb MapPclToBgMeshTbb;
	typedef Step_T3D_CHM_TBB_Task::MapBgMeshToPclTbb MapBgMeshToPclTbb;
	friend class InitPcl;
	friend class MapPclToBgMesh;
	friend class ContactRigidBody;
	friend class UpdateAccelerationAndVelocity;
	friend class CalElemDeAndMapToNode;
	friend class CalNodeDe;
	friend class MapBgMeshToPcl;

	tbb::task_scheduler_init sche_init;
	
	Model_T3D_CHM_mt* pmodel;

	// pcl data
	const double* pcl_m_s;
	const double* pcl_density_s;
	const double* pcl_vol_s;
	const Force* pcl_bf_s;
	const Force* pcl_bf_f;
	const Force* pcl_t;
	const Position* pcl_pos;
	double* pcl_vol;
	MatModel::MaterialModel** pcl_mat_model;

	SortedPclVarArrays spvas[2];
	inline size_t next_spva_id() const noexcept { return (substep_index + 1) & 1; }
	inline size_t prev_spva_id() const noexcept { return substep_index & 1; }

	// elem data
	const ElemNodeIndex* elem_node_id;
	const DShapeFuncABC* elem_N_abc;
	const DShapeFuncD* elem_N_d;
	const double* elem_vol;

	// element calculation data
	double* elem_density_f;
	double* elem_pcl_n;
	double* elem_pcl_m_s;
	double* elem_pcl_m_f;
	StrainInc* elem_de;
	double* elem_p;
	double* elem_n2_miu_div_k_vol;
	Force* elem_seep_force;
	double* elem_m_de_vol_s;
	double* elem_m_de_vol_f;

	// element-node data
	ElemNodeVM* elem_node_vm_s;
	ElemNodeVM* elem_node_vm_f;
	Force* elem_node_force_s;
	Force* elem_node_force_f;

	// node data
	Position* node_pos;
	Acceleration* node_a_s;
	Acceleration* node_a_f;
	Velocity* node_v_s;
	Velocity* node_v_f;
	const NodeHasVBC* node_has_vbc_s;
	const NodeHasVBC* node_has_vbc_f;
	const NodeVBCVec* node_vbc_vec_s;
	const NodeVBCVec* node_vbc_vec_f;
	double* node_am_s;
	double* node_am_f;
	double* node_de_vol_s;
	double* node_de_vol_f;

	double Kf, miu, k;
	// cavitation
	double m_cav, f_cav_end, u_cav_off, u_div_u_cav_lim;
	double* pcl_u_cav;
	double* pcl_is_cavitated;
	double* elem_u_cav;

	ParaUtil::PclSort pcl_sort;
	ParaUtil::NodeElemSort<4> ne_sort;

	size_t* in_pcl_in_elems, * in_prev_pcl_ids;
	const size_t *pcl_in_elems, *prev_pcl_ids;
	const size_t *elem_ids, *node_ids, *node_elem_offs;
	
	InitPcl init_pcl;
	MapPclToBgMesh map_pcl_to_mesh;
	ContactRigidBody cont_rigid_body;
	UpdateAccelerationAndVelocity update_a_and_v;
	CalElemDeAndMapToNode cal_elem_de;
	CalNodeDe cal_node_de;
	MapBgMeshToPcl map_mesh_to_pcl;

	InitPclTbb init_pcl_tbb;
	MapPclToBgMeshTbb map_pcl_to_mesh_tbb;
	MapBgMeshToPclTbb map_mesh_to_pcl_tbb;

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

	// time profiling
	size_t pcl_sort_time;
	size_t ne_sort_time;
	size_t map_pcl_to_mesh_time;
	size_t update_a_and_v_time;
	size_t cal_elem_de_time;
	size_t cal_node_de_time;
	size_t map_mesh_to_pcl_time;

	int init_calculation() override;
	friend int substep_func_T3D_CHM_TBB(void* _self);
	int finalize_calculation() override;
	
public:
	Step_T3D_CHM_TBB(const char* _name);
	~Step_T3D_CHM_TBB();

	friend struct Model_T3D_CHM_mt_hdf5_utilities::ParticleData;
	friend int Model_T3D_CHM_mt_hdf5_utilities::output_pcl_data_to_hdf5_file(
		Model_T3D_CHM_mt& md, Step_T3D_CHM_TBB& stp, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_CHM_mt_hdf5_utilities::time_history_complete_output_to_hdf5_file(
		Model_T3D_CHM_mt& md, Step_T3D_CHM_TBB& stp, ResultFile_hdf5& rf, hid_t frame_grp_id);
	friend int Model_T3D_CHM_mt_hdf5_utilities::load_model_from_hdf5_file(
		Model_T3D_CHM_mt& md, Step_T3D_CHM_TBB& step, const char* hdf5_name,
		const char* th_name, size_t frame_id);
};

#endif