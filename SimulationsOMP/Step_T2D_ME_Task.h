#ifndef __Step_T2D_ME_Task_h__
#define __Step_T2D_ME_Task_h__

#include "tbb/task.h"
#include "SortParticleTask.h"
#include "SortTriMeshNodeTask.h"
#include "Model_T2D_ME_mt.h"

class Step_T2D_ME_TBB;

namespace Step_T2D_ME_Task
{
	struct TaskData
	{
		Model_T2D_ME_mt *pmodel;

		// pcl data
		const double* pcl_m;
		const Model_T2D_ME_mt::Force* pcl_bf;
		const Model_T2D_ME_mt::Force* pcl_t;
		const Model_T2D_ME_mt::Position* pcl_pos;
		double* pcl_vol;
		MatModel::MaterialModel** pcl_mat_model;

		Model_T2D_ME_mt::SortedPclVarArrays spvas[2];

		// elem data
		const Model_T2D_ME_mt::ElemNodeIndex* elem_node_id;
		const Model_T2D_ME_mt::ShapeFuncAB* elem_dN_ab;
		const Model_T2D_ME_mt::ShapeFuncC* elem_dN_c;
		const double* elem_area;

		double* elem_pcl_m;
		double* elem_density;
		Model_T2D_ME_mt::StrainInc* elem_de;
		double* elem_m_de_vol;
		
		// elem node data
		Model_T2D_ME_mt::ElemNodeVM* elem_node_vm;
		Model_T2D_ME_mt::Force* elem_node_force;

		// node data
		Model_T2D_ME_mt::Acceleration* node_a;
		Model_T2D_ME_mt::Velocity* node_v;
		Model_T2D_ME_mt::NodeHasVBC* node_has_vbc;
		double* node_am;
		double* node_de_vol;

#ifdef _DEBUG
		size_t ori_pcl_num;
		size_t elem_num;
		size_t node_num;
		size_t prev_valid_pcl_num;
#endif

		const size_t pcl_num_per_map_pcl_to_mesh_task;
		const size_t node_num_per_update_a_and_v_task;
		const size_t elem_num_per_cal_elem_de_task;
		const size_t node_num_per_cal_node_de_task;
		const size_t pcl_num_per_task_map_mesh_to_pcl;

		// cal data
		size_t sorted_pcl_var_id;
		double dt;
		SortUtils::SortParticleMem pcl_sort_mem;
		SortUtils::SortTriMeshNodeMem node_sort_mem;

		TaskData(
			size_t _pcl_num_per_map_pcl_to_mesh_task,
			size_t _node_num_per_update_a_and_v_task,
			size_t _elem_num_per_cal_elem_de_task,
			size_t _node_num_per_cal_node_de_task,
			size_t _pcl_num_per_task_map_mesh_to_pcl);
		~TaskData() {}
		void set_model(Model_T2D_ME_mt& md) noexcept;
	};
	
	class MapPclToBgMeshTask : public tbb::task
	{
	protected:
		typedef Model_T2D_ME_mt::Force Force;
		typedef Model_T2D_ME_mt::Acceleration Acceleration;
		typedef Model_T2D_ME_mt::Velocity Velocity;
		typedef Model_T2D_ME_mt::Displacement Displacement;
		typedef Model_T2D_ME_mt::Stress Stress;
		typedef Model_T2D_ME_mt::ShapeFunc ShapeFunc;
		typedef Model_T2D_ME_mt::ShapeFuncAB ShapeFuncAB;
		typedef Model_T2D_ME_mt::ShapeFuncC ShapeFuncC;
		typedef Model_T2D_ME_mt::ElemNodeVM ElemNodeVM;
		size_t p_id0, p_id1;
		TaskData& td;
	public:
		MapPclToBgMeshTask(
			size_t start_id,
			size_t end_id,
			TaskData& _td) :
			p_id0(start_id),
			p_id1(end_id),
			td(_td) {}
		~MapPclToBgMeshTask() {}
		tbb::task* execute() override;
	};

	class UpdateAccelerationAndVelocityTask : public tbb::task
	{
	protected:
		typedef Model_T2D_ME_mt::Force Force;
		typedef Model_T2D_ME_mt::ElemNodeVM ElemNodeVM;
		typedef Model_T2D_ME_mt::Acceleration Acceleration;
		typedef Model_T2D_ME_mt::Velocity Velocity;
		typedef Model_T2D_ME_mt::NodeHasVBC NodeHasVBC;
		size_t ve_id0, ve_id1;
		TaskData& td;
	public:
		UpdateAccelerationAndVelocityTask(
			size_t start_id,
			size_t end_id,
			TaskData& _td) :
			ve_id0(start_id),
			ve_id1(end_id),
			td(_td) {}
		~UpdateAccelerationAndVelocityTask() {}
		tbb::task* execute() override;
	};

	class CalElemDeAndMapToNode : public tbb::task
	{
	protected:
		typedef Model_T2D_ME_mt::ElemNodeIndex ElemNodeIndex;
		typedef Model_T2D_ME_mt::Velocity Velocity;
		typedef Model_T2D_ME_mt::ShapeFuncAB ShapeFuncAB;
		typedef Model_T2D_ME_mt::StrainInc StrainInc;

		size_t ve_id0, ve_id1;
		TaskData& td;
	public:
		CalElemDeAndMapToNode(
			size_t start_id,
			size_t end_id,
			TaskData& _td) :
			ve_id0(start_id),
			ve_id1(end_id),
			td(_td) {}
		~CalElemDeAndMapToNode() {}
		tbb::task* execute() override;
	};

	class CalNodeDe : public tbb::task
	{
	protected:
		size_t ve_id0, ve_id1;
		TaskData& td;
	public:
		CalNodeDe(
			size_t start_id,
			size_t end_id,
			TaskData& _td) :
			ve_id0(start_id),
			ve_id1(end_id),
			td(_td) {}
		~CalNodeDe() {}
		tbb::task* execute() override;
	};

	class MapBgMeshToPclTask : public tbb::task
	{
	protected:
		typedef Model_T2D_ME_mt::ElemNodeIndex ElemNodeIndex;
		typedef Model_T2D_ME_mt::Acceleration Acceleration;
		typedef Model_T2D_ME_mt::Velocity Velocity;
		typedef Model_T2D_ME_mt::Displacement Displacement;
		typedef Model_T2D_ME_mt::Position Position;
		typedef Model_T2D_ME_mt::Stress Stress;
		typedef Model_T2D_ME_mt::Strain Strain;
		typedef Model_T2D_ME_mt::StrainInc StrainInc;
		typedef Model_T2D_ME_mt::ShapeFunc ShapeFunc;
		size_t p_id0, p_id1;
		TaskData& td;
	public:
		MapBgMeshToPclTask(
			size_t start_id,
			size_t end_id,
			TaskData& _td) :
			p_id0(start_id),
			p_id1(end_id),
			td(_td) {}
		~MapBgMeshToPclTask() {}
		tbb::task* execute() override;
	};

	class Step_T2D_ME_Task : public tbb::task
	{
	protected:
		TaskData &td;
	public:
		Step_T2D_ME_Task(TaskData& _td) : td(_td) {}
		~Step_T2D_ME_Task() {}
		tbb::task *execute() override;
	};

	class InitPclTask : public tbb::task
	{
	protected:
		typedef Model_T2D_ME_mt::ShapeFunc ShapeFunc;
		typedef Model_T2D_ME_mt::Displacement Displacement;
		typedef Model_T2D_ME_mt::Position Position;
		size_t p_id0, p_id1;
		TaskData& td;
	public:
		InitPclTask(
			size_t start_id,
			size_t end_id,
			TaskData& _td) :
			p_id0(start_id),
			p_id1(end_id),
			td(_td) {}
		~InitPclTask() {}
		tbb::task* execute() override;
	};

	class InitTask : public tbb::task
	{
	protected:
		TaskData& td;
	public:
		InitTask(TaskData& _td) : td(_td) {}
		~InitTask() {}
		tbb::task* execute() override;
	};

	//class ContactWithRigidObejctTask : public tbb::task
	//{
	//protected:

	//public:
	//	ContactWithRigidObejctTask();
	//	~ContactWithRigidObejctTask();
	//	tbb::task* execute() override;
	//};
}

#endif