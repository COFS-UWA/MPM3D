#ifndef __Step_T2D_ME_Task_h__
#define __Step_T2D_ME_Task_h__

#include "tbb/task.h"
#include "SortTask.h"
#include "Model_T2D_ME_mt.h"

class Step_T2D_ME_TBB;

namespace Step_T2D_ME_Task
{
	struct TaskData
	{
		Step_T2D_ME_TBB& stp;
		Model_T2D_ME_mt& md;
		SortUtils::SortMem &pcl_sort_mem;
		SortUtils::SortMem &node_sort_mem;

		// pcl data
		const double* const pcl_m;
		const Model_T2D_ME_mt::Force* const pcl_bf;
		const Model_T2D_ME_mt::Force* const pcl_t;
		const Model_T2D_ME_mt::Position* const pcl_pos;
		double* const pcl_vol;
		MatModel::MaterialModel** const pcl_mat_model;

		Model_T2D_ME_mt::SortedPclVarArrays spvas[2];

		// elem data
		const Model_T2D_ME_mt::ElemNodeIndex* const elem_node_id;
		const Model_T2D_ME_mt::ShapeFuncAB* const elem_dN_ab;
		const Model_T2D_ME_mt::ShapeFuncC* const elem_dN_c;
		const double* const elem_area;

		double* const elem_pcl_m;
		double* const elem_density;
		Model_T2D_ME_mt::StrainInc* const elem_de;
		double* const elem_m_de_vol;

		size_t *const valid_elems;
		size_t* const tmp_valid_elems;

		// elem node data
		Model_T2D_ME_mt::ElemNodeVM* const elem_node_vm;
		Model_T2D_ME_mt::Force* const elem_node_force;

		// node data
		Model_T2D_ME_mt::Acceleration* const node_a;
		Model_T2D_ME_mt::Velocity* const node_v;
		Model_T2D_ME_mt::NodeHasVBC* const node_has_vbc;
		double* const node_am;
		double* const node_de_vol;

		const size_t pcl_num_per_map_pcl_to_mesh_task;
		const size_t node_num_per_update_a_and_v_task;
		const size_t elem_num_per_cal_elem_de_task;
		const size_t node_num_per_cal_node_de_task;
		const size_t pcl_num_per_task_map_mesh_to_pcl;

		const size_t pcl_digit_num;
		const size_t node_digit_num;
#ifdef _DEBUG
		const size_t ori_pcl_num;
		const size_t elem_num;
		const size_t node_num;
#endif
		size_t valid_pcl_num;
		size_t valid_elem_num;
		size_t sorted_pcl_var_id;

		TaskData(Step_T2D_ME_TBB& _stp,
			size_t _pcl_num_per_map_pcl_to_mesh_task,
			size_t _node_num_per_update_a_and_v_task,
			size_t _elem_num_per_cal_elem_de_task,
			size_t _node_num_per_cal_node_de_task,
			size_t _pcl_num_per_task_map_mesh_to_pcl);
		~TaskData() {}
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
		size_t start_id, end_id;
		TaskData& td;
	public:
		MapPclToBgMeshTask(size_t _start_id, size_t _end_id, TaskData &_td);
		~MapPclToBgMeshTask() {}
		tbb::task* execute() override;
	};

	class UpdateAccelerationAndVelocityTask : public tbb::task
	{
	protected:
		size_t start_id, end_id;
		TaskData& td;
	public:
		UpdateAccelerationAndVelocityTask(
			size_t _start_id,
			size_t _end_id,
			TaskData& _td) :
			start_id(_start_id),
			end_id(_end_id),
			td(_td) {}
		~UpdateAccelerationAndVelocityTask() {}
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

	class CalElemDeAndMapToNode : public tbb::task
	{
	protected:
		size_t start_id, end_id;
		TaskData& td;
	public:
		CalElemDeAndMapToNode(
			size_t _start_id,
			size_t _end_id,
			TaskData& _td) :
			start_id(_start_id),
			end_id(_end_id),
			td(_td) {}
		~CalElemDeAndMapToNode() {}
		tbb::task* execute() override;
	};

	class CalNodeDe : public tbb::task
	{
	protected:
		size_t start_id, end_id;
		TaskData& td;
	public:
		CalNodeDe(
			size_t _start_id,
			size_t _end_id,
			TaskData& _td) :
			start_id(_start_id),
			end_id(_end_id),
			td(_td) {}
		~CalNodeDe() {}
		tbb::task* execute() override;
	};

	class MapBgMeshToPclTask : public tbb::task
	{
	protected:
		size_t start_id, end_id;
		TaskData& td;
	public:
		MapBgMeshToPclTask(
			size_t _start_id,
			size_t _end_id,
			TaskData& _td) :
			start_id(_start_id),
			end_id(_end_id),
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
}

#endif