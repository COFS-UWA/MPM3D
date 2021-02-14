#ifndef __Step_T2D_ME_Task_h__
#define __Step_T2D_ME_Task_h__

#include "tbb/task.h"
#include "SortParticleTask.h"
#include "SortTriMeshNodeTask.h"
#include "Model_T2D_ME_mt.h"

class Step_T2D_ME_TBB;

namespace Step_T2D_ME_Task
{
	struct CalData
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
		
		const size_t pcl_num_per_init_pcl_task;
		const size_t pcl_num_per_map_pcl_to_mesh_task;
		const size_t elem_num_per_update_a_and_v_task;
		const size_t elem_num_per_cal_elem_de_task;
		const size_t elem_num_per_cal_node_de_task;
		const size_t pcl_num_per_map_mesh_to_pcl_task;

		size_t thread_num;
		
		// cal data
		SortUtils::SortParticleMem pcl_sort_mem;
		SortUtils::SortTriMeshNodeMem node_sort_mem;
		
#ifdef _DEBUG
		size_t elem_num;
		size_t node_num;
		size_t prev_valid_pcl_num;
		size_t ori_pcl_num;
#endif
		
		size_t init_pcl_task_num;
		size_t map_pcl_to_mesh_task_num;
		size_t update_a_and_v_task_num;
		size_t cal_elem_de_task_num;
		size_t cal_node_de_task_num;
		size_t map_mesh_to_pcl_task_num;
		size_t sorted_pcl_var_id;
		double dt;

		CalData(
			size_t _pcl_num_per_init_pcl_task,
			size_t _pcl_num_per_map_pcl_to_mesh_task,
			size_t _elem_num_per_update_a_and_v_task,
			size_t _elem_num_per_cal_elem_de_task,
			size_t _elem_num_per_cal_node_de_task,
			size_t _pcl_num_per_task_map_mesh_to_pcl);
		~CalData() {}
		void set_model(Model_T2D_ME_mt& md) noexcept;
	};
	
	class InitPcl
	{
	protected:
		typedef Model_T2D_ME_mt::ShapeFunc ShapeFunc;
		typedef Model_T2D_ME_mt::Displacement Displacement;
		typedef Model_T2D_ME_mt::Position Position;
		CalData& cd;
	public:
		InitPcl(CalData& _cd) : cd(_cd) {}
		~InitPcl() {}
		void work(size_t wk_id) const;
	};
	
	class MapPclToBgMesh
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
		CalData& cd;
	public:
		MapPclToBgMesh(CalData &_cd) : cd(_cd) {}
		~MapPclToBgMesh() {}
		void work(size_t wk_id) const;
	};

	class UpdateAccelerationAndVelocity
	{
	protected:
		typedef Model_T2D_ME_mt::Force Force;
		typedef Model_T2D_ME_mt::ElemNodeVM ElemNodeVM;
		typedef Model_T2D_ME_mt::Acceleration Acceleration;
		typedef Model_T2D_ME_mt::Velocity Velocity;
		typedef Model_T2D_ME_mt::NodeHasVBC NodeHasVBC;
		CalData& cd;
	public:
		UpdateAccelerationAndVelocity(CalData& _cd) : cd(_cd) {}
		~UpdateAccelerationAndVelocity() {}
		void work(size_t wk_id) const;
	};

	class CalElemDeAndMapToNode
	{
	protected:
		typedef Model_T2D_ME_mt::ElemNodeIndex ElemNodeIndex;
		typedef Model_T2D_ME_mt::Velocity Velocity;
		typedef Model_T2D_ME_mt::ShapeFuncAB ShapeFuncAB;
		typedef Model_T2D_ME_mt::StrainInc StrainInc;
		CalData& cd;
	public:
		CalElemDeAndMapToNode(CalData &_cd) : cd(_cd) {}
		~CalElemDeAndMapToNode() {}
		void work(size_t wk_id) const;
	};

	class CalNodeDe
	{
	protected:
		CalData& cd;
	public:
		CalNodeDe(CalData& _cd) : cd(_cd) {}
		~CalNodeDe() {}
		void work(size_t wk_id) const;
	};

	class MapBgMeshToPcl
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
		CalData& cd;
	public:
		MapBgMeshToPcl(CalData& _cd) : cd(_cd) {}
		~MapBgMeshToPcl() {}
		void work(size_t wk_id) const;
	};

	//class ContactWithRigidObejctTask : public tbb::task
	//{
	//protected:

	//public:
	//	ContactWithRigidObejctTask();
	//	~ContactWithRigidObejctTask();
	//	void work(size_t wk_id) const;
	//};
}

#endif