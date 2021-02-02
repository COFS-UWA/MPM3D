#ifndef __Step_T2D_ME_TBB_task_h__
#define __Step_T2D_ME_TBB_task_h__

#include "tbb/task.h"
#include "SortTask.h"
#include "Model_T2D_ME_mt.h"

class Step_T2D_ME_TBB;

namespace Step_T2D_ME_TBB_task
{
	struct Range { size_t start_id, end_id; };

	struct CalData_T2D_ME_TBB
	{
		Step_T2D_ME_TBB& stp;
		Model_T2D_ME_mt& md;
		Model_T2D_ME_mt::SortedPclVarArrays& spva0, & spva1;
		SortUtils::SortMem& pcl_sort_mem;

		const double* const pcl_m;
		const Model_T2D_ME_mt::Force* const pcl_bf;
		const Model_T2D_ME_mt::Force* const pcl_t;
		const Model_T2D_ME_mt::Position* const pcl_pos;
		double* const pcl_vol;
		MatModel::MaterialModel** const pcl_mat_model;

		size_t* pcl_index0;
		double* pcl_density0;
		Model_T2D_ME_mt::Velocity* pcl_v0;
		Model_T2D_ME_mt::Displacement* pcl_disp0;
		Model_T2D_ME_mt::ShapeFunc* pcl_N0;
		Model_T2D_ME_mt::Stress* pcl_stress0;
		Model_T2D_ME_mt::Strain* pcl_strain0;
		Model_T2D_ME_mt::Strain* pcl_estrain0;
		Model_T2D_ME_mt::Strain* pcl_pstrain0;

		const size_t* const pcl_index1;
		const double* const pcl_density1;
		const Model_T2D_ME_mt::Velocity* const pcl_v1;
		const Model_T2D_ME_mt::Displacement* const pcl_disp1;
		const Model_T2D_ME_mt::ShapeFunc* const pcl_N1;
		const Model_T2D_ME_mt::Stress* const pcl_stress1;
		const Model_T2D_ME_mt::Strain* const pcl_strain1;
		const Model_T2D_ME_mt::Strain* const pcl_estrain1;
		const Model_T2D_ME_mt::Strain* const pcl_pstrain1;

		const Model_T2D_ME_mt::ElemNodeIndex* const elem_node_id;
		const Model_T2D_ME_mt::ShapeFuncAB* const elem_dN_ab;
		const Model_T2D_ME_mt::ShapeFuncC* const elem_dN_c;
		const double* const elem_area;

		double* const elem_pcl_m;
		double* const elem_density;
		Model_T2D_ME_mt::StrainInc* const elem_de;
		double* const elem_m_de_vol;

		Model_T2D_ME_mt::ElemNodeVM* const elem_node_vm;
		Model_T2D_ME_mt::Force* const elem_node_force;

		Model_T2D_ME_mt::Acceleration* const node_a;
		Model_T2D_ME_mt::Velocity* const node_v;
		Model_T2D_ME_mt::NodeHasVBC* const node_has_vbc;
		double* const node_am;
		double* const node_de_vol;

		size_t **valid_elem_arrays; // thread_num

		CalData_T2D_ME_TBB(Step_T2D_ME_TBB& _stp);
		~CalData_T2D_ME_TBB() {}
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
		const size_t block_id, block_num;
		CalData_T2D_ME_TBB& cd;
	public:
		MapPclToBgMeshTask(
			size_t _block_id,
			size_t _block_num,
			CalData_T2D_ME_TBB &_cd);
		~MapPclToBgMeshTask();
		tbb::task* execute() override;
	};

	class MapBgMeshToPclTask : public tbb::task
	{
	protected:
		

	public:
		MapBgMeshToPclTask(
			size_t _p_id0,
			size_t _p_id1,
			CalData_T2D_ME_TBB& _cd);
		~MapBgMeshToPclTask();
		tbb::task* execute() override;
	};
}

#endif