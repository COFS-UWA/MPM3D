#ifndef __Step_T2D_ME_mt_h__
#define __Step_T2D_ME_mt_h__

#include "Step.h"
#include "Model_T2D_ME_mt.h"

int solve_substep_T2D_ME_mt(void* _self);

class Step_T2D_ME_mt : public Step
{
protected:
	uint32_t thread_num;

	// task division
	uint32_t *pcl_range, *elem_range;
	uint32_t* node_elem_range, *node_range;
	uint32_t* vx_bc_range, * vy_bc_range;

	// for counting sort
	uint32_t* new_to_ori_pcl_map0;
	uint32_t* new_to_ori_pcl_map1;
	uint32_t* pcl_in_elem_array0;
	uint32_t* pcl_in_elem_array1;

	uint32_t *elem_bin;
	uint32_t *elem_bin_range;
	uint32_t *elem_bin_offset;

	void clear_mem();

	int init_calculation() override;
	friend int solve_substep_T2D_ME_mt(void* _self);
	int finalize_calculation() override;

public:
	typedef Model_T2D_ME_mt::PclBodyForce PclBodyForce;
	typedef Model_T2D_ME_mt::PclTraction PclTraction;
	typedef Model_T2D_ME_mt::PclPos PclPos;

	typedef Model_T2D_ME_mt::PclSortedVarArray PclSortedVarArray;
	typedef Model_T2D_ME_mt::PclDisp PclDisp;
	typedef Model_T2D_ME_mt::PclV PclV;
	typedef Model_T2D_ME_mt::PclShapeFunc PclShapeFunc;
	typedef Model_T2D_ME_mt::PclStress PclStress;

	typedef Model_T2D_ME_mt::ElemNodeIndex ElemNodeIndex;
	typedef Model_T2D_ME_mt::ElemShapeFuncAB ElemShapeFuncAB;
	typedef Model_T2D_ME_mt::ElemShapeFuncC ElemShapeFuncC;

	typedef Model_T2D_ME_mt::ElemStrainInc ElemStrainInc;
	typedef Model_T2D_ME_mt::ElemStress ElemStress;

	typedef Model_T2D_ME_mt::ElemNodeVM ElemNodeVM;
	typedef Model_T2D_ME_mt::ElemNodeForce ElemNodeForce;

	Step_T2D_ME_mt(const char* _name);
	~Step_T2D_ME_mt();

	inline void set_thread_num(uint32_t th_num) noexcept { thread_num = th_num; }
};

#endif