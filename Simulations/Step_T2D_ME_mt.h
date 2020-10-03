#ifndef __Step_T2D_ME_mt_h__
#define __Step_T2D_ME_mt_h__

#include "Step.h"
#include "Model_T2D_ME_mt.h"

int solve_substep_T2D_ME_mt(void* _self);

class Step_T2D_ME_mt : public Step
{
protected:
	int init_calculation() override;
	friend int solve_substep_T2D_ME_mt(void* _self);
	int finalize_calculation() override;

public:
	typedef Model_T2D_ME_mt::PclMass PclMass;
	typedef Model_T2D_ME_mt::PclBodyForce PclBodyForce;
	typedef Model_T2D_ME_mt::PclTraction PclTraction;
	typedef Model_T2D_ME_mt::PclPos PclPos;

	typedef Model_T2D_ME_mt::PclSortedVarArray PclSortedVarArray;
	typedef Model_T2D_ME_mt::PclIndex PclIndex;
	typedef Model_T2D_ME_mt::PclDensity PclDensity;
	typedef Model_T2D_ME_mt::PclDisp PclDisp;
	typedef Model_T2D_ME_mt::PclV PclV;
	typedef Model_T2D_ME_mt::PclShapeFunc PclShapeFunc;
	typedef Model_T2D_ME_mt::PclStress PclStress;
	typedef Model_T2D_ME_mt::ElemPclList ElemPclList;

	typedef Model_T2D_ME_mt::ElemNodeIndex ElemNodeIndex;
	typedef Model_T2D_ME_mt::ElemArea ElemArea;
	typedef Model_T2D_ME_mt::ElemShapeFuncAB ElemShapeFuncAB;
	typedef Model_T2D_ME_mt::ElemShapeFuncC ElemShapeFuncC;

	typedef Model_T2D_ME_mt::ElemDensity ElemDensity;
	typedef Model_T2D_ME_mt::ElemStrainInc ElemStrainInc;
	typedef Model_T2D_ME_mt::ElemStress ElemStress;
	typedef Model_T2D_ME_mt::ElemAm ElemAm;
	typedef Model_T2D_ME_mt::ElemAmDeVol ElemAmDeVol;

	typedef Model_T2D_ME_mt::ElemNodeVM ElemNodeVM;
	typedef Model_T2D_ME_mt::ElemNodeForce ElemNodeForce;

	typedef Model_T2D_ME_mt::NodeElemList NodeElemList;
	typedef Model_T2D_ME_mt::NodeA NodeA;
	typedef Model_T2D_ME_mt::NodeV NodeV;
	typedef Model_T2D_ME_mt::NodeAm NodeAm;
	typedef Model_T2D_ME_mt::NodeDeVol NodeDeVol;
	
	Step_T2D_ME_mt(const char* _name);
	~Step_T2D_ME_mt();
};

#endif