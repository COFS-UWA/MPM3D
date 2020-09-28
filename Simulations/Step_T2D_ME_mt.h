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
	Step_T2D_ME_mt(const char* _name);
	~Step_T2D_ME_mt();
};

#endif