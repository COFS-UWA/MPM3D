#ifndef __Step_T2D_CHM_d_h__
#define __Step_T2D_CHM_d_h__

#include "Step.h"
#include "Model_T2D_CHM_d.h"

int solve_substep_T2D_CHM_d(void *_self);

class Step_T2D_CHM_d : public Step
{
public:
	typedef Model_T2D_CHM_d::Element Element;
	typedef Model_T2D_CHM_d::Node Node;
	typedef Model_T2D_CHM_d::SolidParticle SolidParticle;
	typedef Model_T2D_CHM_d::FluidParticle FluidParticle;

protected:
	Model_T2D_CHM_d* model;

	double damping_ratio;

	int init_calculation() override;
	friend int solve_substep_T2D_CHM_d(void *_self);
	int finalize_calculation() override;

public:
	Step_T2D_CHM_d(const char* _name);
	~Step_T2D_CHM_d();

	inline void set_model(Model_T2D_CHM_d &md)
	{
		Step::set_model(md);
		model = &md;
	}

	// restart from previous step
	inline void set_prev_step(Step_T2D_CHM_d &prev_step)
	{
		Step::set_prev_step(prev_step);
		model = prev_step.model;
	}

	inline void set_damping_ratio(double _ratio) { damping_ratio = _ratio; }
};

#endif