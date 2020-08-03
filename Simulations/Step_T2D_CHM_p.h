#ifndef __Step_T2D_CHM_p_h__
#define __Step_T2D_CHM_p_h__

#include "Step.h"
#include "Model_T2D_CHM_p.h"

int solve_substep_T2D_CHM_p(void *_self);

class Step_T2D_CHM_p : public Step
{
public:
	typedef Model_T2D_CHM_p::Particle Particle;
	typedef Model_T2D_CHM_p::Element Element;
	typedef Model_T2D_CHM_p::Node Node;

protected:
	Model_T2D_CHM_p* model;

	int init_calculation() override;
	friend int solve_substep_T2D_CHM_p(void *_self);
	int finalize_calculation() override;

	double damping_ratio;

public:
	Step_T2D_CHM_p(const char* _name);
	~Step_T2D_CHM_p();

	inline void set_model(Model_T2D_CHM_p &md)
	{
		Step::set_model(md);
		model = &md;
	}

	// restart from previous step
	inline void set_prev_step(Step_T2D_CHM_p &prev_step)
	{
		Step::set_prev_step(prev_step);
		model = prev_step.model;
	}

	inline void set_damping_ratio(double _ratio) { damping_ratio = _ratio; }
};

#endif