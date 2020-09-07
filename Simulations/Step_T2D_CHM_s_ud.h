#ifndef __Step_T2D_CHM_s_ud_ud_h__
#define __Step_T2D_CHM_s_ud_ud_h__

#include "Step.h"
#include "Model_T2D_CHM_s.h"

// undrained version of coupled hydro-mechanical code
// no seepage occured

int solve_substep_T2D_CHM_s_ud(void *_self);

class Step_T2D_CHM_s_ud : public Step
{
public:
	typedef Model_T2D_CHM_s::Particle Particle;
	typedef Model_T2D_CHM_s::Element Element;
	typedef Model_T2D_CHM_s::Node Node;

protected:
	Model_T2D_CHM_s* model;

	double damping_ratio;

	int init_calculation() override;
	friend int solve_substep_T2D_CHM_s_ud(void *_self);
	int finalize_calculation() override;

	void apply_rigid_circle(double dtime);

public:
	Step_T2D_CHM_s_ud(const char* _name);
	~Step_T2D_CHM_s_ud();

	inline void set_model(Model_T2D_CHM_s &md)
	{
		Step::set_model(md);
		model = &md;
	}

	// restart from previous step
	inline void set_prev_step(Step_T2D_CHM_s_ud &prev_step)
	{
		Step::set_prev_step(prev_step);
		model = prev_step.model;
	}

	inline void set_damping_ratio(double _ratio) { damping_ratio = _ratio; }
};

#endif