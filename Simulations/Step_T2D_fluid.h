#ifndef __Step_T2D_fluid_H__
#define __Step_T2D_fluid_H__

#include "Step.h"
#include "Model_T2D_fluid.h"

// calculate shear stress when updating pcl variables
int solve_substep_T2D_fluid(void *_self);
// calculate shear stress after mapping velocity
int solve_substep_T2D_fluid2(void *_self);
// mixed integration
int solve_substep_T2D_fluid_MI(void *_self);

// for single object only
class Step_T2D_fluid : public Step
{
public:
	typedef Model_T2D_fluid::Element Element;
	typedef Model_T2D_fluid::Node Node;
	typedef Model_T2D_fluid::Particle Particle;

protected:
	int init_calculation(void) override;

	friend int solve_substep_T2D_fluid(void *_self);
	friend int solve_substep_T2D_fluid2(void *_self);
	friend int solve_substep_T2D_fluid_MI(void *_self);

	int finalize_calculation(void) override;

public:
	Step_T2D_fluid();
	~Step_T2D_fluid();

	inline void set_model(Model_T2D_fluid &md)
	{
		Step::set_model(md);
		model = &md;
	}

	// Restart from previous step
	inline void set_prev_step(Step_T2D_fluid &prev_step)
	{
		Step::set_prev_step(prev_step);
		model = prev_step.model;
	}

	inline void set_dtime(
		double _dt,
		double dt_max_min_raio = 0.1 /*ad hoc number*/,
		double t_tol_r = 0.01
		)
	{
		max_dt = _dt;
		min_dt = max_dt * dt_max_min_raio;
		time_tol_ratio = t_tol_r;
		time_tol = min_dt * t_tol_r;
		dtime = max_dt;
	}

protected:
	Model_T2D_fluid *model;
	double min_dt, max_dt;
};

#endif