#ifndef __Step_T2D_ME_s_H__
#define __Step_T2D_ME_s_H__

#include "Step.h"
#include "Model_T2D_ME_s.h"

int solve_substep_T2D_ME_s(void *_self);

// for single object only
class Step_T2D_ME_s : public Step
{
public:
	typedef Model_T2D_ME_s::Particle Particle;
	typedef Model_T2D_ME_s::Element Element;
	typedef Model_T2D_ME_s::Node Node;

protected:
	int init_calculation(void) override;
	friend int solve_substep_T2D_ME_s(void *_self);
	int finalize_calculation(void) override;

public:
	Step_T2D_ME_s();
	~Step_T2D_ME_s();

	inline void set_model(Model_T2D_ME_s &md)
	{
		Step::set_model(md);
		model = &md;
	}

	// Restart from previous step
	inline void set_prev_step(Step_T2D_ME_s &prev_step)
	{
		Step::set_prev_step(prev_step);
		model = prev_step.model;
	}

	inline void set_dtime(double _dt,
		double dt_max_min_raio = 0.1 /*ad hoc number*/,
		double t_tol_r = 0.01)
	{
		max_dt = _dt;
		min_dt = max_dt * dt_max_min_raio;
		time_tol_ratio = t_tol_r;
		time_tol = min_dt * t_tol_r;
		dtime = max_dt;
	}

	inline void set_damping_ratio(double ra) noexcept { damping_ratio = ra; }
	inline void set_bv_ratio(double ra) noexcept { bv_ratio = ra; }

protected:
	Model_T2D_ME_s *model;
	double min_dt, max_dt;
	double damping_ratio; // local damping
	double bv_ratio; // bulk viscosity damping ratio, should be ??
};

#endif