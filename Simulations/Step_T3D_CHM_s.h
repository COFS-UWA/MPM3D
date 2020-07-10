#ifndef __Step_T3D_CHM_s_h__
#define __Step_T3D_CHM_s_h__

#include "Step.h"
#include "Model_T3D_CHM_s.h"

int solve_substep_T3D_CHM_s(void *_self);

// for single object only
class Step_T3D_CHM_s : public Step
{
public:
	typedef Model_T3D_CHM_s::Node Node;
	typedef Model_T3D_CHM_s::Element Element;
	typedef Model_T3D_CHM_s::Particle Particle;
	
protected:
	int init_calculation() override;
	friend int solve_substep_T3D_CHM_s(void *_self);
	int finalize_calculation() override;

public:
	Step_T3D_CHM_s(const char *_name);
	~Step_T3D_CHM_s();

	inline void set_model(Model_T3D_CHM_s &md)
	{
		Step::set_model(md);
		model = &md;
	}

	// Restart from previous step
	inline void set_prev_step(Step &prev_step)
	{
		Step::set_prev_step(prev_step);
		model = static_cast<Model_T3D_CHM_s *>(&prev_step.get_model());
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
	
protected:
	Model_T3D_CHM_s *model;
	double min_dt, max_dt;
	double damping_ratio; // local damping
};

#endif