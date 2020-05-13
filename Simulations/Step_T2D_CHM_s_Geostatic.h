#ifndef __Step_T2D_CHM_s_SE_Geostatic_H__
#define __Step_T2D_CHM_s_SE_Geostatic_H__

#include <fstream>

#include "Step.h"
#include "Model_T2D_CHM_s.h"

int solve_substep_T2D_CHM_s_SE_Geostatic(void *_self);

// for single object only
class Step_T2D_CHM_s_SE_Geostatic : public Step
{
public:
	typedef Model_T2D_CHM_s::Particle Particle;
	typedef Model_T2D_CHM_s::Element Element;
	typedef Model_T2D_CHM_s::Node Node;

protected:
	int init_calculation(void) override;
	friend int solve_substep_T2D_CHM_s_SE_Geostatic(void *_self);
	int finalize_calculation(void) override;

public:
	Step_T2D_CHM_s_SE_Geostatic();
	~Step_T2D_CHM_s_SE_Geostatic();

	inline void set_model(Model_T2D_CHM_s &md)
	{
		Step::set_model(md);
		model = &md;
	}

	// Restart from previous step
	inline void set_prev_step(Step &prev_step)
	{
		Step::set_prev_step(prev_step);
		model = static_cast<Model_T2D_CHM_s *>(&prev_step.get_model());
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

	inline void set_mass_scale(double _ms_sr, double _mf_sr) noexcept { ms_sr = _ms_sr; mf_sr = _mf_sr; }

protected:
	Model_T2D_CHM_s *model;
	double min_dt, max_dt;
	double ms_sr, mf_sr; // mass scaling ratio

	// convergence criteria
	// unbalanced force
	double init_f_ub;
	bool init_f_ub_is_init;
	double f_ub_ratio;
	// maximum kinematic energy
	double e_kin_max;
	bool e_kin_max_is_init;
	double e_kin_prev;
	double e_kin_ratio;

	double f_ub_ratio_bound;
	double e_kin_ratio_bound;

	// for debugging
public:
	std::fstream out_file;
};

#endif