#ifndef __Step_T2D_CHM_s_SE_h__
#define __Step_T2D_CHM_s_SE_h__

#include "Step.h"
#include "Model_T2D_CHM_s.h"

int solve_substep_T2D_CHM_s_SE(void *_self);

// for single object only
class Step_T2D_CHM_s_SE : public Step
{
public:
	typedef Model_T2D_CHM_s::Particle Particle;
	typedef Model_T2D_CHM_s::Element Element;
	typedef Model_T2D_CHM_s::Node Node;

protected:
	int init_calculation(void) override;
	friend int solve_substep_T2D_CHM_s_SE(void *_self);
	int finalize_calculation(void) override;

public:
	Step_T2D_CHM_s_SE();
	~Step_T2D_CHM_s_SE();

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

	inline void set_damping_ratio(double ra) noexcept { damping_ratio = ra; }
	inline void set_bv_ratio(double ra) noexcept { bv_ratio = ra; }
	inline void set_mass_scale(double _ms_sr, double _mf_sr) noexcept { ms_sr = _ms_sr; mf_sr = _mf_sr; }
	
protected:
	Model_T2D_CHM_s *model;
	double min_dt, max_dt;
	double damping_ratio; // local damping
	double bv_ratio; // bulk viscosity damping ratio
	double ms_sr, mf_sr; // mass scaling ratio
};

#endif