#ifndef __Step_T2D_CHM_p_Geo_h__
#define __Step_T2D_CHM_p_Geo_h__

#include <fstream>

#include "Step.h"
#include "Model_T2D_CHM_p.h"

int solve_substep_T2D_CHM_p_Geo(void *_self);

// for single object only
class Step_T2D_CHM_p_Geo : public Step
{
public:
	typedef Model_T2D_CHM_p::Particle Particle;
	typedef Model_T2D_CHM_p::Element Element;
	typedef Model_T2D_CHM_p::Node Node;

protected:
	Model_T2D_CHM_p* model;

	int init_calculation() override;
	friend int solve_substep_T2D_CHM_p_Geo(void *_self);
	int finalize_calculation() override;

	double damping_ratio;

public:
	Step_T2D_CHM_p_Geo(const char* _name);
	~Step_T2D_CHM_p_Geo();

	inline void set_model(Model_T2D_CHM_p &md)
	{
		Step::set_model(md);
		model = &md;
	}

	// Restart from previous step
	inline void set_prev_step(Step &prev_step)
	{
		Step::set_prev_step(prev_step);
		model = static_cast<Model_T2D_CHM_p *>(&prev_step.get_model());
	}

	inline void set_damping_ratio(double _ratio) noexcept { damping_ratio = _ratio; }

protected:
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