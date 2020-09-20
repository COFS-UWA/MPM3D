#ifndef __Step_T2D_ME_s_Geo_h__
#define __Step_T2D_ME_s_Geo_h__

#include "Step.h"
#include "Model_T2D_ME_s.h"

int solve_substep_T2D_ME_s_Geo(void *_self);
int solve_substep_T2D_ME_s_Geo_avg(void* _self);

// for single object only
class Step_T2D_ME_s_Geo : public Step
{
public:
	typedef Model_T2D_ME_s::Particle Particle;
	typedef Model_T2D_ME_s::Element Element;
	typedef Model_T2D_ME_s::Node Node;

protected:
	Model_T2D_ME_s* model;
	double damping_ratio;

	int init_calculation() override;
	friend int solve_substep_T2D_ME_s_Geo(void *_self);
	friend int solve_substep_T2D_ME_s_Geo_avg(void* _self);
	int finalize_calculation() override;

public:
	Step_T2D_ME_s_Geo(const char *_name);
	~Step_T2D_ME_s_Geo();

	inline void set_model(Model_T2D_ME_s &md)
	{
		Step::set_model(md);
		model = &md;
	}

	// Restart from previous step
	inline void set_prev_step(Step_T2D_ME_s_Geo &prev_step)
	{
		Step::set_prev_step(prev_step);
		model = prev_step.model;
	}

	// convergence criteria
	inline void set_unbalanced_nodal_force_ratio_bound(double bound = 0.0) { f_ub_ratio_bound = bound; }
	inline void set_kinetic_energy_ratio_bound(double bound = 0.0) { e_kin_ratio_bound = bound; }

	inline void set_damping_ratio(double _ratio) noexcept { damping_ratio = _ratio; }

	// unbalanced nodal force
	inline double get_max_nf_ub() const noexcept { return max_f_ub; }
	inline double get_nf_ub() const noexcept { return sqrt(f_ub); }
	inline double get_nf_ub_ratio() const noexcept { return f_ub_ratio; }
	// kinetic energy
	inline double get_kinetic_energy() const noexcept { return 0.5 * e_kin; }
	inline double get_kinetic_energy_ratio() const noexcept { return e_kin_ratio; }

protected:
	// convergence criteria
	bool *node_has_a_or_v_bc;
	// unbalanced force
	double f_ub;
	double max_f_ub;
	double init_f_ub;
	double f_ub_ratio;
	// maximum kinematic energy
	double e_kin;
	double e_kin_max;
	double e_kin_prev;
	double e_kin_ratio;

	double f_ub_ratio_bound;
	double e_kin_ratio_bound;
};

#endif