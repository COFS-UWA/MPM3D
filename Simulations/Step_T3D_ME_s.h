#ifndef __Step_T3D_ME_s_h__
#define __Step_T3D_ME_s_h__

#include "Step.h"
#include "Model_T3D_ME_s.h"

int solve_substep_T3D_ME_s(void *_self);

// for single object only
class Step_T3D_ME_s : public Step
{
public:
	typedef Model_T3D_ME_s::Particle Particle;
	typedef Model_T3D_ME_s::Element Element;
	typedef Model_T3D_ME_s::Node Node;
	typedef Model_T3D_ME_s::SearchingGrid SearchingGrid;
	typedef SearchingGrid::Grid Grid;

protected:
	Model_T3D_ME_s *model;

	// local damping factor
	double damping_ratio;

	int init_calculation() override;
	friend int solve_substep_T3D_ME_s(void *_self);
	int finalize_calculation() override;

	int apply_rb_to_mesh(RigidTetrahedronMesh& rb);
	bool apply_pcl_contact_force(Particle &pcl, RigidTetrahedronMesh& rb);

public:
	explicit Step_T3D_ME_s(const char *_name);
	~Step_T3D_ME_s();

	inline void set_model(Model_T3D_ME_s &md)
	{
		Step::set_model(md);
		model = &md;
	}

	// Restart from previous step
	inline void set_prev_step(Step_T3D_ME_s &prev_step)
	{
		Step::set_prev_step(prev_step);
		model = prev_step.model;
	}

	inline void set_damping_ratio(double _ratio) noexcept { damping_ratio = _ratio; }
};

#endif