#ifndef __Step_T2D_ME_p_Geo_h__
#define __Step_T2D_ME_p_Geo_h__

#include <thread>
#include <mutex>
#include <atomic>

#include "ThreadBarrierFixedNum.h"
#include "Model_T2D_ME_p.h"
#include "Step.h"

int solve_substep_T2D_ME_p_geo(void *_self);

// for single object only
class Step_T2D_ME_p_Geo : public Step
{
public:
	typedef Model_T2D_ME_p::NodeToElem NodeToElem;
	typedef Model_T2D_ME_p::Node Node;
	typedef Model_T2D_ME_p::NodeVarAtElem NodeVarAtElem;
	typedef Model_T2D_ME_p::Element Element;
	typedef Model_T2D_ME_p::Particle Particle;

	explicit Step_T2D_ME_p_Geo(const char *_name);
	~Step_T2D_ME_p_Geo();

	inline void set_model(Model_T2D_ME_p &md)
	{
		Step::set_model(md);
		model = &md;
	}

	inline void set_prev_step(Step_T2D_ME_p_Geo &prev_step)
	{
		Step::set_prev_step(prev_step);
		model = prev_step.model;
	}

	// convergence criteria
	inline void set_unbalanced_nodal_force_ratio_bound(double bound = 0.0) { f_ub_ratio_bound = bound; }
	inline void set_kinetic_energy_ratio_bound(double bound = 0.0) { e_kin_ratio_bound = bound; }

	inline void set_damping_ratio(double _ratio) noexcept { damping_ratio = _ratio; }
	
	inline void set_thread_num(unsigned int num) { thread_num = num; }

	// unbalanced nodal force
	inline double get_nf_ub() const noexcept { return sqrt(f_ub); }
	inline double get_nf_ub_ratio() const noexcept { return f_ub_ratio; }
	// kinetic energy
	inline double get_kinetic_energy() const noexcept { return 0.5 * e_kin; }
	inline double get_kinetic_energy_ratio() const noexcept { return e_kin_ratio; }

protected:
	Model_T2D_ME_p* model;
	double damping_ratio;

	int init_calculation() override;
	friend int solve_substep_T2D_ME_p_geo(void* _self);
	int finalize_calculation() override;

	// convergence criteria
	// unbalanced force
	double f_ub;
	double init_f_ub;
	double f_ub_ratio;
	// maximum kinematic energy
	double e_kin;
	double e_kin_max;
	double e_kin_prev;
	double e_kin_ratio;

	double f_ub_ratio_bound;
	double e_kin_ratio_bound;

	// parallelism
	void cal_thread_func(unsigned int th_id);

	void spawn_all_threads();
	void join_all_threads();
	
	// global data
	unsigned int thread_num;
	std::atomic<bool> not_yet_completed;
	std::atomic<bool> terminate_this_substep;
	//std::atomic_flag not_yet_;
	ThreadBarrierFixedNum step_barrier;
	ThreadBarrierFixedNum cal_barrier;
	std::vector<std::thread> cal_threads;
	// for debug
	std::mutex cout_mutex;

	// actual parallelized calculation
	std::atomic<size_t> cur_pcl_id;
	std::atomic<size_t> cur_elem_id;
	std::atomic<size_t> cur_node_id;
	std::atomic<size_t> cur_ax_bc_id;
	std::atomic<size_t> cur_ay_bc_id;
	std::atomic<size_t> cur_vx_bc_id;
	std::atomic<size_t> cur_vy_bc_id;

	void map_pcl_vars_to_nodes_at_elems(unsigned int th_id);
	void update_node_a_and_v(unsigned int th_id);
	void apply_a_and_v_bcs(unsigned int th_id);
	//void cal_de_at_elem(unsigned int th_id);
	//void map_de_vol_from_elem_to_node(unsigned int th_id);
	void update_pcl_vars(unsigned int th_id);
};

#endif