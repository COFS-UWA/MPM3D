#ifndef __Step_T2D_ME_p_h__
#define __Step_T2D_ME_p_h__

#include <vector>
#include <thread>
#include <mutex>

#include "Step.h"
#include "Model_T2D_ME_p.h"
#include "ThreadBarrierFixedNum.h"

int solve_substep_T2D_ME_p(void* _self);

// parallelism version
class Step_T2D_ME_p : public Step
{
public:
	typedef Model_T2D_ME_p_Internal::NodeToElem NodeToElem;
	typedef Model_T2D_ME_p::Node Node;
	typedef Model_T2D_ME_p_Internal::NodeVarAtElem NodeVarAtElem;
	typedef Model_T2D_ME_p::Element Element;
	typedef Model_T2D_ME_p::Particle Particle;
	
protected:
	Model_T2D_ME_p* model;

	int init_calculation() override;
	friend int solve_substep_T2D_ME_p(void* _self);
	int finalize_calculation() override;

	double damping_ratio;

public:
	Step_T2D_ME_p(const char* _name);
	~Step_T2D_ME_p();

	inline void set_model(Model_T2D_ME_p& md)
	{
		Step::set_model(md);
		model = &md;
	}

	// Restart from previous step
	inline void set_prev_step(Step_T2D_ME_p& prev_step)
	{
		Step::set_prev_step(prev_step);
		model = prev_step.model;
	}

	inline void set_damping_ratio(double _ratio) noexcept { damping_ratio = _ratio; }

protected: // helper structures and functions for parallelism
	unsigned int thread_num; // core number - 1
	std::vector<std::thread> cal_threads;
	ThreadBarrierFixedNum step_barrier;
	ThreadBarrierFixedNum cal_barrier;
	bool not_yet_completed;

	// terminate all other threads
	void join_thread_and_exit();
	// used by thread
	void cal_thread_func(unsigned int th_id);

	struct ThreadElemData
	{
		Particle *pcls;
		inline void init() { pcls = nullptr; }
		inline void add_pcl(Particle &pcl)
		{
			pcl.next_pcl_in_elem = pcls;
			pcls = &pcl;
		}
	};
	ThreadElemData *ted_mem;
	ThreadElemData **pteds;
	void alloc_thread_elem_data(size_t elem_num);
	void clear_thread_elem_data();

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

	void find_pcls_in_which_elems(unsigned int th_id);
	void map_pcl_vars_to_nodes_at_elems(unsigned int th_id);
	void update_node_a_and_v(unsigned int th_id);
	void apply_a_and_v_bcs(unsigned int th_id);
	void cal_de_at_elem(unsigned int th_id);
	void map_de_vol_from_elem_to_node(unsigned int th_id);
	void update_pcl_vars(unsigned int th_id);

public:
	inline void set_thread_num(unsigned int num) { thread_num = num; }
};

#endif