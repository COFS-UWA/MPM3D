#ifndef __Step_T2D_ME_p_h__
#define __Step_T2D_ME_p_h__

#include <vector>
#include <thread>
#include <mutex>

#include "Step.h"
#include "Model_T2D_ME_p.h"
#include "ThreadBarrierFixedNum.h"

int solve_substep_T2D_ME_p(void* _self);
int solve_substep_T2D_ME_p_RigidCircle(void* _self);

// parallelism version
class Step_T2D_ME_p : public Step
{
public:
	typedef Model_T2D_ME_p_Internal::NodeToElem NodeToElem;
	typedef Model_T2D_ME_p::Node Node;
	typedef Model_T2D_ME_p_Internal::NodeVarAtElem NodeVarAtElem;
	typedef Model_T2D_ME_p::Element Element;
	typedef Model_T2D_ME_p::Particle Particle;
	
public:
	Step_T2D_ME_p(const char* _name);
	~Step_T2D_ME_p();

	inline void set_model(Model_T2D_ME_p& md)
	{
		Step::set_model(md);
		model = &md;
	}

	inline void set_prev_step(Step_T2D_ME_p& prev_step)
	{
		Step::set_prev_step(prev_step);
		model = prev_step.model;
	}

	inline void set_damping_ratio(double _ratio) noexcept { damping_ratio = _ratio; }
	
	inline void set_thread_num(unsigned int num) { thread_num = num; }

protected:
	Model_T2D_ME_p* model;
	double damping_ratio;

	int init_calculation() override;
	friend int solve_substep_T2D_ME_p(void* _self);
	friend int solve_substep_T2D_ME_p_RigidCircle(void* _self);
	int finalize_calculation() override;

	// helper structures and functions for parallelism
	class ThreadData;
	void cal_thread_func(unsigned int th_id);
	void cal_thread_func_RigidCircle(unsigned int th_id);

	void spawn_all_threads();
	void join_all_threads();

	// global data
	unsigned int thread_num;
	bool not_yet_completed;
	ThreadBarrierFixedNum step_barrier;
	ThreadBarrierFixedNum cal_barrier;
	std::vector<std::thread> cal_threads;
	// for debug
	std::mutex cout_mutex;

	// thread-wise data
	struct ThreadElemData
	{
		friend Step_T2D_ME_p;
	protected:
#ifdef _DEBUG
		size_t elem_id;
#endif
		// particle stack
		Particle *top_pcl;
		Particle *last_pcl;

	public:
		inline Particle* get_top_pcl() { return top_pcl; }
		inline void init() { top_pcl = nullptr; }
		inline void add_pcl(Particle &pcl)
		{
			if (!top_pcl) // insert the first pcl
				last_pcl = &pcl;
			pcl.next_pcl_in_elem = top_pcl;
			top_pcl = &pcl;
		}
		inline void combine(ThreadElemData& other)
		{
			if (other.top_pcl)
			{
				if (top_pcl)
				{
					other.last_pcl->next_pcl_in_elem = top_pcl;
					top_pcl = other.top_pcl;
				}
				else
				{
					top_pcl = other.top_pcl;
					last_pcl = other.last_pcl;
				}
			}
		}
	};

	struct ThreadData
	{
		RigidCircleForce rcf;
		std::atomic_flag not_combined_yet;
	};
	
	char *thread_data_mem; // main memory
	ThreadData **pth_datas;
	ThreadElemData **pth_elem_datas;

	void alloc_thread_data();
	void clear_thread_data();

	// actual parallelized calculation
	std::atomic<size_t> cur_pcl_id;
	std::atomic<size_t> cur_elem_id;
	std::atomic<size_t> cur_node_id;
	std::atomic<size_t> cur_ax_bc_id;
	std::atomic<size_t> cur_ay_bc_id;
	std::atomic<size_t> cur_vx_bc_id;
	std::atomic<size_t> cur_vy_bc_id;

	void init_cal_vars(unsigned int th_id);
	void find_pcls_in_which_elems(unsigned int th_id);
	void map_pcl_vars_to_nodes_at_elems(unsigned int th_id);
	void cal_contact_force(unsigned int th_id); // only for rigid body
	void update_node_a_and_v(unsigned int th_id);
	void apply_a_and_v_bcs(unsigned int th_id);
	void cal_de_at_elem(unsigned int th_id);
	void map_de_vol_from_elem_to_node(unsigned int th_id);
	void update_pcl_vars(unsigned int th_id);
};

#endif