#ifndef __Step_T2D_ME_p_tbb_h__
#define __Step_T2D_ME_p_tbb_h__

#include <tbb/task_scheduler_init.h>
#include <tbb/blocked_range.h>

#include "FixedSizeMemAllocator.h"
#include "Step_T2D_ME_p.h"

int solve_substep_T2D_ME_p_tbb(void* _self);

class Step_T2D_ME_p_tbb : public Step_T2D_ME_p
{
public:
	struct PclListAtElem
	{
	protected:
		Particle* top_pcl;
		Particle* last_pcl;

	public:
		inline Particle* get_top() const noexcept { return top_pcl; }
		inline void init() noexcept { top_pcl = nullptr; }
		inline void add_pcl(Particle& pcl)
		{
			if (!top_pcl) // the first pcl
				last_pcl = &pcl;
			pcl.next_pcl_in_elem = top_pcl;
			top_pcl = &pcl;
		}
		inline void combine(PclListAtElem& other)
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

protected:
	tbb::task_scheduler_init thread_pool_config;
	FixedSizeMemAllocator<PclListAtElem> plist_mem;

	double init_node_time;
	double find_pcl_in_elem_time;
	double map_pcl_to_node_time;
	double rigid_circle_time;
	double update_node_a_v_time;
	double apply_a_v_bc_time;
	double cal_strain_at_elem_time;
	double map_de_vol_from_elem_to_node_time;
	double update_pcl_time;
	double total_parallel_time;

public:
	explicit Step_T2D_ME_p_tbb(const char* _name);
	~Step_T2D_ME_p_tbb();

	Step_T2D_ME_p_tbb(const Step_T2D_ME_p_tbb& other) = delete;
	Step_T2D_ME_p_tbb &operator=(const Step_T2D_ME_p_tbb& other) = delete;

	int init_calculation() override;
	friend int solve_substep_T2D_ME_p_tbb(void* _self);
	int finalize_calculation() override;
};

#endif