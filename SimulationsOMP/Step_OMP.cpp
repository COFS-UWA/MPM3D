#include "SimulationsOMP_pcp.h"

#include <stdio.h>
#include <omp.h>

#include "Step_OMP.h"

Step_OMP::Step_OMP(
	const char* _name,
	const char* _type,
	CalSubstepFuncOMP _func_omp
	) : Step(_name, _type), thread_num(1),
	cal_substep_func_omp(_func_omp) {}

Step_OMP::~Step_OMP() {}

int Step_OMP::solve()
{
	substep_index = 0;
	current_time = 0.0;
	omp_set_num_threads(thread_num);

	dt = float(dtime);
	stp_time = float(step_time);
	time_tol_f = float(time_tol);
	stp_time_minus_tol = stp_time - time_tol_f;

	init_calculation();
	init_time_history();
	next_output_time_minus_tol = float(next_output_time) - time_tol_f;

	substp_id = uint32_t(substep_index);
	cur_time = float(current_time);

	// new time should <= step time
	new_time = cur_time + dt;
	if (new_time > stp_time)
	{
		dt -= new_time - stp_time;
		new_time = stp_time;
	}

	continue_cal = true;
#pragma omp parallel
	{
		uint32_t my_th_id = uint32_t(omp_get_thread_num());
		do
		{
			do
			{
				(*cal_substep_func_omp)(this, my_th_id, dt, cur_time, substp_id);

//#pragma omp barrier // need sync after each substep
//
//#pragma omp master
//				{
//					++substp_id;
//					cur_time = new_time;
//					output_not_needed = cur_time < next_output_time_minus_tol;
//
//					new_time = cur_time + dt;
//					if (new_time > stp_time)
//					{
//						dt -= new_time - stp_time;
//						new_time = stp_time;
//					}
//				}
//#pragma omp barrier
			
			} while (output_not_needed);

#pragma omp barrier
#pragma omp master
			{
				substep_index = uint32_t(substp_id);
				current_time = double(cur_time);
				current_time_plus_tol = double(cur_time) + time_tol;
				output_time_history();
				next_output_time_minus_tol = float(next_output_time) - time_tol_f;
				step_not_end = cur_time < stp_time_minus_tol;
			}
#pragma omp barrier
		
		} while (step_not_end & continue_cal);
	}

	finalize_calculation();
	finalize_time_history();
	return 0;
}

void Step_OMP::continue_calculation()
{
	++substp_id;
	cur_time = new_time;
	output_not_needed = cur_time < next_output_time_minus_tol;

	new_time = cur_time + dt;
	if (new_time > stp_time)
	{
		dt -= new_time - stp_time;
		new_time = stp_time;
	}
}

void Step_OMP::exit_calculation()
{
	++substp_id;
	cur_time = new_time;
	output_not_needed = false;
	continue_cal = false;
}

int substep_func_omp_default(
	void* _self,
	uint32_t my_th_id,
	float dt,
	float cur_time,
	uint32_t substp_id
	)
{
	Step_OMP& self = *(Step_OMP*)_self;

#pragma omp master
	{
		self.continue_calculation();
	}
#pragma omp barrier

	return 0;
}
