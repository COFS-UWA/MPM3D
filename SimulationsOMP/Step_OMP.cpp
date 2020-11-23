#include "SimulationsOMP_pcp.h"

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

	step_time_minus_tol = step_time - time_tol;

	init_calculation();
	init_time_history();
	next_output_time_minus_tol = next_output_time - time_tol;

	// new time should <= step time
	new_time = current_time + dtime;
	if (new_time > step_time)
	{
		dtime -= new_time - step_time;
		new_time = step_time;
	}

	continue_cal = true;
	cpu_time = 0;
	std::clock_t t0, t1;
#pragma omp parallel
	{
		size_t my_th_id = size_t(omp_get_thread_num());
		do
		{

			do
			{
#pragma omp master
				{
					t0 = std::clock();
				}

				// a barrier is needed at the end of cal_substep_func_omp
				(*cal_substep_func_omp)(this, my_th_id, dtime, current_time, substep_index);			
			
#pragma omp master
				{
					t1 = std::clock();
					cpu_time += t1 - t0;
				}
			} while (output_not_needed);

#pragma omp master
			{
				current_time_plus_tol = current_time + time_tol;
				output_time_history();
				next_output_time_minus_tol = next_output_time - time_tol;
				step_not_end = current_time < step_time_minus_tol;
			}
#pragma omp barrier
		
		} while (step_not_end && continue_cal);
	}

	finalize_calculation();
	finalize_time_history();
	return 0;
}

void Step_OMP::continue_calculation()
{
	++substep_index;
	current_time = new_time;
	output_not_needed = current_time < next_output_time_minus_tol;

	new_time = current_time + dtime;
	if (new_time > step_time)
	{
		dtime -= new_time - step_time;
		new_time = step_time;
	}
}

void Step_OMP::exit_calculation()
{
	++substep_index;
	current_time = new_time;
	output_not_needed = false;
	continue_cal = false;
}

void Step_OMP::abort_calculation()
{
	output_not_needed = false;
	continue_cal = false;
}

int substep_func_omp_default(
	void* _self,
	size_t my_th_id,
	double dt,
	double cur_time,
	size_t substp_id
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
