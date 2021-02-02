#include "SimulationsOMP_pcp.h"

#include "Step_TBB.h"

Step_TBB::Step_TBB(
	const char* _name,
	const char* _type,
	CalSubstepFunc _cal_substep_func
	) : Step(_name, _type, _cal_substep_func),
	thread_num(1) {}

Step_TBB::~Step_TBB() {}

int Step_TBB::solve()
{
	substep_index = 0;
	current_time = 0.0;

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
	cpu_time = std::chrono::microseconds::zero();
	do
	{
		do
		{
			t0 = std::chrono::high_resolution_clock::now();
			cal_substep();
			t1 = std::chrono::high_resolution_clock::now();
			cpu_time += std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0);
		} while (output_not_needed);

		current_time_plus_tol = current_time + time_tol;
		output_time_history();
		next_output_time_minus_tol = next_output_time - time_tol;
		step_not_end = current_time < step_time_minus_tol;

	} while (step_not_end && continue_cal);

	finalize_calculation();
	finalize_time_history();
	return 0;
}

void Step_TBB::continue_calculation()
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

void Step_TBB::exit_calculation()
{
	++substep_index;
	current_time = new_time;
	output_not_needed = false;
	continue_cal = false;
}

void Step_TBB::abort_calculation()
{
	output_not_needed = false;
	continue_cal = false;
}

int cal_substep_func_TBB(
	void* _self,
	size_t my_th_id,
	double dt,
	double cur_time,
	size_t substp_id)
{
	((Step_TBB*)_self)->continue_calculation();
	return 0;
}
