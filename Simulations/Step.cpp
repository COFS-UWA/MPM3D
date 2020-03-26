#include "SimulationCore_pcp.h"

#include "Step.h"

Step::Step(SolveSubstepFunc solve_substep_func, const char *_type) :
	solve_substep(solve_substep_func), type(_type), name(20, '\0'),
	model(nullptr),
	is_first_step(true), start_substep_index(0), substep_num(0),
	step_time(0.0), start_time(0.0), current_time(0.0),
	dtime(0.0), time_tol_ratio(0.01), time_tol(0.0),
	solve_func(&Step::solve_th_only),
	time_history_top(nullptr), model_data_top(nullptr) {}

Step::~Step() {}

int Step::solve_th_only(void)
{
	substep_num = 0;
	current_time = 0.0;
	
	// initialize calculation
	init_calculation();

	// init time history
	for (TimeHistoryOutput *pth = time_history_top; pth; pth = pth->next)
	{
		pth->interval_time = step_time / double(pth->interval_num);
		pth->init_per_step();
		if (pth->need_output_init_state) pth->output();
		pth->next_time = pth->interval_time;
	}

	double step_time_minus_tol = step_time - time_tol;
	double new_current_time;
	do
	{
		new_current_time = current_time + dtime;
		if (new_current_time > step_time)
			dtime -= new_current_time - step_time;

		(*solve_substep)(this);
		++substep_num;
		current_time += dtime;

		output_time_history();

	} while (current_time < step_time_minus_tol);

	// finalize time history
	for (TimeHistoryOutput *pth = time_history_top; pth; pth = pth->next)
		pth->finalize_per_step();

	// finalize calculation
	finalize_calculation();

	return 0;
}

int Step::solve_th_and_md(void)
{
	double next_md_time, next_th_time;
	double new_current_time;
	double dtime_ori = dtime;
	double next_th_time_minus_tol;
	double next_md_time_minus_tol;
	double step_time_minus_tol = step_time - time_tol;

	// initialize calculation
	init_calculation();

	// init next_md_time
	next_md_time = step_time;
	for (ModelDataOutput *pmd = model_data_top; pmd; pmd = pmd->next)
	{
		// get the smallest pmd->time as next_md_time
		if (next_md_time > pmd->current_time)
			next_md_time = pmd->current_time;
	}
	next_md_time_minus_tol = next_md_time - time_tol;

	// init next_th_time
	next_th_time = next_md_time;
	for (TimeHistoryOutput *pth = time_history_top; pth; pth = pth->next)
	{
		// init time history output
		pth->interval_time = step_time / double(pth->interval_num);
		pth->init_per_step();
		pth->next_time = pth->interval_time;
		if (pth->need_output_init_state) pth->output();
		// get the smallest pth->next_time as next_th_time
		if (next_th_time > pth->next_time)
			next_th_time = pth->next_time;
	}
	next_th_time_minus_tol = next_th_time - time_tol;

	current_time = 0.0;
	substep_num = 0;
cal_loop:
	new_current_time = current_time + dtime;
	// see if need to output time history
	if (new_current_time < next_th_time_minus_tol)
	{
		(*solve_substep)(this);
		current_time = new_current_time;
		++substep_num;
		goto cal_loop;
	}
	else if (new_current_time > next_th_time)
	{
		dtime -= new_current_time - next_th_time;
		(*solve_substep)(this);
		++substep_num;
		current_time = next_th_time;
		// recover dtime value
		dtime = dtime_ori;
	}
	else // new_current_time == next_th_time
	{
		(*solve_substep)(this);
		++substep_num;
		current_time = new_current_time;
	}
	// reinit next_th_time
	next_th_time = step_time;
	// output if time history necessary
	for (TimeHistoryOutput *pth = time_history_top; pth; pth = pth->next)
	{
		new_current_time = current_time + time_tol;
		if (pth->next_time <= new_current_time)
		{
			pth->output();
			pth->next_time += pth->interval_time;
			if (pth->next_time < current_time)
				pth->next_time = current_time + pth->interval_time;
		}
		// get the smallest value for next_th_time
		if (next_th_time > pth->next_time)
			next_th_time = pth->next_time;
	}
	next_th_time_minus_tol = next_th_time - time_tol;
	// see if need to output model data
	if (current_time < next_md_time_minus_tol)
		goto cal_loop;
	// reinit next_md_time
	next_md_time = step_time;
	// output model data if necessary
	ModelDataOutput *pmd = model_data_top;
	ModelDataOutput *pmd_prev = nullptr;
	while (pmd)
	{
		if (pmd->current_time - time_tol <= current_time)
		{
			pmd->current_time = get_current_time();
			pmd->total_time = get_total_time();
			pmd->output();
			// remove this model output from list
			pmd = pmd->next;
			if (pmd_prev)
				pmd_prev->next = pmd;
			else
				model_data_top = pmd;
			continue;
		}
		// get the smallest pmd->time as next_md_time
		if (next_md_time > pmd->current_time)
			next_md_time = pmd->current_time;
		pmd_prev = pmd;
		pmd = pmd->next;
	}
	next_md_time_minus_tol = next_md_time - time_tol;
	if (next_th_time > next_md_time)
	{
		next_th_time = next_md_time;
		next_th_time_minus_tol = next_th_time - time_tol;
	}
	// see if this step has completed
	if (current_time < step_time_minus_tol)
		goto cal_loop;

	// finalize time history
	for (TimeHistoryOutput *pth = time_history_top; pth; pth = pth->next)
		pth->finalize_per_step();

	// finalize calculation
	finalize_calculation();

	return 0;
}

int solve_substep_base(void *_self) { return 0; }

void Step::output_all_time_history(void)
{
	for (TimeHistoryOutput *pth = time_history_top; pth; pth = pth->next)
		pth->output();
}

void Step::output_time_history(void)
{
	for (TimeHistoryOutput *pth = time_history_top; pth; pth = pth->next)
	{
		if (pth->next_time - time_tol <= current_time)
		{
			pth->output();
			pth->next_time += pth->interval_time;
			if (pth->next_time < current_time)
				pth->next_time = current_time + pth->interval_time;
		}
	}
}
