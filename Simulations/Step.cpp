#include "Simulations_pcp.h"

#include "Step.h"

Step::Step(const char *_name, const char *_type, CalSubstepFunc _cal_substep_func) :
	name(_name), type(_type),
	is_first_step(true), prev_substep_num(0), start_time(0.0),
	time_tol_ratio(1.0e-3) {}

Step::~Step() {}

int Step::solve()
{
	substep_index = 0;
	current_time = 0.0;

	init_calculation();

	init_time_history();

	double step_time_minus_tol = step_time - time_tol;
	double new_current_time;
	do
	{
		new_current_time = current_time + dtime;
		if (new_current_time > step_time)
			dtime -= new_current_time - step_time;

		cal_substep();
		++substep_index;
		current_time += dtime;

		output_time_history();

	} while (current_time < step_time_minus_tol);

	finalize_time_history();

	finalize_calculation();

	return 0;
}

int solve_substep_base(void *_self) { return 0; }

int Step::init_calculation()
{
	TimeHistoryList th_list_tmp;
	TimeHistory *pth_tmp;
	for (TimeHistory *pth = time_history_list.first(); time_history_list.is_not_end(pth);)
	{
		pth->interval_time = step_time / double(pth->interval_num);
		pth->init_per_step();
		if (pth->need_output_init_state)
			pth->output();
		pth->next_time = pth->interval_time;
		// move to another list
		pth_tmp = pth;
		pth = time_history_list.next(pth);
		time_history_list.del(pth_tmp);
		th_list_tmp.append(pth_tmp);
	}

	// reorder time history according to next_time
	for (TimeHistory *pth = th_list_tmp.first(); th_list_tmp.is_not_end(pth);)
	{
		pth_tmp = pth;
		pth = th_list_tmp.next(pth);
		th_list_tmp.del(pth_tmp);

		TimeHistory *pth = ;
		for (size_t i = 0; i < length; i++)
		{

		}
		time_history_list.insert_before(, pth_tmp);
	}

	return 0;
}

void Step::output_time_history()
{
	double cur_time_tmp = current_time + time_tol;
	for (TimeHistory *pth = time_history_list.first();
		time_history_list.is_not_end(pth);
		pth = time_history_list.next(pth))
	{
		if (cur_time_tmp < pth->next_time)
		{
			break;
		}
		else
		{
			pth->output();
			pth->next_time += pth->interval_time;
			if (pth->next_time < current_time)
				pth->next_time = current_time + pth->interval_time;
			if (next_output_time > pth->next_time)
				next_output_time = pth->next_time;
		}
	}
}

int Step::finalize_calculation()
{
	for (TimeHistory *pth = time_history_list.first();
		time_history_list.is_not_end(pth);
		pth = pth = time_history_list.next(pth))
		pth->finalize_per_step();
	return 0;
}
