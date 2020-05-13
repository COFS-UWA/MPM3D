#include "Simulations_pcp.h"

#include "Step.h"

int cal_substep_base(void *_self)
{
	return 0;
}

Step::Step(const char *_name, const char *_type, CalSubstepFunc _cal_substep_func) :
	name(_name), type(_type),
	cal_substep_func(_cal_substep_func),
	is_first_step(true), prev_substep_num(0), start_time(0.0),
	time_tol_ratio(1.0e-3) {}

Step::~Step() {}

int Step::solve()
{
	substep_index = 0;
	current_time = 0.0;

	init_calculation();
	init_time_history();

	double new_current_time;
	do
	{
		do
		{
			// adjust dtime so that new current time won't > step time
			new_current_time = current_time + dtime;
			if (new_current_time > step_time)
			{
				dtime -= new_current_time - step_time;
				new_current_time = step_time;
			}

			cal_substep();

			++substep_index;
			current_time = new_current_time;
			
			current_time_plus_tol = current_time + time_tol;

		} while (current_time_plus_tol < next_output_time);

		output_time_history();

	} while (current_time_plus_tol < step_time);

	finalize_time_history();
	finalize_calculation();

	return 0;
}

int solve_substep_base(void *_self) { return 0; }

int Step::init_time_history()
{
	if (time_history_list.is_empty())
	{
		next_output_time = step_time;
		return 0;
	}

	for (TimeHistory *pth = time_history_list.first();
		 time_history_list.is_not_end(pth);
		 pth = time_history_list.next(pth))
	{
		pth->interval_time = step_time / double(pth->interval_num);
		pth->init_per_step();
		if (pth->need_output_init_state)
			pth->output();
		pth->next_time = pth->interval_time;
	}

	// sort time history according to "next_time"
	TimeHistoryList th_list_tmp;
	TimeHistory *pth_tmp;
	while ((pth_tmp = get_latest_time_history(time_history_list)) != nullptr)
	{
		th_list_tmp.push(pth_tmp);
	}
	time_history_list.transfer(th_list_tmp);

	// update next output time
	next_output_time = time_history_list.first()->next_time;
	if (next_output_time > step_time)
		next_output_time = step_time;

	return 0;
}

int Step::output_time_history()
{
	if (time_history_list.is_empty())
	{
		next_output_time = step_time;
		return 0;
	}

	TimeHistory *pth = time_history_list.first();
	while (pth->next_time < current_time_plus_tol)
	{
		pth->output();
		pth->next_time += pth->interval_time;
		if (pth->next_time < current_time)
			pth->next_time = current_time + pth->interval_time;
		
		pth = time_history_list.pop();
		insert_time_history_in_ascending_order(time_history_list, pth);

		pth = time_history_list.first();
	}

	// update next output time
	next_output_time = time_history_list.first()->next_time;
	if (next_output_time > step_time)
		next_output_time = step_time;

	return 0;
}

int Step::finalize_time_history()
{
	for (TimeHistory *pth = time_history_list.first();
		time_history_list.is_not_end(pth);
		pth = pth = time_history_list.next(pth))
		pth->finalize_per_step();
	return 0;
}

TimeHistory *Step::get_latest_time_history(TimeHistoryList &th_list)
{
	if (th_list.is_empty())
		return nullptr;

	TimeHistory *th = th_list.first();
	TimeHistory *latest_th = th;
	for (th = th_list.next(th); th_list.is_not_end(th); th = th_list.next(th))
	{
		if (latest_th->next_time < th->next_time)
			latest_th = th;
	}
	th_list.del(latest_th);
	return latest_th;
}

void Step::insert_time_history_in_ascending_order(TimeHistoryList &th_list, TimeHistory *th)
{
	TimeHistory *pth = th_list.first();
	while (true)
	{
		if (th_list.is_end(pth))
		{
			th_list.append(th);
			break;
		}

		if (th->next_time < pth->next_time)
		{
			th_list.insert_before(pth, th);
			break;
		}

		pth = th_list.next(pth);
	}
}
