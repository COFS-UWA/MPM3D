#ifndef __Step_TBB_h__
#define __Step_TBB_h__

#include <chrono>
#include "Step.h"

int cal_substep_func_TBB(void* _self);

class Step_TBB : public Step
{
protected:
	size_t thread_num;
	double new_time;
	double step_time_minus_tol;
	double next_output_time_minus_tol;
	bool output_not_needed, step_not_end;
	bool continue_cal;

	std::chrono::high_resolution_clock::time_point t0, t1;
	std::chrono::microseconds cpu_time;

public:
	Step_TBB(const char* _name, const char* _type = "Step_TBB",
		CalSubstepFunc _cal_substep_func = &cal_substep_func_TBB);
	~Step_TBB();

	inline void set_thread_num(size_t th_num) noexcept { thread_num = th_num; }

	int solve() override;

	// this function need to be put into
	// #pragma omp master
	void continue_calculation();
	void exit_calculation();
	void abort_calculation();

	// in microseconds
	inline long long get_time() const noexcept { return std::chrono::duration_cast<std::chrono::microseconds>(cpu_time).count(); }
};

#endif