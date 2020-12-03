#ifndef __Step_OMP_h__
#define __Step_OMP_h__

#include <chrono>
#include "Step.h"

typedef int (*CalSubstepFuncOMP)(void* _self, size_t my_th_id,
	double dt, double cur_time, size_t substp_id);

int substep_func_omp_default(void *_self, size_t my_th_id, 
	double dt, double cur_time, size_t substp_id);

class Step_OMP : public Step
{
protected:
	size_t thread_num;
	double new_time;
	double step_time_minus_tol;
	double next_output_time_minus_tol;
	bool output_not_needed, step_not_end;
	bool continue_cal;

	CalSubstepFuncOMP cal_substep_func_omp;
	std::chrono::steady_clock::time_point t0, t1;
	std::chrono::nanoseconds cpu_time;

public:
	Step_OMP(const char *_name,	const char *_type = "Step_OMP",
		CalSubstepFuncOMP _func_omp = &substep_func_omp_default);
	~Step_OMP();

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