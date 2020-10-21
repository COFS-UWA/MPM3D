#ifndef __Step_OMP_h__
#define __Step_OMP_h__

#include "Step.h"

typedef int (*CalSubstepFuncOMP)(void* _self, uint32_t my_th_id,
	float dt, float cur_time, uint32_t substp_id);

int substep_func_omp_default(void *_self, uint32_t my_th_id, 
	float dt, float cur_time, uint32_t substp_id);

class Step_OMP : public Step
{
protected:
	uint32_t thread_num;
	float dt;
	float time_tol_f;
	float stp_time;
	float stp_time_minus_tol;
	float next_output_time_minus_tol;
	uint32_t substp_id;
	float cur_time, new_time;
	bool output_not_needed, step_not_end;
	bool continue_cal;

	CalSubstepFuncOMP cal_substep_func_omp;

public:
	Step_OMP(const char *_name,	const char *_type = "Step_OMP",
		CalSubstepFuncOMP _func_omp = &substep_func_omp_default);
	~Step_OMP();

	inline void set_thread_num(uint32_t th_num) noexcept { thread_num = th_num; }

	int solve() override;

	// this function need to be put into
	// #pragma omp master
	void continue_calculation();
	void exit_calculation();
};

#endif