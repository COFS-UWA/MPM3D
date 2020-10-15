#include "SimulationsOMP_pcp.h"

#include <stdio.h>
#include <omp.h>

#include "Step_OMP.h"

Step_OMP::Step_OMP(
	const char* _name,
	const char* _type,
	CalSubstepFuncOMP _func_omp
	) : Step(_name, _type), thread_num(1),
	cal_substep_func_omp(_func_omp)
{

}

Step_OMP::~Step_OMP()
{

}

int Step_OMP::solve()
{
	substep_index = 0;
	current_time = 0.0;

	init_calculation();
	init_time_history();
	
	omp_set_num_threads(thread_num);

#pragma omp parallel
	{
		float time_tol = float(Step::time_tol);
		uint32_t my_th_id = uint32_t(omp_get_thread_num());
		uint32_t substp_id = uint32_t(substep_index);
		float cur_time = float(current_time);
		float dt = float(dtime);
		float stp_time = float(Step::step_time);
		float stp_time_minus_tol = stp_time - time_tol;
		float next_output_time_minus_tol = float(next_output_time) - time_tol;
		do
		{
			do
			{
				// adjust dtime so that new current time won't > step time
				float new_time = cur_time + dt;
				if (new_time > stp_time)
				{
					dt -= new_time - stp_time;
					new_time = stp_time;
				}

				(*cal_substep_func_omp)(this, my_th_id, dt, cur_time, substp_id);
#pragma omp barrier

				++substp_id;
				cur_time = new_time;

			} while (cur_time < next_output_time_minus_tol);

#pragma omp barrier
#pragma omp master
			output_time_history();
			next_output_time_minus_tol = float(next_output_time) - time_tol;
#pragma omp barrier

		} while (cur_time < stp_time_minus_tol);
	
#pragma omp master
		substep_index = uint32_t(substp_id);
		current_time = double(cur_time);
	}

	finalize_calculation();
	finalize_time_history();

	return 0;
}

int substep_func_omp_default(
	void* _self,
	uint32_t my_th_id,
	float dt,
	float cur_time,
	uint32_t substp_id
	)
{
	//printf("th_id %u, cur_time %f, stp id %u\n", my_th_id, cur_time, substp_id);
	return 0;
}
