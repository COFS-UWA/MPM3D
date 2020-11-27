#include "TestsParallel_pcp.h"

#include "Step_OMP.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "test_simulations_omp.h"

namespace
{
	int substep_func_omp_barrier(void* _self, size_t my_th_id,
		double dt, double cur_time, size_t substp_id);

	class Step_Barrier : public Step_OMP
	{
		friend 	int substep_func_omp_barrier(void* _self, size_t my_th_id,
			double dt, double cur_time, size_t substp_id);
	protected:
		int sum;

	public:
		Step_Barrier(const char* _name) : sum(0),
			Step_OMP(_name, "Step_Barrier", &substep_func_omp_barrier) {}
		Step_Barrier::~Step_Barrier() {}
	};

	int substep_func_omp_barrier(void* _self, size_t my_th_id,
		double dt, double cur_time, size_t substp_id)
	{
		Step_Barrier& self = *(Step_Barrier*)_self;

#pragma omp barrier

#pragma omp critical
		++self.sum;

#pragma omp barrier
#pragma omp barrier
#pragma omp barrier

#pragma omp master
		{
			self.sum = 0;
			self.continue_calculation();
		}
#pragma omp barrier
		return 0;
	}
}

// cal the time spend in barrier
void test_omp_barrier_time()
{
	TimeHistory_ConsoleProgressBar out_cpb;

	Step_Barrier step("barrier");
	step.set_step_time(0.5);
	//step.set_step_time(1.0e-4);
	step.set_dtime(1.0e-5);
	step.set_thread_num(6);
	step.add_time_history(out_cpb);
	step.solve();

	std::cout << step.cpu_time_in_ms() << "\n";
	system("pause");
}

