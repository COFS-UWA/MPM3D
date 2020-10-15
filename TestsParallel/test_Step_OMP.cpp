#include "TestsParallel_pcp.h"

#include "Step_OMP.h"
#include "test_simulations_omp.h"

void test_Step_OMP(int argc, char** argv)
{
	Step_OMP step("test");
	step.set_thread_num(4);
	step.set_step_time(1.0);
	step.set_dtime(0.1);
	step.solve();
}
