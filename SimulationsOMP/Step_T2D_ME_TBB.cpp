#include "SimulationsOMP_pcp.h"

#include <fstream>

#include "Step_T2D_ME_TBB_task.h"
#include "Step_T2D_ME_TBB.h"

#define one_third (1.0/3.0)
#define N_min (1.0e-8)
#define Block_Low(th_id, th_num, data_num) ((th_id)*(data_num)/(th_num))

#ifdef _DEBUG
static std::fstream res_file_t2d_me_tbb;
#endif

Step_T2D_ME_TBB::Step_T2D_ME_TBB(const char* _name) :
	Step_TBB(_name, "Step_T2D_ME_TBB", &cal_substep_func_T2D_ME_TBB) {}

Step_T2D_ME_TBB::~Step_T2D_ME_TBB() {}

int Step_T2D_ME_TBB::init_calculation()
{
#ifdef _DEBUG
	res_file_t2d_me_tbb.open("step_t2d_me_tbb.txt", std::ios::out | std::ios::binary);
#endif

	return 0;
}

int Step_T2D_ME_TBB::finalize_calculation()
{
	return 0;
}

int cal_substep_func_T2D_ME_TBB(void* _self)
{
	return 0;
}
