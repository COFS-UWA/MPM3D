#include "SimulationsOMP_pcp.h"

#include <fstream>

#include "tbb/task_scheduler_init.h"
#include "tbb/task.h"

#include "Step_T2D_ME_Task.h"
#include "Step_T2D_ME_TBB.h"

#define one_third (1.0/3.0)
#define N_min (1.0e-8)
#define Block_Low(th_id, th_num, data_num) ((th_id)*(data_num)/(th_num))

#ifdef _DEBUG
static std::fstream res_file_t2d_me_tbb;
#endif

Step_T2D_ME_TBB::Step_T2D_ME_TBB(const char* _name) :
	Step_TBB(_name, "Step_T2D_ME_TBB", &cal_substep_func_T2D_ME_TBB),
	task_data(*this) {}

Step_T2D_ME_TBB::~Step_T2D_ME_TBB() {}

int Step_T2D_ME_TBB::init_calculation()
{
#ifdef _DEBUG
	res_file_t2d_me_tbb.open("step_t2d_me_tbb.txt", std::ios::out | std::ios::binary);
#endif

	Model_T2D_ME_mt& md = *(Model_T2D_ME_mt *)model;
	char *cur_mem = (char *)valid_elem_mem.alloc(md.elem_num * 2);
	valid_elems = (size_t *)cur_mem;
	cur_mem += sizeof(size_t) * md.elem_num;
	tmp_valid_elems = (size_t *)cur_mem;

	tbb::task_scheduler_init init(thread_num);

	return 0;
}

int Step_T2D_ME_TBB::finalize_calculation()
{
	return 0;
}

int cal_substep_func_T2D_ME_TBB(void* _self)
{
	Step_T2D_ME_TBB& stp = *(Step_T2D_ME_TBB*)(_self);
	Step_T2D_ME_Task::TaskData& td = stp.task_data;
	tbb::task::spawn_root_and_wait(
		*new(tbb::task::allocate_root())
			Step_T2D_ME_Task::Step_T2D_ME_Task(td));
	if (td.valid_pcl_num == 0)
		stp.exit_calculation();
	stp.continue_calculation();
	td.sorted_pcl_var_id ^= 0;
	return 0;
}
