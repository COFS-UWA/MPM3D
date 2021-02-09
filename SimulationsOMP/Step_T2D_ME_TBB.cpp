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
	task_data(
		1, // pcl_num_per_map_pcl_to_mesh_task
		1, // node_num_per_update_a_and_v_task
		1, // elem_num_per_cal_elem_de_task
		1, // node_num_per_cal_node_de_task
		1  // pcl_num_per_task_map_mesh_to_pcl
		) {}

Step_T2D_ME_TBB::~Step_T2D_ME_TBB() {}

int Step_T2D_ME_TBB::init_calculation()
{
#ifdef _DEBUG
	res_file_t2d_me_tbb.open("step_t2d_me_tbb.txt", std::ios::out | std::ios::binary);
#endif

	tbb::task_scheduler_init init(thread_num);

	Model_T2D_ME_mt& md = *(Model_T2D_ME_mt*)model;
	task_data.set_model(md);
	task_data.sorted_pcl_var_id = 1;
	auto &pcl_sort_mem = task_data.pcl_sort_mem;
	auto& node_sort_mem = task_data.node_sort_mem;
	const size_t psm_sm_size
		= pcl_sort_mem.get_shared_memory_size(
			thread_num,
			md.pcl_num,
			md.ori_pcl_num);
	const size_t nsm_sm_size
		= node_sort_mem.get_shared_memory_size(
			thread_num,
			md.elem_num,
			md.node_num);
	const size_t sm_size
		= psm_sm_size > nsm_sm_size ?
			psm_sm_size : nsm_sm_size;
	cal_mem.alloc(sm_size);
	pcl_sort_mem.init(
		cal_mem.raw_address(),
		thread_num,
		md.pcl_num,
		md.ori_pcl_num);
	pcl_sort_mem.valid_pcl_num = md.pcl_num;
	node_sort_mem.init(
		cal_mem.raw_address(),
		thread_num,
		md.elem_num,
		md.node_num);

	tbb::task::spawn_root_and_wait(
		*new(tbb::task::allocate_root())
			Step_T2D_ME_Task::InitTask(task_data));

	return 0;
}

int Step_T2D_ME_TBB::finalize_calculation()
{
	Model_T2D_ME_mt& md = *(Model_T2D_ME_mt*)model;
	md.pcl_num = task_data.pcl_sort_mem.valid_pcl_num;
	return 0;
}

int cal_substep_func_T2D_ME_TBB(void* _self)
{
	Step_T2D_ME_TBB& stp = *(Step_T2D_ME_TBB*)(_self);
	Step_T2D_ME_Task::TaskData& td = stp.task_data;
	td.dt = stp.dtime;
	tbb::task::spawn_root_and_wait(
		*new(tbb::task::allocate_root())
			Step_T2D_ME_Task::Step_T2D_ME_Task(td));
	if (td.pcl_sort_mem.valid_pcl_num == 0)
	{
		stp.exit_calculation();
		return 0;
	}
	stp.continue_calculation();
	td.sorted_pcl_var_id ^= 1;
	return 0;
}
