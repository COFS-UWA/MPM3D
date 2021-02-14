#include "SimulationsOMP_pcp.h"

#include <fstream>

#include "tbb/task_scheduler_init.h"
#include "tbb/task.h"

#include "DivideTask.hpp"
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
	cal_data(
		1, // pcl_num_per_init_pcl_task
		1, // pcl_num_per_map_pcl_to_mesh_task
		1, // node_num_per_update_a_and_v_task
		1, // elem_num_per_cal_elem_de_task
		1, // node_num_per_cal_node_de_task
		1  // pcl_num_per_task_map_mesh_to_pcl
		),
	init_pcl(cal_data),
	map_pcl_to_mesh(cal_data),
	update_a_and_v(cal_data),
	cal_elem_de(cal_data),
	cal_node_de(cal_data),
	map_mesh_to_pcl(cal_data) {}

Step_T2D_ME_TBB::~Step_T2D_ME_TBB() {}

int Step_T2D_ME_TBB::init_calculation()
{
#ifdef _DEBUG
	res_file_t2d_me_tbb.open("step_t2d_me_tbb.txt", std::ios::out | std::ios::binary);
#endif

	Model_T2D_ME_mt& md = *(Model_T2D_ME_mt*)model;
	if (md.pcl_num == 0)
		return -1;

	cal_data.set_model(md);
	cal_data.thread_num = thread_num;
	cal_data.sorted_pcl_var_id = 0;

	auto &pcl_sort_mem = cal_data.pcl_sort_mem;
	const size_t psm_sm_size
		= pcl_sort_mem.get_shared_memory_size(
			thread_num,
			md.pcl_num,
			md.ori_pcl_num);
	auto& node_sort_mem = cal_data.node_sort_mem;
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
	node_sort_mem.init(
		cal_mem.raw_address(),
		thread_num,
		md.elem_num,
		md.node_num);
	
	tbb::task_scheduler_init init(thread_num);
	cal_data.init_pcl_task_num = (md.pcl_num
		+ cal_data.pcl_num_per_map_pcl_to_mesh_task - 1)
		/ cal_data.pcl_num_per_map_pcl_to_mesh_task;
	tbb::task::spawn_root_and_wait(
		*new(tbb::task::allocate_root())
			DivideTask<Step_T2D_ME_Task::InitPcl, 8>(
				0, cal_data.init_pcl_task_num, init_pcl));

	return 0;
}

int Step_T2D_ME_TBB::finalize_calculation()
{
	Model_T2D_ME_mt& md = *(Model_T2D_ME_mt*)model;
	md.pcl_num = cal_data.pcl_sort_mem.valid_pcl_num;
	return 0;
}

int cal_substep_func_T2D_ME_TBB(void* _self)
{
	Step_T2D_ME_TBB& self = *(Step_T2D_ME_TBB*)(_self);
	Step_T2D_ME_Task::CalData &cd = self.cal_data;
	cd.sorted_pcl_var_id ^= 1;
	cd.dt = self.dtime;

	// sort pcl id
	SortUtils::SortParticleMem& pcl_sort_mem = cd.pcl_sort_mem;
#ifdef _DEBUG
	cd.prev_valid_pcl_num = pcl_sort_mem.valid_pcl_num;
#endif
	tbb::task::spawn_root_and_wait(
		*new(tbb::task::allocate_root())
			SortUtils::SortParticleTask(pcl_sort_mem));
	
	const size_t valid_pcl_num = pcl_sort_mem.valid_pcl_num;
	if (valid_pcl_num == 0)
	{
		self.exit_calculation();
		return 0;
	}
	cd.map_pcl_to_mesh_task_num = (valid_pcl_num
		+ cd.pcl_num_per_map_pcl_to_mesh_task - 1)
		/ cd.pcl_num_per_map_pcl_to_mesh_task;
	cd.map_mesh_to_pcl_task_num = (valid_pcl_num
		+ cd.pcl_num_per_map_mesh_to_pcl_task - 1)
		/ cd.pcl_num_per_map_mesh_to_pcl_task;

	tbb::task_list tk_list;
	// sort node
	SortUtils::SortTriMeshNodeMem& node_sort_mem = cd.node_sort_mem;
	tk_list.push_back(*new(tbb::task::allocate_root())
		SortUtils::SortTriMeshNodeTask(
#ifdef _DEBUG
			cd.elem_num,
#endif
			pcl_sort_mem.pcl_in_elem,
			valid_pcl_num,
			cd.elem_node_id,
			node_sort_mem));
	// map pcl to bg mesh
	tk_list.push_back(*new(tbb::task::allocate_root())
		DivideTask<Step_T2D_ME_Task::MapPclToBgMesh, 8>(
			0, cd.map_pcl_to_mesh_task_num, self.map_pcl_to_mesh));
	tbb::task::spawn_root_and_wait(tk_list);
	
	const size_t valid_elem_num = node_sort_mem.valid_elem_num;
	const size_t three_valid_elem_num = valid_elem_num * 3;
	cd.update_a_and_v_task_num = (three_valid_elem_num
		+ cd.elem_num_per_update_a_and_v_task - 1)
		/ cd.elem_num_per_update_a_and_v_task;
	cd.cal_elem_de_task_num = (valid_elem_num
		+ cd.elem_num_per_cal_elem_de_task - 1)
		/ cd.elem_num_per_cal_elem_de_task;
	cd.cal_node_de_task_num = (three_valid_elem_num
		+ cd.elem_num_per_cal_node_de_task - 1)
		/ cd.elem_num_per_cal_node_de_task;

	// update nodal a and v
	tbb::task::spawn_root_and_wait(
		*new(tbb::task::allocate_root())
			DivideTask<Step_T2D_ME_Task::UpdateAccelerationAndVelocity, 8>(
				0, cd.update_a_and_v_task_num, self.update_a_and_v));

	// cal element de and map to node
	tbb::task::spawn_root_and_wait(
		*new(tbb::task::allocate_root())
			DivideTask<Step_T2D_ME_Task::CalElemDeAndMapToNode, 8>(
				0, cd.cal_elem_de_task_num, self.cal_elem_de));

	// cal strain increment at node
	tbb::task::spawn_root_and_wait(
		*new(tbb::task::allocate_root())
			DivideTask<Step_T2D_ME_Task::CalNodeDe, 8>(
				0, cd.cal_node_de_task_num, self.cal_node_de));

	// map bg mesh back to pcl
	tbb::task::spawn_root_and_wait(
		*new(tbb::task::allocate_root())
			DivideTask<Step_T2D_ME_Task::MapBgMeshToPcl, 8>(
				0, cd.map_mesh_to_pcl_task_num, self.map_mesh_to_pcl));

	self.continue_calculation();
	return 0;
}
