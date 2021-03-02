#include "SimulationsOMP_pcp.h"

#include <fstream>

#include "tbb/task_scheduler_init.h"
#include "tbb/task.h"

#include "DivideTask.hpp"
#include "MergeTask.hpp"
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

	sche_init.initialize(thread_num);
	
	cal_data.set_model(md);
	cal_data.thread_num = thread_num;
	cal_data.sorted_pcl_var_id = 0;
	cal_data.prev_valid_pcl_num = md.pcl_num;

	init_pcl.cal_task_num(cal_data.prev_valid_pcl_num);
	tbb::task::spawn_root_and_wait(
		*new(tbb::task::allocate_root())
			MergeTask<Step_T2D_ME_Task::InitPcl, size_t, 8>(
				0, init_pcl.get_task_num(), init_pcl,
				cal_data.valid_pcl_num));
	return 0;
}

int Step_T2D_ME_TBB::finalize_calculation()
{
	Model_T2D_ME_mt& md = *(Model_T2D_ME_mt*)model;
	md.pcl_num = cal_data.prev_valid_pcl_num;
	sche_init.terminate();
	return 0;
}

int cal_substep_func_T2D_ME_TBB(void* _self)
{
	Step_T2D_ME_TBB& self = *(Step_T2D_ME_TBB*)(_self);
	Step_T2D_ME_Task::CalData& cd = self.cal_data;
	if (cd.valid_pcl_num == 0)
	{
		self.exit_calculation();
		return 0;
	}

	cd.dt = self.dtime;
	cd.sorted_pcl_var_id ^= 1;

	// sort pcl id
	tbb::task::spawn_root_and_wait(
		*new(tbb::task::allocate_root()) 
			SortParticleTask(cd.pcl_sort_mem, cd.prev_valid_pcl_num));

	tbb::task_list tk_list;
	// sort node
	tk_list.push_back(*new(tbb::task::allocate_root())
		SortTriMeshNodeTask(
			cd.node_sort_mem,
#ifdef _DEBUG
			cd.elem_num,
#endif
			cd.valid_pcl_num,
			cd.pcl_sort_mem.res_keys,
			cd.elem_node_id,
			cd.valid_elem_num));
	// map pcl to bg mesh
	auto &map_pcl_to_mesh = self.map_pcl_to_mesh;
	map_pcl_to_mesh.update();
	tk_list.push_back(*new(tbb::task::allocate_root())
		DivideTask<Step_T2D_ME_Task::MapPclToBgMesh, 8>(
			0, map_pcl_to_mesh.get_task_num(), map_pcl_to_mesh));
	tbb::task::spawn_root_and_wait(tk_list);

	// update nodal a and v
	auto &update_a_and_v = self.update_a_and_v;
	update_a_and_v.update();
	tbb::task::spawn_root_and_wait(
		*new(tbb::task::allocate_root())
			DivideTask<Step_T2D_ME_Task::UpdateAccelerationAndVelocity, 8>(
				0, update_a_and_v.get_task_num(), update_a_and_v));

	// cal element de and map to node
	auto &cal_elem_de = self.cal_elem_de;
	cal_elem_de.update();
	tbb::task::spawn_root_and_wait(
		*new(tbb::task::allocate_root())
			DivideTask<Step_T2D_ME_Task::CalElemDeAndMapToNode, 8>(
				0, cal_elem_de.get_task_num(), cal_elem_de));

	// cal strain increment at node
	auto &cal_node_de = self.cal_node_de;
	cal_node_de.update();
	tbb::task::spawn_root_and_wait(
		*new(tbb::task::allocate_root())
			DivideTask<Step_T2D_ME_Task::CalNodeDe, 8>(
				0, cal_node_de.get_task_num(), cal_node_de));

	// map bg mesh back to pcl
	cd.prev_valid_pcl_num = cd.valid_pcl_num;
	auto& map_mesh_to_pcl = self.map_mesh_to_pcl;
	map_mesh_to_pcl.update();
	tbb::task::spawn_root_and_wait(
		*new(tbb::task::allocate_root())
			MergeTask<Step_T2D_ME_Task::MapBgMeshToPcl, size_t, 8>(
				0, map_mesh_to_pcl.get_task_num(),
				map_mesh_to_pcl, cd.valid_pcl_num));

	self.continue_calculation();
	return 0;
}
