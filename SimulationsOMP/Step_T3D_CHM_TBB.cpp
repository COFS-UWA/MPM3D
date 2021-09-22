#include "SimulationsOMP_pcp.h"

#include <fstream>

#include "ParallelForTask.hpp"
#include "ParallelReduceTask.hpp"
#include "Step_T3D_CHM_TBB.h"

#ifdef _DEBUG
static std::fstream res_file_t3d_me_tbb;
#endif

Step_T3D_CHM_TBB::Step_T3D_CHM_TBB(const char* _name) :
	Step_TBB(_name, "Step_T3D_CHM_TBB", &substep_func_T3D_CHM_TBB),
	init_pcl(cal_data),
	map_pcl_to_mesh(cal_data),
	//cont_rigid_rect(cal_data),
	update_a_and_v(cal_data),
	cal_elem_de(cal_data),
	cal_node_de(cal_data),
	map_mesh_to_pcl(cal_data),
	sche_init(tbb::task_scheduler_init::deferred) {}

Step_T3D_CHM_TBB::~Step_T3D_CHM_TBB() {}

int Step_T3D_CHM_TBB::init_calculation()
{
#ifdef _DEBUG
	res_file_t3d_me_tbb.open("Step_T3D_CHM_tbb.txt", std::ios::out | std::ios::binary);
#endif

	Model_T3D_CHM_mt &md = *(Model_T3D_CHM_mt*)model;
	if (md.pcl_num == 0)
		return -1;

	sche_init.initialize(thread_num);
	
	cal_data.set_model(md);
	cal_data.thread_num = thread_num;
	cal_data.thread_bin_blocks_mem.init(thread_num, 2);
	cal_data.pcl_sort_mem.init(md.pcl_num, md.elem_num, cal_data.thread_bin_blocks_mem);
	cal_data.node_sort_mem.init(md.elem_num, md.node_num, cal_data.thread_bin_blocks_mem);
	cal_data.sorted_pcl_var_id = 0;
	cal_data.prev_valid_pcl_num = md.pcl_num;

	init_pcl.init(thread_num);
	map_pcl_to_mesh.init();
	update_a_and_v.init();
	cal_elem_de.init();
	cal_node_de.init();
	map_mesh_to_pcl.init();
	//if (md.has_rigid_rect())
	//	cont_rigid_rect.init(md);

	init_pcl_res.pcl_num = 0;
	//tbb::parallel_reduce(tbb::blocked_range<size_t>(0, init_pcl.get_task_num(), 1), init_pcl_res);
	ParallelUtils::parallel_reduce(init_pcl, init_pcl_res, init_pcl.get_task_num());
	cal_data.valid_pcl_num = init_pcl_res.pcl_num;

	return 0;
}

int Step_T3D_CHM_TBB::finalize_calculation()
{
	Model_T3D_CHM_mt& md = *(Model_T3D_CHM_mt *)model;
	md.pcl_num = cal_data.prev_valid_pcl_num;
	sche_init.terminate();
	return 0;
}

int substep_func_T3D_CHM_TBB(void* _self)
{
	Step_T3D_CHM_TBB& self = *(Step_T3D_CHM_TBB*)(_self);
	Step_T3D_CHM_Task::CalData& cd = self.cal_data;
	if (cd.valid_pcl_num == 0)
	{
		self.exit_calculation();
		return 0;
	}

	cd.dt = self.dtime;
	cd.sorted_pcl_var_id ^= 1;

	const size_t thread_num = self.thread_num;

	// sort pcl id
	auto& pcl_sort_mem = cd.pcl_sort_mem;
	tbb::task::spawn_root_and_wait(
		*new(tbb::task::allocate_root()) 
			SortParticleTask(pcl_sort_mem, cd.prev_valid_pcl_num));
	pcl_sort_mem.res_keys[cd.prev_valid_pcl_num] = SIZE_MAX;

	// sort node
	auto& node_sort_mem = cd.node_sort_mem;
	//

	// map pcl to bg mesh
	auto &map_pcl_to_mesh = self.map_pcl_to_mesh;
	//tbb::parallel_for(tbb::blocked_range<size_t>(0, map_pcl_to_mesh.get_task_num(), 1), map_pcl_to_mesh);
	ParallelUtils::parallel_for(map_pcl_to_mesh, map_pcl_to_mesh.get_task_num());

	// contact with rigid rect
	Model_T3D_CHM_mt& md = *static_cast<Model_T3D_CHM_mt *>(self.model);
	//if (md.has_rigid_rect())
	//{
	//	auto& cont_rigid_rect = self.cont_rigid_rect;
	//	cont_rigid_rect.update(self.thread_num);
	//	Force2D rr_cf;
	//	tbb::task::spawn_root_and_wait(
	//		*new(tbb::task::allocate_root())
	//			MergeTask<Step_T3D_CHM_Task::ContactRigidRect, Force2D, 8>(
	//				0, cont_rigid_rect.get_task_num(), cont_rigid_rect, rr_cf));
	//	RigidRect& rr = md.get_rigid_rect();
	//	rr.set_cont_force(rr_cf.fx, rr_cf.fy, rr_cf.m);
	//	rr.update_motion(self.dtime);
	//}

	// update nodal a and v
	auto& update_a_and_v = self.update_a_and_v;
	update_a_and_v.update(thread_num);
	//tbb::parallel_for(tbb::blocked_range<size_t>(0, update_a_and_v.get_task_num(), 1), update_a_and_v);
	ParallelUtils::parallel_for(update_a_and_v, update_a_and_v.get_task_num());

	// cal element de and map to node
	auto& cal_elem_de = self.cal_elem_de;
	cal_elem_de.update(thread_num);
	//tbb::parallel_for(tbb::blocked_range<size_t>(0, cal_elem_de.get_task_num(), 1), cal_elem_de);
	ParallelUtils::parallel_for(cal_elem_de, cal_elem_de.get_task_num());

	// cal strain increment at node
	auto &cal_node_de = self.cal_node_de;
	cal_node_de.update(thread_num);
	//tbb::parallel_for(tbb::blocked_range<size_t>(0, cal_node_de.get_task_num(), 1), cal_node_de);
	ParallelUtils::parallel_for(cal_node_de, cal_node_de.get_task_num());

	cd.prev_valid_pcl_num = cd.valid_pcl_num;
	// map bg mesh back to pcl
	auto& map_mesh_to_pcl_res = self.map_mesh_to_pcl_res;
	map_mesh_to_pcl_res.pcl_num = 0;
	auto& map_mesh_to_pcl = self.map_mesh_to_pcl;
	map_mesh_to_pcl.update(thread_num);
	//tbb::parallel_reduce(tbb::blocked_range<size_t>(0, map_mesh_to_pcl.get_task_num(), 1), map_mesh_to_pcl_res);
	ParallelUtils::parallel_reduce(map_mesh_to_pcl, map_mesh_to_pcl_res, map_mesh_to_pcl.get_task_num());
	cd.valid_pcl_num = map_mesh_to_pcl_res.pcl_num;

	self.continue_calculation();
	return 0;
}
