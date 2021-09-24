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
	init_pcl_res(init_pcl),
	init_pcl(cal_data),
	map_pcl_to_mesh(cal_data),
	cont_force(cont_rigid_cylinder),
	cont_rigid_cylinder(cal_data),
	update_a_and_v(cal_data),
	cal_elem_de(cal_data),
	cal_node_de(cal_data),
	map_mesh_to_pcl_res(map_mesh_to_pcl),
	map_mesh_to_pcl(cal_data),
	pcl_ranges(nullptr),
	node_elem_ranges(nullptr),
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
	
	const size_t max_pcl_task_num = ParallelUtils::cal_task_num<
		Step_T3D_CHM_Task::min_pcl_num_per_task,
		Step_T3D_CHM_Task::task_num_per_thread>(thread_num, md.pcl_num);
	if (pcl_ranges) delete[] pcl_ranges;
	pcl_ranges = new Step_T3D_CHM_Task::PclRange[max_pcl_task_num];
	cal_data.pcl_ranges = pcl_ranges;

	const size_t max_ne_task_num = ParallelUtils::cal_task_num<
		Step_T3D_CHM_Task::min_node_elem_num_per_task,
		Step_T3D_CHM_Task::task_num_per_thread>(thread_num, md.elem_num * 4);
	if (node_elem_ranges) delete[] node_elem_ranges;
	node_elem_ranges = new Step_T3D_CHM_Task::NodeElemRange[max_ne_task_num];
	cal_data.node_elem_ranges = node_elem_ranges;
	
	cal_data.sorted_pcl_var_id = 0;
	cal_data.prev_valid_pcl_num = md.pcl_num;

	init_pcl.init(thread_num);
	map_pcl_to_mesh.init();
	update_a_and_v.init();
	cal_elem_de.init();
	cal_node_de.init();
	map_mesh_to_pcl.init();
	if (md.has_rigid_cylinder())
		cont_rigid_cylinder.init(md);

	init_pcl_res.pcl_num = 0;
	//tbb::parallel_reduce(tbb::blocked_range<size_t>(0, init_pcl.get_task_num(), 1), init_pcl_res);
	ParallelUtils::parallel_reduce(init_pcl, init_pcl_res, init_pcl.get_task_num());
	cal_data.valid_pcl_num = init_pcl_res.pcl_num;

	return 0;
}

int Step_T3D_CHM_TBB::finalize_calculation()
{
	if (pcl_ranges) delete[] pcl_ranges;
	if (node_elem_ranges) delete[] node_elem_ranges;
	
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
	const size_t pcl_task_num = ParallelUtils::cal_task_num<
		Step_T3D_CHM_Task::min_pcl_num_per_task,
		Step_T3D_CHM_Task::task_num_per_thread>(
			thread_num, cd.valid_pcl_num);

	// sort pcl id
	auto& pcl_sort_mem = cd.pcl_sort_mem;
	tbb::task::spawn_root_and_wait(
		*new(tbb::task::allocate_root()) 
			SortParticleTask(pcl_sort_mem, cd.prev_valid_pcl_num));
	pcl_sort_mem.res_keys[cd.prev_valid_pcl_num] = SIZE_MAX;

	// sort node
	auto& node_sort_mem = cd.node_sort_mem;
	//...
	//cd.valid_elem_num = node_sort_mem.valid_elem_blocks;

	// map pcl to bg mesh
	auto &map_pcl_to_mesh = self.map_pcl_to_mesh;
	//tbb::parallel_for(tbb::blocked_range<size_t>(0, pcl_task_num, 1), map_pcl_to_mesh);
	ParallelUtils::parallel_for(map_pcl_to_mesh, pcl_task_num);

	// contact
	Model_T3D_CHM_mt& md = *static_cast<Model_T3D_CHM_mt *>(self.model);
	if (md.has_rigid_cylinder())
	{
		auto& cont_rigid_cylinder = self.cont_rigid_cylinder;
		auto& cont_force = self.cont_force;
		// cal contact force
		cont_rigid_cylinder.update();
		cont_force.react_force.reset();
		//tbb::parallel_reduce(tbb::blocked_range<size_t>(0, pcl_task_num, 1), cont_force);
		ParallelUtils::parallel_reduce(cont_rigid_cylinder, cont_force, pcl_task_num);
		// update motion
		RigidCylinder& rcy = md.get_rigid_cylinder();
		rcy.set_cont_force(cont_force.react_force);
		rcy.update_motion(self.dtime);
	}

	// update nodal a and v
	const size_t node_elem_task_num = ParallelUtils::cal_task_num<
		Step_T3D_CHM_Task::min_node_elem_num_per_task,
		Step_T3D_CHM_Task::task_num_per_thread>(
			thread_num, cd.valid_elem_num * 4);
	auto& update_a_and_v = self.update_a_and_v;
	update_a_and_v.update(node_elem_task_num);
	//tbb::parallel_for(tbb::blocked_range<size_t>(0, node_elem_task_num, 1), update_a_and_v);
	ParallelUtils::parallel_for(update_a_and_v, node_elem_task_num);

	// cal element de and map to node
	const size_t elem_task_num = ParallelUtils::cal_task_num<
		Step_T3D_CHM_Task::min_elem_num_per_task,
		Step_T3D_CHM_Task::task_num_per_thread>(
			thread_num, cd.valid_elem_num);
	auto& cal_elem_de = self.cal_elem_de;
	cal_elem_de.update(elem_task_num);
	//tbb::parallel_for(tbb::blocked_range<size_t>(0, elem_task_num, 1), cal_elem_de);
	ParallelUtils::parallel_for(cal_elem_de, elem_task_num);

	// cal strain increment at node
	auto &cal_node_de = self.cal_node_de;
	//tbb::parallel_for(tbb::blocked_range<size_t>(0, node_elem_task_num, 1), cal_node_de);
	ParallelUtils::parallel_for(cal_node_de, node_elem_task_num);

	cd.prev_valid_pcl_num = cd.valid_pcl_num;
	// map bg mesh back to pcl
	auto& map_mesh_to_pcl = self.map_mesh_to_pcl;
	auto& map_mesh_to_pcl_res = self.map_mesh_to_pcl_res;
	map_mesh_to_pcl.update(thread_num);
	map_mesh_to_pcl_res.pcl_num = 0;
	//tbb::parallel_reduce(tbb::blocked_range<size_t>(0, pcl_task_num, 1), map_mesh_to_pcl_res);
	ParallelUtils::parallel_reduce(map_mesh_to_pcl, map_mesh_to_pcl_res, pcl_task_num);
	cd.valid_pcl_num = map_mesh_to_pcl_res.pcl_num;

	self.continue_calculation();
	return 0;
}
