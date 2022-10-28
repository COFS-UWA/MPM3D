#include "SimulationsOMP_pcp.h"

#include <fstream>

#include "ParallelForTask.hpp"
#include "ParallelReduceTask.hpp"
#include "Step_T3D_CHM_up_TBB.h"

static std::fstream res_file_t3d_me_tbb;

Step_T3D_CHM_up_TBB::Step_T3D_CHM_up_TBB(const char* _name) :
	Step_TBB(_name, "Step_T3D_CHM_up_TBB", &cal_substep_func_T3D_CHM_up_TBB),
	sche_init(tbb::task_scheduler_init::deferred),
	init_pcl(*this),
	map_pcl_to_mesh(*this),
	cont_rigid_body(*this),
	update_a_and_v(*this),
	cal_elem_de(*this),
	cal_node_de(*this),
	map_mesh_to_pcl(*this) {}

Step_T3D_CHM_up_TBB::~Step_T3D_CHM_up_TBB() {}

int Step_T3D_CHM_up_TBB::init_calculation()
{
	res_file_t3d_me_tbb.open("step_t3d_me_tbb.csv", std::ios::out | std::ios::binary);
	res_file_t3d_me_tbb << "substep_id, pcl_sort, ne_sort, map_pcl_to_mesh, "
		"update_a_and_v, cal_elem_de, cal_node_de, map_mesh_to_pcl, "
		"map_mesh_to_pcl0, map_mesh_to_pcl1\n";

	Model_T3D_CHM_up_mt &md = *(Model_T3D_CHM_up_mt*)model;
	if (md.pcl_num == 0)
		return -1;

	sche_init.initialize(thread_num);
	
	pmodel = &md;

	size_t e_num = md.get_elem_num();
	for (size_t e_id = 0; e_id < e_num; ++e_id)
		md.elem_has_pcls[e_id] = SIZE_MAX;

	auto& spva0 = spvas[0];
	const auto& md_spva0 = md.sorted_pcl_var_arrays[0];
	spva0.pcl_index = md_spva0.pcl_index;
	spva0.pcl_n = md_spva0.pcl_n;
	spva0.pcl_density_f = md_spva0.pcl_density_f;
	spva0.pcl_v_s = md_spva0.pcl_v_s;
	spva0.pcl_u = md_spva0.pcl_u;
	spva0.pcl_stress = md_spva0.pcl_stress;
	spva0.pcl_p = md_spva0.pcl_p;
	spva0.pcl_strain = md_spva0.pcl_strain;
	spva0.pcl_estrain = md_spva0.pcl_estrain;
	spva0.pcl_pstrain = md_spva0.pcl_pstrain;
	spva0.pcl_N = md_spva0.pcl_N;

	auto& spva1 = spvas[1];
	const auto& md_spva1 = md.sorted_pcl_var_arrays[1];
	spva1.pcl_index = md_spva1.pcl_index;
	spva1.pcl_n = md_spva1.pcl_n;
	spva1.pcl_density_f = md_spva1.pcl_density_f;
	spva1.pcl_v_s = md_spva1.pcl_v_s;
	spva1.pcl_u = md_spva1.pcl_u;
	spva1.pcl_stress = md_spva1.pcl_stress;
	spva1.pcl_p = md_spva1.pcl_p;
	spva1.pcl_strain = md_spva1.pcl_strain;
	spva1.pcl_estrain = md_spva1.pcl_estrain;
	spva1.pcl_pstrain = md_spva1.pcl_pstrain;
	spva1.pcl_N = md_spva1.pcl_N;

	pcl_sort.init(md.ori_pcl_num, thread_num);
	ne_sort.init(md.elem_num, md.node_num, md.ori_pcl_num, thread_num,
		pcl_sort.out_pcl_in_elems(), (size_t*)md.elem_node_id);
	
	in_pcl_in_elems = pcl_sort.in_pcl_in_elems();
	in_prev_pcl_ids = pcl_sort.in_prev_pcl_ids();
	pcl_in_elems = pcl_sort.out_pcl_in_elems();
	prev_pcl_ids = pcl_sort.out_prev_pcl_ids();
	elem_ids = ne_sort.elem_ids();
	node_ids = ne_sort.node_ids();
	node_elem_offs = ne_sort.node_elem_offs();
		
	prev_valid_pcl_num = md.pcl_num;
#ifdef _DEBUG
	ori_pcl_num = md.ori_pcl_num;
	elem_num = md.elem_num;
	node_num = md.node_num;
#endif

	init_pcl.init(thread_num);
	ParaUtil::parallel_reduce(init_pcl, init_pcl_res, init_pcl.get_task_num());	
	valid_pcl_num = init_pcl_res.pcl_num;
	
	map_pcl_to_mesh.init();
	cont_rigid_body.init(init_pcl_res.max_pcl_vol, is_first_step);
	update_a_and_v.init();
	cal_elem_de.init();
	cal_node_de.init();
	map_mesh_to_pcl.init();
	map_pcl_to_mesh_res.react_force.reset();

	return 0;
}

int Step_T3D_CHM_up_TBB::finalize_calculation()
{
	Model_T3D_CHM_up_mt& md = *(Model_T3D_CHM_up_mt *)model;
	md.pcl_num = prev_valid_pcl_num;

	const auto& spva0 = spvas[prev_spva_id()];
	auto& md_spva0 = md.sorted_pcl_var_arrays[0];
	md_spva0.pcl_index = spva0.pcl_index;
	md_spva0.pcl_n = spva0.pcl_n;
	md_spva0.pcl_density_f = spva0.pcl_density_f;
	md_spva0.pcl_v_s = spva0.pcl_v_s;
	md_spva0.pcl_u = spva0.pcl_u;
	md_spva0.pcl_stress = spva0.pcl_stress;
	md_spva0.pcl_p = spva0.pcl_p;
	md_spva0.pcl_strain = spva0.pcl_strain;
	md_spva0.pcl_estrain = spva0.pcl_estrain;
	md_spva0.pcl_pstrain = spva0.pcl_pstrain;
	md_spva0.pcl_N = spva0.pcl_N;

	const auto& spva1 = spvas[next_spva_id()];
	auto& md_spva1 = md.sorted_pcl_var_arrays[1];
	md_spva1.pcl_index = spva1.pcl_index;
	md_spva1.pcl_n = spva1.pcl_n;
	md_spva1.pcl_density_f = spva1.pcl_density_f;
	md_spva1.pcl_v_s = spva1.pcl_v_s;
	md_spva1.pcl_u = spva1.pcl_u;
	md_spva1.pcl_stress = spva1.pcl_stress;
	md_spva1.pcl_p = spva1.pcl_p;
	md_spva1.pcl_strain = spva1.pcl_strain;
	md_spva1.pcl_estrain = spva1.pcl_estrain;
	md_spva1.pcl_pstrain = spva1.pcl_pstrain;
	md_spva1.pcl_N = spva1.pcl_N;

	sche_init.terminate();
	return 0;
}

int cal_substep_func_T3D_CHM_up_TBB(void* _self)
{
	Step_T3D_CHM_up_TBB& self = *(Step_T3D_CHM_up_TBB*)(_self);
	if (self.valid_pcl_num == 0)
	{
		self.exit_calculation();
		return 0;
	}

	// sort pcl id
	self.pcl_sort.sort(self.prev_valid_pcl_num);

	// sort node
	self.ne_sort.sort(self.valid_pcl_num);
	self.valid_elem_num = self.ne_sort.elem_num();

	size_t task_num;

	// map pcl to bg mesh
	task_num = ParaUtil::cal_task_num<
		Step_T3D_CHM_up_TBB_Task::min_pcl_num_per_task,
		Step_T3D_CHM_up_TBB_Task::map_pcl_to_mesh_task_num_per_thread>(
			self.thread_num, self.valid_pcl_num);
	self.map_pcl_to_mesh.update(task_num);
	self.cont_rigid_body.update();
	ParaUtil::parallel_reduce(self.map_pcl_to_mesh, self.map_pcl_to_mesh_res, task_num);
	
	// update rb motion
	Model_T3D_CHM_up_mt& md = *static_cast<Model_T3D_CHM_up_mt*>(self.model);
	if (md.has_t3d_rigid_mesh())
	{
		RigidObjectByT3DMesh &rmesh = md.get_t3d_rigid_mesh();
		rmesh.set_cont_force(self.react_force);
		rmesh.update_motion(self.dtime);
	}
	
	// update nodal a and v
	task_num = ParaUtil::cal_task_num<
		Step_T3D_CHM_up_TBB_Task::min_node_elem_num_per_task,
		Step_T3D_CHM_up_TBB_Task::update_node_av_task_num_per_thread>(
			self.thread_num, self.valid_elem_num * 4);
	self.update_a_and_v.update(task_num);
	ParaUtil::parallel_for(self.update_a_and_v, task_num);

	// cal element de and map to node
	task_num = ParaUtil::cal_task_num<
		Step_T3D_CHM_up_TBB_Task::min_elem_num_per_task,
		Step_T3D_CHM_up_TBB_Task::cal_elem_de_task_num_per_thread>(
			self.thread_num, self.valid_elem_num);
	self.cal_elem_de.update(task_num);
	ParaUtil::parallel_for(self.cal_elem_de, task_num);

	// cal strain increment at node
	task_num = ParaUtil::cal_task_num<
		Step_T3D_CHM_up_TBB_Task::min_node_elem_num_per_task,
		Step_T3D_CHM_up_TBB_Task::cal_node_de_task_num_per_thread>(
			self.thread_num, self.valid_elem_num * 4);
	self.cal_node_de.update(task_num);
	ParaUtil::parallel_for(self.cal_node_de, task_num);

	// map bg mesh back to pcl
#ifdef _DEBUG
	self.prev_valid_pcl_num_tmp = self.prev_valid_pcl_num;
#endif // _DEBUG
	self.prev_valid_pcl_num = self.valid_pcl_num;
	task_num = ParaUtil::cal_task_num<
		Step_T3D_CHM_up_TBB_Task::min_pcl_num_per_task,
		Step_T3D_CHM_up_TBB_Task::map_mesh_to_pcl_task_num_per_thread>(
			self.thread_num, self.prev_valid_pcl_num);
	self.map_mesh_to_pcl.update(task_num);
	ParaUtil::parallel_reduce(self.map_mesh_to_pcl, self.map_mesh_to_pcl_res, task_num);
	
	self.continue_calculation();
	return 0;
}
