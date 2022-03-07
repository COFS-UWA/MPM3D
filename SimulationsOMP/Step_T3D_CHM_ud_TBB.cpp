#include "SimulationsOMP_pcp.h"

#include <chrono>
#include <fstream>

#include "ParallelForTask.hpp"
#include "ParallelReduceTask.hpp"
#include "Step_T3D_CHM_ud_TBB.h"

#define TIMING

static std::fstream res_file_t3d_chm_ud_tbb;

Step_T3D_CHM_ud_TBB::Step_T3D_CHM_ud_TBB(const char* _name) :
	Step_TBB(_name, "Step_T3D_CHM_ud_TBB", &substep_func_T3D_CHM_ud_TBB),
	sche_init(tbb::task_scheduler_init::deferred),
	init_pcl(*this),
	map_pcl_to_mesh(*this),
	cont_rigid_body(*this),
	update_a_and_v(*this),
	cal_elem_de(*this),
	cal_node_de(*this),
	map_mesh_to_pcl(*this),
	init_pcl_tbb(init_pcl),
	map_pcl_to_mesh_tbb(map_pcl_to_mesh),
	map_mesh_to_pcl_tbb(map_mesh_to_pcl) {}

Step_T3D_CHM_ud_TBB::~Step_T3D_CHM_ud_TBB() {}

int Step_T3D_CHM_ud_TBB::init_calculation()
{
#ifdef _DEBUG
	res_file_t3d_chm_ud_tbb.open("Step_T3D_CHM_ud_TBB.csv", std::ios::out | std::ios::binary);
	res_file_t3d_chm_ud_tbb << "pcl_sort, ne_sort, map_pcl_to_mesh, "
		"update_a_and_v, cal_elem_de, cal_node_de, map_mesh_to_pcl\n";
#endif

	sche_init.initialize(thread_num);

	Model_T3D_CHM_mt &md = *(Model_T3D_CHM_mt*)model;
	if (md.pcl_num == 0)
		return -1;
	
	pmodel = &md;

	pcl_m_s = md.pcl_m_s;
	pcl_density_s = md.pcl_density_s;
	pcl_vol_s = md.pcl_vol_s;
	pcl_bf_s = md.pcl_bf_s;
	pcl_bf_f = md.pcl_bf_f;
	pcl_t = md.pcl_t;
	pcl_pos = md.pcl_pos;
	pcl_vol = md.pcl_vol;
	pcl_mat_model = md.pcl_mat_model;

	auto& spva0 = spvas[0];
	const auto& md_spva0 = md.sorted_pcl_var_arrays[0];
	spva0.pcl_index = md_spva0.pcl_index;
	spva0.pcl_n = md_spva0.pcl_n;
	spva0.pcl_density_f = md_spva0.pcl_density_f;
	spva0.pcl_v_s = md_spva0.pcl_v_s;
	spva0.pcl_v_f = md_spva0.pcl_v_f;
	spva0.pcl_u_s = md_spva0.pcl_u_s;
	spva0.pcl_u_f = md_spva0.pcl_u_f;
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
	spva1.pcl_v_f = md_spva1.pcl_v_f;
	spva1.pcl_u_s = md_spva1.pcl_u_s;
	spva1.pcl_u_f = md_spva1.pcl_u_f;
	spva1.pcl_stress = md_spva1.pcl_stress;
	spva1.pcl_p = md_spva1.pcl_p;
	spva1.pcl_strain = md_spva1.pcl_strain;
	spva1.pcl_estrain = md_spva1.pcl_estrain;
	spva1.pcl_pstrain = md_spva1.pcl_pstrain;
	spva1.pcl_N = md_spva1.pcl_N;

	elem_node_id = md.elem_node_id;
	elem_N_abc = md.elem_N_abc;
	elem_N_d = md.elem_N_d;
	elem_vol = md.elem_vol;
	elem_density_f = md.elem_density_f;
	elem_pcl_n = md.elem_pcl_n;
	elem_pcl_m_s = md.elem_pcl_m_s;
	elem_pcl_m_f = md.elem_pcl_m_f;
	elem_de = md.elem_de;
	elem_p = md.elem_p;
	elem_n2_miu_div_k_vol = md.elem_n2_miu_div_k_vol;
	elem_seep_force = md.elem_seep_force;
	elem_m_de_vol_s = md.elem_m_de_vol_s;
	elem_m_de_vol_f = md.elem_m_de_vol_f;

	// element-node data
	elem_node_vm_s = md.elem_node_vm_s;
	elem_node_vm_f = md.elem_node_vm_f;
	elem_node_force_s = md.elem_node_force_s;
	elem_node_force_f = md.elem_node_force_f;

	// node data
	node_pos = md.node_pos;
	node_a_s = md.node_a_s;
	node_a_f = md.node_a_f;
	node_v_s = md.node_v_s;
	node_v_f = md.node_v_f;
	node_has_vbc_s = md.node_has_vbc_s;
	node_has_vbc_f = md.node_has_vbc_f;
	node_am_s = md.node_am_s;
	node_am_f = md.node_am_f;
	node_de_vol_s = md.node_de_vol_s;
	node_de_vol_f = md.node_de_vol_f;

	Kf = md.Kf;
	miu = md.miu;
	k = md.k;

	pcl_sort.init(md.ori_pcl_num, thread_num);
	ne_sort.init(md.elem_num, md.node_num, md.ori_pcl_num, thread_num,
		pcl_sort.out_pcl_in_elems(), (size_t *)elem_node_id);

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
	//init_pcl_tbb.reset();
	//tbb::parallel_reduce(tbb::blocked_range<size_t>(0, init_pcl.get_task_num(), 1), init_pcl_res);
	//valid_pcl_num = init_pcl_tbb.pcl_num;

	valid_pcl_num = init_pcl_res.pcl_num;

	map_pcl_to_mesh.init();
	cont_rigid_body.init(init_pcl_res.max_pcl_vol);
	update_a_and_v.init();
	cal_elem_de.init();
	cal_node_de.init();
	map_mesh_to_pcl.init();

	pcl_sort_time = 0;
	ne_sort_time = 0;
	map_pcl_to_mesh_time = 0;
	update_a_and_v_time = 0;
	cal_elem_de_time = 0;
	cal_node_de_time = 0;
	map_mesh_to_pcl_time = 0;

	return 0;
}

int Step_T3D_CHM_ud_TBB::finalize_calculation()
{
	Model_T3D_CHM_mt& md = *(Model_T3D_CHM_mt *)model;
	md.pcl_num = prev_valid_pcl_num;

	const auto& spva0 = spvas[prev_spva_id()];
	auto& md_spva0 = md.sorted_pcl_var_arrays[0];
	md_spva0.pcl_index = spva0.pcl_index;
	md_spva0.pcl_n = spva0.pcl_n;
	md_spva0.pcl_density_f = spva0.pcl_density_f;
	md_spva0.pcl_v_s = spva0.pcl_v_s;
	md_spva0.pcl_v_f = spva0.pcl_v_f;
	md_spva0.pcl_u_s = spva0.pcl_u_s;
	md_spva0.pcl_u_f = spva0.pcl_u_f;
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
	md_spva1.pcl_v_f = spva1.pcl_v_f;
	md_spva1.pcl_u_s = spva1.pcl_u_s;
	md_spva1.pcl_u_f = spva1.pcl_u_f;
	md_spva1.pcl_stress = spva1.pcl_stress;
	md_spva1.pcl_p = spva1.pcl_p;
	md_spva1.pcl_strain = spva1.pcl_strain;
	md_spva1.pcl_estrain = spva1.pcl_estrain;
	md_spva1.pcl_pstrain = spva1.pcl_pstrain;
	md_spva1.pcl_N = spva1.pcl_N;

	sche_init.terminate();
	return 0;
}

int substep_func_T3D_CHM_ud_TBB(void* _self)
{
	Step_T3D_CHM_ud_TBB& self = *(Step_T3D_CHM_ud_TBB*)(_self);
	if (self.valid_pcl_num == 0)
	{
		self.exit_calculation();
		return 0;
	}
	
	// timing
	std::chrono::high_resolution_clock::time_point t0, t1;

	// sort pcl id
	t0 = std::chrono::high_resolution_clock::now();
	self.pcl_sort.sort(self.prev_valid_pcl_num);
	t1 = std::chrono::high_resolution_clock::now();
	self.pcl_sort_time += (t1 - t0).count();

	// sort node
	t0 = std::chrono::high_resolution_clock::now();
	self.ne_sort.sort(self.valid_pcl_num);
	t1 = std::chrono::high_resolution_clock::now();
	self.ne_sort_time += (t1 - t0).count();
	self.valid_elem_num = self.ne_sort.elem_num();

	// map pcl to bg mesh
	size_t task_num = ParaUtil::cal_task_num<
		Step_T3D_CHM_ud_TBB_Task::min_pcl_num_per_task,
		Step_T3D_CHM_ud_TBB_Task::map_pcl_to_mesh_task_num_per_thread>(
			self.thread_num, self.valid_pcl_num);
	self.map_pcl_to_mesh.update(task_num);
	self.cont_rigid_body.update();
	t0 = std::chrono::high_resolution_clock::now();
	ParaUtil::parallel_reduce(self.map_pcl_to_mesh, self.map_pcl_to_mesh_res, task_num);
	t1 = std::chrono::high_resolution_clock::now();
	self.map_pcl_to_mesh_time += (t1 - t0).count();
	//self.map_pcl_to_mesh_tbb.reset();
	//tbb::parallel_reduce(tbb::blocked_range<size_t>(0, pcl_task_num, 1), self.map_pcl_to_mesh_tbb);
	//self.react_force = self.map_pcl_to_mesh_tbb.res.react_force;

	// update rigid body motion
	Model_T3D_CHM_mt& md = *static_cast<Model_T3D_CHM_mt *>(self.model);
	if (md.has_rigid_cylinder())
	{
		RigidCylinder& rcy = md.get_rigid_cylinder();
		rcy.set_cont_force(self.react_force);
		rcy.update_motion(self.dtime);
	}
	if (md.has_t3d_rigid_mesh())
	{
		RigidObjectByT3DMesh &rmesh = md.get_t3d_rigid_mesh();
		rmesh.set_cont_force(self.react_force);
		rmesh.update_motion(self.dtime);
	}

	// update nodal a and v
	task_num = ParaUtil::cal_task_num<
		Step_T3D_CHM_ud_TBB_Task::min_node_elem_num_per_task,
		Step_T3D_CHM_ud_TBB_Task::update_node_av_task_num_per_thread>(
			self.thread_num, self.valid_elem_num * 4);
	self.update_a_and_v.update(task_num);
	t0 = std::chrono::high_resolution_clock::now();
	ParaUtil::parallel_for(self.update_a_and_v, task_num);
	//tbb::parallel_for(tbb::blocked_range<size_t>(0, node_elem_task_num, 1), update_a_and_v);
	t1 = std::chrono::high_resolution_clock::now();
	self.update_a_and_v_time += (t1 - t0).count();

	// cal element de and map to node
	task_num = ParaUtil::cal_task_num<
		Step_T3D_CHM_ud_TBB_Task::min_elem_num_per_task,
		Step_T3D_CHM_ud_TBB_Task::cal_elem_de_task_num_per_thread>(
			self.thread_num, self.valid_elem_num);
	self.cal_elem_de.update(task_num);
	t0 = std::chrono::high_resolution_clock::now();
	ParaUtil::parallel_for(self.cal_elem_de, task_num);
	//tbb::parallel_for(tbb::blocked_range<size_t>(0, elem_task_num, 1), self.cal_elem_de);
	t1 = std::chrono::high_resolution_clock::now();
	self.cal_elem_de_time += (t1 - t0).count();

	// cal strain increment at node
	task_num = ParaUtil::cal_task_num<
		Step_T3D_CHM_ud_TBB_Task::min_elem_num_per_task,
		Step_T3D_CHM_ud_TBB_Task::cal_node_de_task_num_per_thread>(
			self.thread_num, self.valid_elem_num);
	t0 = std::chrono::high_resolution_clock::now();
	self.cal_node_de.update(task_num);
	ParaUtil::parallel_for(self.cal_node_de, task_num);
	//tbb::parallel_for(tbb::blocked_range<size_t>(0, node_elem_task_num, 1), self.cal_node_de);
	t1 = std::chrono::high_resolution_clock::now();
	self.cal_node_de_time += (t1 - t0).count();

	// map bg mesh back to pcl
#ifdef _DEBUG
	self.prev_valid_pcl_num_tmp = self.prev_valid_pcl_num;
#endif
	self.prev_valid_pcl_num = self.valid_pcl_num;
	task_num = ParaUtil::cal_task_num<
		Step_T3D_CHM_ud_TBB_Task::min_pcl_num_per_task,
		Step_T3D_CHM_ud_TBB_Task::map_mesh_to_pcl_task_num_per_thread>(
			self.thread_num, self.prev_valid_pcl_num);
	self.map_mesh_to_pcl.update(task_num);
	t0 = std::chrono::high_resolution_clock::now();
	ParaUtil::parallel_reduce(self.map_mesh_to_pcl, self.map_mesh_to_pcl_res, task_num);
	//self.map_mesh_to_pcl_tbb.reset();
	//tbb::parallel_reduce(tbb::blocked_range<size_t>(0, pcl_task_num, 1), self.map_mesh_to_pcl_res);
	//self.valid_pcl_num = self.map_mesh_to_pcl_tbb.res.pcl_num;
	t1 = std::chrono::high_resolution_clock::now();
	self.map_mesh_to_pcl_time += (t1 - t0).count();

//	if (self.substep_index % 100 == 99)
//	{
//#ifdef TIMING
//		res_file_t3d_chm_ud_tbb << self.pcl_sort_time << ", "
//			<< self.ne_sort_time << ", "
//			<< self.map_pcl_to_mesh_time << ", "
//			<< self.update_a_and_v_time << ", "
//			<< self.cal_elem_de_time << ", "
//			<< self.cal_node_de_time << ", "
//			<< self.map_mesh_to_pcl_time << "\n";
//#endif
//		self.pcl_sort_time = 0;
//		self.ne_sort_time = 0;
//		self.map_pcl_to_mesh_time = 0;
//		self.update_a_and_v_time = 0;
//		self.cal_elem_de_time = 0;
//		self.cal_node_de_time = 0;
//		self.map_mesh_to_pcl_time = 0;
//	}

	self.continue_calculation();
	return 0;
}
