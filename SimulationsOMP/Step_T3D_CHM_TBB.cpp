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
	sche_init(tbb::task_scheduler_init::deferred),
	init_pcl_res(init_pcl),
	init_pcl(*this),
	map_pcl_to_mesh(*this),
	cont_force(cont_rigid_cylinder),
	cont_rigid_cylinder(*this),
	update_a_and_v(*this),
	cal_elem_de(*this),
	cal_node_de(*this),
	map_mesh_to_pcl_res(map_mesh_to_pcl),
	map_mesh_to_pcl(*this),
	pcl_ranges(nullptr),
	node_elem_ranges(nullptr) {}

Step_T3D_CHM_TBB::~Step_T3D_CHM_TBB() {}

int Step_T3D_CHM_TBB::init_calculation()
{
#ifdef _DEBUG
	res_file_t3d_me_tbb.open("Step_T3D_CHM_tbb.txt", std::ios::out | std::ios::binary);
#endif

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

#ifdef _DEBUG
	ori_pcl_num = md.ori_pcl_num;
	elem_num = md.elem_num;
	node_num = md.node_num;
#endif
	
	sche_init.initialize(thread_num);

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
	
	const size_t max_pcl_task_num = ParaUtil::cal_task_num<
		Step_T3D_CHM_TBB_Task::min_pcl_num_per_task,
		Step_T3D_CHM_TBB_Task::task_num_per_thread>(
			thread_num, md.ori_pcl_num);
	const size_t max_ne_task_num = ParaUtil::cal_task_num<
		Step_T3D_CHM_TBB_Task::min_node_elem_num_per_task,
		Step_T3D_CHM_TBB_Task::task_num_per_thread>(
			thread_num, md.elem_num * 4);
	pcl_ranges = (PclRange *)range_mem.alloc(
		  max_pcl_task_num * sizeof(PclRange)
		+ max_ne_task_num * sizeof(NodeElemRange));
	node_elem_ranges = (NodeElemRange *)(((char *)pcl_ranges) + max_pcl_task_num * sizeof(PclRange));

	prev_valid_pcl_num = md.pcl_num;

	init_pcl.init(thread_num);
	map_pcl_to_mesh.init();
	if (md.has_rigid_cylinder())
		cont_rigid_cylinder.init(md);
	update_a_and_v.init();
	cal_elem_de.init();
	cal_node_de.init();
	map_mesh_to_pcl.init();

	init_pcl_res.pcl_num = 0;
	//tbb::parallel_reduce(tbb::blocked_range<size_t>(0, init_pcl.get_task_num(), 1), init_pcl_res);
	ParaUtil::parallel_reduce(init_pcl, init_pcl_res, init_pcl.get_task_num());
	valid_pcl_num = init_pcl_res.pcl_num;

	return 0;
}

int Step_T3D_CHM_TBB::finalize_calculation()
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

int substep_func_T3D_CHM_TBB(void* _self)
{
	Step_T3D_CHM_TBB& self = *(Step_T3D_CHM_TBB*)(_self);
	if (self.valid_pcl_num == 0)
	{
		self.exit_calculation();
		return 0;
	}
	
	// sort pcl id
	self.pcl_sort.sort(self.valid_pcl_num);
	
	// sort node
	self.ne_sort.sort(self.valid_pcl_num);
	self.valid_elem_num = self.ne_sort.elem_num();

	// map pcl to bg mesh
	const size_t pcl_task_num = ParaUtil::cal_task_num<
		Step_T3D_CHM_TBB_Task::min_pcl_num_per_task,
		Step_T3D_CHM_TBB_Task::task_num_per_thread>(
			self.thread_num, self.valid_pcl_num);
	self.map_pcl_to_mesh.update(pcl_task_num);
	//tbb::parallel_for(tbb::blocked_range<size_t>(0, pcl_task_num, 1), self.map_pcl_to_mesh);
	ParaUtil::parallel_for(self.map_pcl_to_mesh, pcl_task_num);

	// contact
	Model_T3D_CHM_mt& md = *static_cast<Model_T3D_CHM_mt *>(self.model);
	if (md.has_rigid_cylinder())
	{
		// cal contact force
		self.cont_force.react_force.reset();
		self.cont_rigid_cylinder.update();
		//tbb::parallel_reduce(tbb::blocked_range<size_t>(0, pcl_task_num, 1), cont_force);
		ParaUtil::parallel_reduce(self.cont_rigid_cylinder, self.cont_force, pcl_task_num);
		
		// update motion
		RigidCylinder& rcy = md.get_rigid_cylinder();
		rcy.set_cont_force(self.cont_force.react_force);
		rcy.update_motion(self.dtime);
	}

	// update nodal a and v
	const size_t node_elem_task_num = ParaUtil::cal_task_num<
		Step_T3D_CHM_TBB_Task::min_node_elem_num_per_task,
		Step_T3D_CHM_TBB_Task::task_num_per_thread>(
			self.thread_num, self.valid_elem_num * 4);
	self.update_a_and_v.update(node_elem_task_num);
	//tbb::parallel_for(tbb::blocked_range<size_t>(0, node_elem_task_num, 1), update_a_and_v);
	ParaUtil::parallel_for(self.update_a_and_v, node_elem_task_num);

	// cal element de and map to node
	const size_t elem_task_num = ParaUtil::cal_task_num<
		Step_T3D_CHM_TBB_Task::min_elem_num_per_task,
		Step_T3D_CHM_TBB_Task::task_num_per_thread>(
			self.thread_num, self.valid_elem_num);
	auto& cal_elem_de = self.cal_elem_de;
	cal_elem_de.update(elem_task_num);
	//tbb::parallel_for(tbb::blocked_range<size_t>(0, elem_task_num, 1), self.cal_elem_de);
	ParaUtil::parallel_for(self.cal_elem_de, elem_task_num);

	// cal strain increment at node
	//tbb::parallel_for(tbb::blocked_range<size_t>(0, node_elem_task_num, 1), self.cal_node_de);
	ParaUtil::parallel_for(self.cal_node_de, node_elem_task_num);

	// map bg mesh back to pcl
	self.prev_valid_pcl_num = self.valid_pcl_num;
	self.map_mesh_to_pcl.update(pcl_task_num);
	self.map_mesh_to_pcl_res.pcl_num = 0;
	//tbb::parallel_reduce(tbb::blocked_range<size_t>(0, pcl_task_num, 1), self.map_mesh_to_pcl_res);
	ParaUtil::parallel_reduce(self.map_mesh_to_pcl, self.map_mesh_to_pcl_res, pcl_task_num);
	self.valid_pcl_num = self.map_mesh_to_pcl_res.pcl_num;

	self.continue_calculation();
	return 0;
}
