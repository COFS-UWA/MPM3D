#include "Tests_pcp.h"

#include "ItemArray.hpp"
#include "utils.h"
#include "Model_T2D_ME_p.h"
#include "Step_T2D_ME_p.h"
#include "ModelData_T2D_ME_p.h"
#include "TimeHistory_T2D_ME_p_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "QtApp_Prep_T2D_ME_p.h"

#include "test_simulations.h"

void test_t2d_me_p_rigid_circle_penetration(int argc, char** argv)
{
	Model_T2D_ME_p model;
	model.load_mesh_from_hdf5("..\\..\\Asset\\rect_test_rc_mesh.h5");
	model.init_search_grid(0.5, 0.5);

	ParticleGenerator2D<Model_T2D_ME_p> pcl_generator;
	pcl_generator.generate_pcls_in_grid_layout(Rect(0.0, 20.0, 0.0, 15.0), 0.2, 0.2);
	model.init_pcls(pcl_generator, 20.0);

	size_t pcl_num = model.get_pcl_num();
	Model_T2D_ME_p::Particle* pcls = model.get_pcls();
	MatModel::LinearElasticity* mms = model.add_LinearElasticity(pcl_num);
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		Model_T2D_ME_p::Particle& pcl = pcls[pcl_id];
		MatModel::LinearElasticity& mm = mms[pcl_id];
		mm.set_param(100.0, 0.0);
		pcl.set_mat_model(mm);
	}

	//model.init_rigid_circle(150.0, 10.0, 17.6, 2.5);
	//model.set_rigid_circle_velocity(0.0, -1.0, 0.0);

	// vx bc
	IndexArray vx_bc_pt_array(50);
	find_2d_nodes_on_x_line(model, vx_bc_pt_array, 0.0);
	find_2d_nodes_on_x_line(model, vx_bc_pt_array, 20.0, false);
	size_t* vx_bc_n_id = vx_bc_pt_array.get_mem();
	model.init_vxs(vx_bc_pt_array.get_num());
	size_t vx_num = model.get_vx_num();
	VelocityBC* vxs = model.get_vxs();
	for (size_t v_id = 0; v_id < vx_num; ++v_id)
	{
		VelocityBC& vbc = vxs[v_id];
		vbc.node_id = vx_bc_n_id[v_id];
		vbc.v = 0.0;
	}

	// vy bc
	IndexArray vy_bc_pt_array(50);
	find_2d_nodes_on_y_line(model, vy_bc_pt_array, 0.0);
	size_t* vy_bc_n_id = vy_bc_pt_array.get_mem();
	model.init_vys(vy_bc_pt_array.get_num());
	size_t vy_num = model.get_vy_num();
	VelocityBC* vys = model.get_vys();
	for (size_t v_id = 0; v_id < vy_num; ++v_id)
	{
		VelocityBC& vbc = vys[v_id];
		vbc.node_id = vy_bc_n_id[v_id];
		vbc.v = 0.0;
	}

	//QtApp_Prep_T2D_ME_p md_disp(argc, argv);
	//md_disp.set_win_size(900, 900);
	//md_disp.set_model(model);
	//md_disp.set_pts_from_node_id(vx_bc_pt_array.get_mem(), vx_bc_pt_array.get_num(), 0.2);
	////md_disp.set_pts_from_node_id(vy_bc_pt_array.get_mem(), vy_bc_pt_array.get_num(), 0.2);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_me_p_rigid_circle_penetration.h5");

	ModelData_T2D_ME_p md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T2D_ME_p_complete out1("penetration");
	out1.set_interval_num(100);
	out1.set_res_file(res_file_hdf5);
	out1.set_output_init_state();

	TimeHistory_ConsoleProgressBar out_pb;

	Step_T2D_ME_p step("step1");
	step.set_model(model);
	step.set_step_time(1.0); // 3.0
	step.set_dtime(5.0e-6);
	step.set_thread_num(6);
	step.add_time_history(out1);
	step.add_time_history(out_pb);
	step.solve();
}

#include "QtApp_Posp_T2D_ME_p.h"
#include "test_model_view.h"

void test_t2d_me_p_rigid_circle_penetration_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t2d_me_p_rigid_circle_penetration.h5");

	QtApp_Posp_T2D_ME_p app(argc, argv, QtApp_Posp_T2D_ME_p::Animation);
	app.set_win_size(900, 900);
	app.set_fld_range(0.0, 1.0);
	app.set_res_file(rf, "penetration", "s22");
	app.set_ani_time(5.0);
	//app.set_png_name("t2d_me_p_rigid_circle_penetration");
	//app.set_gif_name("t2d_me_p_rigid_circle_penetration");
	app.start();
}