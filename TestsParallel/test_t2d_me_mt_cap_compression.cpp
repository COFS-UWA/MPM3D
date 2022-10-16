#include "TestsParallel_pcp.h"

#include "test_parallel_utils.h"
#include "TriangleMesh.h"
#include "ParticleGenerator2D.hpp"
#include "Model_T2D_ME_mt.h"
#include "Model_T2D_ME_mt_hdf5_utilities.h"
#include "ModelData_T2D_ME_mt.h"
#include "Step_T2D_ME_mt.h"
#include "TimeHistory_T2D_ME_mt_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "QtApp_Prep_T2D_ME_mt.h"
#include "test_simulations_omp.h"

void test_t2d_me_mt_cap_compression(int argc, char** argv)
{
	TriangleMesh tri_mesh;
	tri_mesh.load_mesh_from_hdf5("../../Asset/rect_mesh.h5");
	//tri_mesh.load_mesh_from_hdf5("../../Asset/rect_mesh_005.h5");
	//tri_mesh.load_mesh_from_hdf5("../../Asset/rect_mesh_0025.h5");

	Model_T2D_ME_mt model;
	model.init_mesh(tri_mesh);
	model.init_search_grid(tri_mesh, 0.05, 0.05);
	
	ParticleGenerator2D<TriangleMesh> pcl_generator;
	pcl_generator.generate_pcls_in_grid_layout(Rect(0.0, 0.2, 0.0, 1.0), 0.02, 0.02);
	//pcl_generator.generate_pcls_in_grid_layout(Rect(0.0, 0.2, 0.0, 1.0), 0.015, 0.015);
	//pcl_generator.generate_pcls_in_grid_layout(Rect(0.0, 0.2, 0.0, 1.0), 0.01, 0.01);
	//pcl_generator.generate_pcls_in_grid_layout(Rect(0.0, 0.2, 0.0, 1.0), 0.005, 0.005);
	model.init_pcls(pcl_generator, 10.0);
	MatModel::MaterialModel** mms = model.get_mat_models();
	MatModel::LinearElasticity* les = model.add_LinearElasticity(model.get_pcl_num());
	for (uint32_t p_id = 0; p_id < model.get_pcl_num(); ++p_id)
	{
		mms[p_id] = les;
		les->set_param(1000.0, 0.0);
		les = model.following_LinearElasticity(les);
	}

	// rigid rect
	model.init_rigid_rect(0.1, 1.03, 0.3, 0.06, 1.0);
	model.set_rigid_rect_velocity(0.0, -0.01, 0.0);
	// rigid body
	////model.init_rb(2.0, "../../Asset/rect_mesh_w006_1by5.h5", -0.05, 1.0, 0.0, 0.04, 0.04);
	//model.init_rb(2.0, "C:\\MyData\\Work\\Contact algorithm 3D\\pap2_md\\rect_cap_pp2.h5", -0.05, 1.0, 0.0, 0.04, 0.04);
	//model.set_rb_velocity(0.0, -0.01, 0.0);
	// 5000.0, 10000.0, 20000.0, 50000.0
	model.set_contact_param(20000.0, 20000.0, 0.2, 3.0);

	// vx bc
	IndexArray vx_bc_pt_array(50);
	find_2d_nodes_on_x_line(model, vx_bc_pt_array, 0.0);
	find_2d_nodes_on_x_line(model, vx_bc_pt_array, 0.2, false);
	model.init_fixed_vx_bc(vx_bc_pt_array.get_num(), vx_bc_pt_array.get_mem());

	// vy bc
	IndexArray vy_bc_pt_array(50);
	find_2d_nodes_on_y_line(model, vy_bc_pt_array, 0.0);
	model.init_fixed_vy_bc(vy_bc_pt_array.get_num(), vy_bc_pt_array.get_mem());

	//QtApp_Prep_T2D_ME_mt md_disp(argc, argv);
	//md_disp.set_win_size(600, 950);
	//md_disp.set_model(model);
	//md_disp.set_bg_color(1.0f, 1.0f, 1.0f);
	////md_disp.set_pts_from_node_id(vx_bc_pt_array.get_mem(), vx_bc_pt_array.get_num(), 0.01);
	////md_disp.set_pts_from_node_id(vy_bc_pt_array.get_mem(), vy_bc_pt_array.get_num(), 0.01);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_me_mt_cap_compression.h5");

	ModelData_T2D_ME_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T2D_ME_mt_complete out("loading");
	out.set_interval_num(100);
	out.set_output_init_state();
	out.set_output_final_state();
	out.set_res_file(res_file_hdf5);
	TimeHistory_ConsoleProgressBar out_pb;

	Step_T2D_ME_mt step("step1");
	step.set_model(model);
	step.set_step_time(5.0);
	//step.set_step_time(5.0e-5);
	step.set_dtime(1.0e-5);
	//step.set_thread_num(4);
	step.add_time_history(out);
	step.add_time_history(out_pb);
	step.solve();
}

#include "QtApp_Posp_T2D_ME_mt.h"
#include "test_model_view_omp.h"

void test_t2d_me_mt_cap_compression_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t2d_me_mt_cap_compression.h5");

	//QtApp_Posp_T2D_ME_mt app(argc, argv, QtApp_Posp_T2D_ME_mt::Animation);
	//app.set_win_size(1200, 800);
	//app.set_ani_time(5.0);
	//app.set_res_file(rf, "loading", Hdf5Field::s22);
	//app.set_bg_color(1.0, 1.0, 1.0);
	//app.set_rb_color(0.92941, 0.49, 0.19216);
	//app.set_mesh_color(0.75, 0.75, 0.75);
	//app.set_color_map_fld_range(-50.0, 0.0);
	//app.set_color_map_geometry(1.0f, 0.45f, 0.5f);
	//app.set_color_map_char_color(0.0, 0.0, 0.0);
	////app.set_png_name("t2d_me_mt_cap_compression");
	//app.set_gif_name("t2d_me_mt_cap_compression");
	//app.start();

	QtApp_Posp_T2D_ME_mt app(argc, argv, QtApp_Posp_T2D_ME_mt::SingleFrame);
	app.set_win_size(1200, 800);
	app.set_res_file(rf, "loading", 100, Hdf5Field::s22);
	app.set_mono_color_pcl();
	app.set_pcl_color(0.26667, 0.44706, 0.76863);
	app.set_bg_color(1.0, 1.0, 1.0);
	app.set_rb_color(0.92941, 0.49, 0.19216);
	app.set_mesh_color(0.75, 0.75, 0.75);
	//app.set_png_name("t2d_me_mt_cap_compression");
	app.start();
}
