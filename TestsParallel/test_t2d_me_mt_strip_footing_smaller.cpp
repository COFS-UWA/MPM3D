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

void test_t2d_me_mt_strip_footing_smaller(int argc, char** argv)
{
	TriangleMesh tri_mesh;
	tri_mesh.load_mesh_from_hdf5("../../Asset/smaller_rect_mesh.h5");

	Model_T2D_ME_mt model;
	model.init_mesh(tri_mesh);
	model.init_search_grid(tri_mesh, 0.25, 0.25);

	ParticleGenerator2D<TriangleMesh> pcl_generator;
	pcl_generator.generate_pcls_in_grid_layout(Rect(0.0, 20.0, 0.0, 15.0), 0.2, 0.2);
	model.init_pcls(pcl_generator, 20.0);
	MatModel::MaterialModel** mms = model.get_mat_models();
	MatModel::VonMises* vms = model.add_VonMises(model.get_pcl_num());
	for (size_t p_id = 0; p_id < model.get_pcl_num(); ++p_id)
	{
		vms[p_id].set_param(4000.0, 0.3, 10.0);
		mms[p_id] = &vms[p_id];
	}
	//MatModel::LinearElasticity *les = model.add_LinearElasticity(model.get_pcl_num());
	//for (size_t p_id = 0; p_id < model.get_pcl_num(); ++p_id)
	//{
	//	les[p_id].set_param(4000.0, 0.0);
	//	mms[p_id] = &les[p_id];
	//}

	model.init_rigid_rect(200.0, 10.0, 15.25, 4.0, 0.5);
	model.set_rigid_rect_velocity(0.0, -0.01, 0.0);

	// vx bc
	IndexArray vx_bc_pt_array(50);
	find_2d_nodes_on_x_line(model, vx_bc_pt_array, 0.0);
	find_2d_nodes_on_x_line(model, vx_bc_pt_array, 20.0, false);
	model.init_fixed_vx_bc(vx_bc_pt_array.get_num(), vx_bc_pt_array.get_mem());

	// vy bc
	IndexArray vy_bc_pt_array(50);
	find_2d_nodes_on_y_line(model, vy_bc_pt_array, 0.0);
	model.init_fixed_vy_bc(vy_bc_pt_array.get_num(), vy_bc_pt_array.get_mem());

	//QtApp_Prep_T2D_ME_mt md_disp(argc, argv);
	//md_disp.set_win_size(900, 900);
	//md_disp.set_model(model);
	//md_disp.set_pts_from_node_id(vx_bc_pt_array.get_mem(), vx_bc_pt_array.get_num(), 0.1);
	////md_disp.set_pts_from_node_id(vy_bc_pt_array.get_mem(), vy_bc_pt_array.get_num(), 0.1);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_me_mt_strip_footing_smaller.h5");

	ModelData_T2D_ME_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T2D_ME_mt_complete out1("loading");
	out1.set_output_init_state();
	out1.set_interval_num(100);
	out1.set_res_file(res_file_hdf5);
	TimeHistory_ConsoleProgressBar out_pb;

	Step_T2D_ME_mt step("step1");
	step.set_model(model);
	step.set_step_time(20.0);
	step.set_dtime(1.0e-5);
	step.set_thread_num(5);
	step.add_time_history(out1);
	step.add_time_history(out_pb);
	step.solve();
}

#include "QtApp_Posp_T2D_ME_mt.h"
#include "test_model_view_omp.h"

void test_t2d_me_mt_strip_footing_smaller_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t2d_me_mt_strip_footing_smaller.h5");

	QtApp_Posp_T2D_ME_mt app(argc, argv, QtApp_Posp_T2D_ME_mt::Animation);
	app.set_win_size(900, 900);
	app.set_ani_time(5.0);
	app.set_res_file(rf, "loading", Hdf5Field::s22);
	app.set_color_map_fld_range(-10.0, 0.0);
	//app.set_png_name("t2d_me_mt_strip_footing_smaller");
	//app.set_gif_name("t2d_me_mt_strip_footing_smaller");
	app.start();
}
