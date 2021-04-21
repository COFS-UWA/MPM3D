#include "TestsParallel_pcp.h"

#include <iomanip>

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

void test_t2d_me_mt_strip_footing(int argc, char** argv)
{
	TriangleMesh tri_mesh;
	//tri_mesh.load_mesh_from_hdf5("../../Asset/strip_footing_soil_mesh.h5");
	//model.init_search_grid(tri_mesh, 0.1, 0.1);
	tri_mesh.load_mesh_from_hdf5("../../Asset/strip_footing_soil_mesh_denser.h5");
	tri_mesh.init_search_grid(0.0333333, 0.0333333);

	Model_T2D_ME_mt model;
	model.init_mesh(tri_mesh);
	//model.init_search_grid(tri_mesh, 0.1, 0.1);
	model.init_search_grid(tri_mesh, 0.0333333, 0.0333333);

	ParticleGenerator2D<TriangleMesh> pcl_generator;
	//pcl_generator.generate_pcls_in_grid_layout(Rect(-5.0, 5.0, -3.5, 0.0), 0.025, 0.025);
	//pcl_generator.generate_pcls_in_grid_layout(Rect(-5.0, 5.0, -5.0, -3.5), 0.07, 0.07);
	pcl_generator.generate_pcls_in_grid_layout(Rect(-5.0, 5.0, -3.5, 0.0), 0.01, 0.01);
	pcl_generator.generate_pcls_in_grid_layout(Rect(-5.0, 5.0, -5.0, -3.5), 0.03, 0.03);
	pcl_generator.adjust_pcl_size_to_fit_elems(tri_mesh);
	model.init_pcls(pcl_generator, 20.0);
	std::cout << model.get_pcl_num() << "\n";
	MatModel::MaterialModel** mms = model.get_mat_models();
	//MatModel::VonMises* vms = model.add_VonMises(model.get_pcl_num());
	//for (szie_t p_id = 0; p_id < model.get_pcl_num(); ++p_id)
	//{
	//	vms[p_id].set_param(4000.0, 0.3, 10.0);
	//	mms[p_id] = &vms[p_id];
	//}
	MatModel::Tresca* tes = model.add_Tresca(model.get_pcl_num());
	for (size_t p_id = 0; p_id < model.get_pcl_num(); ++p_id)
	{
		tes[p_id].set_param(4000.0, 0.3, 5.0);
		mms[p_id] = &tes[p_id];
	}

	model.init_rigid_rect(0.0, 0.1, 1.0, 0.2, 1.0);
	model.set_rigid_rect_velocity(0.0, -0.01, 0.0);
	model.set_contact_param(20000.0, 20000.0, 0.2, 1.5);
	//model.set_frictional_contact_between_pcl_and_rect();
	//model.set_rough_contact_between_pcl_and_rect();

	// vx bc
	IndexArray vx_bc_pt_array(50);
	find_2d_nodes_on_x_line(model, vx_bc_pt_array, -5.0);
	find_2d_nodes_on_x_line(model, vx_bc_pt_array,  5.0, false);
	model.init_fixed_vx_bc(vx_bc_pt_array.get_num(), vx_bc_pt_array.get_mem());

	// vy bc
	IndexArray vy_bc_pt_array(50);
	find_2d_nodes_on_y_line(model, vy_bc_pt_array, -5.0);
	model.init_fixed_vy_bc(vy_bc_pt_array.get_num(), vy_bc_pt_array.get_mem());

	//QtApp_Prep_T2D_ME_mt md_disp(argc, argv);
	//md_disp.set_win_size(1500, 950);
	//md_disp.set_model(model);
	////md_disp.set_display_range(-1.0, 1.0, -1.5, -0.5);
	////md_disp.set_pts_from_node_id(vx_bc_pt_array.get_mem(), vx_bc_pt_array.get_num(), 0.02);
	//md_disp.set_pts_from_node_id(vy_bc_pt_array.get_mem(), vy_bc_pt_array.get_num(), 0.02);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_me_mt_strip_footing.h5");

	ModelData_T2D_ME_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T2D_ME_mt_complete out1("loading");
	out1.set_output_init_state();
	out1.set_interval_num(100);
	out1.set_res_file(res_file_hdf5);
	TimeHistory_ConsoleProgressBar out_pb;
	out_pb.set_interval_num(2000);

	Step_T2D_ME_mt step("step1");
	step.set_model(model);
	step.set_step_time(10.0); // 10.0
	step.set_dtime(5.0e-6);
	step.set_thread_num(20);
	step.add_time_history(out1);
	step.add_time_history(out_pb);
	step.solve();
}

#include "QtApp_Posp_T2D_ME_mt.h"
#include "test_model_view_omp.h"

void test_t2d_me_mt_strip_footing_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t2d_me_mt_strip_footing.h5");

	QtApp_Posp_T2D_ME_mt app(argc, argv, QtApp_Posp_T2D_ME_mt::Animation);
	app.set_win_size(1600, 1000);
	app.set_ani_time(8.0);
	//app.set_display_range(-2.5, 2.5, -2.4, 0.1);
	app.set_res_file(rf, "loading", Hdf5Field::s22);
	app.set_color_map_fld_range(-20.0, 0.0);
	//app.set_res_file(rf, "loading", Hdf5Field::plastic_mises_strain_2d);
	//app.set_color_map_fld_range(0.0, 0.01);
	//app.set_res_file(rf, "loading", Hdf5Field::vx);
	//app.set_color_map_fld_range(-0.003, 0.003);
	//app.set_png_name("t2d_me_mt_strip_footing");
	app.set_gif_name("t2d_me_mt_strip_footing");
	app.start();
}
