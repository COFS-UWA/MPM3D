#include "TestsParallel_pcp.h"

#include "ParticleGenerator3D.hpp"
#include "Model_T3D_ME_mt.h"
#include "ModelData_T3D_ME_mt.h"
#include "QtApp_Prep_T3D_ME_mt.h"
#include "Step_T3D_ME_TBB.h"
#include "TimeHistory_T3D_ME_tbb_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "test_parallel_utils.h"
#include "test_simulations_omp.h"

void test_t3d_me_tbb_piezofoundation(int argc, char** argv)
{
	Model_T3D_ME_mt model;
	Step_T3D_ME_TBB step("step2");
	Model_T3D_ME_mt_hdf5_utilities::load_me_mt_model_from_hdf5_file(
		model, step, "t3d_me_mt_piezofoundation_geo.h5", "geostatic", 21); // 21
	std::cout << "Load model completed.\n";

	QtApp_Prep_T3D_ME_mt md_disp(argc, argv);
	md_disp.set_model(model);
	md_disp.set_win_size(1200, 950);
	md_disp.set_view_dir(-90.0f, -20.0f);
	md_disp.set_light_dir(-60.0f, 15.0f);
	md_disp.set_display_bg_mesh(false);
	//md_disp.set_view_dist_scale(0.6);
	md_disp.set_pts_from_vx_bc(0.05);
	//md_disp.set_pts_from_vy_bc(0.05);
	//md_disp.set_pts_from_vz_bc(0.05);
	md_disp.start();
	return;

	model.set_rigid_cylinder_velocity(0.0, 0.0, -0.5);
	constexpr double sml_pcl_size = 0.03125;
	model.set_contact_param(20000.0 / (sml_pcl_size * sml_pcl_size),
		20000.0 / (sml_pcl_size * sml_pcl_size), 0.1, 5.0);
	
	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_me_tbb_piezofoundation.h5");

	ModelData_T3D_ME_mt md;
	md.output_model(model, res_file_hdf5);
	std::cout << "Output model completed.\n";

	TimeHistory_T3D_ME_TBB_complete out1("penetration");
	out1.set_interval_num(100);
	out1.set_output_init_state();
	out1.set_output_final_state();
	out1.set_res_file(res_file_hdf5);
	TimeHistory_ConsoleProgressBar out_cpb;
	out_cpb.set_interval_num(2000);

	std::cout << "Start solving...\n";
	step.set_thread_num(22);
	step.set_step_time(0.9); // 0.3 when v=-1.5 0.6=0.45D
	//step.set_thread_num(4);
	//step.set_step_time(6.0e-6);
	step.set_dtime(2.0e-6);
	step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();
}

#include "QtApp_Posp_T3D_ME_mt_Div.h"
#include "test_model_view_omp.h"

void test_t3d_me_tbb_piezofoundation_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t3d_me_tbb_piezofoundation.h5");

	// Single frame
	//QtApp_Posp_T3D_ME_mt_Div<PlaneDivisionSet> app(argc, argv, QtApp_Posp_T3D_ME_mt_Div<PlaneDivisionSet>::SingleFrame);
	//app.set_res_file(rf, "geostatic", 21, Hdf5Field::s33);
	//app.set_color_map_fld_range(-8.0e4, 0.0);
	//app.set_color_map_geometry(1.2f, 0.4f, 0.45f);
	//app.set_win_size(1600, 950);
	//app.set_view_dir(-90.0f, 20.0f);
	//app.set_light_dir(-135.0f, 20.0f);
	//app.set_light_dist_scale(1.0f);
	//app.set_fog_coef(0.02f);
	//app.set_view_dist_scale(0.85);
	//app.set_display_bg_mesh(false);
	//app.start();

	//// Animation
	//QtApp_Posp_T3D_ME_mt_Div<PlaneDivisionSet> app(argc, argv,
	//	QtApp_Posp_T3D_ME_mt_Div<PlaneDivisionSet>::Animation);
	//app.get_div_set().set_by_normal_and_point(0.0, 1.0, 0.0, 0.0, 0.06, 0.0);
	//app.set_ani_time(5.0);
	//app.set_win_size(1600, 950);
	//app.set_view_dir(-90.0f, 10.0f);
	//app.set_light_dir(-135.0f, 20.0f);
	//app.set_light_dist_scale(1.0f);
	//app.set_fog_coef(0.02f);
	//app.set_view_dist_scale(0.85);
	//app.set_display_bg_mesh(false);
	////app.set_mono_color_pcl(true);
	//// s33
	////app.set_res_file(rf, "penetration", Hdf5Field::s33);
	////app.set_color_map_fld_range(-8.0e4, 0.0);
	//// shear stress
	////app.set_res_file(rf, "penetration", Hdf5Field::max_shear_stress);
	////app.set_color_map_fld_range(0.0, 5.0);
	//// mises strain
	////app.set_res_file(rf, "penetration", Hdf5Field::plastic_mises_strain_2d);
	////app.set_color_map_fld_range(0.0, 0.08);
	//// mat_e
	//app.set_res_file(rf, "penetration", Hdf5Field::mat_e);
	//app.set_color_map_fld_range(0.6, 0.8);
	//// mat_s11
	////app.set_res_file(rf, "penetration", Hdf5Field::mat_s11);
	////app.set_color_map_fld_range(-4.0e4, 0.0);
	//// mat_s33
	////app.set_res_file(rf, "penetration", Hdf5Field::mat_s33);
	////app.set_color_map_fld_range(-1.0e5, 0.0);
	////
	//app.set_color_map_geometry(1.2f, 0.4f, 0.45f);
	////app.set_png_name("t3d_me_mt_piezofoundation");
	//app.set_gif_name("t3d_me_mt_piezofoundation");

	QtApp_Posp_T3D_ME_mt_Div<PlaneDivisionSet> app(argc, argv,
		QtApp_Posp_T3D_ME_mt_Div<PlaneDivisionSet>::SingleFrame);
	app.get_div_set().set_by_normal_and_point(0.0, 1.0, 0.0, 0.0, 0.06, 0.0);
	app.set_win_size(1600, 950);
	app.set_view_dir(-90.0f, 10.0f);
	app.set_light_dir(-135.0f, 20.0f);
	app.set_light_dist_scale(1.0f);
	app.set_fog_coef(0.02f);
	app.set_view_dist_scale(0.85);
	app.set_display_bg_mesh(false);
	app.set_color_map_geometry(1.2f, 0.4f, 0.45f);
	// mat_e
	app.set_res_file(rf, "penetration", 101, Hdf5Field::mat_e);
	app.set_color_map_fld_range(0.6, 0.8);
	//
	app.set_png_name("t3d_me_tbb_piezofoundation");

	app.start();
}
