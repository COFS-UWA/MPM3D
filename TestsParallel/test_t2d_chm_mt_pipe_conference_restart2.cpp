#include "TestsParallel_pcp.h"

#include "Model_T2D_CHM_mt.h"
#include "Model_T2D_CHM_mt_hdf5_utilities.h"
#include "ModelData_T2D_CHM_mt.h"
#include "Step_T2D_CHM_mt.h"
#include "TimeHistory_T2D_CHM_mt_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "QtApp_Prep_T2D_CHM_mt.h"
#include "test_simulations_omp.h"

void test_t2d_chm_mt_pipe_conference_restart2(int argc, char** argv)
{
	Model_T2D_CHM_mt model;
	Step_T2D_CHM_mt step("step2");
	
	using Model_T2D_CHM_mt_hdf5_utilities::load_chm_mt_model_from_hdf5_file;
	load_chm_mt_model_from_hdf5_file(
		model, step,
		"t2d_chm_s_pipe_conference_restart1.h5",
		"penetration",
		464
		);

	//size_t pt_id_array[] = { 54337 };
	//QtApp_Prep_T2D_CHM_mt md_disp(argc, argv);
	//md_disp.set_win_size(900, 900);
	//md_disp.set_model(model);
	////md_disp.set_pts_from_node_id(left_right_bc_pt_array.get_mem(), left_right_bc_pt_array.get_num(), 0.05);
	////md_disp.set_pts_from_node_id(bottom_bc_pt_array.get_mem(), bottom_bc_pt_array.get_num(), 0.05);
	////md_disp.set_pts_from_pcl_id(mid_tbc_pt_array.get_mem(), mid_tbc_pt_array.get_num(), 0.01);
	////md_disp.set_pts_from_pcl_id(left_right_tbc_pt_array.get_mem(), left_right_tbc_pt_array.get_num(), 0.015);
	//md_disp.set_pts_from_pcl_id(pt_id_array, sizeof(pt_id_array) / sizeof(size_t), 0.015);
	//// all
	////md_disp.set_display_range(-3.6, 3.6, -5.1, 1.1);
	//// left
	////md_disp.set_display_range(-3.8, -2.2, -1.0, 1.0);
	//// middle
	//md_disp.set_display_range(-2.2, -1.8, -0.75, 0.25);
	//// right
	////md_disp.set_display_range(2.2, 3.8, -1.0, 1.0);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_chm_mt_pipe_conference_restart2.h5");

	ModelData_T2D_CHM_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T2D_CHM_mt_complete out("penetration");
	out.set_output_init_state();
	out.set_interval_num(500);
	out.set_res_file(res_file_hdf5);
	TimeHistory_ConsoleProgressBar out_pb;

	step.set_model(model);
	step.set_step_time(0.01);
	step.set_dtime(2.0e-6);
	step.add_time_history(out);
	step.add_time_history(out_pb);
	step.solve();
}

void test_t2d_chm_mt_pipe_conference_restart3(int argc, char** argv)
{
	Model_T2D_CHM_mt model;
	Step_T2D_CHM_mt step("step3");

	using Model_T2D_CHM_mt_hdf5_utilities::load_chm_mt_model_from_hdf5_file;
	load_chm_mt_model_from_hdf5_file(
		model,
		step,
		"t2d_chm_mt_pipe_conference_restart2.h5",
		"penetration",
		423
		);
	
	//size_t pt_id_array[] = { 409 }; // { 95, 8, 409 };
	//QtApp_Prep_T2D_CHM_mt md_disp(argc, argv);
	//md_disp.set_win_size(900, 900);
	//md_disp.set_model(model);
	////md_disp.set_pts_from_node_id(left_right_bc_pt_array.get_mem(), left_right_bc_pt_array.get_num(), 0.05);
	////md_disp.set_pts_from_node_id(bottom_bc_pt_array.get_mem(), bottom_bc_pt_array.get_num(), 0.05);
	////md_disp.set_pts_from_pcl_id(mid_tbc_pt_array.get_mem(), mid_tbc_pt_array.get_num(), 0.01);
	////md_disp.set_pts_from_pcl_id(left_right_tbc_pt_array.get_mem(), left_right_tbc_pt_array.get_num(), 0.015);
	//md_disp.set_pts_from_node_id(pt_id_array, sizeof(pt_id_array) / sizeof(size_t), 0.015);
	//// all
	////md_disp.set_display_range(-3.6, 3.6, -5.1, 1.1);
	//// left
	////md_disp.set_display_range(-3.8, -2.2, -1.0, 1.0);
	//// middle
	//md_disp.set_display_range(-2.2, -1.8, -0.75, 0.25);
	//// right
	////md_disp.set_display_range(2.2, 3.8, -1.0, 1.0);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_chm_mt_pipe_conference_restart3.h5");

	ModelData_T2D_CHM_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T2D_CHM_mt_complete out("penetration");
	out.set_output_init_state();
	out.set_interval_num(100);
	out.set_res_file(res_file_hdf5);
	TimeHistory_ConsoleProgressBar out_pb;

	step.set_model(model);
	step.set_step_time(2.0e-5);
	step.set_dtime(2.0e-6);
	step.add_time_history(out);
	//step.add_time_history(out_pb);
	step.solve();
}

#include "QtApp_Posp_T2D_CHM_mt.h"
#include "test_model_view_omp.h"

void test_t2d_chm_mt_pipe_conference_restart2_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t2d_chm_mt_pipe_conference_restart2.h5");

	QtApp_Posp_T2D_CHM_mt app(argc, argv, QtApp_Posp_T2D_CHM_mt::Animation);
	app.set_ani_time(5.0);
	app.set_res_file(rf, "penetration", Hdf5Field::s22);
	app.set_win_size(900, 900);
	app.set_display_range(-3.6, 3.6, -5.1, 0.6);
	app.set_color_map_fld_range(-11.0, -9.0);
	app.set_color_map_geometry(0.7f, 0.45f, 0.5f);
	//app.set_png_name("t2d_chm_s_pipe_conference_restart2");
	//app.set_gif_name("t2d_chm_s_pipe_conference_restart2");
	app.start();
}