#include "Tests_pcp.h"

#include "ItemArray.hpp"
#include "utils.h"
#include "ModifiedCamClay.h"
#include "Model_T2D_CHM_s.h"
#include "Model_T2D_CHM_s_hdf5_utilities.h"
#include "ModelData_T2D_CHM_s.h"
#include "Step_T2D_CHM_s.h"
#include "TimeHistory_T2D_CHM_s_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "QtApp_Prep_T2D_CHM_s.h"

#include "test_simulations.h"

void test_t2d_chm_s_pipe_conference_restart2(int argc, char** argv)
{
	Model_T2D_CHM_s model;

	using Model_T2D_CHM_s_hdf5_utilities::load_CHM_s_model_from_hdf5_file;
	load_CHM_s_model_from_hdf5_file(
		model,
		"t2d_chm_s_pipe_conference_restart1.h5",
		"penetration",
		501
		);

	//QtApp_Prep_T2D_CHM_s md_disp(argc, argv);
	//md_disp.set_win_size(900, 900);
	//md_disp.set_model(model);
	////md_disp.set_pts_from_node_id(left_right_bc_pt_array.get_mem(), left_right_bc_pt_array.get_num(), 0.05);
	////md_disp.set_pts_from_node_id(bottom_bc_pt_array.get_mem(), bottom_bc_pt_array.get_num(), 0.05);
	////md_disp.set_pts_from_pcl_id(mid_tbc_pt_array.get_mem(), mid_tbc_pt_array.get_num(), 0.01);
	////md_disp.set_pts_from_pcl_id(left_right_tbc_pt_array.get_mem(), left_right_tbc_pt_array.get_num(), 0.015);
	//// all
	//md_disp.set_display_range(-3.6, 3.6, -5.1, 1.1);
	//// left
	////md_disp.set_display_range(-3.8, -2.2, -1.0, 1.0);
	//// middle
	////md_disp.set_display_range(-1.5, 1.5, -0.75, 0.25);
	//// right
	////md_disp.set_display_range(2.2, 3.8, -1.0, 1.0);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_chm_s_pipe_conference_restart2.h5");

	ModelData_T2D_CHM_s md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T2D_CHM_s_complete out("penetration");
	out.set_res_file(res_file_hdf5);
	out.set_output_init_state();
	TimeHistory_ConsoleProgressBar out_pb;

	Step_T2D_CHM_s step("step2");
	step.set_model(model);
	step.set_step_time(5.0); // 1.0s
	step.set_dtime(2.0e-6); // 2.0e-7
	out.set_interval_num(500);
	step.add_time_history(out);
	step.add_time_history(out_pb);
	step.solve();
}

#include "QtApp_Posp_T2D_CHM_s.h"
#include "test_model_view.h"

void test_t2d_chm_s_pipe_conference_restart2_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t2d_chm_s_pipe_conference_restart2.h5");

	QtApp_Posp_T2D_CHM_s app(argc, argv, QtApp_Posp_T2D_CHM_s::Animation);
	app.set_win_size(900, 900);
	app.set_ani_time(5.0);
	app.set_res_file(rf, "penetration", Hdf5Field::s22);
	app.set_display_range(-3.6, 3.6, -5.1, 0.6);
	app.set_color_map_fld_range(-11.0, -9.0);
	app.set_color_map_geometry(0.7f, 0.45f, 0.5f);
	//app.set_png_name("t2d_chm_s_pipe_conference_restart2");
	//app.set_gif_name("t2d_chm_s_pipe_conference_restart2");
	app.start();
}