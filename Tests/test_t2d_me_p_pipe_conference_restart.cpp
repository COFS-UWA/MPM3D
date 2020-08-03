#include "Tests_pcp.h"

#include "ItemArray.hpp"
#include "utils.h"
#include "ModifiedCamClay.h"
#include "Model_T2D_ME_p.h"
#include "Model_T2D_ME_p_hdf5_utilities.h"
#include "ModelData_T2D_ME_p.h"
#include "Step_T2D_ME_p.h"
#include "TimeHistory_T2D_ME_p_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "QtApp_Prep_T2D_ME_p.h"
#include "test_simulations.h"

void test_t2d_me_p_pipe_conference_restart(int argc, char** argv)
{
	Model_T2D_ME_p model;

	using Model_T2D_ME_p_hdf5_utilities::load_me_s_model_from_hdf5_file;
	load_me_s_model_from_hdf5_file(
		model,
		"t2d_me_p_pipe_conference_geo.h5",
		"geostatic",
		11
		);

	model.set_rigid_circle_velocity(0.0, -0.05, 0.0);

	//IndexArray left_right_bc_pt_array;
	//left_right_bc_pt_array.reserve(model.get_vx_num());
	//for (size_t v_id = 0; v_id < model.get_vx_num(); ++v_id)
	//	left_right_bc_pt_array.add(model.get_vxs()[v_id].node_id);

	//IndexArray bottom_bc_pt_array;
	//bottom_bc_pt_array.reserve(model.get_vy_num());
	//for (size_t v_id = 0; v_id < model.get_vy_num(); ++v_id)
	//	bottom_bc_pt_array.add(model.get_vys()[v_id].node_id);

	//IndexArray mid_tbc_pt_array, left_right_tbc_pt_array;
	//mid_tbc_pt_array.reserve(100);
	//left_right_tbc_pt_array.reserve(100);
	//for (size_t t_id = 0; t_id < model.get_ty_num(); ++t_id)
	//{
	//	if (model.get_tys()[t_id].t > -500)
	//		mid_tbc_pt_array.add(model.get_tys()[t_id].pcl_id);
	//	else
	//		left_right_tbc_pt_array.add(model.get_tys()[t_id].pcl_id);
	//}

	//QtApp_Prep_T2D_ME_p md_disp(argc, argv);
	//md_disp.set_win_size(900, 900);
	//md_disp.set_model(model);
	////md_disp.set_pts_from_node_id(left_right_bc_pt_array.get_mem(), left_right_bc_pt_array.get_num(), 0.05);
	////md_disp.set_pts_from_node_id(bottom_bc_pt_array.get_mem(), bottom_bc_pt_array.get_num(), 0.05);
	//md_disp.set_pts_from_pcl_id(mid_tbc_pt_array.get_mem(), mid_tbc_pt_array.get_num(), 0.01);
	////md_disp.set_pts_from_pcl_id(left_right_tbc_pt_array.get_mem(), left_right_tbc_pt_array.get_num(), 0.015);
	//// all
	////md_disp.set_display_range(-3.6, 3.6, -5.1, 1.1);
	//// left
	////md_disp.set_display_range(-3.8, -2.2, -1.0, 1.0);
	//// middle
	////md_disp.set_display_range(-1.5, 1.5, -0.75, 0.25);
	//// right
	//md_disp.set_display_range(2.2, 3.8, -1.0, 1.0);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_me_p_pipe_conference_restart.h5");

	ModelData_T2D_ME_p md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T2D_ME_p_complete out("penetration");
	out.set_res_file(res_file_hdf5);
	out.set_output_init_state();
	out.set_interval_num(500);
	TimeHistory_ConsoleProgressBar out_pb;

	Step_T2D_ME_p step("step1");
	step.set_model(model);
	step.set_step_time(5.0e-5);
	step.set_dtime(2.0e-6);
	step.set_thread_num(1);
	step.add_time_history(out);
	step.add_time_history(out_pb);
	step.solve();
}

#include "QtApp_Posp_T2D_ME_p.h"
#include "test_model_view.h"

void test_t2d_me_p_pipe_conference_restart_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t2d_me_p_pipe_conference_restart.h5");

	QtApp_Posp_T2D_ME_p app(argc, argv, QtApp_Posp_T2D_ME_p::Animation);
	app.set_win_size(900, 900);
	app.set_res_file(rf, "penetration", "s22");
	app.set_ani_time(5.0);
	app.set_display_range(-3.6, 3.6, -5.1, 0.6);
	app.set_fld_range(-20010.0, -19990.0);
	app.set_color_map_pos(0.6, 0.45, 0.5);
	//app.set_png_name("t2d_me_p_pipe_conference_restart");
	//app.set_gif_name("t2d_me_p_pipe_conference_restart");
	app.start();
}