#include "Tests_pcp.h"

#include "ItemArray.hpp"
#include "utils.h"
#include "ModifiedCamClay.h"
#include "Model_T2D_ME_s.h"
#include "Model_T2D_ME_s_hdf5_utilities.h"
#include "ModelData_T2D_ME_s.h"
#include "Step_T2D_ME_s.h"
#include "TimeHistory_T2D_ME_s_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "QtApp_Prep_T2D_ME_s.h"

#include "test_simulations.h"

void test_t2d_me_s_pipe_conference_restart2(int argc, char** argv)
{
	Model_T2D_ME_s model;

	using Model_T2D_ME_s_hdf5_utilities::load_me_s_model_from_hdf5_file;
	load_me_s_model_from_hdf5_file(
		model,
		"t2d_me_s_pipe_conference_restart.h5",
		"penetration",
		500
		);

	for (size_t pcl_id = 0; pcl_id < model.get_pcl_num(); ++pcl_id)
	{
		Model_T2D_ME_s::Particle& pcl = model.get_pcls()[pcl_id];
		static_cast<MatModel::UndrainedModifiedCamClay*>(pcl.mm)->set_Kw(2.0e7);
	}
	
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

	//QtApp_Prep_T2D_ME_s md_disp(argc, argv);
	//md_disp.set_win_size(900, 900);
	//md_disp.set_model(model);
	//md_disp.set_pts_from_node_id(left_right_bc_pt_array.get_mem(), left_right_bc_pt_array.get_num(), 0.05);
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
	res_file_hdf5.create("t2d_me_s_pipe_conference_restart2.h5");

	ModelData_T2D_ME_s md("md1");
	md.set_model(model);
	md.set_res_file(res_file_hdf5);
	md.output();

	TimeHistory_T2D_ME_s_complete out("penetration");
	out.set_res_file(res_file_hdf5);
	out.set_interval_num(500);
	out.set_output_init_state();
	TimeHistory_ConsoleProgressBar out_pb;

	Step_T2D_ME_s step("step2");
	step.set_model(model);
	step.set_step_time(6.0);
	step.set_dtime(2.0e-6);
	step.set_damping_ratio(0.15);
	step.add_time_history(out);
	step.add_time_history(out_pb);
	step.solve();
}

#include "QtApp_Posp_T2D_ME_s.h"
#include "test_model_view.h"

void test_t2d_me_s_pipe_conference_restart_result2(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t2d_me_s_pipe_conference_restart2.h5");

	//// single frame
	//QtApp_Posp_T2D_ME_s app(argc, argv);
	//app.set_win_size(900, 900);
	//app.set_res_file(rf, "penetration", 100, "s22");
	//app.set_display_range(-3.6, 3.6, -5.1, 1.1);
	//app.set_fld_range(-20100.0, -19900.0);
	//app.set_color_map_pos(0.8, 0.65, 0.3);
	////app.set_png_name("t2d_me_s_pipe_conference_restart2");
	//app.start();

	// animation
	QtApp_Posp_T2D_ME_s app(argc, argv, QtApp_Posp_T2D_ME_s::Animation);
	app.set_win_size(900, 900);
	app.set_res_file(rf, "penetration", "s22");
	app.set_ani_time(5.0);
	app.set_display_range(-3.6, 3.6, -5.1, 0.6);
	app.set_fld_range(-30000.0, -10000.0);
	app.set_color_map_pos(0.8, 0.65, 0.3);
	//app.set_png_name("t2d_me_s_pipe_conference_restart2");
	app.set_gif_name("t2d_me_s_pipe_conference_restart2");
	app.start();
}