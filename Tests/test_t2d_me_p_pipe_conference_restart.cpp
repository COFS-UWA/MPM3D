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

	//DisplayModel_T2D disp_model;
	//disp_model.init_win();
	//disp_model.init_model(model);
	//disp_model.init_rigid_circle(model.get_rigid_circle());
	//disp_model.init_points(pt_array.get_mem(), pt_array.get_num() / 3);
	//// all
	//disp_model.display(-3.6, 3.6, -5.1, 1.1);
	//// left
	////disp_model.display(-3.8, -2.2, -1.0, 1.0);
	//// middle
	////disp_model.display(2.3, 2.7, -0.25, 0.25);
	//// right
	////disp_model.display(2.2, 3.8, -1.0, 1.0);
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