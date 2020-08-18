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

void test_t2d_me_s_geostatic_restart(int argc, char** argv)
{
	Model_T2D_ME_s model;

	using Model_T2D_ME_s_hdf5_utilities::load_me_s_model_from_hdf5_file;
	load_me_s_model_from_hdf5_file(
		model,
		"t2d_me_s_geostatic.h5",
		"geostatic",
		101
		);

	for (size_t pcl_id = 0; pcl_id < model.get_pcl_num(); ++pcl_id)
	{
		Model_T2D_ME_s::Particle& pcl = model.get_pcls()[pcl_id];
		static_cast<MatModel::UndrainedModifiedCamClay*>(pcl.mm)->set_Kw(1.0e8);
	}

	IndexArray vx_bc_pt_array;
	for (size_t n_id = 0; n_id < model.get_vx_num(); ++n_id)
		vx_bc_pt_array.add(model.get_vxs()[n_id].node_id);

	IndexArray vy_bc_pt_array;
	for (size_t n_id = 0; n_id < model.get_vy_num(); ++n_id)
		vy_bc_pt_array.add(model.get_vys()[n_id].node_id);

	IndexArray tbc_pt_array(100);
	for (size_t t_id = 0; t_id < model.get_ty_num(); ++t_id)
		tbc_pt_array.add(model.get_tys()[t_id].pcl_id);

	//QtApp_Prep_T2D_ME_s md_disp(argc, argv);
	//md_disp.set_win_size(900, 900);
	//md_disp.set_model(model);
	////md_disp.set_pts_from_node_id(vx_bc_pt_array.get_mem(), vx_bc_pt_array.get_num(), 0.01);
	////md_disp.set_pts_from_node_id(vy_bc_pt_array.get_mem(), vy_bc_pt_array.get_num(), 0.01);
	//md_disp.set_pts_from_pcl_id(tbc_pt_array.get_mem(), tbc_pt_array.get_num(), 0.01);
	//md_disp.start();
	//return;

	TimeHistory_ConsoleProgressBar out3;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_me_s_geostatic_restart.h5");

	ModelData_T2D_ME_s md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T2D_ME_s_complete out1("restart");
	out1.set_res_file(res_file_hdf5);
	out1.set_interval_num(100);
	out1.set_output_init_state();

	TimeHistory_ConsoleProgressBar out_pb;

	Step_T2D_ME_s step("step1");
	step.set_model(model);
	step.set_step_time(1.0);
	step.set_dtime(1.0e-6);
	step.add_time_history(out1);
	step.add_time_history(out_pb);
	step.solve();
}

#include "QtApp_Posp_T2D_ME_s.h"
#include "test_model_view.h"

void test_t2d_me_s_geostatic_restart_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t2d_me_s_geostatic_restart.h5");

	// single frame
	//QtApp_Posp_T2D_ME_s app(argc, argv);
	//app.set_win_size(900, 900);
	//app.set_res_file(rf, "geostatic", 0, "s22");
	//app.set_fld_range(-11.0, -9.0);
	//app.set_color_map_pos(0.7, 0.45, 0.5); // color map legend
	////app.set_png_name("t2d_me_s_geostatic_restart");
	//app.start();

	// animation
	QtApp_Posp_T2D_ME_s app(argc, argv, QtApp_Posp_T2D_ME_s::Animation);
	app.set_win_size(900, 900);
	app.set_res_file(rf, "restart", "s22");
	app.set_ani_time(5.0);
	app.set_fld_range(-20100.0, -19900.0);
	app.set_color_map_pos(0.7, 0.45, 0.5); // color map legend
	//app.set_png_name("t2d_me_s_geostatic_restart");
	app.set_gif_name("t2d_me_s_geostatic_restart");
	app.start();
}
