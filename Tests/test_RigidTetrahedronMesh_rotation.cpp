#include "Tests_pcp.h"

#include <iostream>
#include "ItemArray.hpp"
#include "Model_T3D_ME_s.h"
#include "Step_T3D_ME_s.h"
#include "ModelData_T3D_ME_s.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "TimeHistory_T3D_ME_s_complete.h"
#include "QtApp_Prep_T3D_ME_s.h"

#include "utils.h"
#include "test_simulations.h"

void test_RigidTetrahedronMesh_rotation(int argc, char** argv)
{
	Model_T3D_ME_s model;
	model.load_mesh_from_hdf5("../../Asset/brick_mesh_1.00_1x1x1.h5");

	model.init_rb(1.0, "../../Asset/ball_mesh_r1.00.h5", 0.0, 0.0, 2.0);
	RigidTetrahedronMesh& rb = model.get_rb();
	rb.set_vx_bc(0.0);
	rb.set_vy_bc(0.0);
	rb.set_vz_bc(0.0);
	rb.add_f_ext(1.0, 0.0, 0.0, 0.0, 0.0, 3.0);

	QtApp_Prep_T3D_ME_s md_disp(argc, argv);
	md_disp.set_win_size(900, 900);
	md_disp.set_view_dir(30.0, 30.0);
	//md_disp.set_view_dist_scale(4.0);
	md_disp.set_light_dir(90.0, 30.0);
	//md_disp.set_rb_display_mode(QtRigidTetrahedronMeshGLObject::LineFrame);
	md_disp.set_model(model);
	md_disp.start();
	return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_me_s_ball_rotation.h5");

	ModelData_T3D_ME_s md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T3D_ME_s_complete out1("rotation");
	out1.set_res_file(res_file_hdf5);
	out1.set_output_init_state();
	out1.set_interval_num(50);
	TimeHistory_ConsoleProgressBar out_cpb;

	Step_T3D_ME_s step("step1");
	step.set_model(model);
	step.set_step_time(10.0);
	step.set_dtime(1.0e-5);
	step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();
}

#include "QtApp_Posp_T3D_ME_s.h"
#include "test_model_view.h"

void test_RigidTetrahedronMesh_rotation_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t3d_me_s_ball_rotation.h5");

	//QtApp_Posp_T3D_ME_s app(argc, argv);
	//app.set_res_file(rf, "rotation", 0, Hdf5Field::s33);
	QtApp_Posp_T3D_ME_s app(argc, argv, QtApp_Posp_T3D_ME_s::Animation);
	app.set_res_file(rf, "rotation", Hdf5Field::s33);
	app.set_ani_time(5.0);
	app.set_win_size(900, 900);
	app.set_view_dir(30.0f, 30.0f);
	//app.set_view_dist_scale(3.0);
	app.set_light_dir(90.0f, 30.0f);
	//app.set_color_map_fld_range(-0.01, 0.0);
	//app.set_color_map_geometry(0.7f, 0.45f, 0.5f);
	//app.set_rb_display_mode(QtRigidTetrahedronMeshGLObject::LineFrame);
	//app.set_png_name("t3d_me_s_ball_rotation");
	app.set_gif_name("t3d_me_s_ball_rotation");
	app.start();
}