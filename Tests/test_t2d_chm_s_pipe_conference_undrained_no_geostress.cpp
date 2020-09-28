#include "Tests_pcp.h"

#include "ItemArray.hpp"
#include "utils.h"
#include "ModifiedCamClay.h"
#include "Model_T2D_CHM_s.h"
#include "ModelData_T2D_CHM_s.h"
#include "Step_T2D_CHM_s_ud.h"
#include "TimeHistory_T2D_CHM_s_ud_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "QtApp_Prep_T2D_ME_s.h"

#include "test_simulations.h"

void test_t2d_chm_s_pipe_conference_undrained_no_geostress(int argc, char** argv)
{
	Model_T2D_CHM_s model;
	model.load_mesh_from_hdf5("../../Asset/rect_pipe_conference_mesh.h5");
	model.init_search_grid(0.05, 0.05);

	ParticleGenerator2D<Model_T2D_CHM_s> pcl_generator;
	pcl_generator.generate_pcls_in_grid_layout(Rect(-3.5, 3.5, -3.5, 0.0), 0.04, 0.04);
	pcl_generator.generate_pcls_in_grid_layout(Rect(-3.5, 3.5, -5.0, -3.5), 0.04, 0.04);
	pcl_generator.replace_with_pcls_in_grid_layout(Rect(-2.5, 2.5, -3.5, 0.0), 0.02, 0.02);
	pcl_generator.adjust_pcl_size_to_fit_elems(model);
	model.init_pcls(pcl_generator, 0.6, 2650.0, 1000.0, 2.0e6, 1.0e-11, 1.0e-3);

	size_t pcl_num = model.get_pcl_num();
	Model_T2D_CHM_s::Particle* pcls = model.get_pcls();
	// mcc
	MatModel::ModifiedCamClay* mms = model.add_ModifiedCamClay(pcl_num);
	double K = 1.0 - sin(23.5 / 180.0 * 3.14159165359);
	double ini_stress[6] = { -12025.0, -20000.0, -12025.0, 0.0, 0.0, 0.0 };
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		Model_T2D_CHM_s::Particle& pcl = pcls[p_id];
		MatModel::ModifiedCamClay& mm = mms[p_id];
		mm.set_param_NC(0.3, 0.044, 0.205, 23.5, 3.6677, ini_stress);
		pcl.set_mat_model(mm);
	}

	model.init_rigid_circle(1.0e5, 1.0e3, 0.0, 0.5, 0.5);
	model.set_rigid_circle_velocity(0.0, -0.05, 0.0);

	IndexArray left_right_bc_pt_array(50);
	find_2d_nodes_on_x_line(model, left_right_bc_pt_array, -3.5);
	find_2d_nodes_on_x_line(model, left_right_bc_pt_array, 3.5, false);
	size_t* left_right_bc_n_id = left_right_bc_pt_array.get_mem();
	model.init_vsxs(left_right_bc_pt_array.get_num());
	size_t vsx_num = model.get_vsx_num();
	VelocityBC* vsxs = model.get_vsxs();
	for (size_t v_id = 0; v_id < vsx_num; ++v_id)
	{
		VelocityBC& vbc = vsxs[v_id];
		vbc.node_id = left_right_bc_n_id[v_id];
		vbc.v = 0.0;
	}

	IndexArray bottom_bc_pt_array(50);
	find_2d_nodes_on_y_line(model, bottom_bc_pt_array, -5.0);
	size_t* bottom_bc_n_ids = bottom_bc_pt_array.get_mem();
	model.init_vsys(bottom_bc_pt_array.get_num());
	size_t vsy_num = model.get_vsy_num();
	VelocityBC* vsys = model.get_vsys();
	for (size_t v_id = 0; v_id < vsy_num; ++v_id)
	{
		VelocityBC& vbc = vsys[v_id];
		vbc.node_id = bottom_bc_n_ids[v_id];
		vbc.v = 0.0;
	}

	//QtApp_Prep_T2D_ME_s md_disp(argc, argv);
	//md_disp.set_win_size(900, 900);
	//md_disp.set_model(model);
	////md_disp.set_pts_from_node_id(left_right_bc_pt_array.get_mem(), left_right_bc_pt_array.get_num(), 0.05);
	////md_disp.set_pts_from_node_id(bottom_bc_pt_array.get_mem(), bottom_bc_pt_array.get_num(), 0.05);
	//// all
	////md_disp.set_display_range(-3.6, 3.6, -5.1, 1.1);
	//// left
	//md_disp.set_display_range(-3.8, -2.2, -1.0, 1.0);
	//// middle
	////md_disp.set_display_range(-1.5, 1.5, -0.75, 0.25);
	//// right
	////md_disp.set_display_range(2.2, 3.8, -1.0, 1.0);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_chm_s_pipe_conference_undrained_no_geostress.h5");

	ModelData_T2D_CHM_s md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T2D_CHM_s_ud_complete out("penetration");
	out.set_res_file(res_file_hdf5);
	out.set_interval_num(100);
	out.set_output_init_state();
	TimeHistory_ConsoleProgressBar out_pb;

	Step_T2D_CHM_s_ud step("step1");
	step.set_model(model);
	step.set_step_time(5.0);
	step.set_dtime(2.0e-6);
	step.add_time_history(out);
	step.add_time_history(out_pb);
	step.solve();
}

#include "QtApp_Posp_T2D_CHM_s.h"
#include "test_model_view.h"

void test_t2d_chm_s_pipe_conference_undrained_no_geostress_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t2d_chm_s_pipe_conference_undrained_no_geostress.h5");

	// single frame
	//QtApp_Posp_T2D_CHM_s app(argc, argv);
	//app.set_win_size(900, 900);
	//app.set_res_file(rf, "geostatic", 100, Hdf5Field::s22);
	//app.set_display_range(-3.6, 3.6, -5.1, 1.1);
	//app.set_color_map_fld_range(-20100.0, -19900.0);
	//app.set_color_map_geometry(0.8, 0.65, 0.3);
	////app.set_png_name("t2d_me_s_pipe_conference_geo");
	//app.start();

	// animation
	QtApp_Posp_T2D_CHM_s app(argc, argv, QtApp_Posp_T2D_CHM_s::Animation);
	app.set_win_size(900, 900);
	app.set_res_file(rf, "penetration", Hdf5Field::s22);
	app.set_ani_time(5.0);
	//app.set_res_file(rf, "penetration", Hdf5Field::s22);
	//app.set_res_file(rf, "penetration", Hdf5Field::p);
	app.set_res_file(rf, "penetration", Hdf5Field::mises_strain_2d);
	app.set_ani_time(5.0);
	//app.set_display_range(-3.6, 3.6, -5.1, 0.6);
	app.set_display_range(-1.5, 1.5, -2.0, 0.5);
	//app.set_color_map_fld_range(-10000.0, 10000.0);
	app.set_color_map_fld_range(0.0, 2.0);
	app.set_color_map_geometry(0.8, 0.65, 0.3);
	//app.set_png_name("t2d_chm_s_pipe_conference_undrained_no_geostress");
	//app.set_gif_name("t2d_chm_s_pipe_conference_undrained_no_geostress");
	app.start();
}