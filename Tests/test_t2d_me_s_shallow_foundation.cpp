#include "Tests_pcp.h"

#include "ItemArray.hpp"
#include "utils.h"
#include "Model_T2D_ME_s.h"
#include "Step_T2D_ME_s.h"
#include "ModelData_T2D_ME_s.h"
#include "TimeHistory_T2D_ME_s_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "QtApp_Prep_T2D_ME_s.h"

#include "test_simulations.h"

void test_t2d_me_s_shallow_foundation(int argc, char** argv)
{
	Model_T2D_ME_s model;
	model.load_mesh_from_hdf5("../../Asset/rect_pipe_conference_mesh.h5");
	model.init_search_grid(0.05, 0.05);

	ParticleGenerator2D<Model_T2D_ME_s> pcl_generator;
	pcl_generator.generate_pcls_in_grid_layout(Rect(-3.5, 3.5, -3.5, 0.0), 0.04, 0.04);
	pcl_generator.generate_pcls_in_grid_layout(Rect(-3.5, 3.5, -5.0, -3.5), 0.04, 0.04);
	pcl_generator.replace_with_pcls_in_grid_layout(Rect(-2.5, 2.5, -3.5, 0.0), 0.02, 0.02);
	pcl_generator.adjust_pcl_size_to_fit_elems(model);
	model.init_pcls(pcl_generator, 10.0);

	size_t pcl_num = model.get_pcl_num();
	Model_T2D_ME_s::Particle* pcls = model.get_pcls();
	MatModel::VonMises* vms = model.add_VonMises(pcl_num);
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		Model_T2D_ME_s::Particle& pcl = pcls[pcl_id];
		MatModel::VonMises &mm = vms[pcl_id];
		mm.set_param(1000.0, 0.2, 5.0);
		pcl.set_mat_model(mm);
	}

	model.init_rigid_rect(200.0, 0.0, 0.15, 1.0, 0.3);
	model.set_rigid_rect_velocity(0.0, -0.03, 0.0);

	IndexArray left_right_bc_pt_array(50);
	find_2d_nodes_on_x_line(model, left_right_bc_pt_array, -3.5);
	find_2d_nodes_on_x_line(model, left_right_bc_pt_array, 3.5, false);
	size_t* left_right_bc_n_id = left_right_bc_pt_array.get_mem();
	model.init_vxs(left_right_bc_pt_array.get_num());
	size_t vsx_num = model.get_vx_num();
	VelocityBC* vsxs = model.get_vxs();
	for (size_t v_id = 0; v_id < vsx_num; ++v_id)
	{
		VelocityBC& vbc = vsxs[v_id];
		vbc.node_id = left_right_bc_n_id[v_id];
		vbc.v = 0.0;
	}

	IndexArray bottom_bc_pt_array(50);
	find_2d_nodes_on_y_line(model, bottom_bc_pt_array, -5.0);
	size_t* bottom_bc_n_ids = bottom_bc_pt_array.get_mem();
	model.init_vys(bottom_bc_pt_array.get_num());
	size_t vsy_num = model.get_vy_num();
	VelocityBC* vsys = model.get_vys();
	for (size_t v_id = 0; v_id < vsy_num; ++v_id)
	{
		VelocityBC& vbc = vsys[v_id];
		vbc.node_id = bottom_bc_n_ids[v_id];
		vbc.v = 0.0;
	}

	//QtApp_Prep_T2D_ME_s md_disp(argc, argv);
	//md_disp.set_win_size(900, 900);
	//md_disp.set_model(model);
	//md_disp.set_pts_from_node_id(left_right_bc_pt_array.get_mem(), left_right_bc_pt_array.get_num(), 0.05);
	////md_disp.set_pts_from_node_id(bottom_bc_pt_array.get_mem(), bottom_bc_pt_array.get_num(), 0.05);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_me_s_shallow_foundation.h5");

	ModelData_T2D_ME_s md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T2D_ME_s_complete out1("loading");
	out1.set_output_init_state();
	out1.set_interval_num(100);
	out1.set_res_file(res_file_hdf5);

	TimeHistory_ConsoleProgressBar out_pb;

	Step_T2D_ME_s step("step1");
	step.set_model(model);
	step.set_step_time(5.0);
	step.set_dtime(2.0e-6);
	step.add_time_history(out1);
	step.add_time_history(out_pb);
	step.solve();
}

#include "QtApp_Posp_T2D_ME_s.h"
#include "test_model_view.h"

void test_t2d_me_s_shallow_foundation_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t2d_me_s_shallow_foundation.h5");

	QtApp_Posp_T2D_ME_s app(argc, argv, QtApp_Posp_T2D_ME_s::Animation);
	app.set_ani_time(5.0);
	app.set_res_file(rf, "loading", Hdf5Field::mises_strain_2d);
	app.set_win_size(900, 900);
	app.set_color_map_fld_range(0.0, 1.0e-3);
	//app.set_png_name("t2d_me_s_shallow_foundation");
	//app.set_gif_name("t2d_me_s_shallow_foundation");
	app.start();
}