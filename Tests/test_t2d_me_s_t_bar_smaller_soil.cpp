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

void test_t2d_me_s_t_bar_smaller_soil(int argc, char** argv)
{
	Model_T2D_ME_s model;
	model.load_mesh_from_hdf5("../../Asset/smaller_rect_mesh.h5");
	model.init_search_grid(0.25, 0.25);
	
	ParticleGenerator2D<Model_T2D_ME_s> pcl_generator;
	pcl_generator.generate_pcls_in_grid_layout(Rect(0.0, 20.0, 0.0, 15.0), 0.2, 0.2);
		
	// elasticity
	//model.init_pcls(pcl_generator, 0.4, 20.0, 10.0, 1000.0, 1.0e-5, 1.0);
	//size_t pcl_num = model.get_pcl_num();
	//Model_T2D_ME_s::Particle* pcls = model.get_pcls();
	//MatModel::LinearElasticity* mms = model.add_LinearElasticity(pcl_num);
	//for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	//{
	//	Model_T2D_ME_s::Particle& pcl = pcls[pcl_id];
	//	MatModel::LinearElasticity& mm = mms[pcl_id];
	//	mm.set_param(100.0, 0.0);
	//	pcl.set_mat_model(mm);
	//}
	//model.init_rigid_circle(100.0, 10.0, 10.0, 17.5, 2.5);
	// mcc
	model.init_pcls(pcl_generator, 2650.0);
	size_t pcl_num = model.get_pcl_num();
	Model_T2D_ME_s::Particle* pcls = model.get_pcls();
	MatModel::ModifiedCamClay* mms = model.add_ModifiedCamClay(pcl_num);
	double ini_stress[6] = { -24050.0, -40000.0, -24050.0, 0.0, 0.0, 0.0 };
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		Model_T2D_ME_s::Particle &pcl = pcls[p_id];
		MatModel::ModifiedCamClay &mm = mms[p_id];
		mm.set_param_OC(0.3, 0.044, 0.205, 23.5, 3.6677, ini_stress, 39610.0);
		pcl.set_mat_model(mm);
	}
	model.init_rigid_circle(1.0e5, 1.0e3, 10.0, 17.41, 2.5);

	model.set_rigid_circle_velocity(0.0, -0.25, 0.0);

	// vx
	IndexArray vx_bc_pt_array(50);
	find_2d_nodes_on_x_line(model, vx_bc_pt_array, 0.0);
	find_2d_nodes_on_x_line(model, vx_bc_pt_array, 20.0, false);
	size_t* vx_bc_n_id = vx_bc_pt_array.get_mem();
	model.init_vxs(vx_bc_pt_array.get_num());
	size_t vsx_num = model.get_vx_num();
	VelocityBC* vsxs = model.get_vxs();
	for (size_t v_id = 0; v_id < vsx_num; ++v_id)
	{
		VelocityBC& vbc = vsxs[v_id];
		vbc.node_id = vx_bc_n_id[v_id];
		vbc.v = 0.0;
	}

	// vy
	IndexArray vy_bc_pt_array(50);
	find_2d_nodes_on_y_line(model, vy_bc_pt_array, 0.0);
	size_t* vy_bc_n_id = vy_bc_pt_array.get_mem();
	model.init_vys(vy_bc_pt_array.get_num());
	size_t vsy_num = model.get_vy_num();
	VelocityBC* vsys = model.get_vys();
	for (size_t v_id = 0; v_id < vsy_num; ++v_id)
	{
		VelocityBC& vbc = vsys[v_id];
		vbc.node_id = vy_bc_n_id[v_id];
		vbc.v = 0.0;
	}

	//QtApp_Prep_T2D_ME_s md_disp(argc, argv);
	//md_disp.set_win_size(900, 900);
	//md_disp.set_model(model);
	////md_disp.set_pts_from_node_id(vx_bc_pt_array.get_mem(), vx_bc_pt_array.get_num(), 0.1);
	////md_disp.set_pts_from_node_id(vy_bc_pt_array.get_mem(), vy_bc_pt_array.get_num(), 0.1);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file;
	res_file.create("t2d_me_s_t_bar_smaller_soil.h5");

	ModelData_T2D_ME_s md("md1");
	md.output_model(model, res_file);

	TimeHistory_T2D_ME_s_complete out("penetration");
	out.set_output_init_state();
	out.set_interval_num(100);
	out.set_res_file(res_file);
	
	TimeHistory_ConsoleProgressBar out_pb;

	Step_T2D_ME_s step("step1");
	step.set_model(model);
	step.set_step_time(6.0);
	step.set_dtime(3.0e-6);
	step.add_time_history(out);
	step.add_time_history(out_pb);
	step.solve();
}

#include "QtApp_Posp_T2D_ME_s.h"
#include "test_model_view.h"

void test_t2d_me_s_t_bar_smaller_soil_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t2d_me_s_t_bar_smaller_soil.h5");

	QtApp_Posp_T2D_ME_s app(argc, argv, QtApp_Posp_T2D_ME_s::Animation);
	app.set_ani_time(5.0);
	app.set_res_file(rf, "penetration", Hdf5Field::s22);
	app.set_win_size(900, 900);
	app.set_color_map_fld_range(0.0, 10.0);
	app.set_color_map_geometry(0.6, 0.5, 0.4);
	//app.set_png_name("t2d_me_s_t_bar_smaller_soil");
	//app.set_gif_name("t2d_me_s_t_bar_smaller_soil");
	app.start();
}