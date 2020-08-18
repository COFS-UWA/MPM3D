#include "Tests_pcp.h"

#include "ItemArray.hpp"
#include "utils.h"
#include "Model_T2D_ME_s.h"
#include "ModelData_T2D_ME_s.h"
#include "Step_T2D_ME_s_Geo.h"
#include "TimeHistory_T2D_ME_s_Geo_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "QtApp_Prep_T2D_ME_s.h"

#include "test_simulations.h"

void test_t2d_me_s_geostatic(int argc, char** argv)
{
	Model_T2D_ME_s model;
	model.load_mesh_from_hdf5("..\\..\\Asset\\rect_mesh.h5");
	model.init_search_grid(0.05, 0.05);

	ParticleGenerator2D<Model_T2D_ME_s> pcl_generator;
	pcl_generator.generate_pcls_in_grid_layout(Rect(0.0, 0.2, 0.0, 1.0), 0.02, 0.02);
	pcl_generator.adjust_pcl_size_to_fit_elems(model);
	model.init_pcls(pcl_generator, 10.0);
	//model.sum_vol_for_all_elements();

	size_t pcl_num = model.get_pcl_num();
	Model_T2D_ME_s::Particle* pcls = model.get_pcls();
	// linear elasticity
	//double K = 0.0;
	//MatModel::LinearElasticity* mms = model.add_LinearElasticity(pcl_num);
	//for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	//{
	//	Model_T2D_ME_s::Particle& pcl = pcls[pcl_id];
	//	MatModel::LinearElasticity& mm = mms[pcl_id];
	//	mm.set_param(100000.0, 0.0);
	//	pcl.set_mat_model(mm);
	//	pcl.s22 = -20000.0;
	//	pcl.s11 = K * pcl.s22;
	//	pcl.s12 = 0.0;
	//}
	// mcc
	//MatModel::ModifiedCamClay* mms = model.add_ModifiedCamClay(pcl_num);
	//double K = 1.0 - sin(23.5 / 180.0 * 3.14159165359);
	//double ini_stress[6] = { -12025.0, -20000.0, -12025.0, 0.0, 0.0, 0.0 };
	//for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	//{
	//	Model_T2D_ME_s::Particle& pcl = pcls[p_id];
	//	pcl.s11 = ini_stress[0];
	//	pcl.s22 = ini_stress[1];
	//	pcl.s12 = 0.0;
	//	MatModel::ModifiedCamClay& mm = mms[p_id];
	//	mm.set_param_NC(0.3, 0.044, 0.205, 23.5, 3.6677, ini_stress);
	//  pcl.set_mat_model(mm);
	//}
	// undrained mcc
	MatModel::UndrainedModifiedCamClay* mms = model.add_UndrainedModifiedCamClay(pcl_num);
	double K = 1.0 - sin(23.5 / 180.0 * 3.14159165359);
	double ini_stress[6] = { -12025.0, -20000.0, -12025.0, 0.0, 0.0, 0.0 };
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		Model_T2D_ME_s::Particle& pcl = pcls[p_id];
		pcl.s11 = ini_stress[0];
		pcl.s22 = ini_stress[1];
		pcl.s12 = 0.0;
		MatModel::UndrainedModifiedCamClay& mm = mms[p_id];
		mm.set_param_NC(0.3, 0.044, 0.205, 23.5, 3.6677, ini_stress, 0.0, 0.0);
		pcl.set_mat_model(mm);
	}

	// vx bc
	IndexArray vx_bc_pt_array(50);
	find_2d_nodes_on_x_line(model, vx_bc_pt_array, 0.0);
	find_2d_nodes_on_x_line(model, vx_bc_pt_array, 0.2, false);
	size_t* vx_bc_n_id = vx_bc_pt_array.get_mem();
	model.init_vxs(vx_bc_pt_array.get_num());
	size_t vx_num = model.get_vx_num();
	VelocityBC* vxs = model.get_vxs();
	for (size_t v_id = 0; v_id < vx_num; ++v_id)
	{
		VelocityBC& vbc = vxs[v_id];
		vbc.node_id = vx_bc_n_id[v_id];
		vbc.v = 0.0;
	}

	// vy bc
	IndexArray vy_bc_pt_array(50);
	find_2d_nodes_on_y_line(model, vy_bc_pt_array, 0.0);
	size_t* vy_bc_n_id = vy_bc_pt_array.get_mem();
	model.init_vys(vy_bc_pt_array.get_num());
	size_t vy_num = model.get_vy_num();
	VelocityBC* vys = model.get_vys();
	for (size_t v_id = 0; v_id < vy_num; ++v_id)
	{
		VelocityBC& vbc = vys[v_id];
		vbc.node_id = vy_bc_n_id[v_id];
		vbc.v = 0.0;
	}

	// traction bc
	IndexArray tbc_pt_array(50);
	find_2d_pcls(model, tbc_pt_array, Rect(0.0, 0.2, 0.987, 1.0));
	size_t* tbc_pcl_id = tbc_pt_array.get_mem();
	model.init_tys(tbc_pt_array.get_num());
	size_t ty_num = model.get_ty_num();
	TractionBCAtPcl* tys = model.get_tys();
	for (size_t t_id = 0; t_id < ty_num; ++t_id)
	{
		TractionBCAtPcl& tbc = tys[t_id];
		tbc.pcl_id = tbc_pcl_id[t_id];
		tbc.t = 0.02 * -20000.0;
	}

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
	res_file_hdf5.create("t2d_me_s_geostatic.h5");

	ModelData_T2D_ME_s md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T2D_ME_s_Geo_complete out1("geostatic");
	out1.set_res_file(res_file_hdf5);
	out1.set_interval_num(100);
	out1.set_output_init_state();
	TimeHistory_ConsoleProgressBar out_pb;

	Step_T2D_ME_s_Geo step("step1");
	step.set_model(model);
	step.set_step_time(1.0);
	step.set_dtime(2.5e-6);
	//step.set_damping_ratio(0.1);
	step.add_time_history(out1);
	step.add_time_history(out_pb);
	step.solve();
}

#include "QtApp_Posp_T2D_ME_s.h"
#include "test_model_view.h"

void test_t2d_me_s_geostatic_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t2d_me_s_geostatic.h5");

	// single frame
	//QtApp_Posp_T2D_ME_s app(argc, argv);
	//app.set_win_size(900, 900);
	//app.set_res_file(rf, "geostatic", 0, "s22");
	//app.set_fld_range(-11.0, -9.0);
	//app.set_color_map_pos(0.7, 0.45, 0.5); // color map legend
	////app.set_png_name("t2d_me_s_geostatic");
	//app.start();

	// animation
	QtApp_Posp_T2D_ME_s app(argc, argv, QtApp_Posp_T2D_ME_s::Animation);
	app.set_win_size(900, 900);
	app.set_res_file(rf, "geostatic", "s22");
	app.set_ani_time(5.0);
	app.set_fld_range(-20100.0, -19900.0);
	app.set_color_map_pos(0.7, 0.45, 0.5); // color map legend
	//app.set_png_name("t2d_me_s_geostatic");
	app.set_gif_name("t2d_me_s_geostatic");
	app.start();
}
