#include "Tests_pcp.h"

#include "ItemArray.hpp"
#include "utils.h"
#include "ModifiedCamClay.h"
#include "Model_T2D_CHM_s.h"
#include "ModelData_T2D_CHM_s.h"
#include "Model_T2D_CHM_s_hdf5_utilities.h"
#include "Step_T2D_CHM_s_Geo.h"
#include "TimeHistory_T2D_CHM_s_Geo_complete.h"
#include "Step_T2D_CHM_s.h"
#include "TimeHistory_T2D_CHM_s_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"

#include "QtApp_Prep_T2D_CHM_s.h"

#include "test_simulations.h"

void test_t2d_chm_s_test_rigid_circle(int argc, char** argv)
{
	Model_T2D_CHM_s model;
	model.load_mesh_from_hdf5("..\\..\\Asset\\rect_test_rc_mesh.h5");
	model.init_search_grid(0.5, 0.5);

	ParticleGenerator2D<Model_T2D_CHM_s> pcl_generator;
	pcl_generator.generate_pcls_in_grid_layout(Rect(0.0, 20.0, 0.0, 15.0), 0.2, 0.2);
	model.init_pcls(pcl_generator, 0.6, 2650.0, 1000.0, 5.0e6, 5.0e-12, 1.0e-3);
	size_t pcl_num = model.get_pcl_num();
	Model_T2D_CHM_s::Particle* pcls = model.get_pcls();
	std::cout << "pcl num: " << pcl_num << "\n";

	// elasticity
	//MatModel::LinearElasticity* mms = model.add_LinearElasticity(pcl_num);
	//for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	//{
	//	Model_T2D_CHM_s::Particle& pcl = pcls[pcl_id];
	//	MatModel::LinearElasticity& mm = mms[pcl_id];
	//	mm.set_param(2.0e6, 0.3);
	//	pcl.set_mat_model(mm);
	//}
	// mcc
	MatModel::ModifiedCamClay *mms = model.add_ModifiedCamClay(pcl_num);
	double ini_stress[6] = { -40000.0, -24050.0, -24050.0, 0.0, 0.0, 0.0 };
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		Model_T2D_CHM_s::Particle& pcl = pcls[p_id];
		MatModel::ModifiedCamClay &mcc = mms[p_id];
		mcc.set_param_OC(0.3, 0.044, 0.205, 23.5, 3.6677, ini_stress, 39610.0);
		pcl.set_mat_model(mcc);
	}

	model.init_rigid_circle(1.0e5, 1.0e5, 10.0, 17.5, 2.5);
	model.set_rigid_circle_velocity(0.0, -0.25, 0.0);

	IndexArray vx_bc_pt_array(50);
	find_2d_nodes_on_x_line(model, vx_bc_pt_array, 0.0);
	find_2d_nodes_on_x_line(model, vx_bc_pt_array, 20.0, false);
	size_t* vx_bc_n_id = vx_bc_pt_array.get_mem();
	model.init_vsxs(vx_bc_pt_array.get_num());
	size_t vx_num = model.get_vsx_num();
	VelocityBC* vxs = model.get_vsxs();
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
	model.init_vsys(vy_bc_pt_array.get_num());
	size_t vy_num = model.get_vsy_num();
	VelocityBC* vys = model.get_vsys();
	for (size_t v_id = 0; v_id < vy_num; ++v_id)
	{
		VelocityBC& vbc = vys[v_id];
		vbc.node_id = vy_bc_n_id[v_id];
		vbc.v = 0.0;
	}

	QtApp_Prep_T2D_CHM_s md_disp(argc, argv);
	md_disp.set_win_size(900, 900);
	md_disp.set_model(model);
	md_disp.set_pts_from_node_id(vx_bc_pt_array.get_mem(), vx_bc_pt_array.get_num(), 0.01);
	//md_disp.set_pts_from_node_id(vy_bc_pt_array.get_mem(), vy_bc_pt_array.get_num(), 0.01);
	//md_disp.set_pts_from_pcl_id(tbc_pt_array.get_mem(), tbc_pt_array.get_num(), 0.01);
	md_disp.start();
	return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_chm_s_test_rigid_circle.h5");

	ModelData_T2D_CHM_s md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T2D_CHM_s_complete out1("circle_penetration");
	out1.set_res_file(res_file_hdf5);
	out1.set_interval_num(100);
	out1.set_output_init_state();

	TimeHistory_ConsoleProgressBar out_pb;

	Step_T2D_CHM_s step("rigid_circle");
	step.set_model(model);
	step.set_step_time(3.0);
	step.set_dtime(5.0e-6);
	step.add_time_history(out1);
	step.add_time_history(out_pb);
	step.solve();
}

#include "QtApp_Posp_T2D_CHM_s.h"
#include "test_model_view.h"

void test_t2d_chm_s_test_rigid_circle_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t2d_chm_s_test_rigid_circle.h5");

	QtApp_Posp_T2D_CHM_s app(argc, argv);
	app.set_win_size(900, 900);
	app.set_fld_range(0.0, 1.0);
	app.set_res_file(rf, "circle_penetration", 0, "y");
	//app.set_png_name("t2d_chm_s_test_rigid_circle");
	app.start();
}