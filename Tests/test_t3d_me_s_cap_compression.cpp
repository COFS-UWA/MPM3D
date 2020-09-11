#include "Tests_pcp.h"

#include "ItemArray.hpp"
#include "ParticleGenerator3D.hpp"
#include "LinearElasticity.h"
#include "Model_T3D_ME_s.h"
#include "Step_T3D_ME_s.h"
#include "ModelData_T3D_ME_s.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "TimeHistory_T3D_ME_s_complete.h"
#include "QtApp_Prep_T3D_ME_s.h"

#include "utils.h"
#include "test_simulations.h"

void test_t3d_me_s_cap_compression(int argc, char **argv)
{
	Model_T3D_ME_s model;
	model.load_mesh_from_hdf5("../../Asset/brick_mesh_1.00_2x2x10.h5");
	model.init_search_grid(0.05, 0.05, 0.05);

	ParticleGenerator3D<Model_T3D_ME_s> pg;
	pg.generate_pcls_grid(Cube(0.0, 0.2, 0.0, 0.2, 0.0, 1.0), 0.025, 0.025, 0.025);
	model.init_pcls(pg, 10.0);
	size_t pcl_num = model.get_pcl_num();
	Model_T3D_ME_s::Particle *pcls = model.get_pcls();
	MatModel::LinearElasticity *mms = model.add_LinearElasticity(pcl_num);
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		Model_T3D_ME_s::Particle &pcl = pcls[pcl_id];
		MatModel::LinearElasticity &mm = mms[pcl_id];
		mm.set_param(100.0, 0.0);
		pcl.set_mat_model(mm);
	}

	model.init_rb(1.0, "../../Asset/square_cap_mesh.h5", -0.05, -0.05, 1.0);
	model.set_contact_params(10.0, 10.0, 10.0);
	RigidTetrahedronMesh& rb = model.get_rb();
	rb.init_bg_grids(0.05, 0.07); // 0.06, 0.075
	rb.set_dist_max(0.05); // must > particle radius
	rb.set_v_bc(0.0, 0.0, -0.05);

	IndexArray vx_bc_pt_array(100);
	find_3d_nodes_on_x_plane(model, vx_bc_pt_array, 0.0);
	find_3d_nodes_on_x_plane(model, vx_bc_pt_array, 0.2, false);
	size_t *vx_bc_n_id = vx_bc_pt_array.get_mem();
	model.init_vxs(vx_bc_pt_array.get_num());
	for (size_t v_id = 0; v_id < model.vx_num; ++v_id)
	{
		VelocityBC &vbc = model.vxs[v_id];
		vbc.node_id = vx_bc_n_id[v_id];
		vbc.v = 0.0;
	}

	IndexArray vy_bc_pt_array(100);
	find_3d_nodes_on_y_plane(model, vy_bc_pt_array, 0.0);
	find_3d_nodes_on_y_plane(model, vy_bc_pt_array, 0.2, false);
	size_t *vy_bc_n_id = vy_bc_pt_array.get_mem();
	model.init_vys(vy_bc_pt_array.get_num());
	for (size_t v_id = 0; v_id < model.vy_num; ++v_id)
	{
		VelocityBC &vbc = model.vys[v_id];
		vbc.node_id = vy_bc_n_id[v_id];
		vbc.v = 0.0;
	}
	
	IndexArray vz_bc_pt_array(100);
	find_3d_nodes_on_z_plane(model, vz_bc_pt_array, 0.0);
	size_t *vz_bc_n_id = vz_bc_pt_array.get_mem();
	model.init_vzs(vz_bc_pt_array.get_num());
	for (size_t v_id = 0; v_id < model.vz_num; ++v_id)
	{
		VelocityBC &vbc = model.vzs[v_id];
		vbc.node_id = vz_bc_n_id[v_id];
		vbc.v = 0.0;
	}

	//QtApp_Prep_T3D_ME_s md_disp(argc, argv);
	//md_disp.set_win_size(900, 900);
	////md_disp.set_view_dir(30.0, 30.0);
	//md_disp.set_view_dir(0.0, 0.0);
	//md_disp.set_light_dir(20.0, 30.0);
	//md_disp.set_view_dist_scale(1.5);
	//md_disp.set_model(model);
	////md_disp.set_pts_from_node_id(vx_bc_pt_array.get_mem(), vx_bc_pt_array.get_num(), 0.01);
	////md_disp.set_pts_from_node_id(vy_bc_pt_array.get_mem(), vy_bc_pt_array.get_num(), 0.01);
	////md_disp.set_pts_from_node_id(vz_bc_pt_array.get_mem(), vz_bc_pt_array.get_num(), 0.01);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_me_s_cap_compression.h5");

	ModelData_T3D_ME_s md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T3D_ME_s_complete out1("compression");
	out1.set_res_file(res_file_hdf5);
	out1.set_output_init_state();
	out1.set_interval_num(50);
	TimeHistory_ConsoleProgressBar out_cpb;

	Step_T3D_ME_s step("step1");
	step.set_model(model);
	step.set_step_time(1.0);
	step.set_dtime(1.0e-5);
	step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();
}

#include "QtApp_Posp_T3D_ME_s.h"
#include "test_model_view.h"

void test_t3d_me_s_cap_compression_result(int argc, char **argv)
{
	ResultFile_hdf5 rf;
	rf.open("t3d_me_s_cap_compression.h5");
	
	//QtApp_Posp_T3D_ME_s app(argc, argv, QtApp_Posp_T3D_ME_s::SingleFrame);
	//app.set_res_file(rf, "compression", 2, Hdf5Field::z);
	QtApp_Posp_T3D_ME_s app(argc, argv, QtApp_Posp_T3D_ME_s::Animation);
	app.set_res_file(rf, "compression", Hdf5Field::s33);
	app.set_ani_time(5.0);
	app.set_win_size(900, 900);
	app.set_view_dir(30.0f, 0.0f);
	app.set_light_dir(35.0f, -50.0f);
	app.set_view_dist_scale(1.5);
	app.set_color_map_fld_range(-6.0, 0.0);
	app.set_color_map_geometry(0.7f, 0.45f, 0.5f);
	//app.set_png_name("t3d_me_s_cap_compression");
	app.set_gif_name("t3d_me_s_cap_compression");
	app.start();
}
