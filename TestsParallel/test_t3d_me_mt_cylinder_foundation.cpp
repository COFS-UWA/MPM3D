#include "TestsParallel_pcp.h"

#include "ParticleGenerator3D.hpp"
#include "LinearElasticity.h"
#include "Model_T3D_ME_mt.h"
#include "Step_T3D_ME_mt.h"
#include "ModelData_T3D_ME_mt.h"
#include "TimeHistory_T3D_ME_mt_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "QtApp_Prep_T3D_ME_mt_Div.h"
#include "test_parallel_utils.h"
#include "test_simulations_omp.h"

void test_t3d_me_mt_cylinder_foundation_create_model(int argc, char** argv)
{
	TetrahedronMesh teh_mesh;
	teh_mesh.load_mesh_from_hdf5("../../Asset/cube_block.h5");
	teh_mesh.init_search_grid(0.05, 0.05, 0.05);

	Model_T3D_ME_mt model;
	model.init_mesh(teh_mesh);
	model.init_search_grid(teh_mesh);

	ParticleGenerator3D<TetrahedronMesh> pcl_generator;
	pcl_generator.generate_pcls_grid(Cube(0.0, 7.0, 0.0, 7.0, 0.0, 4.0), 0.05, 0.05, 0.05);
	pcl_generator.adjust_pcl_size_to_fit_elems(teh_mesh);
	model.init_pcls(pcl_generator, 10.0);
	size_t pcl_num = model.get_pcl_num();
	std::cout << "pcl_num: " << pcl_num << "\n";
	MatModel::MaterialModel** mms = model.get_mat_models();
	MatModel::LinearElasticity* les = model.add_LinearElasticity(pcl_num);
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		MatModel::LinearElasticity& le = les[pcl_id];
		le.set_param(1000.0, 0.0);
		mms[pcl_id] = &le;
	}

	model.init_rigid_cylinder(3.5, 3.5, 4.1, 0.2, 0.5);
	model.set_rigid_cylinder_velocity(0.0, 0.0, -0.1);
	model.set_contact_param(20000.0, 20000.0, 0.1);

	IndexArray vx_bc_pt_array(100);
	find_3d_nodes_on_x_plane(model, vx_bc_pt_array, 0.0);
	find_3d_nodes_on_x_plane(model, vx_bc_pt_array, 7.0, false);
	model.init_fixed_vx_bc(vx_bc_pt_array.get_num(), vx_bc_pt_array.get_mem());

	IndexArray vy_bc_pt_array(100);
	find_3d_nodes_on_y_plane(model, vy_bc_pt_array, 0.0);
	find_3d_nodes_on_y_plane(model, vy_bc_pt_array, 7.0, false);
	model.init_fixed_vy_bc(vy_bc_pt_array.get_num(), vy_bc_pt_array.get_mem());

	IndexArray vz_bc_pt_array(100);
	find_3d_nodes_on_z_plane(model, vz_bc_pt_array, 0.0);
	model.init_fixed_vz_bc(vz_bc_pt_array.get_num(), vz_bc_pt_array.get_mem());

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_me_mt_cylinder_foundation_model.h5");

	ModelData_T3D_ME_mt md;
	md.output_model(model, res_file_hdf5);
}

#include "Model_T3D_ME_mt_hdf5_utilities.h"

void test_t3d_me_mt_cylinder_foundation(int argc, char** argv)
{
	Model_T3D_ME_mt model;
	Model_T3D_ME_mt_hdf5_utilities::load_me_mt_model_from_hdf5_file(model, "t3d_me_mt_cylinder_foundation_model.h5");

	QtApp_Prep_T3D_ME_mt_Div<TwoPlaneDivisionSet> md_disp(argc, argv);
	auto &div_set = md_disp.get_div_set();
	div_set.seta().set_by_normal_and_point(0.0, 1.0, 0.0, 3.5, 3.5, 0.0);
	div_set.setb().set_by_normal_and_point(1.0, 0.0, 0.0, 3.5, 3.5, 0.0);
	md_disp.set_win_size(1200, 950);
	md_disp.set_view_dir(30.0f, 30.0f);
	md_disp.set_light_dir(90.0f, 30.0f);
	md_disp.set_model(model);
	md_disp.set_display_bg_mesh(false);
	//md_disp.set_pts_from_node_id(vx_bc_pt_array.get_mem(), vx_bc_pt_array.get_num(), 0.01);
	//md_disp.set_pts_from_node_id(vy_bc_pt_array.get_mem(), vy_bc_pt_array.get_num(), 0.01);
	//md_disp.set_pts_from_node_id(vz_bc_pt_array.get_mem(), vz_bc_pt_array.get_num(), 0.01);
	md_disp.start();
	return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_me_mt_cylinder_foundation.h5");

	ModelData_T3D_ME_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T3D_ME_mt_complete out1("penetration");
	out1.set_res_file(res_file_hdf5);
	out1.set_output_init_state();
	out1.set_interval_num(50);
	TimeHistory_ConsoleProgressBar out_cpb;

	Step_T3D_ME_mt step("step1");
	step.set_model(model);
	step.set_step_time(1.0);
	//step.set_step_time(1.0e-4);
	step.set_dtime(1.0e-5);
	step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();
}

#include "QtApp_Posp_T3D_ME_mt.h"
#include "test_model_view_omp.h"

void test_t3d_me_mt_cylinder_foundation_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t3d_me_mt_cylinder_foundation.h5");

	//QtApp_Posp_T3D_ME_mt app(argc, argv, QtApp_Posp_T3D_ME_mt::SingleFrame);
	//app.set_res_file(rf, "compression", 2, Hdf5Field::z);
	QtApp_Posp_T3D_ME_mt app(argc, argv, QtApp_Posp_T3D_ME_mt::Animation);
	app.set_res_file(rf, "penetration", Hdf5Field::s33);
	app.set_ani_time(5.0);
	app.set_win_size(1200, 950);
	app.set_view_dir(30.0f, 0.0f);
	app.set_light_dir(35.0f, -50.0f);
	app.set_view_dist_scale(1.5);
	app.set_color_map_fld_range(-100.0, 0.0);
	app.set_color_map_geometry(0.7f, 0.45f, 0.5f);
	//app.set_png_name("t3d_me_mt_cylinder_foundation");
	//app.set_gif_name("t3d_me_mt_cylinder_foundation");
	app.start();
}