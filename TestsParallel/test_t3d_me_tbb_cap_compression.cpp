#include "TestsParallel_pcp.h"

#include "ParticleGenerator3D.hpp"
#include "LinearElasticity.h"
#include "Model_T3D_ME_mt.h"
#include "ModelData_T3D_ME_mt.h"
#include "Step_T3D_ME_TBB.h"
#include "TimeHistory_T3D_ME_TBB_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "QtApp_Prep_T3D_ME_mt_Div.h"
#include "test_parallel_utils.h"
#include "test_simulations_omp.h"

void test_t3d_me_tbb_cap_compression(int argc, char **argv)
{
	TetrahedronMesh teh_mesh;
	teh_mesh.load_mesh_from_hdf5("../../Asset/brick_mesh_1.00_2x2x10.h5");
	teh_mesh.init_search_grid(0.05, 0.05, 0.05);

	ParticleGenerator3D<TetrahedronMesh> pcl_generator;
	pcl_generator.generate_pcls_grid(Cube(0.0, 0.2, 0.0, 0.2, 0.0, 1.0), 0.025, 0.025, 0.025);
	pcl_generator.adjust_pcl_size_to_fit_elems(teh_mesh);

	Model_T3D_ME_mt model;
	model.init_mesh(teh_mesh);
	model.init_search_grid(teh_mesh);
	model.init_pcls(pcl_generator, 10.0);
	const size_t pcl_num = model.get_pcl_num();
	MatModel::MaterialModel** mms = model.get_mat_models();
	// Linear elasticity
	MatModel::LinearElasticity* les = model.add_LinearElasticity(pcl_num);
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		les->set_param(1000.0, 0.0);
		mms[pcl_id] = les;
		les = model.following_LinearElasticity(les);
	}

	model.init_rigid_cylinder(0.1, 0.1, 1.025, 0.05, 0.2);
	model.set_rigid_cylinder_velocity(0.0, 0.0, -0.02);
	//model.set_cylinder_vz_bc_ramp_up_time(1.0);

	//model.init_t3d_rigid_mesh(1.0, "../../Asset/cone_pap2.h5",
	//	0.1, 0.1, 1.0, 90.0, 0.0, 0.0, 0.05, 0.05, 0.05); // D = 3 m
	//model.set_t3d_rigid_mesh_velocity(0.0, 0.0, -0.02);

	model.set_contact_param(20.0 / (0.025*0.025), 20.0 / (0.025*0.025), 0.1, 0.1);

	IndexArray vx_bc_pt_array(100);
	find_3d_nodes_on_x_plane(model, vx_bc_pt_array, 0.0);
	find_3d_nodes_on_x_plane(model, vx_bc_pt_array, 0.2, false);
	model.init_fixed_vx_bc(vx_bc_pt_array.get_num(), vx_bc_pt_array.get_mem());

	IndexArray vy_bc_pt_array(100);
	find_3d_nodes_on_y_plane(model, vy_bc_pt_array, 0.0);
	find_3d_nodes_on_y_plane(model, vy_bc_pt_array, 0.2, false);
	model.init_fixed_vy_bc(vy_bc_pt_array.get_num(), vy_bc_pt_array.get_mem());

	IndexArray vz_bc_pt_array(100);
	find_3d_nodes_on_z_plane(model, vz_bc_pt_array, 0.0);
	model.init_fixed_vz_bc(vz_bc_pt_array.get_num(), vz_bc_pt_array.get_mem());

	//QtApp_Prep_T3D_ME_mt_Div<> md_disp(argc, argv);
	////QtApp_Prep_T3D_ME_mt_Div<BoxDivisionSet> md_disp(argc, argv);
	////md_disp.get_div_set().set_param(0.0, 0.1, 0.0, 0.1, 0.0, 0.18);
	////QtApp_Prep_T3D_ME_mt_Div<PlaneDivisionSet> md_disp(argc, argv);
	////md_disp.get_div_set().set_param(0.0, 0.0, -1.0, 0.45);
	//md_disp.set_win_size(1200, 950);
	//md_disp.set_view_dir(30.0f, 0.0f);
	//md_disp.set_light_dir(30.0f, 20.0f);
	//md_disp.set_bg_color(QVector3D(1.0f, 1.0f, 1.0f));
	////md_disp.set_view_dist_scale(0.5);
	//md_disp.set_model(model);
	////md_disp.set_pts_from_node_id(vx_bc_pt_array.get_mem(), vx_bc_pt_array.get_num(), 0.01);
	////md_disp.set_pts_from_node_id(vy_bc_pt_array.get_mem(), vy_bc_pt_array.get_num(), 0.01);
	////md_disp.set_pts_from_node_id(vz_bc_pt_array.get_mem(), vz_bc_pt_array.get_num(), 0.01);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_me_tbb_cap_compression.h5");

	ModelData_T3D_ME_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T3D_ME_TBB_complete out1("compression");
	out1.set_res_file(res_file_hdf5);
	out1.set_output_init_state();
	out1.set_interval_num(50);
	TimeHistory_ConsoleProgressBar out_cpb;

	Step_T3D_ME_TBB step("step1");
	step.set_model(model);
	step.set_step_time(2.5);
	//step.set_step_time(1.0e-5);
	step.set_dtime(1.0e-5);
	step.set_thread_num(4);
	step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();
}

#include "QtApp_Posp_T3D_ME_mt.h"
#include "test_model_view_omp.h"

void test_t3d_me_tbb_cap_compression_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t3d_me_tbb_cap_compression.h5");

	QtApp_Posp_T3D_ME_mt app(argc, argv, QtApp_Posp_T3D_ME_mt::SingleFrame);
	app.set_res_file(rf, "compression", 50, Hdf5Field::s33);
	//QtApp_Posp_T3D_ME_mt app(argc, argv, QtApp_Posp_T3D_ME_mt::Animation);
	app.set_win_size(1200, 950);
	app.set_ani_time(7.5);
	app.set_view_dir(30.0f, 0.0f);
	app.set_light_dir(30.0f, 20.0f);
	//app.set_view_dist_scale(1.1);
	app.set_bg_color(1.0f, 1.0f, 1.0f);
	app.set_color_map_char_color(0.0f, 0.0f, 0.0f);
	app.set_color_map_geometry(0.85f, 0.45f, 0.5f);
	// s33
	app.set_res_file(rf, "compression", Hdf5Field::s33);
	app.set_color_map_fld_range(-66.0, 0.0);
	// shear stress
	//app.set_res_file(rf, "compression", Hdf5Field::max_shear_stress);
	//app.set_color_map_fld_range(0.0, 30.0);
	// shear strain
	//app.set_res_file(rf, "compression", Hdf5Field::mises_strain_3d);
	//app.set_color_map_fld_range(0.0, 0.1);
	//
	//app.set_png_name("t3d_me_tbb_cap_compression");
	//app.set_gif_name("t3d_me_tbb_cap_compression");
	app.start();
}
