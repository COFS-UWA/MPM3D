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

void test_t3d_me_mt_cap_compression(int argc, char **argv)
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
	//model.init_pcls(pcl_generator, 1800.0);
	const size_t pcl_num = model.get_pcl_num();
	MatModel::MaterialModel** mms = model.get_mat_models();
	// Linear elasticity
	//MatModel::LinearElasticity* les = model.add_LinearElasticity(pcl_num);
	//for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	//{
	//	MatModel::LinearElasticity& le = les[pcl_id];
	//	le.set_param(1000.0, 0.0);
	//	mms[pcl_id] = &le;
	//}
	// Tresca
	//MatModel::LinearElasticity* les = model.add_LinearElasticity(pcl_num);
	//for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	//{
	//	MatModel::LinearElasticity& le = les[pcl_id];
	//	le.set_param(1000.0, 0.0);
	//	mms[pcl_id] = &le;
	//}
	// Mohr Coulomb
	const double ini_stress[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	MatModel::MohrCoulombWrapper* mcs = model.add_MohrCoulombWrapper(pcl_num);
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		mcs->set_param(ini_stress, 30.0, 0.0, 100.0, 10000.0, 0.3);
		model.add_mat_model(pcl_id, *mcs, sizeof(MatModel::MohrCoulombWrapper));
		mcs = model.following_MohrCoulombWrapper(mcs);
	}
	// Stb hypoplasticity
	//const double ini_stress[6] = { -100.0e3, -100.0e3, -100.0e3, 0.0, 0.0, 0.0 };
	//MatModel::SandHypoplasticityWrapper* shps = model.add_SandHypoplasticityWrapper(pcl_num);
	//for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	//{
	//	MatModel::SandHypoplasticityWrapper& shp = shps[pcl_id];
	//	shp.set_param(
	//		ini_stress, 0.55,
	//		30.0, 1354.0e6, 0.34,
	//		0.18, 1.27,
	//		0.49, 0.76, 0.86);
	//	mms[pcl_id] = &shp;
	//}
	// Stb hypoplasticity
	//const double ini_stress[6] = { -100.0e3, -100.0e3, -100.0e3, 0.0, 0.0, 0.0 };
	//MatModel::SandHypoplasticityStbWrapper* shps = model.add_SandHypoplasticityStbWrapper(pcl_num);
	//for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	//{
	//	MatModel::SandHypoplasticityStbWrapper& shp = shps[pcl_id];
	//	shp.set_param(
	//		ini_stress, 0.55, //0.75,// 0.55,
	//		30.0, 1354.0e6, 0.34,
	//		0.18, 1.27,
	//		0.49, 0.76, 0.86,
	//		//1.5, 43.0, 100, // drained
	//		1.5, 43.0, 180, //100, //180,
	//		200.0, 0.2);
	//	mms[pcl_id] = &shp;
	//}

	model.init_rigid_cylinder(0.0, 0.0, 1.025, 0.05, 0.3);
	model.set_rigid_cylinder_velocity(0.0, 0.0, -0.1);
	const double Kct = 100.0 / (0.025 * 0.025); 
	//const double Kct = 1.0e5 / (0.025 * 0.025); // hypo 
	model.set_contact_param(Kct, Kct, 0.1, 0.2);

	//model.init_rigid_cone(0.0, 0.5, 0.5, 0.34641, 0.3, 0.2);
	//model.set_rigid_cone_velocity(0.0, 0.0, 0.0);

	IndexArray vx_bc_pt_array(100);
	find_3d_nodes_on_x_plane(model, vx_bc_pt_array, 0.0);
	model.init_fixed_vx_bc(vx_bc_pt_array.get_num(), vx_bc_pt_array.get_mem());

	IndexArray vy_bc_pt_array(100);
	find_3d_nodes_on_y_plane(model, vy_bc_pt_array, 0.0);
	model.init_fixed_vy_bc(vy_bc_pt_array.get_num(), vy_bc_pt_array.get_mem());

	IndexArray vz_bc_pt_array(100);
	find_3d_nodes_on_z_plane(model, vz_bc_pt_array, 0.0);
	model.init_fixed_vz_bc(vz_bc_pt_array.get_num(), vz_bc_pt_array.get_mem());

	std::cout << "pcl_num: " << model.get_pcl_num() << ",\n"
		"elem_num: " << model.get_elem_num() << ",\n"
		"node_num: " << model.get_node_num() << ",\n";

	//QtApp_Prep_T3D_ME_mt_Div<> md_disp(argc, argv);
	////QtApp_Prep_T3D_ME_mt_Div<BoxDivisionSet> md_disp(argc, argv);
	////md_disp.get_div_set().set_param(0.0, 0.1, 0.0, 0.1, 0.0, 0.18);
	////QtApp_Prep_T3D_ME_mt_Div<PlaneDivisionSet> md_disp(argc, argv);
	////md_disp.get_div_set().set_param(0.0, 0.0, -1.0, 0.45);
	//md_disp.set_win_size(1200, 950);
	//md_disp.set_view_dir(30.0f, 25.0f); // (-130.0f, -5.0f)
	//md_disp.set_light_dir(40.0f, 30.0f); // (-140.0f, -10.0f)
	//md_disp.set_view_dist_scale(1.0);
	//md_disp.set_model(model);
	////md_disp.set_pts_from_node_id(vx_bc_pt_array.get_mem(), vx_bc_pt_array.get_num(), 0.01);
	////md_disp.set_pts_from_node_id(vy_bc_pt_array.get_mem(), vy_bc_pt_array.get_num(), 0.01);
	////md_disp.set_pts_from_node_id(vz_bc_pt_array.get_mem(), vz_bc_pt_array.get_num(), 0.01);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_me_mt_cap_compression.h5");

	ModelData_T3D_ME_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T3D_ME_mt_complete out1("compression");
	//out1.set_output_init_state();
	out1.set_interval_num(200);
	out1.set_res_file(res_file_hdf5);
	TimeHistory_ConsoleProgressBar out_cpb;

	Step_T3D_ME_mt step("step1");
	step.set_model(model);
	step.set_thread_num(4);
	step.set_step_time(0.5);
	//step.set_step_time(5.0e-4);
	step.set_dtime(5.0e-6);
	step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();
}

void test_t3d_me_mt_cap_compression_restart(int argc, char** argv)
{
	Model_T3D_ME_mt model;
	Model_T3D_ME_mt_hdf5_utilities::load_me_mt_model_from_hdf5_file(model, "t3d_me_mt_cap_compression.h5");

	//model.clear_rigid_cylinder();

	//QtApp_Prep_T3D_ME_mt_Div<> md_disp(argc, argv);
	////QtApp_Prep_T3D_ME_mt_Div<TwoPlaneDivisionSet> md_disp(argc, argv);
	////auto& div_set = md_disp.get_div_set();
	////div_set.seta().set_by_normal_and_point(0.0, 1.0, 0.0, 3.5, 3.5, 0.0);
	////div_set.setb().set_by_normal_and_point(1.0, 0.0, 0.0, 3.5, 3.5, 0.0);
	//md_disp.set_win_size(1200, 950);
	//md_disp.set_view_dir(30.0f, 30.0f);
	//md_disp.set_light_dir(30.0f, 30.0f);
	//md_disp.set_model(model);
	////md_disp.set_pts_from_vx_bc(0.01);
	//md_disp.set_pts_from_vy_bc(0.01);
	////md_disp.set_pts_from_vz_bc(0.01);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_me_mt_cap_compression_restart.h5");

	ModelData_T3D_ME_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T3D_ME_mt_complete out1("compression");
	out1.set_res_file(res_file_hdf5);
	out1.set_output_init_state();
	out1.set_output_final_state();
	out1.set_interval_num(100);
	TimeHistory_ConsoleProgressBar out_cpb;

	Step_T3D_ME_mt step("step1");
	step.set_model(model);
	step.set_step_time(0.5);
	//step.set_step_time(5.0e-5);
	step.set_dtime(1.0e-5);
	step.set_thread_num(6);
	step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();
}

void test_t3d_me_mt_cap_compression_restart2(int argc, char **argv)
{
	Model_T3D_ME_mt model;
	Step_T3D_ME_mt step("step1");
	Model_T3D_ME_mt_hdf5_utilities::load_me_mt_model_from_hdf5_file(
		model,
		step,
		"t3d_me_mt_cap_compression.h5",
		"compression",
		0
		);

	//QtApp_Prep_T3D_ME_mt_Div<> md_disp(argc, argv);
	////QtApp_Prep_T3D_ME_mt_Div<TwoPlaneDivisionSet> md_disp(argc, argv);
	////auto& div_set = md_disp.get_div_set();
	////div_set.seta().set_by_normal_and_point(0.0, 1.0, 0.0, 3.5, 3.5, 0.0);
	////div_set.setb().set_by_normal_and_point(1.0, 0.0, 0.0, 3.5, 3.5, 0.0);
	//md_disp.set_win_size(1200, 950);
	//md_disp.set_view_dir(30.0f, 30.0f);
	//md_disp.set_light_dir(30.0f, 30.0f);
	//md_disp.set_model(model);
	////md_disp.set_pts_from_vx_bc(0.01);
	//md_disp.set_pts_from_vy_bc(0.01);
	////md_disp.set_pts_from_vz_bc(0.01);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_me_mt_cap_compression_restart2.h5");

	ModelData_T3D_ME_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T3D_ME_mt_complete out1("compression");
	out1.set_res_file(res_file_hdf5);
	out1.set_output_init_state();
	out1.set_output_final_state();
	out1.set_interval_num(100);
	TimeHistory_ConsoleProgressBar out_cpb;

	step.set_step_time(0.5);
	//step.set_step_time(5.0e-5);
	step.set_dtime(1.0e-5);
	step.set_thread_num(6);
	step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();
}

#include "QtApp_Posp_T3D_ME_mt.h"
#include "test_model_view_omp.h"

void test_t3d_me_mt_cap_compression_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t3d_me_mt_cap_compression.h5");
	//rf.open("t3d_me_mt_cap_compression_restart.h5");
	//rf.open("t3d_me_mt_cap_compression_restart2.h5");

	//QtApp_Posp_T3D_ME_mt app(argc, argv, QtApp_Posp_T3D_ME_mt::SingleFrame);
	//app.set_win_size(1200, 950);
	//app.set_view_dir(30.0f, 0.0f);
	//app.set_light_dir(30.0f, 20.0f);
	//app.set_color_map_geometry(0.85f, 0.45f, 0.5f);
	//app.set_update_rb_pos();
	//// s33
	//app.set_res_file(rf, "compression", 19, Hdf5Field::s33);
	//app.set_color_map_fld_range(-250.0e3, 0.0);
	//// e
	////app.set_res_file(rf, "compression", 119, Hdf5Field::mat_e);
	////app.set_color_map_fld_range(0.67, 0.75);
	////app.start();
	
	QtApp_Posp_T3D_ME_mt app(argc, argv, QtApp_Posp_T3D_ME_mt::Animation);
	app.set_win_size(1200, 950);
	app.set_ani_time(5.0);
	app.set_view_dir(30.0f, 0.0f);
	app.set_light_dir(30.0f, 20.0f);
	//app.set_view_dist_scale(1.1);
	app.set_color_map_geometry(0.95f, 0.45f, 0.5f);
	// s33
	app.set_res_file(rf, "compression", Hdf5Field::s33);
	app.set_color_map_fld_range(-100.0, 0.0);
	// shear stress
	//app.set_res_file(rf, "compression", Hdf5Field::max_shear_stress);
	//app.set_color_map_fld_range(0.0, 30.0);
	//
	//app.set_png_name("t3d_me_mt_cap_compression");
	//app.set_gif_name("t3d_me_mt_cap_compression");
	app.start();
}

#include "QtApp_Posp_T3D_ME_mt_Div.h"

void test_t3d_me_mt_cap_compression_result_div(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t3d_me_mt_cap_compression.h5");
	//rf.open("t3d_me_mt_cap_compression_restart.h5");
	//rf.open("t3d_me_mt_cap_compression_restart2.h5");

	//QtApp_Posp_T3D_ME_mt_Div<PlaneDivisionSet> app(argc, argv, QtApp_Posp_T3D_ME_mt_Div<PlaneDivisionSet>::SingleFrame);
	//app.set_res_file(rf, "compression", 2, Hdf5Field::z);
	QtApp_Posp_T3D_ME_mt_Div<PlaneDivisionSet> app(argc, argv, QtApp_Posp_T3D_ME_mt_Div<PlaneDivisionSet>::Animation);
	app.set_res_file(rf, "compression", Hdf5Field::s33);
	app.get_div_set().set_by_normal_and_point(0.0, 0.0, -1.0, 0.0, 0.0, 0.4);
	app.set_win_size(1200, 950);
	app.set_ani_time(5.0);
	app.set_view_dir(30.0f, 30.0f);
	app.set_light_dir(35.0f, 30.0f);
	app.set_view_dist_scale(1.1);
	app.set_color_map_fld_range(-50.0, 0.0);
	app.set_color_map_geometry(0.7f, 0.45f, 0.5f);
	//app.set_png_name("t3d_me_mt_cap_compression_div");
	//app.set_gif_name("t3d_me_mt_cap_compression_div");
	app.start();
}
