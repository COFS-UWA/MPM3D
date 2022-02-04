#include "TestsParallel_pcp.h"

#include "ParticleGenerator3D.hpp"
#include "Model_T3D_ME_mt.h"
#include "ModelData_T3D_ME_mt.h"
#include "QtApp_Prep_T3D_ME_mt.h"
#include "Step_T3D_ME_mt.h"
#include "Step_T3D_ME_TBB.h"
#include "TimeHistory_T3D_ME_mt_complete.h"
#include "TimeHistory_T3D_ME_TBB_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "test_parallel_utils.h"
#include "test_simulations_omp.h"

#include "ivt_utils.h"

void test_t3d_me_tbb_piezofoundation_sim_mat_model(int argc, char** argv)
{
	TetrahedronMesh teh_mesh;
	teh_mesh.load_mesh_from_hdf5("../../Asset/piezofoundation_soil_quarter.h5");
	teh_mesh.init_search_grid(0.1, 0.1, 0.1);
	std::cout << "node_num: " << teh_mesh.get_node_num() << "\n"
		<< "elem_num: " << teh_mesh.get_elem_num() << "\n";

	constexpr double pzf_radius = 1.0;
	constexpr double depth = 7.0 * pzf_radius;
	constexpr double coarse_depth = 3.8 * pzf_radius;
	constexpr double width = 5.5 * pzf_radius;
	constexpr double coarse_width = 3.3 * pzf_radius;
	constexpr double sml_pcl_size = 0.03125;
	constexpr double lgr_pcl_size = 0.0625;
	ParticleGenerator3D<TetrahedronMesh> pcl_generator;
	// dense area
	pcl_generator.generate_pcls_grid(
		Cube(0.0, coarse_width,
			0.0, coarse_width,
			-coarse_depth, 0.0),
		sml_pcl_size, sml_pcl_size, sml_pcl_size);
	// peripheral coarse area
	// 4 edges
	pcl_generator.generate_pcls_grid(
		Cube(0.0, coarse_width,
			coarse_width, width,
			-coarse_depth, 0.0),
		lgr_pcl_size, lgr_pcl_size, lgr_pcl_size);
	pcl_generator.generate_pcls_grid(
		Cube(coarse_width, width,
			0.0, coarse_width,
			-coarse_depth, 0.0),
		lgr_pcl_size, lgr_pcl_size, lgr_pcl_size);
	// 4 corners
	pcl_generator.generate_pcls_grid(
		Cube(coarse_width, width, coarse_width, width,
			-coarse_depth, 0.0),
		lgr_pcl_size, lgr_pcl_size, lgr_pcl_size);
	// bottom
	pcl_generator.generate_pcls_grid(
		Cube(0.0, width, 0.0, width, -depth, -coarse_depth),
		lgr_pcl_size, lgr_pcl_size, lgr_pcl_size);
	//
	pcl_generator.adjust_pcl_size_to_fit_elems(teh_mesh);
	std::cout << "pcl_num: " << pcl_generator.get_num() << "\n";

	Model_T3D_ME_mt model;
	model.init_mesh(teh_mesh);
	model.init_search_grid(teh_mesh);
	teh_mesh.clear();
	model.init_pcls(pcl_generator, 2000.0);
	pcl_generator.clear();

	const size_t pcl_num = model.get_pcl_num();
	MatModel::MaterialModel** mms = model.get_mat_models();
	// Linear elasticity
	MatModel::LinearElasticity* les = model.add_LinearElasticity(pcl_num);
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		MatModel::LinearElasticity& le = *les;
		le.set_param(1.0e6, 0.2);
		mms[pcl_id] = les;
		les = model.following_LinearElasticity(les);
	}
	// Tresca
	//MatModel::Tresca *trcs = model.add_Tresca(pcl_num);
	//for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	//{
	//	MatModel::Tresca &trc = trcs[pcl_id];
	//	trc.set_param(1.0e6, 0.2, 5.0e3);
	//	mms[pcl_id] = &trc;
	//}
	
	model.init_rigid_cylinder(0.0, 0.0, 1.125, 2.25, 1.0, 2000.0);
	model.set_rigid_cylinder_velocity(0.0, 0.0, -1.5);
	model.set_contact_param(1.0e5 / (sml_pcl_size * sml_pcl_size),
		1.0e3 / (sml_pcl_size * sml_pcl_size), 0.1, 5.0);

	IndexArray vx_bc_pt_array(100);
	find_3d_nodes_on_x_plane(model, vx_bc_pt_array, 0.0);
	find_3d_nodes_on_x_plane(model, vx_bc_pt_array, width, false);
	model.init_fixed_vx_bc(vx_bc_pt_array.get_num(), vx_bc_pt_array.get_mem());

	IndexArray vy_bc_pt_array(100);
	find_3d_nodes_on_y_plane(model, vy_bc_pt_array, 0.0);
	find_3d_nodes_on_y_plane(model, vy_bc_pt_array, width, false);
	model.init_fixed_vy_bc(vy_bc_pt_array.get_num(), vy_bc_pt_array.get_mem());

	IndexArray vz_bc_pt_array(100);
	find_3d_nodes_on_z_plane(model, vz_bc_pt_array, -depth);
	model.init_fixed_vz_bc(vz_bc_pt_array.get_num(), vz_bc_pt_array.get_mem());

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_me_tbb_piezofoundation_sim_mat_model.h5");

	ModelData_T3D_ME_mt md;
	md.output_model(model, res_file_hdf5);

	QtApp_Prep_T3D_ME_mt md_disp(argc, argv);
	md_disp.set_model(model);
	md_disp.set_win_size(1200, 950);
	md_disp.set_view_dir(-70.0f, 30.0f);
	md_disp.set_light_dir(-60.0f, 15.0f);
	md_disp.set_display_bg_mesh(false);
	//md_disp.set_pts_from_vx_bc(0.05);
	md_disp.set_pts_from_vy_bc(0.05);
	//md_disp.set_pts_from_vz_bc(0.05);
	md_disp.start();
}

// ===============================================================
// =================== mt version and restart ====================
void test_t3d_me_mt_piezofoundation_sim_mat(int argc, char** argv)
{
	Model_T3D_ME_mt model;
	Model_T3D_ME_mt_hdf5_utilities::load_me_mt_model_from_hdf5_file(
		model, "t3d_me_tbb_piezofoundation_sim_mat_model.h5");
	std::cout << "Completed loading model.\n";

	model.set_rigid_cylinder_velocity(0.0, 0.0, 0.0);

	//QtApp_Prep_T3D_ME_mt md_disp(argc, argv);
	//md_disp.set_model(model);
	//md_disp.set_win_size(1200, 950);
	//md_disp.set_view_dir(-70.0f, 30.0f);
	//md_disp.set_light_dir(-60.0f, 15.0f);
	//md_disp.set_display_bg_mesh(false);
	////md_disp.set_pts_from_vx_bc(0.05);
	//md_disp.set_pts_from_vy_bc(0.05);
	////md_disp.set_pts_from_vz_bc(0.05);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_me_mt_piezofoundation_sim_mat.h5");

	//ModelData_T3D_ME_mt md;
	//md.output_model(model, res_file_hdf5);

	//TimeHistory_T3D_ME_mt_complete out1("penetration");
	//out1.set_interval_num(2);
	//out1.set_output_init_state();
	//out1.set_output_final_state();
	//out1.set_res_file(res_file_hdf5);

	TimeHistory_ConsoleProgressBar out_cpb;
	out_cpb.set_interval_num(100);

	std::cout << "Start solving...\n";
	Step_T3D_ME_mt step("step2");
	step.set_model(model);
	step.set_thread_num(8);
	step.set_step_time(2.0e-4);
	step.set_dtime(2.0e-6);
	//step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();
}

void test_t3d_me_mt_piezofoundation_sim_mat_restart(int argc, char** argv)
{
	Model_T3D_ME_mt model;
	Step_T3D_ME_mt step("step2");
	Model_T3D_ME_mt_hdf5_utilities::load_me_mt_model_from_hdf5_file(
		model, step, "t3d_me_mt_piezofoundation_sim_mat.h5", "penetration", 3);
	std::cout << "Completed loading model.\n";

	model.set_rigid_cylinder_velocity(0.0, 0.0, 0.0);

	//QtApp_Prep_T3D_ME_mt md_disp(argc, argv);
	//md_disp.set_model(model);
	//md_disp.set_win_size(1200, 950);
	//md_disp.set_view_dir(-70.0f, 30.0f);
	//md_disp.set_light_dir(-60.0f, 15.0f);
	//md_disp.set_display_bg_mesh(false);
	////md_disp.set_pts_from_vx_bc(0.05);
	//md_disp.set_pts_from_vy_bc(0.05);
	////md_disp.set_pts_from_vz_bc(0.05);
	//md_disp.start();
	//return;

	//ResultFile_hdf5 res_file_hdf5;
	//res_file_hdf5.create("t3d_me_tbb_piezofoundation_sim_mat.h5");

	//ModelData_T3D_ME_mt md;
	//md.output_model(model, res_file_hdf5);

	//TimeHistory_T3D_ME_TBB_complete out1("penetration");
	//out1.set_interval_num(100);
	//out1.set_output_init_state();
	//out1.set_output_final_state();
	//out1.set_res_file(res_file_hdf5);
	TimeHistory_ConsoleProgressBar out_cpb;
	out_cpb.set_interval_num(100);

	std::cout << "Start solving...\n";
	step.set_thread_num(8);
	step.set_step_time(2.0e-4);
	step.set_dtime(2.0e-6);
	//step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();

	//Step_T3D_ME_mt step("step2");
	//step.set_model(model);
	//step.set_thread_num(3);
	//step.set_step_time(2.0e-3);
	//step.set_dtime(2.0e-6);
	////step.add_time_history(out1);
	//step.add_time_history(out_cpb);
	//step.solve();
}

// ===============================================================
// =================== tbb version and restart ===================
void test_t3d_me_tbb_piezofoundation_sim_mat(int argc, char** argv)
{
	IVT_PAUSE;

	Model_T3D_ME_mt model;
	Model_T3D_ME_mt_hdf5_utilities::load_me_mt_model_from_hdf5_file(
		model, "t3d_me_tbb_piezofoundation_sim_mat_model.h5");
	std::cout << "Completed loading model.\n";

	model.set_rigid_cylinder_velocity(0.0, 0.0, 0.0);

	//QtApp_Prep_T3D_ME_mt md_disp(argc, argv);
	//md_disp.set_model(model);
	//md_disp.set_win_size(1200, 950);
	//md_disp.set_view_dir(-70.0f, 30.0f);
	//md_disp.set_light_dir(-60.0f, 15.0f);
	//md_disp.set_display_bg_mesh(false);
	////md_disp.set_pts_from_vx_bc(0.05);
	//md_disp.set_pts_from_vy_bc(0.05);
	////md_disp.set_pts_from_vz_bc(0.05);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_me_tbb_piezofoundation_sim_mat.h5");

	ModelData_T3D_ME_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T3D_ME_TBB_complete out1("penetration");
	out1.set_interval_num(2);
	out1.set_output_init_state();
	out1.set_output_final_state();
	out1.set_res_file(res_file_hdf5);
	TimeHistory_ConsoleProgressBar out_cpb;
	out_cpb.set_interval_num(100);

	std::cout << "Start solving...\n";
	Step_T3D_ME_TBB step("step2");
	step.set_model(model);
	step.set_thread_num(12);
	step.set_step_time(2.0e-4);
	step.set_dtime(2.0e-6);
	//step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();
}

void test_t3d_me_tbb_piezofoundation_sim_mat_restart(int argc, char** argv)
{
	Model_T3D_ME_mt model;
	Step_T3D_ME_TBB step("step2");
	Model_T3D_ME_mt_hdf5_utilities::load_me_mt_model_from_hdf5_file(
		model, step, "t3d_me_tbb_piezofoundation_sim_mat.h5", "penetration", 2);
	std::cout << "Completed loading model.\n";

	model.set_rigid_cylinder_velocity(0.0, 0.0, 0.0);

	//QtApp_Prep_T3D_ME_mt md_disp(argc, argv);
	//md_disp.set_model(model);
	//md_disp.set_win_size(1200, 950);
	//md_disp.set_view_dir(-70.0f, 30.0f);
	//md_disp.set_light_dir(-60.0f, 15.0f);
	//md_disp.set_display_bg_mesh(false);
	////md_disp.set_pts_from_vx_bc(0.05);
	//md_disp.set_pts_from_vy_bc(0.05);
	////md_disp.set_pts_from_vz_bc(0.05);
	//md_disp.start();
	//return;

	//ResultFile_hdf5 res_file_hdf5;
	//res_file_hdf5.create("t3d_me_tbb_piezofoundation_sim_mat.h5");

	//ModelData_T3D_ME_mt md;
	//md.output_model(model, res_file_hdf5);

	//TimeHistory_T3D_ME_TBB_complete out1("penetration");
	//out1.set_interval_num(10);
	//out1.set_output_init_state();
	//out1.set_output_final_state();
	//out1.set_res_file(res_file_hdf5);
	
	TimeHistory_ConsoleProgressBar out_cpb;
	out_cpb.set_interval_num(100);

	std::cout << "Start solving...\n";
	step.set_thread_num(4);
	step.set_step_time(2.0e-4);
	step.set_dtime(2.0e-6);
	//step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();
}

#include "QtApp_Posp_T3D_ME_mt_Div.h"
#include "test_model_view_omp.h"

void test_t3d_me_tbb_piezofoundation_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t3d_me_tbb_piezofoundation.h5");

	// Single frame
	//QtApp_Posp_T3D_ME_mt_Div<PlaneDivisionSet> app(argc, argv, QtApp_Posp_T3D_ME_mt_Div<PlaneDivisionSet>::SingleFrame);
	//app.set_res_file(rf, "geostatic", 21, Hdf5Field::s33);
	//app.set_color_map_fld_range(-8.0e4, 0.0);
	//app.set_color_map_geometry(1.2f, 0.4f, 0.45f);
	//app.set_win_size(1600, 950);
	//app.set_view_dir(-90.0f, 20.0f);
	//app.set_light_dir(-135.0f, 20.0f);
	//app.set_light_dist_scale(1.0f);
	//app.set_fog_coef(0.02f);
	//app.set_view_dist_scale(0.85);
	//app.set_display_bg_mesh(false);
	//app.start();

	//// Animation
	//QtApp_Posp_T3D_ME_mt_Div<PlaneDivisionSet> app(argc, argv,
	//	QtApp_Posp_T3D_ME_mt_Div<PlaneDivisionSet>::Animation);
	//app.get_div_set().set_by_normal_and_point(0.0, 1.0, 0.0, 0.0, 0.06, 0.0);
	//app.set_ani_time(5.0);
	//app.set_win_size(1600, 950);
	//app.set_view_dir(-90.0f, 10.0f);
	//app.set_light_dir(-135.0f, 20.0f);
	//app.set_light_dist_scale(1.0f);
	//app.set_fog_coef(0.02f);
	//app.set_view_dist_scale(0.85);
	//app.set_display_bg_mesh(false);
	////app.set_mono_color_pcl(true);
	//// s33
	////app.set_res_file(rf, "penetration", Hdf5Field::s33);
	////app.set_color_map_fld_range(-8.0e4, 0.0);
	//// shear stress
	////app.set_res_file(rf, "penetration", Hdf5Field::max_shear_stress);
	////app.set_color_map_fld_range(0.0, 5.0);
	//// mises strain
	////app.set_res_file(rf, "penetration", Hdf5Field::plastic_mises_strain_2d);
	////app.set_color_map_fld_range(0.0, 0.08);
	//// mat_e
	//app.set_res_file(rf, "penetration", Hdf5Field::mat_e);
	//app.set_color_map_fld_range(0.6, 0.8);
	//// mat_s11
	////app.set_res_file(rf, "penetration", Hdf5Field::mat_s11);
	////app.set_color_map_fld_range(-4.0e4, 0.0);
	//// mat_s33
	////app.set_res_file(rf, "penetration", Hdf5Field::mat_s33);
	////app.set_color_map_fld_range(-1.0e5, 0.0);
	////
	//app.set_color_map_geometry(1.2f, 0.4f, 0.45f);
	////app.set_png_name("t3d_me_mt_piezofoundation");
	//app.set_gif_name("t3d_me_mt_piezofoundation");

	QtApp_Posp_T3D_ME_mt_Div<PlaneDivisionSet> app(argc, argv,
		QtApp_Posp_T3D_ME_mt_Div<PlaneDivisionSet>::SingleFrame);
	app.get_div_set().set_by_normal_and_point(0.0, 1.0, 0.0, 0.0, 0.06, 0.0);
	app.set_win_size(1600, 950);
	app.set_view_dir(-90.0f, 10.0f);
	app.set_light_dir(-135.0f, 20.0f);
	app.set_light_dist_scale(1.0f);
	app.set_fog_coef(0.02f);
	app.set_view_dist_scale(0.85);
	app.set_display_bg_mesh(false);
	app.set_color_map_geometry(1.2f, 0.4f, 0.45f);
	// mat_e
	app.set_res_file(rf, "penetration", 101, Hdf5Field::mat_e);
	app.set_color_map_fld_range(0.6, 0.8);
	//
	app.set_png_name("t3d_me_tbb_piezofoundation");

	app.start();
}
