#include "TestsParallel_pcp.h"

#include "ParticleGenerator3D.hpp"
#include "Model_T3D_CHM_mt.h"
#include "Step_T3D_CHM_mt.h"
#include "Step_T3D_CHM_mt_Geo.h"
#include "ModelData_T3D_CHM_mt.h"
#include "TimeHistory_T3D_CHM_mt_complete.h"
#include "TimeHistory_T3D_CHM_mt_Geo_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "QtApp_Prep_T3D_CHM_mt.h"
#include "QtApp_Prep_T3D_CHM_mt_Div.h"
#include "test_parallel_utils.h"
#include "test_simulations_omp.h"

void test_t3d_chm_mt_spudcan_model(int argc, char** argv)
{
	TetrahedronMesh teh_mesh;
	teh_mesh.load_mesh_from_hdf5("../../Asset/spudcan_soil_quarter2.h5");
	teh_mesh.init_search_grid(0.08, 0.08, 0.08);
	std::cout << "node_num: " << teh_mesh.get_node_num() << "\n"
			  << "elem_num: " << teh_mesh.get_elem_num() << "\n";

	constexpr double depth = 5.25;
	constexpr double coarse_depth = 3.45;
	constexpr double width = 4.05;
	constexpr double coarse_width = 3.45;
	constexpr double sml_pcl_size = 0.04;
	constexpr double lgr_pcl_size = 0.08;
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
	
	constexpr double e0 = 0.55;
	constexpr double den_grain = 2670.0;
	constexpr double den_dry = den_grain / (e0 + 1.0);
	constexpr double den_sat = den_grain / (e0 + 1.0) + 1000 * e0 / (e0 + 1.0);
	constexpr double den_float = den_sat - 1000.0;
	constexpr double stress_depth_limit = -0.01;
	constexpr double void_depth_limit = -0.05;
	Model_T3D_CHM_mt model;
	model.init_mesh(teh_mesh);
	model.init_search_grid(teh_mesh);
	teh_mesh.clear();
	model.init_pcls(pcl_generator, e0 / (1.0+e0), den_grain, 1000.0, 2.0e7, 1.5e-7, 0.7);
	pcl_generator.clear();

	const size_t pcl_num = model.get_pcl_num();
	MatModel::MaterialModel** mms = model.get_mat_models();
	const double K0 = 1.0 - sin(30.0 / 180.0 * 3.14159265359);
	double ini_stress[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	// Sand hypoplasticity
	constexpr double fric_ang = 30.0;
	constexpr double hs = 1354.0e6;
	constexpr double n = 0.34;
	constexpr double alpha = 0.18;
	constexpr double beta = 1.27;
	constexpr double ed0 = 0.49;
	constexpr double ec0 = 0.76;
	constexpr double ei0 = 0.86;
	//
	constexpr double Ig = 200.0;
	constexpr double niu = 0.2;
	//
	constexpr double N = 1.5;
	constexpr double chi = 43.0;
	constexpr double H = 100.0;
	MatModel::SandHypoplasticityStbWrapper* shps = model.add_SandHypoplasticityStbWrapper(pcl_num);
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		double pcl_z = model.get_pcl_pos()[pcl_id].z;
		auto& pcl_s = model.get_pcl_stress0()[pcl_id];
		pcl_s.s33 = pcl_z * 9.81 * den_float;
		pcl_s.s22 = K0 * pcl_s.s33;
		pcl_s.s11 = pcl_s.s22;
		if (pcl_z > stress_depth_limit) // shallow depth
			pcl_z = stress_depth_limit;
		ini_stress[2] = pcl_z * 9.81 * den_float;
		ini_stress[0] = K0 * ini_stress[2];
		ini_stress[1] = ini_stress[0];
		MatModel::SandHypoplasticityStbWrapper& shp = shps[pcl_id];
		mms[pcl_id] = &shp;
		if (pcl_z > void_depth_limit)
			shp.set_param(
				ini_stress, 0.76, // critical state
				fric_ang, hs, n,
				alpha, beta,
				ed0, ec0, ei0,
				N, chi, H,
				Ig, niu);
		else // normal
			shp.set_param(
				ini_stress, e0,
				fric_ang, hs, n,
				alpha, beta,
				ed0, ec0, ei0,
				N, chi, H,
				Ig, niu);
	}

	model.init_t3d_rigid_mesh(1.0, "../../Asset/spudcan_model_flat_tip.h5",
		0.0, 0.0, 0.0, 90.0, 0.0, 0.0, 0.3, 0.3, 0.3);
	model.set_t3d_rigid_mesh_velocity(0.0, 0.0, -0.15);
	constexpr double K_cont = 5.0e4 / (sml_pcl_size * sml_pcl_size);
	model.set_contact_param(K_cont, K_cont, 0.1, 5.0, K_cont/50.0, K_cont/50.0);

	// gravity force, float unit weight
	IndexArray bfz_pcl_array(pcl_num);
	MemoryUtils::ItemArray<double> bfz_array(pcl_num);
	double bfz = -9.81 * den_float / den_dry;
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		bfz_pcl_array.add(pcl_id);
		bfz_array.add(bfz);
	}
	model.init_bfz_ss(pcl_num, bfz_pcl_array.get_mem(), bfz_array.get_mem());
	
	IndexArray vx_bc_pt_array(100);
	find_3d_nodes_on_x_plane(model, vx_bc_pt_array, 0.0);
	model.init_fixed_vx_f_bc(vx_bc_pt_array.get_num(), vx_bc_pt_array.get_mem());
	find_3d_nodes_on_x_plane(model, vx_bc_pt_array, width, false);
	model.init_fixed_vx_s_bc(vx_bc_pt_array.get_num(), vx_bc_pt_array.get_mem());
	
	IndexArray vy_bc_pt_array(100);
	find_3d_nodes_on_y_plane(model, vy_bc_pt_array, 0.0);
	model.init_fixed_vy_f_bc(vy_bc_pt_array.get_num(), vy_bc_pt_array.get_mem());
	find_3d_nodes_on_y_plane(model, vy_bc_pt_array, width, false);
	model.init_fixed_vy_s_bc(vy_bc_pt_array.get_num(), vy_bc_pt_array.get_mem());

	IndexArray vz_bc_pt_array(100);
	find_3d_nodes_on_z_plane(model, vz_bc_pt_array, -depth);
	model.init_fixed_vz_s_bc(vz_bc_pt_array.get_num(), vz_bc_pt_array.get_mem());

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_chm_mt_spudcan_model.h5");

	ModelData_T3D_CHM_mt md;
	md.output_model(model, res_file_hdf5);

	QtApp_Prep_T3D_CHM_mt md_disp(argc, argv);
	md_disp.set_model(model);
	md_disp.set_win_size(1200, 950);
	md_disp.set_view_dir(-100.0f, 5.0f);
	md_disp.set_light_dir(-90.0f, -5.0f);
	md_disp.set_display_bg_mesh(false);
	md_disp.set_view_dist_scale(0.7);
	//md_disp.set_pts_from_vx_bc_s(0.025);
	//md_disp.set_pts_from_vy_bc_s(0.025);
	md_disp.set_pts_from_vz_bc_s(0.025);
	//md_disp.set_pts_from_vx_bc_f(0.025);
	//md_disp.set_pts_from_vy_bc_f(0.025);
	//md_disp.set_pts_from_vz_bc_f(0.025);
	md_disp.start();
}

void test_t3d_chm_mt_spudcan_geostatic(int argc, char** argv)
{
	Model_T3D_CHM_mt model;
	Model_T3D_CHM_mt_hdf5_utilities::load_model_from_hdf5_file(
		model, "t3d_chm_mt_spudcan_model.h5");

	//QtApp_Prep_T3D_CHM_mt_Div<EmptyDivisionSet> md_disp(argc, argv);
	////QtApp_Prep_T3D_CHM_mt_Div<PlaneDivisionSet> md_disp(argc, argv);
	////md_disp.get_div_set().set_by_normal_and_point(1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	//md_disp.set_model(model);
	//md_disp.set_win_size(1200, 950);
	//md_disp.set_view_dir(0.0f, 5.0f);
	//md_disp.set_light_dir(10.0f, 5.0f);
	//md_disp.set_display_bg_mesh(false);
	//md_disp.set_view_dist_scale(0.6);
	//md_disp.set_pts_from_vx_s_bc(0.05);
	////md_disp.set_pts_from_vy_s_bc(0.05);
	////md_disp.set_pts_from_vz_s_bc(0.05);
	////md_disp.set_pts_from_vx_f_bc(0.05);
	////md_disp.set_pts_from_vy_f_bc(0.05);
	////md_disp.set_pts_from_vz_f_bc(0.05);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_chm_mt_spudcan_geo.h5");

	ModelData_T3D_CHM_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T3D_CHM_mt_Geo_complete out1("geostatic");
	out1.set_interval_num(20);
	out1.set_output_init_state();
	out1.set_output_final_state();
	out1.set_res_file(res_file_hdf5);
	TimeHistory_ConsoleProgressBar out_cpb;
	out_cpb.set_interval_num(2000);

	Step_T3D_CHM_mt_Geo step("step1");
	step.set_model(model);
	step.set_thread_num(22);
	step.set_step_time(1.0); // 1.0
	//step.set_thread_num(5);
	//step.set_step_time(1.0e-5);
	step.set_dtime(1.0e-5);
	step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();
}

void test_t3d_chm_mt_spudcan(int argc, char** argv)
{
	Model_T3D_CHM_mt model;
	Step_T3D_CHM_mt step("step2");
	Model_T3D_CHM_mt_hdf5_utilities::load_model_from_hdf5_file(
		model, step, "t3d_chm_mt_spudcan_geo.h5", "geostatic", 21);
	
	// modified velocity
	model.set_t3d_rigid_mesh_velocity(0.0, 0.0, -0.15);
	// modified contact stiffness
	//constexpr double sml_pcl_size = 0.04;
	//constexpr double K_cont = 5.0e4 / (sml_pcl_size * sml_pcl_size);
	//model.set_contact_param(K_cont, K_cont, 0.1, 5.0, K_cont / 50.0, K_cont / 50.0);
	// modified permeability
	model.set_k(1.0e-9);

	//QtApp_Prep_T3D_CHM_mt_Div<EmptyDivisionSet> md_disp(argc, argv);
	////QtApp_Prep_T3D_CHM_mt_Div<PlaneDivisionSet> md_disp(argc, argv);
	////md_disp.get_div_set().set_by_normal_and_point(1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	//md_disp.set_model(model);
	//md_disp.set_win_size(1200, 950);
	//md_disp.set_view_dir(0.0f, 5.0f);
	//md_disp.set_light_dir(10.0f, 5.0f);
	//md_disp.set_display_bg_mesh(false);
	//md_disp.set_view_dist_scale(0.6);
	//md_disp.set_pts_from_vx_s_bc(0.05);
	////md_disp.set_pts_from_vy_s_bc(0.05);
	////md_disp.set_pts_from_vz_s_bc(0.05);
	////md_disp.set_pts_from_vx_f_bc(0.05);
	////md_disp.set_pts_from_vy_f_bc(0.05);
	////md_disp.set_pts_from_vz_f_bc(0.05);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_chm_mt_spudcan.h5");

	ModelData_T3D_CHM_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T3D_CHM_mt_complete out1("penetration");
	out1.set_interval_num(100);
	out1.set_output_init_state();
	out1.set_output_final_state();
	out1.set_res_file(res_file_hdf5);
	TimeHistory_ConsoleProgressBar out_cpb;
	out_cpb.set_interval_num(2000);

	step.set_model(model);
	step.set_thread_num(22);
	step.set_step_time(3.0); // 3.0
	//step.set_thread_num(3);
	//step.set_step_time(1.0e-5);
	step.set_dtime(5.0e-6);
	step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();
}

#include "QtApp_Posp_T3D_CHM_mt.h"
#include "QtApp_Posp_T3D_CHM_mt_Div.h"
#include "test_model_view_omp.h"

void test_t3d_chm_mt_spudcan_geo_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t3d_chm_mt_spudcan_geo.h5");

	//QtApp_Posp_T3D_CHM_mt_Div<PlaneDivisionSet> app(argc, argv, QtApp_Posp_T3D_CHM_mt_Div<PlaneDivisionSet>::SingleFrame);
	////app.set_res_file(rf, "penetration", 50, Hdf5Field::s33);
	////app.set_res_file(rf, "penetration", 50, Hdf5Field::max_shear_stress);
	//app.set_res_file(rf, "penetration", 50, Hdf5Field::plastic_mises_strain_2d);
	QtApp_Posp_T3D_CHM_mt_Div<PlaneDivisionSet> app(argc, argv,
		QtApp_Posp_T3D_CHM_mt_Div<PlaneDivisionSet>::Animation);
	app.get_div_set().set_by_normal_and_point(0.0, 1.0, 0.0, 0.0, 1.0, 0.0);
	//QtApp_Posp_T3D_CHM_mt app(argc, argv, QtApp_Posp_T3D_CHM_mt::Animation);
	app.set_ani_time(5.0);
	app.set_win_size(1200, 800);
	app.set_view_dir(-90.0f, 10.0f);
	app.set_fog_coef(0.02f);
	app.set_light_dir(-135.0f, 20.0f);
	app.set_light_dist_scale(1.0f);
	app.set_view_dist_scale(0.7f);
	app.set_display_bg_mesh(false);
	// s33
	app.set_res_file(rf, "geostatic", Hdf5Field::s33);
	app.set_color_map_fld_range(-56000.0, 0.0);
	// shear stress
	//app.set_res_file(rf, "geostatic", Hdf5Field::max_shear_stress);
	//app.set_color_map_fld_range(0.0, 5000.0);
	// plastic mises strain
	//app.set_res_file(rf, "geostatic", Hdf5Field::plastic_mises_strain_2d);
	//app.set_color_map_fld_range(0.0, 0.35);
	// p
	//app.set_res_file(rf, "geostatic", Hdf5Field::p);
	//app.set_color_map_fld_range(-1000.0, 1000.0);
	//
	app.set_color_map_geometry(1.2f, 0.4f, 0.45f);
	app.set_png_name("t3d_chm_mt_spudcan_geo");
	//app.set_gif_name("t3d_chm_mt_spudcan_geo");
	app.start();
}

void test_t3d_chm_mt_spudcan_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t3d_chm_mt_spudcan.h5");

	//QtApp_Posp_T3D_CHM_mt_Div<PlaneDivisionSet> app(argc, argv, QtApp_Posp_T3D_CHM_mt_Div<PlaneDivisionSet>::SingleFrame);
	////app.set_res_file(rf, "penetration", 50, Hdf5Field::s33);
	////app.set_res_file(rf, "penetration", 50, Hdf5Field::max_shear_stress);
	//app.set_res_file(rf, "penetration", 50, Hdf5Field::plastic_mises_strain_2d);
	QtApp_Posp_T3D_CHM_mt_Div<PlaneDivisionSet> app(argc, argv,
		QtApp_Posp_T3D_CHM_mt_Div<PlaneDivisionSet>::Animation);
	app.get_div_set().set_by_normal_and_point(0.0, 1.0, 0.0, 0.0, 1.0, 0.0);
	//QtApp_Posp_T3D_CHM_mt app(argc, argv, QtApp_Posp_T3D_CHM_mt::Animation);
	app.set_ani_time(5.0);
	app.set_win_size(1200, 800);
	app.set_view_dir(-90.0f, 10.0f);
	app.set_fog_coef(0.02f);
	app.set_light_dir(-135.0f, 20.0f);
	app.set_light_dist_scale(1.0f);
	app.set_view_dist_scale(0.7f);
	app.set_display_bg_mesh(false);
	// s33
	//app.set_res_file(rf, "penetration", Hdf5Field::s33);
	//app.set_color_map_fld_range(-56000.0, 0.0);
	// shear stress
	//app.set_res_file(rf, "penetration", Hdf5Field::max_shear_stress);
	//app.set_color_map_fld_range(0.0, 5000.0);
	// plastic mises strain
	//app.set_res_file(rf, "penetration", Hdf5Field::plastic_mises_strain_2d);
	//app.set_color_map_fld_range(0.0, 0.35);
	// p
	app.set_res_file(rf, "penetration", Hdf5Field::p);
	app.set_color_map_fld_range(-3000.0, 3000.0);
	//
	app.set_color_map_geometry(1.2f, 0.4f, 0.45f);
	//app.set_png_name("t3d_chm_mt_spudcan");
	app.set_gif_name("t3d_chm_mt_spudcan");
	app.start();
}
