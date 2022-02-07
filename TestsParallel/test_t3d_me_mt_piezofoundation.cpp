#include "TestsParallel_pcp.h"

#include "ParticleGenerator3D.hpp"
#include "Model_T3D_ME_mt.h"
#include "ModelData_T3D_ME_mt.h"
#include "Step_T3D_ME_mt_Geo.h"
#include "TimeHistory_T3D_ME_mt_Geo_complete.h"
#include "Step_T3D_ME_mt.h"
#include "TimeHistory_T3D_ME_mt_complete.h"
#include "Step_T3D_ME_TBB.h"
#include "TimeHistory_T3D_ME_TBB_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "QtApp_Prep_T3D_ME_mt.h"
#include "test_parallel_utils.h"
#include "test_simulations_omp.h"

void test_t3d_me_mt_piezofoundation_model(int argc, char** argv)
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

	constexpr double e0 = 0.66;
	constexpr double den_grain = 2670.0;
	constexpr double den_sat = den_grain / (e0 + 1.0) + 1000 * e0 / (e0 + 1.0);
	constexpr double den_float = den_sat - 1000.0;
	constexpr double stress_depth_limit = -0.02;
	Model_T3D_ME_mt model;
	model.init_mesh(teh_mesh);
	model.init_search_grid(teh_mesh);
	teh_mesh.clear();
	model.init_pcls(pcl_generator, den_sat);
	pcl_generator.clear();

	const size_t pcl_num = model.get_pcl_num();
	MatModel::MaterialModel** mms = model.get_mat_models();
	const double K0 = 1.0 - sin(30.0 / 180.0 * 3.14159265359);
	double ini_stress[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	// Mohr Coulomb
	//MatModel::MohrCoulombWrapper *mcs = model.add_MohrCoulombWrapper(pcl_num);
	//for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	//{
	//	//const double pcl_z = -1.0; //debug
	//	//const double pcl_z = model.get_pcl_pos()[pcl_id].z - 1.0;
	//	double pcl_z = model.get_pcl_pos()[pcl_id].z;
	//	auto &pcl_s = model.get_pcl_stress0()[pcl_id];
	//	pcl_s.s33 = pcl_z * 9.81 * den_float;
	//	pcl_s.s22 = K0 * pcl_s.s33;
	//	pcl_s.s11 = pcl_s.s22;
	//	ini_stress[2] = pcl_s.s33;
	//	ini_stress[0] = pcl_s.s22;
	//	ini_stress[1] = pcl_s.s11;
	//	MatModel::MohrCoulombWrapper &mc = mcs[pcl_id];
	//	mc.set_param(ini_stress, 30.0, 0.0, 5.0, 1.0e6, 0.15);
	//	mms[pcl_id] = &mc;
	//}
	// Sand hypoplasticity
	MatModel::SandHypoplasticityStbWrapper* shps = model.add_SandHypoplasticityStbWrapper(pcl_num);
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		mms[pcl_id] = shps;
		
		double pcl_z = model.get_pcl_pos()[pcl_id].z;
		auto &pcl_s = model.get_pcl_stress0()[pcl_id];
		pcl_s.s33 = pcl_z * 9.81 * den_float;
		pcl_s.s22 = K0 * pcl_s.s33;
		pcl_s.s11 = pcl_s.s22;
		if (pcl_z > stress_depth_limit) // shallow depth
			pcl_z = stress_depth_limit;
		ini_stress[2] = pcl_z * 9.81 * den_float;
		ini_stress[0] = K0 * ini_stress[2];
		ini_stress[1] = ini_stress[0];
		shps->set_param(
			ini_stress, e0,
			30.0, 1354.0e6, 0.34,
			0.18, 1.27,
			0.49, 0.76, 0.86,
			1.5, 43.0, 100.0, // 43.0, 45.0
			200.0, 0.2);
		shps = model.following_SandHypoplasticityStbWrapper(shps);
	}

	// gravity force, float unit weight
	IndexArray bfz_pcl_array(pcl_num);
	MemoryUtils::ItemArray<double> bfz_array(pcl_num);
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		double bfz = -9.81 * den_float / den_sat;
		bfz_pcl_array.add(pcl_id);
		bfz_array.add(bfz);
	}
	model.init_bfzs(pcl_num, bfz_pcl_array.get_mem(), bfz_array.get_mem());

	model.init_rigid_cylinder(0.0, 0.0, 1.125, 2.25, 1.0, 2000.0);
	// velocity bc
	model.set_rigid_cylinder_velocity(0.0, 0.0, -1.5); //1.5
	// free motion
	//model.set_rigid_cylinder_ext_force(0.0, 0.0, -10.0 * model.get_rigid_cylinder_m());
	//model.fix_rigid_cylinder_vx();
	//model.fix_rigid_cylinder_vy();
	//model.set_rigid_cylinder_init_velocity(0.0, 0.0, -2.0);
	//
	model.set_contact_param(400.0 / (sml_pcl_size * sml_pcl_size),
		400.0 / (sml_pcl_size * sml_pcl_size), 0.1, 5.0);

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
	res_file_hdf5.create("t3d_me_mt_piezofoundation_model.h5");

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

void test_t3d_me_mt_piezofoundation_geo(int argc, char** argv)
{
	Model_T3D_ME_mt model;
	Model_T3D_ME_mt_hdf5_utilities::load_me_mt_model_from_hdf5_file(
		model, "t3d_me_mt_piezofoundation_model.h5");

	//QtApp_Prep_T3D_ME_mt md_disp(argc, argv);
	//md_disp.set_model(model);
	//md_disp.set_win_size(1200, 950);
	//md_disp.set_view_dir(-90.0f, 30.0f);
	//md_disp.set_light_dir(-60.0f, 15.0f);
	//md_disp.set_display_bg_mesh(false);
	////md_disp.set_pts_from_vx_bc(0.05);
	//md_disp.set_pts_from_vy_bc(0.05);
	////md_disp.set_pts_from_vz_bc(0.05);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_me_mt_piezofoundation_geo.h5");

	ModelData_T3D_ME_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T3D_ME_mt_Geo_complete out1("geostatic");
	out1.set_interval_num(20);
	out1.set_output_init_state();
	out1.set_output_final_state();
	out1.set_res_file(res_file_hdf5);
	TimeHistory_ConsoleProgressBar out_cpb;
	out_cpb.set_interval_num(2000);

	Step_T3D_ME_mt_Geo step("step1");
	step.set_model(model);
	step.set_thread_num(22);
	step.set_step_time(0.5);
	//step.set_thread_num(4);
	//step.set_step_time(5.0e-5);
	step.set_dtime(2.5e-5);
	step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();
}

void test_t3d_me_mt_piezofoundation(int argc, char** argv)
{
	Model_T3D_ME_mt model;
	//Step_T3D_ME_mt step("step2");
	Step_T3D_ME_TBB step("step2");
	Model_T3D_ME_mt_hdf5_utilities::load_me_mt_model_from_hdf5_file(
		model, step, "t3d_me_mt_piezofoundation_geo_N43.h5", "geostatic", 21); // 21
	
	std::cout << "Load model completed.\n";
	
	//QtApp_Prep_T3D_ME_mt md_disp(argc, argv);
	//md_disp.set_model(model);
	//md_disp.set_win_size(1200, 950);
	//md_disp.set_view_dir(-90.0f, -20.0f);
	//md_disp.set_light_dir(-60.0f, 15.0f);
	//md_disp.set_display_bg_mesh(false);
	////md_disp.set_view_dist_scale(0.6);
	//md_disp.set_pts_from_vx_bc(0.05);
	////md_disp.set_pts_from_vy_bc(0.05);
	////md_disp.set_pts_from_vz_bc(0.05);
	//md_disp.start();
	//return;

	model.set_rigid_cylinder_velocity(0.0, 0.0, -0.5);
	constexpr double sml_pcl_size = 0.03125;
	model.set_contact_param(20000.0 / (sml_pcl_size * sml_pcl_size),
		20000.0 / (sml_pcl_size * sml_pcl_size), 0.1, 5.0);
	
	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_me_mt_piezofoundation.h5");

	ModelData_T3D_ME_mt md;
	md.output_model(model, res_file_hdf5);
	std::cout << "Output model completed.\n";

	TimeHistory_ConsoleProgressBar out_cpb;
	out_cpb.set_interval_num(2000);
	
	// omp
	//TimeHistory_T3D_ME_mt_complete out1("penetration");
	//out1.set_interval_num(100);
	//out1.set_output_init_state();
	//out1.set_output_final_state();
	//out1.set_res_file(res_file_hdf5);
	// tbb
	TimeHistory_T3D_ME_TBB_complete out1("penetration");
	out1.set_interval_num(150);
	out1.set_output_init_state();
	out1.set_output_final_state();
	out1.set_res_file(res_file_hdf5);

	std::cout << "Start solving...\n";
	step.set_thread_num(24);
	step.set_step_time(1.4); // 0.9
	step.set_dtime(2.0e-6); 
	step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();
}

void test_t3d_me_mt_piezofoundation2(int argc, char** argv)
{
	Model_T3D_ME_mt model;
	//Step_T3D_ME_mt step("step3");
	Step_T3D_ME_TBB step("step3");
	Model_T3D_ME_mt_hdf5_utilities::load_me_mt_model_from_hdf5_file(
		model, step, "t3d_me_mt_piezofoundation.h5", "penetration", 101); // 21
	std::cout << "Load model completed.\n";

	//QtApp_Prep_T3D_ME_mt md_disp(argc, argv);
	//md_disp.set_model(model);
	//md_disp.set_win_size(1200, 950);
	//md_disp.set_view_dir(-90.0f, -20.0f);
	//md_disp.set_light_dir(-60.0f, 15.0f);
	//md_disp.set_display_bg_mesh(false);
	////md_disp.set_view_dist_scale(0.6);
	//md_disp.set_pts_from_vx_bc(0.05);
	////md_disp.set_pts_from_vy_bc(0.05);
	////md_disp.set_pts_from_vz_bc(0.05);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_me_mt_piezofoundation2.h5");

	ModelData_T3D_ME_mt md;
	md.output_model(model, res_file_hdf5);
	std::cout << "Output model completed.\n";

	TimeHistory_ConsoleProgressBar out_cpb;
	out_cpb.set_interval_num(2000);

	//TimeHistory_T3D_ME_mt_complete out1("penetration");
	TimeHistory_T3D_ME_TBB_complete out1("penetration");
	out1.set_interval_num(100);
	out1.set_output_init_state();
	out1.set_output_final_state();
	out1.set_res_file(res_file_hdf5);

	std::cout << "Start solving...\n";
	step.set_thread_num(24);
	step.set_step_time(0.5);
	step.set_dtime(1.0e-6); // 2.0e-6
	step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();
}

#include "QtApp_Posp_T3D_ME_mt_Div.h"
#include "test_model_view_omp.h"

void test_t3d_me_mt_piezofoundation_geo_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t3d_me_mt_piezofoundation_geo.h5");

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

	// Animation
	QtApp_Posp_T3D_ME_mt_Div<PlaneDivisionSet> app(argc, argv,
		QtApp_Posp_T3D_ME_mt_Div<PlaneDivisionSet>::Animation);
	app.get_div_set().set_by_normal_and_point(0.0, 1.0, 0.0, 0.0, 0.06, 0.0);
	app.set_ani_time(5.0);
	app.set_win_size(1600, 950);
	app.set_view_dir(-90.0f, 10.0f);
	app.set_light_dir(-135.0f, 20.0f);
	app.set_light_dist_scale(1.0f);
	app.set_fog_coef(0.02f);
	app.set_view_dist_scale(0.85);
	app.set_display_bg_mesh(false);
	//app.set_mono_color_pcl(true);
	// s33
	app.set_res_file(rf, "geostatic", Hdf5Field::s33);
	app.set_color_map_fld_range(-8.0e4, 0.0);
	// shear stress
	//app.set_res_file(rf, "geostatic", Hdf5Field::max_shear_stress);
	//app.set_color_map_fld_range(0.0, 5.0);
	// mises strain
	//app.set_res_file(rf, "geostatic", Hdf5Field::plastic_mises_strain_2d);
	//app.set_color_map_fld_range(0.0, 0.01);
	// mat_e
	//app.set_res_file(rf, "geostatic", Hdf5Field::mat_e);
	//app.set_color_map_fld_range(0.5, 0.8);
	// mat_s11
	//app.set_res_file(rf, "geostatic", Hdf5Field::mat_s11);
	//app.set_color_map_fld_range(-4.0e4, 0.0);
	// mat_s33
	//app.set_res_file(rf, "geostatic", Hdf5Field::mat_s33);
	//app.set_color_map_fld_range(-8.0e4, 0.0);
	//
	app.set_color_map_geometry(1.2f, 0.4f, 0.45f);
	//app.set_png_name("t3d_me_mt_piezofoundation_geo");
	app.set_gif_name("t3d_me_mt_piezofoundation_geo");
	app.start();
}

void test_t3d_me_mt_piezofoundation_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t3d_me_mt_piezofoundation.h5");
	//rf.open("t3d_me_mt_piezofoundation2.h5");

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

	// Animation
	QtApp_Posp_T3D_ME_mt_Div<PlaneDivisionSet> app(argc, argv,
		QtApp_Posp_T3D_ME_mt_Div<PlaneDivisionSet>::Animation);
	app.get_div_set().set_by_normal_and_point(0.0, 1.0, 0.0, 0.0, 0.06, 0.0);
	app.set_ani_time(5.0);
	app.set_win_size(1600, 950);
	app.set_view_dir(-90.0f, 10.0f);
	app.set_light_dir(-135.0f, 20.0f);
	app.set_light_dist_scale(1.0f);
	app.set_fog_coef(0.02f);
	app.set_view_dist_scale(0.85);
	app.set_display_bg_mesh(false);
	//app.set_mono_color_pcl(true);
	// s33
	//app.set_res_file(rf, "penetration", Hdf5Field::s33);
	//app.set_color_map_fld_range(-8.0e4, 0.0);
	// shear stress
	//app.set_res_file(rf, "penetration", Hdf5Field::max_shear_stress);
	//app.set_color_map_fld_range(0.0, 5.0);
	// mises strain
	//app.set_res_file(rf, "penetration", Hdf5Field::plastic_mises_strain_2d);
	//app.set_color_map_fld_range(0.0, 0.08);
	// mat_e
	//app.set_res_file(rf, "penetration", Hdf5Field::mat_e);
	//app.set_color_map_fld_range(0.6, 0.8);
	// mat_s11
	//app.set_res_file(rf, "penetration", Hdf5Field::mat_s11);
	//app.set_color_map_fld_range(-4.0e4, 0.0);
	// mat_s33
	app.set_res_file(rf, "penetration", Hdf5Field::mat_s33);
	app.set_color_map_fld_range(-1.0e5, 0.0);
	//
	app.set_color_map_geometry(1.2f, 0.4f, 0.45f);
	//app.set_png_name("t3d_me_mt_piezofoundation");
	app.set_gif_name("t3d_me_mt_piezofoundation");

	//QtApp_Posp_T3D_ME_mt_Div<PlaneDivisionSet> app(argc, argv,
	//	QtApp_Posp_T3D_ME_mt_Div<PlaneDivisionSet>::SingleFrame);
	//app.get_div_set().set_by_normal_and_point(0.0, 1.0, 0.0, 0.0, 0.06, 0.0);
	//app.set_win_size(1600, 950);
	//app.set_view_dir(-90.0f, 10.0f);
	//app.set_light_dir(-135.0f, 20.0f);
	//app.set_light_dist_scale(1.0f);
	//app.set_fog_coef(0.02f);
	//app.set_view_dist_scale(0.85);
	//app.set_display_bg_mesh(false);
	//app.set_color_map_geometry(1.2f, 0.4f, 0.45f);
	//// mat_e
	//app.set_res_file(rf, "penetration", 101, Hdf5Field::mat_e);
	//app.set_color_map_fld_range(0.6, 0.8);
	////
	////app.set_png_name("t3d_me_mt_piezofoundation2");

	app.start();
}
