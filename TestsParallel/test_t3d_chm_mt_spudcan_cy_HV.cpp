#include "TestsParallel_pcp.h"

#include "ParticleGenerator3D.hpp"
#include "Model_T3D_CHM_mt.h"
#include "ModelData_T3D_CHM_mt.h"
#include "Step_T3D_CHM_mt_Geo.h"
#include "TimeHistory_T3D_CHM_mt_Geo_complete.h"
#include "Step_T3D_CHM_mt.h"
#include "Step_T3D_CHM_ud_mt.h"
#include "TimeHistory_T3D_CHM_mt_complete.h"
#include "Step_T3D_CHM_TBB.h"
#include "TimeHistory_T3D_CHM_TBB_complete.h"
#include "Step_T3D_CHM_ud_TBB.h"
#include "TimeHistory_T3D_CHM_ud_TBB_complete.h"

#include "Model_T3D_CHM_up_mt.h"
#include "ModelData_T3D_CHM_up_mt.h"
#include "Step_T3D_CHM_up_TBB.h"
#include "TimeHistory_T3D_CHM_up_TBB_complete.h"

#include "TimeHistory_ConsoleProgressBar.h"
#include "QtApp_Prep_T3D_CHM_mt.h"
#include "QtApp_Prep_T3D_CHM_mt_Div.h"
#include "QtApp_Prep_T3D_CHM_up_mt_Div.h"
#include "test_parallel_utils.h"
#include "test_simulations_omp.h"

#define min_prin_stress 1000.0

void test_t3d_chm_mt_spudcan_cy_HV_model(int argc, char** argv)
{
	constexpr double footing_radius = 1.5;

	constexpr double cy_radius = 8.0 * footing_radius; // 6.0, 8.0
	constexpr double cy_coarse_radius = 4.5 * footing_radius; // 4.5
	constexpr double cy_top = 0.5 * footing_radius;
	constexpr double cy_depth = 8.0 * footing_radius; // 7.0, 8.0
	constexpr double cy_coarse_depth = 4.5 * footing_radius; // 4.5
	constexpr double cy_len = cy_top + cy_depth;
	constexpr double cy_len2 = cy_top + cy_coarse_depth;
	constexpr double dense_elem_size = 0.16666667 * footing_radius;
	constexpr double coarse_elem_size = 0.33333333 * footing_radius;
	constexpr double sml_pcl_size = dense_elem_size * 0.18;
	constexpr double lgr_pcl_size = coarse_elem_size * 0.18;

	TetrahedronMesh teh_mesh;
	//teh_mesh.load_mesh_from_hdf5("../../Asset/spudcan_soil_half_cylinder_3_5R.h5");
	//teh_mesh.init_search_grid(0.35, 0.35, 0.35);
	teh_mesh.load_mesh_from_hdf5("../../Asset/spudcan_soil_half_cylinder_8R.h5");
	teh_mesh.init_search_grid(0.2, 0.2, 0.2);
	std::cout << "node_num: " << teh_mesh.get_node_num() << "\n"
			  << "elem_num: " << teh_mesh.get_elem_num() << "\n";

	ParticleGenerator3D<TetrahedronMesh> pcl_generator;
	// dense area
	pcl_generator.generate_pcls_in_cylinder(
		0.0, 0.0, -cy_coarse_depth,
		0.0, 0.0, 1.0,
		0.0, cy_coarse_radius,
		-180.0, 0.0,
		cy_coarse_depth,
		sml_pcl_size, sml_pcl_size, sml_pcl_size);
	// coarse area
	pcl_generator.generate_pcls_in_cylinder(
		0.0, 0.0, -cy_coarse_depth,
		0.0, 0.0, 1.0,
		cy_coarse_radius, cy_radius,
		-180.0, 0.0,
		cy_coarse_depth,
		lgr_pcl_size, lgr_pcl_size, lgr_pcl_size);
	pcl_generator.generate_pcls_in_cylinder(
		0.0, 0.0, -cy_depth,
		0.0, 0.0, 1.0,
		0.0, cy_radius,
		-180.0, 0.0,
		cy_depth - cy_coarse_depth,
		lgr_pcl_size, lgr_pcl_size, lgr_pcl_size);
	//
	pcl_generator.adjust_pcl_size_to_fit_elems(teh_mesh);
	std::cout << "pcl_num: " << pcl_generator.get_num() << "\n";
	
	constexpr double e0 = 0.55;
	constexpr double den_grain = 2670.0;
	constexpr double den_dry = den_grain / (e0 + 1.0);
	constexpr double den_sat = den_grain / (e0 + 1.0) + 1000.0 * e0 / (e0 + 1.0);
	constexpr double den_float = den_sat - 1000.0;
	constexpr double n0 = e0 / (1.0 + e0);
	//Model_T3D_CHM_mt model;
	Model_T3D_CHM_up_mt model;
	model.init_mesh(teh_mesh);
	model.init_search_grid(teh_mesh);
	teh_mesh.clear();
	model.init_pcls(pcl_generator, n0, den_grain, 1000.0, 2.0e7, 0.0, 1.0e-3);
	pcl_generator.clear();

	const size_t pcl_num = model.get_pcl_num();
	MatModel::MaterialModel** mms = model.get_mat_models();
	constexpr double fric_ang = 30.0; // 30.02298846
	const double K0 = 1.0 - sin(fric_ang / 180.0 * 3.14159265359);
	double ini_stress[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	// Norsand
	constexpr double gamma = 0.875; // 0.875
	constexpr double lambda = 0.0058;
	constexpr double N = 0.3;
	constexpr double chi = 2.5;
	constexpr double H = 230.0; // 200.0
	constexpr double niu = 0.2;
	constexpr double Ig = 230.0; // 200.0
	MatModel::NorsandWrapper *ns = model.add_NorsandWrapper(pcl_num);
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		mms[pcl_id] = ns;
		double pcl_z = model.get_pcl_pos()[pcl_id].z;
		auto& pcl_s = model.get_pcl_stress0()[pcl_id];
		pcl_s.s33 = pcl_z * 9.81 * den_float;
		pcl_s.s22 = K0 * pcl_s.s33;
		pcl_s.s11 = pcl_s.s22;
		ini_stress[2] = pcl_s.s33 < -min_prin_stress ? pcl_s.s33 : -min_prin_stress;
		ini_stress[1] = pcl_s.s22 < -min_prin_stress ? pcl_s.s22 : -min_prin_stress;
		ini_stress[0] = pcl_s.s11 < -min_prin_stress ? pcl_s.s11 : -min_prin_stress;
		ns->set_param(
			ini_stress, e0, 
			fric_ang, gamma, lambda,
			N, chi, H, Ig, niu);
		ns->set_min_prin_s(min_prin_stress);
		ns = model.following_NorsandWrapper(ns);
	}

	model.init_t3d_rigid_mesh(1.0, "../../Asset/spudcan_model_flat_tip.h5",
		0.0, -1.0, 0.308, 65.0, 0.0, 0.0, 0.3, 0.3, 0.3);
	constexpr double K_cont = 1.0e6 / (sml_pcl_size * sml_pcl_size);
	model.set_contact_param(K_cont, K_cont, 0.36, 5.0);
	model.set_frictional_contact();

	model.set_t3d_rigid_mesh_velocity(0.0, 0.2, -0.05);
	model.set_k(1.0e-8);

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
	
	// x face bcs
	IndexArray vx_bc_pt_array(500);
	find_3d_nodes_on_x_plane(model, vx_bc_pt_array, 0.0);
	//model.init_fixed_vx_s_bc(vx_bc_pt_array.get_num(), vx_bc_pt_array.get_mem());
	//model.init_fixed_vx_f_bc(vx_bc_pt_array.get_num(), vx_bc_pt_array.get_mem());
	model.init_fixed_vx_bc(vx_bc_pt_array.get_num(), vx_bc_pt_array.get_mem());

	// y face bcs
	//IndexArray vy_bc_pt_array(500);
	//find_3d_nodes_on_y_plane(model, vy_bc_pt_array, 0.0);
	//model.init_fixed_vy_s_bc(vy_bc_pt_array.get_num(), vy_bc_pt_array.get_mem());
	//model.init_fixed_vy_f_bc(vy_bc_pt_array.get_num(), vy_bc_pt_array.get_mem());

	// arc bcs
	IndexArray arc_bc_pt_array(1000);
	//find_3d_nodes_on_cylinder(model, arc_bc_pt_array,
	//	Vector3D(0.0, 0.0, -cy_coarse_depth), Vector3D(0.0, 0.0, 1.0), cy_coarse_radius, cy_len2);
	find_3d_nodes_on_cylinder(model, arc_bc_pt_array,
		Vector3D(0.0, 0.0, -cy_depth), Vector3D(0.0, 0.0, 1.0), cy_radius, cy_len);
	const auto* node_pos = model.get_node_pos();
	for (size_t n_id = 0; n_id < arc_bc_pt_array.get_num(); n_id++)
	{
		const size_t nid = arc_bc_pt_array[n_id];
		const auto& np = node_pos[nid];
		model.set_vbc_vec(nid, np.x, np.y, 0.0);
	}

	// bottom bcs
	IndexArray vz_bc_pt_array(100);
	//find_3d_nodes_on_z_plane(model, vz_bc_pt_array, -cy_coarse_depth);
	find_3d_nodes_on_z_plane(model, vz_bc_pt_array, -cy_depth);
	model.init_fixed_vz_bc(vz_bc_pt_array.get_num(), vz_bc_pt_array.get_mem());

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_chm_mt_spudcan_cy_HV_model.h5");

	ModelData_T3D_CHM_up_mt md;
	md.output_model(model, res_file_hdf5);

	QtApp_Prep_T3D_CHM_up_mt_Div<PlaneDivisionSet> md_disp(argc, argv);
	//md_disp.get_div_set().set_by_normal_and_point(-1.0, 0.0, 0.0, -0.5, 0.0, 0.0);
	md_disp.set_model(model);
	md_disp.set_win_size(1200, 950);
	md_disp.set_view_dir(70.0f, -30.0f);
	md_disp.set_light_dir(80.0f, -25.0f);
	md_disp.set_display_bg_mesh(true);
	md_disp.set_view_dist_scale(0.6);
	//md_disp.set_pts_from_vx_bc(0.04);
	md_disp.set_pts_from_vz_bc(0.04);
	//md_disp.set_pts_from_vec_bc(0.04);
	//md_disp.set_pts_from_vx_bc_f(0.04);
	//md_disp.set_pts_from_vz_bc_f(0.04);
	md_disp.start();
}

void test_t3d_chm_mt_spudcan_cy_HV_geostatic(int argc, char** argv)
{
	Model_T3D_CHM_up_mt model;
	Model_T3D_CHM_up_mt_hdf5_utilities::load_model_from_hdf5_file(
		model, "t3d_chm_mt_spudcan_cy_HV_model.h5");

	// set tension cut-off surface
	//const size_t pcl_num = model.get_pcl_num();
	//MatModel::MaterialModel** mms = model.get_mat_models();
	//for (size_t p_id = 0; p_id < pcl_num; p_id++)
	//	((MatModel::NorsandWrapper*)mms[p_id])->set_min_prin_s(min_prin_stress);
	
	// contact
	constexpr double K_cont = 2.0e9; // 1e8
	//model.set_contact_param(K_cont, K_cont, 0.2, 5.0, K_cont / 50.0, K_cont / 50.0);
	model.set_contact_param(K_cont, K_cont, 0.2, 5.0);
	model.set_frictional_contact();

	// modified velocity and permeability
	model.set_t3d_rigid_mesh_velocity(0.0, 12.0, -3.0); // 0.4, -0.1
	model.set_k(1.5e-11);
	
	// -100.0e3 - 0m; -200.0e3 - 10m, -500.0e3 - 40m
	//model.set_cavitation(100.0, -1600.0e3, 0.05, 0.0, 0.0, 10.0e3);

	//QtApp_Prep_T3D_CHM_up_mt md_disp(argc, argv);
	////QtApp_Prep_T3D_CHM_mt_Div<PlaneDivisionSet> md_disp(argc, argv);
	////md_disp.get_div_set().set_by_normal_and_point(1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	//md_disp.set_model(model);
	//md_disp.set_win_size(1200, 950);
	//md_disp.set_view_dir(0.0f, 5.0f);
	//md_disp.set_light_dir(10.0f, 5.0f);
	//md_disp.set_display_bg_mesh(false);
	//md_disp.set_view_dist_scale(1.2);
	////md_disp.set_pts_from_vx_s_bc(0.05);
	////md_disp.set_pts_from_vy_s_bc(0.05);
	////md_disp.set_pts_from_vz_s_bc(0.05);
	////md_disp.set_pts_from_vx_f_bc(0.05);
	////md_disp.set_pts_from_vy_f_bc(0.05);
	////md_disp.set_pts_from_vz_f_bc(0.05);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_chm_mt_spudcan_cy_HV_geo.h5");

	ModelData_T3D_CHM_up_mt md;
	md.output_model(model, res_file_hdf5);

	//TimeHistory_T3D_CHM_mt_Geo_complete out1("geostatic");
	TimeHistory_T3D_CHM_up_TBB_complete out1("geostatic");
	out1.set_interval_num(30);
	out1.set_output_init_state();
	out1.set_output_final_state();
	out1.set_res_file(res_file_hdf5);
	TimeHistory_ConsoleProgressBar out_cpb;
	out_cpb.set_interval_num(2000);

	//Step_T3D_CHM_ud_mt step("step1");
	Step_T3D_CHM_up_TBB step("step1");
	step.set_model(model);
	step.set_thread_num(27);
	step.set_step_time(0.16); // 2.5
	step.set_dtime(1.0e-5); // 1.0e-5
	step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();
}

void test_t3d_chm_mt_spudcan_cy_HV(int argc, char** argv)
{
	Model_T3D_CHM_up_mt model;
#ifndef Undrained
	//Step_T3D_CHM_mt step("step2");
	Step_T3D_CHM_up_TBB step("step2");
#else
	Step_T3D_CHM_ud_TBB step("step2");
#endif
	Model_T3D_CHM_up_mt_hdf5_utilities::load_model_from_hdf5_file(
		model, step, "t3d_chm_mt_spudcan_cy_HV_geo.h5", "geostatic", 31);

	// contact
	constexpr double K_cont = 2.0e9;
	model.set_contact_param(K_cont, K_cont, 0.2, 5.0);
	model.set_frictional_contact();

	// modified velocity and permeability
	model.set_t3d_rigid_mesh_velocity(0.0, 0.4, -0.1);

	model.set_k(2.0e-8);
	
	//QtApp_Prep_T3D_CHM_mt_Div<EmptyDivisionSet> md_disp(argc, argv);
	////QtApp_Prep_T3D_CHM_mt_Div<PlaneDivisionSet> md_disp(argc, argv);
	////md_disp.get_div_set().set_by_normal_and_point(1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	//md_disp.set_model(model);
	//md_disp.set_win_size(1200, 950);
	//md_disp.set_view_dir(-150.0f, 10.0f);
	//md_disp.set_light_dir(-130.0f, -15.0f);
	//md_disp.set_display_bg_mesh(false);
	//md_disp.set_view_dist_scale(0.75);
	////md_disp.set_pts_from_vx_s_bc(0.05);
	////md_disp.set_pts_from_vy_s_bc(0.05);
	////md_disp.set_pts_from_vz_s_bc(0.05);
	////md_disp.set_pts_from_vec_s_bc(0.05);
	//md_disp.set_pts_from_vx_f_bc(0.05);
	//md_disp.set_pts_from_vy_f_bc(0.05);
	//md_disp.set_pts_from_vz_f_bc(0.05);
	//md_disp.set_pts_from_vec_f_bc(0.05);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_chm_mt_spudcan_cy_HV.h5");

	ModelData_T3D_CHM_up_mt md;
	md.output_model(model, res_file_hdf5);

#ifndef Undrained
	//TimeHistory_T3D_CHM_mt_complete out1("penetration");
	TimeHistory_T3D_CHM_up_TBB_complete out1("penetration");
#else
	TimeHistory_T3D_CHM_ud_TBB_complete out1("penetration");
#endif
	out1.set_interval_num(30);
	out1.set_output_init_state();
	out1.set_output_final_state();
	out1.set_res_file(res_file_hdf5);
	TimeHistory_ConsoleProgressBar out_cpb;
	out_cpb.set_interval_num(2000);

	step.set_model(model);
	step.set_thread_num(27);
	step.set_step_time(0.75); // 2.25
	step.set_dtime(1.0e-5); // 1.0e-5
	step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();
}

#include "QtApp_Posp_T3D_CHM_mt.h"
#include "QtApp_Posp_T3D_CHM_mt_Div.h"
#include "QtApp_Posp_T3D_CHM_up_mt_Div.h"
#include "test_model_view_omp.h"

void test_t3d_chm_mt_spudcan_cy_HV_geo_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	//rf.open("E:\\HV_spudcan\\t3d_chm_mt_spudcan_cy_HV_geo_2e-8.h5");
	//rf.open("E:\\HV_spudcan\\t3d_chm_mt_spudcan_cy_HV_geo_3e-11.h5");
	//rf.open("E:\\HV_spudcan\\t3d_chm_mt_spudcan_cy_HV_geo_0.h5");
	//rf.open("E:\\HV_spudcan\\t3d_chm_mt_spudcan_cy_HV_geo_dyn_inf_04.h5");
	//rf.open("E:\\HV_spudcan\\t3d_chm_mt_spudcan_cy_HV_geo_dyn_inf_1.h5");
	//rf.open("E:\\HV_spudcan\\t3d_chm_mt_spudcan_cy_HV_geo_dyn_inf_2.h5");
	//rf.open("E:\\HV_spudcan\\t3d_chm_mt_spudcan_cy_HV_geo_dyn_inf_3.h5");
	//rf.open("F:\\t3d_chm_mt_spudcan_cy_HV_geo_dyn_60_3.h5");
	//rf.open("E:\\HV_spudcan\\t3d_chm_mt_spudcan_cy_HV_geo_dyn_400_3.h5");
	//rf.open("F:\\t3d_chm_mt_spudcan_cy_HV_geo_dyn_3000_3.h5");

	// cavitation
	rf.open("F:\\t3d_chm_mt_spudcan_cy_HV_geo_cav-1000.h5");

	QtApp_Posp_T3D_CHM_mt_Div<PlaneDivisionSet> app(argc, argv, QtApp_Posp_T3D_CHM_mt_Div<PlaneDivisionSet>::SingleFrame);
	//app.set_res_file(rf, "geostatic", 28, Hdf5Field::s33);
	//app.set_color_map_fld_range(-220000.0, 0.0);
	//app.set_res_file(rf, "geostatic", 28, Hdf5Field::p);
	//app.set_color_map_fld_range(-200000.0, 20000.0);
	app.set_res_file(rf, "geostatic", 28, Hdf5Field::mat_e);
	app.set_color_map_fld_range(0.5, 0.665);
	//app.set_res_file(rf, "geostatic", 28, Hdf5Field::is_cavitated);
	//app.set_color_map_fld_range(0.0, 1.0);

	//QtApp_Posp_T3D_CHM_up_mt_Div<PlaneDivisionSet> app(argc, argv, QtApp_Posp_T3D_CHM_up_mt_Div<PlaneDivisionSet>::Animation);
	//app.get_div_set().set_by_normal_and_point(-1.0, 0.0, 0.0, -0.5, 0.0, 0.0);
	app.set_ani_time(8.0);
	app.set_win_size(1200, 800);
	app.set_view_dir(0.0f, 0.0f);
	app.set_light_dir(10.0f, 5.0f);
	app.set_fog_coef(0.02f);
	app.move_view_pos(0.0, 0.0, 2.0);
	app.set_view_dist_scale(0.58f);
	app.set_display_bg_mesh(false);
	app.set_color_map_geometry(3.0f, 0.4f, 0.45f);
	// s33
	//app.set_res_file(rf, "geostatic", Hdf5Field::s33);
	//app.set_color_map_fld_range(-110000.0, 0.0);
	// pore
	//app.set_res_file(rf, "geostatic", Hdf5Field::p);
	//app.set_color_map_fld_range(-200000.0, 200000.0);
	// void ratio
	//app.set_res_file(rf, "geostatic", Hdf5Field::mat_e);
	//app.set_color_map_fld_range(0.5, 0.65);
	//
	app.set_png_name("t3d_chm_mt_spudcan_cy_HV_geo");
	//app.set_gif_name("t3d_chm_mt_spudcan_cy_HV_geo");
	app.start();
}

void test_t3d_chm_mt_spudcan_cy_HV_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	//rf.open("E:\\HV_spudcan\\t3d_chm_mt_spudcan_cy_HV_2e-8.h5");
	//rf.open("E:\\HV_spudcan\\t3d_chm_mt_spudcan_cy_HV_3e-11.h5");
	//rf.open("E:\\HV_spudcan\\t3d_chm_mt_spudcan_cy_HV_0.h5");
	//rf.open("E:\\HV_spudcan\\t3d_chm_mt_spudcan_cy_HV_2e-10.h5");
	rf.open("E:\\HV_spudcan\\t3d_chm_mt_spudcan_cy_HV_4e-12.h5");

	QtApp_Posp_T3D_CHM_mt_Div<PlaneDivisionSet> app(argc, argv, QtApp_Posp_T3D_CHM_mt_Div<PlaneDivisionSet>::SingleFrame);
	app.set_res_file(rf, "penetration", 27, Hdf5Field::s33);
	app.set_color_map_fld_range(-220000.0, 0.0);
	//app.set_res_file(rf, "penetration", 27, Hdf5Field::p);
	//app.set_color_map_fld_range(-200000.0, 20000.0);
	//app.set_res_file(rf, "penetration", 27, Hdf5Field::mat_e);
	//app.set_color_map_fld_range(0.5, 0.665);

	//QtApp_Posp_T3D_CHM_mt_Div<PlaneDivisionSet> app(argc, argv,
	//	QtApp_Posp_T3D_CHM_mt_Div<PlaneDivisionSet>::Animation);
	//QtApp_Posp_T3D_CHM_mt_Div<PlaneDivisionSet> app(argc, argv,
	//	QtApp_Posp_T3D_CHM_mt_Div<PlaneDivisionSet>::SingleFrame);
	//app.get_div_set().set_by_normal_and_point(-1.0, 0.0, 0.0, -0.5, 0.0, 0.0);
	app.set_ani_time(5.0);
	app.set_win_size(1200, 800);
	app.set_view_dir(0.0f, 0.0f);
	app.set_light_dir(10.0f, 5.0f);
	app.set_fog_coef(0.02f);
	app.set_view_dist_scale(0.58f);
	app.set_display_bg_mesh(false);
	app.move_view_pos(0.0, 0.0, 2.0);
	app.set_color_map_geometry(3.0f, 0.4f, 0.45f);
	// showing legend
	//app.move_view_pos(0.0, 0.0, 100.0);
	//app.set_color_map_geometry(0.1f, 0.1f, 0.8f);
	// s33
	//app.set_res_file(rf, "penetration", Hdf5Field::s33);
	//app.set_color_map_fld_range(-11000.0, 0.0);
	// pore
	//app.set_res_file(rf, "penetration", Hdf5Field::p);
	//app.set_color_map_fld_range(-200000.0, 200000.0);
	// void ratio
	//app.set_res_file(rf, "penetration", Hdf5Field::mat_e);
	//app.set_color_map_fld_range(0.5, 0.65);
	//
	app.set_png_name("t3d_chm_mt_spudcan_cy_HV");
	//app.set_gif_name("t3d_chm_mt_spudcan_cy_HV");
	app.start();
}
