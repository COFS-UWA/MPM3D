#include "TestsParallel_pcp.h"

#include "ParticleGenerator3D.hpp"
#include "Model_T3D_ME_mt.h"
#include "ModelData_T3D_ME_mt.h"
#include "Step_T3D_ME_TBB.h"
#include "TimeHistory_T3D_ME_TBB_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "QtApp_Prep_T3D_ME_mt.h"
#include "test_parallel_utils.h"
#include "test_simulations_omp.h"

void test_t3d_me_mt_spudcan_cy_Hossain_2006_model(int argc, char** argv)
{
	constexpr double footing_radius = 1.5;

	constexpr double cy_radius = 8.0 * footing_radius; // 6.0, 8.0
	constexpr double cy_coarse_radius = 3.5 * footing_radius;
	constexpr double cy_top = 1.0 * footing_radius;
	constexpr double cy_depth = 10.0 * footing_radius; // 7.0, 8.0
	constexpr double cy_coarse_depth = 5.5 * footing_radius;
	constexpr double cy_len = cy_top + cy_depth;
	constexpr double dense_elem_size = 0.16 * footing_radius;
	constexpr double coarse_elem_size = 0.30 * footing_radius;
	constexpr double sml_pcl_size = dense_elem_size * 0.2;
	constexpr double lgr_pcl_size = coarse_elem_size * 0.2;

	TetrahedronMesh teh_mesh;
	teh_mesh.load_mesh_from_hdf5("../../Asset/spudcan_soil_quarter_Hossain_4D.h5");
	teh_mesh.init_search_grid(0.1, 0.1, 0.1);
	std::cout << "node_num: " << teh_mesh.get_node_num() << "\n"
			  << "elem_num: " << teh_mesh.get_elem_num() << "\n";

	ParticleGenerator3D<TetrahedronMesh> pcl_generator;
	// dense area
	pcl_generator.generate_pcls_in_cylinder(
		0.0, 0.0, -cy_coarse_depth,
		0.0, 0.0, 1.0,
		0.0, cy_coarse_radius,
		0.0, 90.0,
		cy_coarse_depth,
		sml_pcl_size, sml_pcl_size, sml_pcl_size);
	// coarse area
	pcl_generator.generate_pcls_in_cylinder(
		0.0, 0.0, -cy_coarse_depth,
		0.0, 0.0, 1.0,
		cy_coarse_radius, cy_radius,
		0.0, 90.0,
		cy_coarse_depth,
		lgr_pcl_size, lgr_pcl_size, lgr_pcl_size);
	pcl_generator.generate_pcls_in_cylinder(
		0.0, 0.0, -cy_depth,
		0.0, 0.0, 1.0,
		0.0, cy_radius,
		0.0, 90.0,
		cy_depth - cy_coarse_depth,
		lgr_pcl_size, lgr_pcl_size, lgr_pcl_size);
	//
	pcl_generator.adjust_pcl_size_to_fit_elems(teh_mesh);
	std::cout << "pcl_num: " << pcl_generator.get_num() << "\n";
	
	constexpr double den_sat = 1700.0;
	Model_T3D_ME_mt model;
	model.init_mesh(teh_mesh);
	model.init_search_grid(teh_mesh);
	teh_mesh.clear();
	model.init_pcls(pcl_generator, den_sat);
	pcl_generator.clear();

	const size_t pcl_num = model.get_pcl_num();
	MatModel::MaterialModel** mms = model.get_mat_models();
	double ini_stress[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	constexpr double K0 = 0.47 / (1.0 - 0.47);
	// Tresca
	MatModel::Tresca* tcs = model.add_Tresca(pcl_num);
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		double pcl_z = model.get_pcl_pos()[pcl_id].z;
		auto& pcl_s = model.get_pcl_stress0()[pcl_id];
		pcl_s.s33 = pcl_z * 9.81 * den_sat;
		pcl_s.s22 = K0 * pcl_s.s33;
		pcl_s.s11 = pcl_s.s22;
		ini_stress[2] = pcl_s.s33;
		ini_stress[1] = pcl_s.s22;
		ini_stress[0] = pcl_s.s11;
		// 
		mms[pcl_id] = tcs;
		const double su_clay = 5.0e3 + 6.66667e3 * (-pcl_z);
		const double E_clay = 200.0 * su_clay;
		tcs->set_param(E_clay, 0.47, su_clay, ini_stress);
		mms[pcl_id] = tcs;
		tcs = model.following_Tresca(tcs);
	}

	model.init_t3d_rigid_mesh(1.0, "../../Asset/spudcan_model_Hossain_2006.h5",
		0.0, 0.0, 0.0, 90.0, 0.0, 0.0, 0.3, 0.3, 0.3);
	model.set_t3d_rigid_mesh_velocity(0.0, 0.0, -0.25);
	constexpr double K_cont = 1.0e7 / (sml_pcl_size * sml_pcl_size);
	model.set_contact_param(K_cont, K_cont, 0.36, 5.0);
	model.set_smooth_contact_between_pcl_and_rect();

	// gravity force, float unit weight
	IndexArray bfz_pcl_array(pcl_num);
	MemoryUtils::ItemArray<double> bfz_array(pcl_num);
	double bfz = -9.81;
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		bfz_pcl_array.add(pcl_id);
		bfz_array.add(bfz);
	}
	model.init_bfzs(pcl_num, bfz_pcl_array.get_mem(), bfz_array.get_mem());
	
	// x face bcs
	IndexArray vx_bc_pt_array(500);
	find_3d_nodes_on_x_plane(model, vx_bc_pt_array, 0.0);
	model.init_fixed_vx_bc(vx_bc_pt_array.get_num(), vx_bc_pt_array.get_mem());

	// y face bcs
	IndexArray vy_bc_pt_array(500);
	find_3d_nodes_on_y_plane(model, vy_bc_pt_array, 0.0);
	model.init_fixed_vy_bc(vy_bc_pt_array.get_num(), vy_bc_pt_array.get_mem());

	// arc bcs
	IndexArray arc_bc_pt_array(1000);
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
	find_3d_nodes_on_z_plane(model, vz_bc_pt_array, -cy_depth);
	model.init_fixed_vz_bc(vz_bc_pt_array.get_num(), vz_bc_pt_array.get_mem());

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_me_mt_spudcan_cy_model.h5");

	ModelData_T3D_ME_mt md;
	md.output_model(model, res_file_hdf5);

	QtApp_Prep_T3D_ME_mt md_disp(argc, argv);
	md_disp.set_model(model);
	md_disp.set_win_size(1200, 950);
	md_disp.set_view_dir(-80.0f, 20.0f);
	md_disp.set_light_dir(-70.0f, 15.0f);
	md_disp.set_display_bg_mesh(false);
	md_disp.set_view_dist_scale(0.9);
	//md_disp.set_pts_from_vx_bc(0.04);
	//md_disp.set_pts_from_vy_bc(0.04);
	//md_disp.set_pts_from_vz_bc(0.04);
	//md_disp.set_pts_from_vec_bc(0.04);
	md_disp.start();
}

void test_t3d_me_mt_spudcan_cy_Hossain_2006(int argc, char** argv)
{
	Model_T3D_ME_mt model;
	Model_T3D_ME_mt_hdf5_utilities::load_me_mt_model_from_hdf5_file(
		model, "t3d_me_mt_spudcan_cy_model.h5");
	
	constexpr double footing_radius = 1.5;
	constexpr double dense_elem_size = 0.16 * footing_radius;
	constexpr double sml_pcl_size = dense_elem_size * 0.2;
	constexpr double K_cont = 1.0e7 / (sml_pcl_size * sml_pcl_size);
	model.set_contact_param(K_cont, K_cont, 0.2, 5.0);
	model.set_smooth_contact_between_pcl_and_rect();

	model.set_t3d_rigid_mesh_velocity(0.0, 0.0, -0.3);
	
	//QtApp_Prep_T3D_ME_mt md_disp(argc, argv);
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
	res_file_hdf5.create("t3d_me_mt_spudcan_cy.h5");

	ModelData_T3D_ME_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T3D_ME_TBB_complete out1("penetration");
	out1.set_interval_num(100);
	out1.set_output_init_state();
	out1.set_output_final_state();
	out1.set_res_file(res_file_hdf5);
	TimeHistory_ConsoleProgressBar out_cpb;
	out_cpb.set_interval_num(2000);

	Step_T3D_ME_TBB step("step2");
	step.set_model(model);
	step.set_thread_num(31);
	step.set_step_time(8.0);
	step.set_dtime(5.0e-5); // 5.0e-6
	step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();
}

void test_t3d_me_mt_spudcan_cy_Hossain_2006_restart(int argc, char** argv)
{
	Model_T3D_ME_mt model;
	Step_T3D_ME_TBB step("step2");
	Model_T3D_ME_mt_hdf5_utilities::load_me_mt_model_from_hdf5_file(
		model, step, "t3d_me_mt_spudcan_cy.h5", "penetration", 101);

	constexpr double footing_radius = 1.5;
	constexpr double dense_elem_size = 0.16 * footing_radius;
	constexpr double sml_pcl_size = dense_elem_size * 0.2;
	constexpr double K_cont = 1.0e7 / (sml_pcl_size * sml_pcl_size);
	model.set_contact_param(K_cont, K_cont, 0.2, 5.0);
	model.set_smooth_contact_between_pcl_and_rect();
	model.set_t3d_rigid_mesh_velocity(0.0, 0.0, -0.3);

	//QtApp_Prep_T3D_ME_mt md_disp(argc, argv);
	//md_disp.set_model(model);
	//md_disp.set_win_size(1200, 950);
	//md_disp.set_view_dir(-150.0f, 10.0f);
	//md_disp.set_light_dir(-130.0f, -15.0f);
	//md_disp.set_display_bg_mesh(false);
	//md_disp.set_view_dist_scale(0.75);
	//md_disp.set_pts_from_vx_bc(0.05);
	////md_disp.set_pts_from_vy_bc(0.05);
	////md_disp.set_pts_from_vz_bc(0.05);
	////md_disp.set_pts_from_vec_bc(0.05);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_me_mt_spudcan_cy2.h5");

	ModelData_T3D_ME_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T3D_ME_TBB_complete out1("penetration");
	out1.set_interval_num(100);
	out1.set_output_init_state();
	out1.set_output_final_state();
	out1.set_res_file(res_file_hdf5);
	TimeHistory_ConsoleProgressBar out_cpb;
	out_cpb.set_interval_num(2000);

	step.set_model(model);
	step.set_thread_num(31);
	step.set_step_time(6.0);
	step.set_dtime(5.0e-5);
	step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();
}

#include "QtApp_Posp_T3D_ME_mt_Div.h"
#include "test_model_view_omp.h"

void test_t3d_me_mt_spudcan_cy_Hossain_2006_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t3d_me_mt_spudcan_cy.h5");

	//QtApp_Posp_T3D_ME_mt_Div<PlaneDivisionSet> app(argc, argv, QtApp_Posp_T3D_ME_mt_Div<PlaneDivisionSet>::SingleFrame);
	//app.get_div_set().set_by_normal_and_point(0.0, 1.0, 0.0, 0.0, 0.1, 0.0);
	//app.set_res_file(rf, "penetration", 75, Hdf5Field::s33);
	//app.set_res_file(rf, "penetration", 50, Hdf5Field::max_shear_stress);
	//app.set_res_file(rf, "penetration", 50, Hdf5Field::plastic_mises_strain_2d);
	QtApp_Posp_T3D_ME_mt_Div<PlaneDivisionSet> app(argc, argv,
		QtApp_Posp_T3D_ME_mt_Div<PlaneDivisionSet>::Animation);
	app.get_div_set().set_by_normal_and_point(0.0, 1.0, 0.0, 0.0, 0.1, 0.0);
	app.set_ani_time(10.0);
	app.set_win_size(1600, 900);
	app.set_view_dir(-90.0f, 5.0f);
	app.set_fog_coef(0.02f);
	app.set_light_dir(-90.0f, 5.0f);
	app.set_light_dist_scale(1.0f);
	app.move_view_pos(0.0, 0.0, 4.0);
	app.set_view_dist_scale(0.5f);
	app.set_display_bg_mesh(false);
	//app.set_update_rb_pos();
	//app.set_bg_color(QVector3D(1.0, 1.0, 1.0));
	//app.set_color_map_char_color(0.0, 0.0, 0.0);
	app.set_color_map_geometry(3.5f, 0.4f, 0.45f);
	// s33
	app.set_res_file(rf, "penetration", Hdf5Field::s33);
	app.set_color_map_fld_range(-100.0e3, 0.0);
	//
	app.set_color_map_geometry(1.5f, 0.4f, 0.45f);
	//app.set_png_name("t3d_me_mt_spudcan_cy");
	app.set_gif_name("t3d_me_mt_spudcan_cy");
	app.start();
}
