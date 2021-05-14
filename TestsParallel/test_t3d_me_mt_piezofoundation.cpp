#include "TestsParallel_pcp.h"

#include "ParticleGenerator3D.hpp"
#include "Model_T3D_ME_mt.h"
#include "Step_T3D_ME_mt.h"
#include "ModelData_T3D_ME_mt.h"
#include "TimeHistory_T3D_ME_mt_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "QtApp_Prep_T3D_ME_mt.h"
#include "test_parallel_utils.h"
#include "test_simulations_omp.h"

void test_t3d_me_mt_piezofoundation_model(int argc, char** argv)
{
	TetrahedronMesh teh_mesh;
	teh_mesh.load_mesh_from_hdf5("../../Asset/piezofoundation_soil_quarter.h5");
	teh_mesh.init_search_grid(2.0e-3, 2.0e-3, 2.0e-3);
	std::cout << "node_num: " << teh_mesh.get_node_num() << "\n"
		<< "elem_num: " << teh_mesh.get_elem_num() << "\n";

	constexpr double depth = 0.14;
	constexpr double coarse_depth = 0.076;
	constexpr double width = 0.11;
	constexpr double coarse_width = 0.066;
	constexpr double sml_pcl_size = 0.625e-3;
	constexpr double lgr_pcl_size = 1.250e-3;
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

	constexpr double e0 = 0.63;
	constexpr double den_grain = 2670.0;
	const double den_dry = den_grain / (e0 + 1.0);
	const double den_float = den_dry - 1000.0;
	Model_T3D_ME_mt model;
	model.init_mesh(teh_mesh);
	model.init_search_grid(teh_mesh);
	teh_mesh.clear();
	model.init_pcls(pcl_generator, den_dry);
	pcl_generator.clear();

	const size_t pcl_num = model.get_pcl_num();
	MatModel::MaterialModel** mms = model.get_mat_models();
	// Sand hypoplasticity
	const double K0 = 1.0 - sin(30.0);
	double ini_stress[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	MatModel::SandHypoplasticityWrapper* shps = model.add_SandHypoplasticityWrapper(pcl_num);
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		const double pcl_z = model.get_pcl_pos()[pcl_id].z;
		ini_stress[2] = pcl_z * den_float * 50.0 * 9.8;
		ini_stress[0] = K0 * ini_stress[2];
		ini_stress[1] = ini_stress[0];
		MatModel::SandHypoplasticityWrapper& shp = shps[pcl_id];
		shp.set_param(ini_stress, e0,
			30.0, 1354.0e6, 0.34,
			0.49, 0.76, 0.86,
			0.18, 1.27);
		mms[pcl_id] = &shp;
	}

	model.init_rigid_cylinder(0.0, 0.0, 0.0225, 0.045, 0.02);
	model.set_rigid_cylinder_velocity(0.0, 0.0, -0.01);
	model.set_contact_param(20000.0 / (0.15 * 0.15), 20000.0 / (0.15 * 0.15), 0.1, 5.0);

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
	md_disp.set_view_dir(-90.0f, -30.0f);
	md_disp.set_light_dir(-60.0f, 15.0f);
	md_disp.set_display_bg_mesh(false);
	//md_disp.set_view_dist_scale(0.7);
	//md_disp.set_pts_from_vx_bc(9.0e-4);
	md_disp.set_pts_from_vy_bc(9.0e-4);
	//md_disp.set_pts_from_vz_bc(9.0e-4);
	md_disp.start();
}

void test_t3d_me_mt_piezofoundation(int argc, char** argv)
{
	Model_T3D_ME_mt model;
	Model_T3D_ME_mt_hdf5_utilities::load_me_mt_model_from_hdf5_file(
		model, "t3d_me_mt_piezofoundation_model.h5");

	//QtApp_Prep_T3D_ME_mt md_disp(argc, argv);
	//md_disp.set_model(model);
	//md_disp.set_win_size(1200, 950);
	//md_disp.set_view_dir(-90.0f, -30.0f);
	//md_disp.set_light_dir(-60.0f, 15.0f);
	//md_disp.set_display_bg_mesh(false);
	////md_disp.set_view_dist_scale(0.6);
	////md_disp.set_pts_from_vx_bc(5.0e-4);
	////md_disp.set_pts_from_vy_bc(5.0e-4);
	//md_disp.set_pts_from_vz_bc(5.0e-4);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_me_mt_piezofoundation.h5");

	ModelData_T3D_ME_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T3D_ME_mt_complete out1("penetration");
	out1.set_interval_num(100);
	out1.set_output_init_state();
	out1.set_output_final_state();
	out1.set_res_file(res_file_hdf5);
	TimeHistory_ConsoleProgressBar out_cpb;
	out_cpb.set_interval_num(2000);

	Step_T3D_ME_mt step("step1");
	step.set_model(model);
	//step.set_thread_num(22);
	//step.set_step_time(2.0); // 0.5D
	step.set_thread_num(5);
	step.set_step_time(4.0e-3);
	step.set_dtime(5.0e-6);
	//step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();
}

#include "QtApp_Posp_T3D_ME_mt_Div.h"
#include "test_model_view_omp.h"

void test_t3d_me_mt_piezofoundation_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t3d_me_mt_piezofoundation.h5");

	// Single frame
	//QtApp_Posp_T3D_ME_mt_Div<PlaneDivisionSet> app(argc, argv, QtApp_Posp_T3D_ME_mt_Div<PlaneDivisionSet>::SingleFrame);
	////app.set_res_file(rf, "penetration", 50, Hdf5Field::s33);
	////app.set_res_file(rf, "penetration", 50, Hdf5Field::max_shear_stress);
	//app.set_res_file(rf, "penetration", 50, Hdf5Field::plastic_mises_strain_2d);
	
	// Animation
	QtApp_Posp_T3D_ME_mt_Div<PlaneDivisionSet> app(argc, argv,
		QtApp_Posp_T3D_ME_mt_Div<PlaneDivisionSet>::Animation);
	app.get_div_set().set_by_normal_and_point(0.0, 1.0, 0.0, 0.0, 1.0, 0.0);
	app.set_ani_time(5.0);
	app.set_win_size(1200, 800);
	app.set_view_dir(-90.0f, 10.0f);
	app.set_fog_coef(0.02f);
	app.set_light_dir(-135.0f, 20.0f);
	app.set_light_dist_scale(1.0f);
	app.set_view_dist_scale(0.7f);
	app.set_display_bg_mesh(false);
	//app.set_res_file(rf, "penetration", Hdf5Field::s33);
	//app.set_color_map_fld_range(-20.0, 0.0);
	//app.set_res_file(rf, "penetration", Hdf5Field::max_shear_stress);
	//app.set_color_map_fld_range(0.0, 5.0);
	app.set_res_file(rf, "penetration", Hdf5Field::plastic_mises_strain_2d);
	app.set_color_map_fld_range(0.0, 0.1);
	app.set_color_map_geometry(1.2f, 0.4f, 0.45f);
	//app.set_png_name("t3d_me_mt_piezofoundation");
	//app.set_gif_name("t3d_me_mt_piezofoundation");
	app.start();
}
