#include "TestsParallel_pcp.h"

#include "ParticleGenerator3D.hpp"
#include "Model_T3D_CHM_mt.h"
#include "Step_T3D_CHM_mt.h"
#include "ModelData_T3D_CHM_mt.h"
#include "TimeHistory_T3D_CHM_mt_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "QtApp_Prep_T3D_CHM_mt.h"
#include "test_parallel_utils.h"
#include "test_simulations_omp.h"

void test_t3d_chm_mt_spudcan_sand_hypo_model(int argc, char** argv)
{
	TetrahedronMesh teh_mesh;
	teh_mesh.load_mesh_from_hdf5("../../Asset/spudcan_soil_quarter.h5");
	teh_mesh.init_search_grid(0.1, 0.1, 0.1);
	std::cout << "node_num: " << teh_mesh.get_node_num() << "\n"
			  << "elem_num: " << teh_mesh.get_elem_num() << "\n";

	constexpr double depth = 20.0;
	constexpr double coarse_depth = 11.5;
	constexpr double width = 17.5;
	constexpr double coarse_width = 11.5;
	constexpr double sml_pcl_size = 0.15;
	constexpr double lgr_pcl_size = 0.3;
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
	//pcl_generator.adjust_pcl_size_to_fit_elems(teh_mesh);
	std::cout << "pcl_num: " << pcl_generator.get_num() << "\n";
	
	Model_T3D_CHM_mt model;
	model.init_mesh(teh_mesh);
	model.init_search_grid(teh_mesh);
	teh_mesh.clear();
	model.init_pcls(pcl_generator, 0.3, 20.0, 10.0, 1.0e4, 1.0e-5, 1.0);
	pcl_generator.clear();

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
	//MatModel::Tresca *tcs = model.add_Tresca(pcl_num);
	//for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	//{
	//	MatModel::Tresca &tc = tcs[pcl_id];
	//	tc.set_param(4000.0, 0.3, 5.0);
	//	mms[pcl_id] = &tc;
	//}
	// Sand hypoplasticity
	const double ini_stress[6] = { -100.0, -100.0, -100.0, 0.0, 0.0, 0.0 };
	MatModel::SandHypoplasticityWrapper *shps = model.add_SandHypoplasticityWrapper(pcl_num);
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		MatModel::SandHypoplasticityWrapper &shp = shps[pcl_id];
		shp.set_param(ini_stress, 0.817,
			33.1, 4.0e9, 0.27,
			0.677, 1.054, 1.212,
			0.14, 2.5);
		mms[pcl_id] = &shp;
	}

	model.init_t3d_rigid_mesh(1.0, "../../Asset/spudcan_model.h5",
		0.0, 0.0, 0.0, 90.0, 0.0, 0.0, 0.3, 0.3, 0.3);
	model.set_t3d_rigid_mesh_velocity(0.0, 0.0, -1.0);
	model.set_contact_param(20000.0/(0.15*0.15), 20000.0/(0.15*0.15),
		0.1, 5.0, 20.0/(0.15*0.15), 20.0/(0.15*0.15));

	IndexArray vx_bc_pt_array(100);
	find_3d_nodes_on_x_plane(model, vx_bc_pt_array, 0.0);
	find_3d_nodes_on_x_plane(model, vx_bc_pt_array, width, false);
	model.init_fixed_vx_s_bc(vx_bc_pt_array.get_num(), vx_bc_pt_array.get_mem());
	model.init_fixed_vx_f_bc(vx_bc_pt_array.get_num(), vx_bc_pt_array.get_mem());
	
	IndexArray vy_bc_pt_array(100);
	find_3d_nodes_on_y_plane(model, vy_bc_pt_array, 0.0);
	find_3d_nodes_on_y_plane(model, vy_bc_pt_array, width, false);
	model.init_fixed_vy_s_bc(vy_bc_pt_array.get_num(), vy_bc_pt_array.get_mem());
	model.init_fixed_vy_f_bc(vy_bc_pt_array.get_num(), vy_bc_pt_array.get_mem());

	IndexArray vz_bc_pt_array(100);
	find_3d_nodes_on_z_plane(model, vz_bc_pt_array, -depth);
	model.init_fixed_vz_s_bc(vz_bc_pt_array.get_num(), vz_bc_pt_array.get_mem());
	model.init_fixed_vz_f_bc(vz_bc_pt_array.get_num(), vz_bc_pt_array.get_mem());

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_chm_mt_spudcan_sand_hypo_model.h5");

	ModelData_T3D_CHM_mt md;
	md.output_model(model, res_file_hdf5);

	QtApp_Prep_T3D_CHM_mt md_disp(argc, argv);
	md_disp.set_model(model);
	md_disp.set_win_size(1200, 950);
	md_disp.set_view_dir(-90.0f, 30.0f);
	md_disp.set_light_dir(-60.0f, 15.0f);
	md_disp.set_display_bg_mesh(false);
	md_disp.set_view_dist_scale(0.7);
	//md_disp.set_pts_from_vx_bc_s(0.1);
	//md_disp.set_pts_from_vy_bc_s(0.1);
	//md_disp.set_pts_from_vz_bc_s(0.1);
	//md_disp.set_pts_from_vx_bc_f(0.1);
	//md_disp.set_pts_from_vy_bc_f(0.1);
	md_disp.set_pts_from_vz_bc_f(0.1);
	md_disp.start();
}

void test_t3d_chm_mt_spudcan_sand_hypo(int argc, char** argv)
{
	Model_T3D_CHM_mt model;
	Model_T3D_CHM_mt_hdf5_utilities::load_model_from_hdf5_file(
		model, "t3d_chm_mt_spudcan_sand_hypo_model.h5");

	//QtApp_Prep_T3D_CHM_mt_Div<PlaneDivisionSet> md_disp(argc, argv);
	//md_disp.get_div_set().set_by_normal_and_point(1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	//md_disp.set_model(model);
	//md_disp.set_win_size(1200, 950);
	//md_disp.set_view_dir(0.0f, 5.0f);
	//md_disp.set_light_dir(10.0f, 5.0f);
	//md_disp.set_display_bg_mesh(false);
	//md_disp.set_view_dist_scale(0.6);
	////md_disp.set_pts_from_vx_bc(0.2);
	////md_disp.set_pts_from_vy_bc(0.2);
	////md_disp.set_pts_from_vz_bc(0.2);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_chm_mt_spudcan_sand_hypo.h5");

	ModelData_T3D_CHM_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T3D_CHM_mt_complete out1("penetration");
	out1.set_interval_num(100);
	out1.set_output_init_state();
	out1.set_output_final_state();
	out1.set_res_file(res_file_hdf5);
	TimeHistory_ConsoleProgressBar out_cpb;
	out_cpb.set_interval_num(2000);

	Step_T3D_CHM_mt step("step1");
	step.set_model(model);
	step.set_thread_num(22);
	step.set_step_time(3.0);
	//step.set_step_time(1.0e-4);
	step.set_dtime(5.0e-6);
	step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();
}

#include "QtApp_Posp_T3D_CHM_mt.h"
#include "QtApp_Posp_T3D_CHM_mt_Div.h"
#include "test_model_view_omp.h"

void test_t3d_chm_mt_spudcan_sand_hypo_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t3d_me_mt_spudcan_sand_hypo.h5");

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
	//app.set_res_file(rf, "penetration", Hdf5Field::s33);
	//app.set_color_map_fld_range(-20.0, 0.0);
	//app.set_res_file(rf, "penetration", Hdf5Field::max_shear_stress);
	//app.set_color_map_fld_range(0.0, 5.0);
	app.set_res_file(rf, "penetration", Hdf5Field::plastic_mises_strain_2d);
	app.set_color_map_fld_range(0.0, 0.1);
	app.set_color_map_geometry(1.2f, 0.4f, 0.45f);
	//app.set_png_name("t3d_me_mt_spudcan_sand_hypo");
	app.set_gif_name("t3d_me_mt_spudcan_sand_hypo");
	app.start();
}
