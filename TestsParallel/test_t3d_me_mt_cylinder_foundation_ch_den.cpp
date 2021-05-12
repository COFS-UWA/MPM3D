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

void test_t3d_me_mt_cylinder_foundation_create_model_ch_den(int argc, char** argv)
{
	TetrahedronMesh teh_mesh;
	teh_mesh.load_mesh_from_hdf5("../../Asset/cube_block_ch_den.h5");
	teh_mesh.init_search_grid(0.05, 0.05, 0.05);
	std::cout << "node_num: " << teh_mesh.get_node_num() << "\n";
	std::cout << "elem_num: " << teh_mesh.get_elem_num() << "\n";

	constexpr double depth = 4.0;
	constexpr double coarse_depth = 3.0;
	constexpr double width = 7.0;
	constexpr double coarse_width = 5.0;
	ParticleGenerator3D<TetrahedronMesh> pcl_generator;
	// dense area
	pcl_generator.generate_pcls_grid(
		Cube(-0.5*coarse_width, 0.5*coarse_width,
			 -0.5*coarse_width, 0.5*coarse_width,
			 -coarse_depth, 0.0),
		0.024, 0.024, 0.024);
	// peripheral coarse area
	// 4 edges
	pcl_generator.generate_pcls_grid(
		Cube(-0.5*width, -0.5*coarse_width,
			 -0.5*coarse_width, 0.5*coarse_width,
			-coarse_depth, 0.0),
			 0.05, 0.05, 0.05);
	pcl_generator.generate_pcls_grid(
		Cube(0.5*coarse_width, 0.5*width,
			-0.5*coarse_width, 0.5*coarse_width,
			-coarse_depth, 0.0),
		0.05, 0.05, 0.05);
	pcl_generator.generate_pcls_grid(
		Cube(-0.5*coarse_width, 0.5*coarse_width,
			 -0.5*width, -0.5*coarse_width,
			 -coarse_depth, 0.0),
		0.05, 0.05, 0.05);
	pcl_generator.generate_pcls_grid(
		Cube(-0.5*coarse_width, 0.5*coarse_width,
			  0.5*coarse_width, 0.5*width,
			 -coarse_depth, 0.0),
		0.05, 0.05, 0.05);
	// 4 corners
	pcl_generator.generate_pcls_grid(
		Cube(-0.5 * width, -0.5 * coarse_width,
			 -0.5 * width, -0.5 * coarse_width,
			 -coarse_depth, 0.0),
			 0.05, 0.05, 0.05);
	pcl_generator.generate_pcls_grid(
		Cube(0.5 * coarse_width, 0.5 * width,
			-0.5 * width, -0.5 * coarse_width,
			-coarse_depth, 0.0),
			 0.05, 0.05, 0.05);
	pcl_generator.generate_pcls_grid(
		Cube(0.5 * coarse_width, 0.5*width,
			 0.5 * coarse_width, 0.5*width,
			 -coarse_depth, 0.0),
			 0.05, 0.05, 0.05);
	pcl_generator.generate_pcls_grid(
		Cube(-0.5 * width, -0.5 * coarse_width,
			  0.5 * coarse_width, 0.5 * width,
			 -coarse_depth, 0.0),
			 0.05, 0.05, 0.05);
	// bottom
	pcl_generator.generate_pcls_grid(
		Cube(-0.5*width, 0.5*width,
			 -0.5*width, 0.5*width,
			 -depth, -coarse_depth),
		0.05, 0.05, 0.05);
	pcl_generator.adjust_pcl_size_to_fit_elems(teh_mesh);
	std::cout << "pcl_num: " << pcl_generator.get_num() << "\n";

	Model_T3D_ME_mt model;
	model.init_mesh(teh_mesh);
	model.init_search_grid(teh_mesh);
	teh_mesh.clear();
	model.init_pcls(pcl_generator, 10.0);
	pcl_generator.clear();

	const size_t pcl_num = model.get_pcl_num();
	MatModel::MaterialModel** mms = model.get_mat_models();
	//MatModel::LinearElasticity* les = model.add_LinearElasticity(pcl_num);
	//for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	//{
	//	MatModel::LinearElasticity& le = les[pcl_id];
	//	le.set_param(1000.0, 0.0);
	//	mms[pcl_id] = &le;
	//}
	MatModel::Tresca *tcs = model.add_Tresca(pcl_num);
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		MatModel::Tresca &tc = tcs[pcl_id];
		tc.set_param(4000.0, 0.3, 5.0);
		mms[pcl_id] = &tc;
	}

	model.init_rigid_cylinder(0.0, 0.0, 0.1, 0.2, 0.5);
	model.set_rigid_cylinder_velocity(0.0, 0.0, -0.1);
	model.set_contact_param(20000.0, 20000.0, 0.1, 5.0);

	IndexArray vx_bc_pt_array(100);
	find_3d_nodes_on_x_plane(model, vx_bc_pt_array, -0.5*width);
	find_3d_nodes_on_x_plane(model, vx_bc_pt_array, 0.5*width, false);
	model.init_fixed_vx_bc(vx_bc_pt_array.get_num(), vx_bc_pt_array.get_mem());

	IndexArray vy_bc_pt_array(100);
	find_3d_nodes_on_y_plane(model, vy_bc_pt_array, -0.5*width);
	find_3d_nodes_on_y_plane(model, vy_bc_pt_array, 0.5*width, false);
	model.init_fixed_vy_bc(vy_bc_pt_array.get_num(), vy_bc_pt_array.get_mem());

	IndexArray vz_bc_pt_array(100);
	find_3d_nodes_on_z_plane(model, vz_bc_pt_array, -depth);
	model.init_fixed_vz_bc(vz_bc_pt_array.get_num(), vz_bc_pt_array.get_mem());

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_me_mt_cylinder_foundation_model_ch_den.h5");

	ModelData_T3D_ME_mt md;
	md.output_model(model, res_file_hdf5);
}

#include "QtApp_Prep_T3D_ME_mt.h"

void test_t3d_me_mt_cylinder_foundation_ch_den(int argc, char** argv)
{
	Model_T3D_ME_mt model;
	Model_T3D_ME_mt_hdf5_utilities::load_me_mt_model_from_hdf5_file(model,
		"t3d_me_mt_cylinder_foundation_model_ch_den.h5");

	model.set_rigid_cylinder_velocity(0.0, 0.0, -0.02);

	////QtApp_Prep_T3D_ME_mt md_disp(argc, argv);
	//QtApp_Prep_T3D_ME_mt_Div<PlaneDivisionSet> md_disp(argc, argv);
	//md_disp.get_div_set().set_by_normal_and_point(1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	//md_disp.set_model(model);
	//md_disp.set_win_size(1200, 950);
	//md_disp.set_view_dir(0.0f, 10.0f);
	//md_disp.set_light_dir(10.0f, 10.0f);
	//md_disp.set_view_dist_scale(0.8);
	//md_disp.set_display_bg_mesh(false);
	//md_disp.set_pts_from_vx_bc(0.02);
	////md_disp.set_pts_from_vy_bc(0.02);
	////md_disp.set_pts_from_vz_bc(0.02);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_me_mt_cylinder_foundation_ch_den.h5");

	ModelData_T3D_ME_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T3D_ME_mt_complete out1("penetration");
	out1.set_interval_num(50);
	out1.set_output_init_state();
	out1.set_res_file(res_file_hdf5);
	TimeHistory_ConsoleProgressBar out_cpb;
	out_cpb.set_interval_num(5000);

	Step_T3D_ME_mt step("step1");
	step.set_model(model);
	step.set_thread_num(30);
	//step.set_step_time(0.35); // 0.035
	step.set_step_time(1.5); // 1.75
	//step.set_step_time(1.0e-5);
	step.set_dtime(5.0e-6);
	step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();
}

void test_t3d_me_mt_cylinder_foundation_restart_ch_den(int argc, char** argv)
{
	Model_T3D_ME_mt model;
	Step_T3D_ME_mt step("step2");
	Model_T3D_ME_mt_hdf5_utilities::load_me_mt_model_from_hdf5_file(
		model,
		step,
		"t3d_me_mt_cylinder_foundation_ch_den.h5",
		"penetration",
		50
		);

	//QtApp_Prep_T3D_ME_mt_Div<TwoPlaneDivisionSet> md_disp(argc, argv);
	//auto &div_set = md_disp.get_div_set();
	//div_set.seta().set_by_normal_and_point(0.0, 1.0, 0.0, 3.5, 3.5, 0.0);
	//div_set.setb().set_by_normal_and_point(1.0, 0.0, 0.0, 3.5, 3.5, 0.0);
	//md_disp.set_model(model);
	//md_disp.set_win_size(1200, 950);
	//md_disp.set_view_dir(30.0f, 30.0f);
	//md_disp.set_light_dir(0.0f, 30.0f);
	//md_disp.set_display_bg_mesh(false);
	////md_disp.set_pts_from_node_id(vx_bc_pt_array.get_mem(), vx_bc_pt_array.get_num(), 0.01);
	////md_disp.set_pts_from_node_id(vy_bc_pt_array.get_mem(), vy_bc_pt_array.get_num(), 0.01);
	////md_disp.set_pts_from_node_id(vz_bc_pt_array.get_mem(), vz_bc_pt_array.get_num(), 0.01);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_me_mt_cylinder_foundation_ch_den2.h5");

	ModelData_T3D_ME_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T3D_ME_mt_complete out1("penetration");
	out1.set_res_file(res_file_hdf5);
	out1.set_output_init_state();
	out1.set_interval_num(50);
	TimeHistory_ConsoleProgressBar out_cpb;

	step.set_thread_num(6);
	step.set_step_time(0.00034);
	step.set_dtime(1.0e-5);
	step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();
}

void test_t3d_me_mt_cylinder_foundation_restart_ch_den2(int argc, char** argv)
{
	Model_T3D_ME_mt model;
	Step_T3D_ME_mt step("step3");
	Model_T3D_ME_mt_hdf5_utilities::load_me_mt_model_from_hdf5_file(
		model,
		step,
		"t3d_me_mt_cylinder_foundation_ch_den2.h5",
		"penetration",
		34
		);

	//QtApp_Prep_T3D_ME_mt_Div<TwoPlaneDivisionSet> md_disp(argc, argv);
	//auto &div_set = md_disp.get_div_set();
	//div_set.seta().set_by_normal_and_point(0.0, 1.0, 0.0, 3.5, 3.5, 0.0);
	//div_set.setb().set_by_normal_and_point(1.0, 0.0, 0.0, 3.5, 3.5, 0.0);
	//md_disp.set_model(model);
	//md_disp.set_win_size(1200, 950);
	//md_disp.set_view_dir(30.0f, 30.0f);
	//md_disp.set_light_dir(0.0f, 30.0f);
	//md_disp.set_display_bg_mesh(false);
	////md_disp.set_pts_from_node_id(vx_bc_pt_array.get_mem(), vx_bc_pt_array.get_num(), 0.01);
	////md_disp.set_pts_from_node_id(vy_bc_pt_array.get_mem(), vy_bc_pt_array.get_num(), 0.01);
	////md_disp.set_pts_from_node_id(vz_bc_pt_array.get_mem(), vz_bc_pt_array.get_num(), 0.01);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_me_mt_cylinder_foundation_ch_den3.h5");

	ModelData_T3D_ME_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T3D_ME_mt_complete out1("penetration");
	out1.set_res_file(res_file_hdf5);
	out1.set_output_init_state();
	out1.set_interval_num(50);
	TimeHistory_ConsoleProgressBar out_cpb;

	//step.set_thread_num(6);
	step.set_step_time(0.00005);
	step.set_dtime(1.0e-5);
	step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();
}

#include "QtApp_Posp_T3D_ME_mt_Div.h"
#include "test_model_view_omp.h"

void test_t3d_me_mt_cylinder_foundation_result_ch_den2(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t3d_me_mt_cylinder_foundation_ch_den_v0.02.h5");

	//QtApp_Posp_T3D_ME_mt_Div<PlaneDivisionSet> app(argc, argv, QtApp_Posp_T3D_ME_mt_Div<PlaneDivisionSet>::SingleFrame);
	////app.set_res_file(rf, "penetration", 50, Hdf5Field::s33);
	////app.set_res_file(rf, "penetration", 50, Hdf5Field::max_shear_stress);
	//app.set_res_file(rf, "penetration", 50, Hdf5Field::plastic_mises_strain_2d);
	QtApp_Posp_T3D_ME_mt_Div<PlaneDivisionSet> app(argc, argv, QtApp_Posp_T3D_ME_mt_Div<PlaneDivisionSet>::Animation);
	app.get_div_set().set_by_normal_and_point(1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	app.set_win_size(1200, 800);
	app.set_ani_time(5.0);
	app.set_view_dir(0.0f, 10.0f);
	app.set_light_dir(0.0f, 10.0f);
	app.set_view_dist_scale(0.8);
	app.set_display_bg_mesh(false);
	//app.set_res_file(rf, "penetration", Hdf5Field::s33);
	//app.set_color_map_fld_range(-20.0, 0.0);
	//app.set_res_file(rf, "penetration", Hdf5Field::max_shear_stress);
	//app.set_color_map_fld_range(0.0, 5.0);
	app.set_res_file(rf, "penetration", Hdf5Field::plastic_mises_strain_2d);
	app.set_color_map_fld_range(0.0, 0.03);
	app.set_color_map_geometry(1.2f, 0.4f, 0.45f);
	//app.set_png_name("t3d_me_mt_cylinder_foundation_ch_den");
	app.set_gif_name("t3d_me_mt_cylinder_foundation_ch_den");
	app.start();
}
