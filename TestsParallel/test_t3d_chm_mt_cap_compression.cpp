#include "TestsParallel_pcp.h"

#include "ParticleGenerator3D.hpp"
#include "LinearElasticity.h"
#include "Model_T3D_CHM_mt.h"
#include "ModelData_T3D_CHM_mt.h"
#include "Step_T3D_CHM_ud_mt_subiter.h"
#include "TimeHistory_T3D_CHM_ud_mt_subiter_complete.h"
#include "TimeHistory_T3D_CHM_ud_mt_subiter_ratio.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "QtApp_Prep_T3D_CHM_mt_Div.h"
#include "test_parallel_utils.h"
#include "test_simulations_omp.h"

void test_t3d_chm_mt_cap_compression(int argc, char **argv)
{
	TetrahedronMesh teh_mesh;
	teh_mesh.load_mesh_from_hdf5("../../Asset/brick_mesh_1.00_2x2x10.h5");
	teh_mesh.init_search_grid(0.05, 0.05, 0.05);

	ParticleGenerator3D<TetrahedronMesh> pcl_generator;
	pcl_generator.generate_pcls_grid(Cube(0.0, 0.2, 0.0, 0.2, 0.0, 1.0), 0.025, 0.025, 0.025);
	pcl_generator.adjust_pcl_size_to_fit_elems(teh_mesh);

	constexpr double e0 = 0.55;
	constexpr double den_grain = 2670.0;
	Model_T3D_CHM_mt model;
	model.init_mesh(teh_mesh);
	model.init_search_grid(teh_mesh);
	//model.init_pcls(pcl_generator, e0 / (1.0 + e0), 2.0, 1.0, 2.0e4, 5.0e-9, 1.0); // elastic
	//model.init_pcls(pcl_generator, e0 / (1.0 + e0), den_grain, 1000.0, 3.6e8, 5.0e-9, 1.0); // hypo, undrained
	model.init_pcls(pcl_generator, e0 / (1.0 + e0), den_grain, 1000.0, 0.0, 5.0e-9, 1.0); // hypo, drained
	const size_t pcl_num = model.get_pcl_num();
	// Linear elasticity
	//MatModel::LinearElasticity* les = model.add_LinearElasticity(pcl_num);
	//// Linear elasticity
	//for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	//{
	//	MatModel::LinearElasticity& le = les[pcl_id];
	//	le.set_param(1000.0, 0.0);
	//	model.add_mat_model(pcl_id, le, sizeof(MatModel::LinearElasticity));
	//}
	// Stb hypoplasticity
	const double ini_stress[6] = { -100.0e3, -100.0e3, -100.0e3, 0.0, 0.0, 0.0 };
	MatModel::SandHypoplasticityStbWrapper* shps = model.add_SandHypoplasticityStbWrapper(pcl_num);
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		MatModel::SandHypoplasticityStbWrapper& shp = shps[pcl_id];
		shp.set_param(
			ini_stress, 0.55,
			30.0, 1354.0e6, 0.34,
			0.18, 1.27,
			0.49, 0.76, 0.86,
			2.0, 350.0, 300.0, // 1.0, 180.0, 200.0 // higher N, lower peak
			200.0, 0.2);
	 	model.add_mat_model(pcl_id, shp, sizeof(MatModel::SandHypoplasticityStbWrapper));
	}

	model.init_rigid_cylinder(0.1, 0.1, 1.025, 0.05, 0.2);
	model.set_rigid_cylinder_velocity(0.0, 0.0, -0.05);
	//const double Kct = 20.0 / (0.025 * 0.025); // elastic 
	const double Kct = 1.0e5 / (0.025 * 0.025); // hypo 
	model.set_contact_param(Kct, Kct, 0.1, 0.2, Kct / 10.0, Kct / 10.0);

	IndexArray vx_bc_pt_array(100);
	find_3d_nodes_on_x_plane(model, vx_bc_pt_array, 0.0);
	//find_3d_nodes_on_x_plane(model, vx_bc_pt_array, 0.2, false);
	model.init_fixed_vx_s_bc(vx_bc_pt_array.get_num(), vx_bc_pt_array.get_mem());
	model.init_fixed_vx_f_bc(vx_bc_pt_array.get_num(), vx_bc_pt_array.get_mem());

	IndexArray vy_bc_pt_array(100);
	find_3d_nodes_on_y_plane(model, vy_bc_pt_array, 0.0);
	//find_3d_nodes_on_y_plane(model, vy_bc_pt_array, 0.2, false);
	model.init_fixed_vy_s_bc(vy_bc_pt_array.get_num(), vy_bc_pt_array.get_mem());
	model.init_fixed_vy_f_bc(vy_bc_pt_array.get_num(), vy_bc_pt_array.get_mem());

	IndexArray vz_bc_pt_array(100);
	find_3d_nodes_on_z_plane(model, vz_bc_pt_array, 0.0);
	model.init_fixed_vz_s_bc(vz_bc_pt_array.get_num(), vz_bc_pt_array.get_mem());
	model.init_fixed_vz_f_bc(vz_bc_pt_array.get_num(), vz_bc_pt_array.get_mem());

	//QtApp_Prep_T3D_CHM_mt_Div<> md_disp(argc, argv);
	////QtApp_Prep_T3D_ME_mt_Div<BoxDivisionSet> md_disp(argc, argv);
	////md_disp.get_div_set().set_param(0.0, 0.1, 0.0, 0.1, 0.0, 0.18);
	////QtApp_Prep_T3D_ME_mt_Div<PlaneDivisionSet> md_disp(argc, argv);
	////md_disp.get_div_set().set_param(0.0, 0.0, -1.0, 0.45);
	//md_disp.set_win_size(1200, 950);
	//md_disp.set_view_dir(200.0f, 30.0f);
	//md_disp.set_light_dir(200.0f, 30.0f);
	////md_disp.set_view_dist_scale(0.5);
	//md_disp.set_model(model);
	////md_disp.set_pts_from_node_id(vx_bc_pt_array.get_mem(), vx_bc_pt_array.get_num(), 0.01);
	////md_disp.set_pts_from_node_id(vy_bc_pt_array.get_mem(), vy_bc_pt_array.get_num(), 0.01);
	////md_disp.set_pts_from_node_id(vz_bc_pt_array.get_mem(), vz_bc_pt_array.get_num(), 0.01);
	//size_t disp_p_id = 200;
	//md_disp.set_pts_from_pcl_id(&disp_p_id, 1, 0.01);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_chm_mt_cap_compression.h5");

	ModelData_T3D_CHM_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_ConsoleProgressBar out_cpb;
	out_cpb.set_interval_num(500);
	TimeHistory_T3D_CHM_ud_mt_subiter_complete out1("compression");
	out1.set_res_file(res_file_hdf5);
	out1.set_interval_num(100);
	out1.set_output_init_state();
	TimeHistory_T3D_CHM_ud_mt_subiter_ratio out2("compression");
	out2.set_res_file("t3d_chm_mt_cap_compression_ratio.csv");
	out2.set_interval_num(200);
	out2.set_output_init_state();

	Step_T3D_CHM_ud_mt_subiter step("step1");
	step.set_model(model);
	step.set_thread_num(5);
	step.set_step_time(1.0); // 1.0
	//step.set_mass_factor(0.1); // debug
	//step.set_step_time(1.0e-4); // debug
	step.set_max_subiter_num(0); // debug
	step.set_dtime(5.0e-6);
	step.set_pdt(5.0e-6);
	step.add_time_history(out_cpb);
	step.add_time_history(out1);
	step.add_time_history(out2);
	step.solve();
}

#include "QtApp_Posp_T3D_CHM_mt.h"
#include "test_model_view_omp.h"

void test_t3d_chm_mt_cap_compression_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t3d_chm_mt_cap_compression.h5");

	QtApp_Posp_T3D_CHM_mt app(argc, argv, QtApp_Posp_T3D_CHM_mt::Animation);
	app.set_win_size(1200, 950);
	app.set_ani_time(5.0);
	app.set_view_dir(30.0f, 0.0f);
	app.set_light_dir(30.0f, 20.0f);
	//app.set_view_dist_scale(1.1);
	app.set_color_map_geometry(0.85f, 0.45f, 0.5f);
	//app.set_png_name("t3d_chm_mt_cap_compression");
	app.set_gif_name("t3d_chm_mt_cap_compression");
	// s33
	//app.set_res_file(rf, "compression", Hdf5Field::s33);
	////app.set_color_map_fld_range(-50.0, 0.0); // elastic
	//app.set_color_map_fld_range(-500.0e3, 0.0); // hypo
	// 	p
	app.set_res_file(rf, "compression", Hdf5Field::p);
	app.set_color_map_fld_range(-50.0e3, 50.0e3); // hypo
	// shear stress
	//app.set_res_file(rf, "compression", Hdf5Field::max_shear_stress);
	//app.set_color_map_fld_range(0.0, 30.0);
	//
	app.start();
}
