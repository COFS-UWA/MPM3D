#include "TestsParallel_pcp.h"

#include "ParticleGenerator3D.hpp"
#include "LinearElasticity.h"
#include "Model_T3D_CHM_up_mt.h"
#include "Step_T3D_CHM_up_TBB.h"
#include "ModelData_T3D_CHM_up_mt.h"
#include "TimeHistory_T3D_CHM_up_TBB_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "QtApp_Prep_T3D_CHM_up_mt.h"
#include "test_parallel_utils.h"
#include "test_simulations_omp.h"

void test_t3d_chm_mt_1d_consolidation_up(int argc, char **argv)
{
	TetrahedronMesh teh_mesh;
	teh_mesh.load_mesh_from_hdf5("../../Asset/brick_mesh_1.20_2x2x12.h5");
	teh_mesh.init_search_grid(0.05, 0.05, 0.05);

	Model_T3D_CHM_up_mt model;
	model.init_mesh(teh_mesh);
	model.init_search_grid(teh_mesh);

	ParticleGenerator3D<TetrahedronMesh> pcl_generator;
	pcl_generator.generate_pcls_grid(Cube(0.0, 0.2, 0.0, 0.2, 0.0, 1.0), 0.025, 0.025, 0.025);
	//pcl_generator.generate_pcls_grid(Cube(0.0, 0.2, 0.0, 0.2, 0.0, 0.95), 0.025, 0.025, 0.025);
	pcl_generator.adjust_pcl_size_to_fit_elems(teh_mesh);
	model.init_pcls(pcl_generator, 0.4, 20.0, 10.0, 40000.0, 1.0e-4, 1.0);
	//model.init_pcls(pcl_generator, 0.4, 20.0, 10.0, 40000.0, 0.0, 1.0);
	const size_t pcl_num = model.get_pcl_num();
	MatModel::MaterialModel** mms = model.get_mat_models();
	MatModel::LinearElasticity* les = model.add_LinearElasticity(pcl_num);
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		mms[pcl_id] = les;
		les->set_param(1000.0, 0.0);
		les = model.following_LinearElasticity(les);
	}

	model.init_t3d_rigid_mesh(1.0, "../../Asset/cylinder_cap.h5",
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3, 0.3, 0.3);
	model.set_t3d_rigid_mesh_velocity(0.0, 0.0, -0.01);
	constexpr double K_cont = 1.04;
	model.set_contact_param(K_cont, K_cont, 0.2, 5.0);
	
	//IndexArray tbc_pcl_array(100);
	//find_3d_pcls(model, tbc_pcl_array, Cube(0.0, 0.2, 0.0, 0.2, 1.0 - 0.013, 1.0));
	////find_3d_pcls(model, tbc_pcl_array, Cube(0.0, 0.2, 0.0, 0.2, 0.95 - 0.013, 1.0));
	//MemoryUtils::ItemArray<double> tzs_mem(tbc_pcl_array.get_num());
	//double tz_mag = 0.025 * 0.025 * -1.0;
	//for (size_t t_id = 0; t_id < tbc_pcl_array.get_num(); ++t_id)
	//	tzs_mem.add(tz_mag);
	//model.init_tzs(tbc_pcl_array.get_num(), tbc_pcl_array.get_mem(), tzs_mem.get_mem());

	//IndexArray drained_bc_pt_array(100);
	//find_3d_nodes_on_z_plane(model, drained_bc_pt_array, 1.0);
	//model.init_drained_bc(drained_bc_pt_array.get_num(), drained_bc_pt_array.get_mem());

	IndexArray vx_bc_pt_array(100);
	find_3d_nodes_on_x_plane(model, vx_bc_pt_array, 0.0);
	find_3d_nodes_on_x_plane(model, vx_bc_pt_array, 0.2, false);
	model.init_fixed_vx_bc(vx_bc_pt_array.get_num(), vx_bc_pt_array.get_mem());

	IndexArray vy_bc_pt_array(100);
	find_3d_nodes_on_y_plane(model, vy_bc_pt_array, 0.0);
	find_3d_nodes_on_y_plane(model, vy_bc_pt_array, 0.2, false);
	model.init_fixed_vy_bc(vy_bc_pt_array.get_num(), vy_bc_pt_array.get_mem());

	IndexArray vz_bc_pt_array(100);
	find_3d_nodes_on_z_plane(model, vz_bc_pt_array, 0.0);
	model.init_fixed_vz_bc(vz_bc_pt_array.get_num(), vz_bc_pt_array.get_mem());

	//QtApp_Prep_T3D_CHM_up_mt md_disp(argc, argv);
	//md_disp.set_win_size(900, 900);
	//md_disp.set_view_dir(30.0, 30.0);
	//md_disp.set_light_dir(30.0, 15.0);
	//md_disp.set_model(model);
	//md_disp.set_view_dist_scale(1.1);
	////md_disp.set_pts_from_drained_bc(0.01);
	////md_disp.set_pts_from_vx_bc(0.01);
	////md_disp.set_pts_from_vy_bc(0.01);
	////md_disp.set_pts_from_vz_bc(0.01);
	////md_disp.set_pts_from_pcl_id(tbc_pcl_array.get_mem(), tbc_pcl_array.get_num(), 0.012);
	////size_t ids[] = { 90, 91, 92, 93, 94, 95, 96, 97, 98 };
	//size_t ids[] = { 82, 84, 86, 88 };
	//md_disp.set_pts_from_node_id(ids, sizeof(ids)/sizeof(ids[0]), 0.01);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_chm_mt_1d_consolidation_up.h5");

	ModelData_T3D_CHM_up_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T3D_CHM_up_TBB_complete out1("consolidation");
	out1.set_interval_num(100);
	out1.set_output_init_state();
	out1.set_res_file(res_file_hdf5);
	TimeHistory_ConsoleProgressBar out_cpb;
	out_cpb.set_interval_num(1000);

	Step_T3D_CHM_up_TBB step("step1");
	step.set_model(model);
	//step.set_step_time(10.0);
	step.set_step_time(3.0e-5);
	step.set_dtime(1.0e-5);
	step.set_thread_num(1);
	step.add_time_history(out1);
	//step.add_time_history(out_cpb);
	step.solve();
}

void test_t3d_chm_mt_1d_consolidation_up_restart(int argc, char** argv)
{
	Model_T3D_CHM_up_mt model;
	Model_T3D_CHM_up_mt_hdf5_utilities::load_model_from_hdf5_file(
		model, "t3d_chm_mt_1d_consolidation.h5");

	//QtApp_Prep_T3D_CHM_mt md_disp(argc, argv);
	//md_disp.set_win_size(900, 900);
	//md_disp.set_view_dir(30.0, -30.0);
	//md_disp.set_light_dir(90.0, -15.0);
	//md_disp.set_model(model);
	////md_disp.set_pts_from_vx_bc_s(0.01);
	////md_disp.set_pts_from_vy_bc_s(0.01);
	////md_disp.set_pts_from_vz_bc_s(0.01);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_chm_mt_1d_consolidation_up2.h5");

	ModelData_T3D_CHM_up_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T3D_CHM_up_TBB_complete out1("consolidation");
	out1.set_interval_num(100);
	out1.set_output_init_state();
	out1.set_output_final_state();
	out1.set_res_file(res_file_hdf5);
	TimeHistory_ConsoleProgressBar out_cpb;
	out_cpb.set_interval_num(2000);

	Step_T3D_CHM_up_TBB step("step1");
	step.set_model(model);
	//step.set_step_time(10.0);
	step.set_step_time(1.0e-5);
	step.set_dtime(5.0e-6);
	step.set_thread_num(5);
	step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();
}

#include "QtApp_Posp_T3D_CHM_up_mt.h"
#include "test_model_view_omp.h"

void test_t3d_chm_mt_1d_consolidation_up_result(int argc, char **argv)
{
	ResultFile_hdf5 rf;
	rf.open("t3d_chm_mt_1d_consolidation_up.h5");

	//QtApp_Posp_T3D_CHM_up_mt app(argc, argv, QtApp_Posp_T3D_CHM_up_mt::SingleFrame);
	//app.set_win_size(900, 900);
	//app.set_view_dir(30.0f, 30.0f);
	//app.set_view_dist_scale(1.1);
	//app.set_light_dir(90.0f, 30.0f);
	//app.set_color_map_geometry(0.7f, 0.45f, 0.5f);
	////app.set_png_name("t3d_chm_mt_1d_consolidation");
	////app.set_gif_name("t3d_chm_mt_1d_consolidation");
	//// p
	//app.set_res_file(rf, "consolidation", 0, Hdf5Field::p);
	//app.set_color_map_fld_range(0.0, 1.0);
	////
	//app.start();

	QtApp_Posp_T3D_CHM_up_mt app(argc, argv, QtApp_Posp_T3D_CHM_up_mt::Animation);
	app.set_ani_time(5.0);
	app.set_win_size(900, 900);
	app.set_view_dir(30.0f, 30.0f);
	app.set_light_dir(90.0f, 30.0f);
	app.set_color_map_geometry(0.7f, 0.45f, 0.5f);
	//app.set_png_name("t3d_chm_mt_1d_consolidation");
	//app.set_gif_name("t3d_chm_mt_1d_consolidation");
	// p
	app.set_res_file(rf, "consolidation", Hdf5Field::p);
	app.set_color_map_fld_range(0.0, 2.0);
	//
	app.start();
}
