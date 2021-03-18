#include "TestsParallel_pcp.h"

#include "ParticleGenerator3D.hpp"
#include "LinearElasticity.h"
#include "Model_T3D_CHM_mt.h"
#include "ModelData_T3D_CHM_mt.h"
#include "Step_T3D_CHM_TBB.h"
#include "TimeHistory_T3D_CHM_TBB_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "QtApp_Prep_T3D_CHM_mt.h"
#include "test_parallel_utils.h"
#include "test_simulations_omp.h"

void test_t3d_chm_tbb_1d_consolidation(int argc, char **argv)
{
	TetrahedronMesh teh_mesh;
	teh_mesh.load_mesh_from_hdf5("../../Asset/brick_mesh_1.00_2x2x10.h5");
	teh_mesh.init_search_grid(0.05, 0.05, 0.05);

	Model_T3D_CHM_mt model;
	model.init_mesh(teh_mesh);
	model.init_search_grid(teh_mesh);

	ParticleGenerator3D<TetrahedronMesh> pcl_generator;
	pcl_generator.generate_pcls_grid(Cube(0.0, 0.2, 0.0, 0.2, 0.0, 1.0), 0.025, 0.025, 0.025);
	model.init_pcls(pcl_generator, 0.3, 20.0, 10.0, 10000.0, 1.0e-4, 1.0);
	const size_t pcl_num = model.get_pcl_num();
	MatModel::MaterialModel** mms = model.get_mat_models();
	MatModel::LinearElasticity* les = model.add_LinearElasticity(pcl_num);
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		MatModel::LinearElasticity& le = les[pcl_id];
		le.set_param(1000.0, 0.0);
		mms[pcl_id] = &le;
	}

	IndexArray tbc_pcl_array(100);
	find_3d_pcls(model, tbc_pcl_array, Cube(0.0, 0.2, 0.0, 0.2, 1.0 - 0.013, 1.0));
	MemoryUtils::ItemArray<double> tzs_mem(tbc_pcl_array.get_num());
	double tz_mag = 0.025 * 0.025 * -10.0;
	for (size_t t_id = 0; t_id < tbc_pcl_array.get_num(); ++t_id)
		tzs_mem.add(tz_mag);
	model.init_tzs(tbc_pcl_array.get_num(), tbc_pcl_array.get_mem(), tzs_mem.get_mem());

	IndexArray vx_bc_pt_array(100);
	find_3d_nodes_on_x_plane(model, vx_bc_pt_array, 0.0);
	find_3d_nodes_on_x_plane(model, vx_bc_pt_array, 0.2, false);
	model.init_fixed_vx_s_bc(vx_bc_pt_array.get_num(), vx_bc_pt_array.get_mem());
	model.init_fixed_vx_f_bc(vx_bc_pt_array.get_num(), vx_bc_pt_array.get_mem());

	IndexArray vy_bc_pt_array(100);
	find_3d_nodes_on_y_plane(model, vy_bc_pt_array, 0.0);
	find_3d_nodes_on_y_plane(model, vy_bc_pt_array, 0.2, false);
	model.init_fixed_vy_s_bc(vy_bc_pt_array.get_num(), vy_bc_pt_array.get_mem());
	model.init_fixed_vy_f_bc(vy_bc_pt_array.get_num(), vy_bc_pt_array.get_mem());

	IndexArray vz_bc_pt_array(100);
	find_3d_nodes_on_z_plane(model, vz_bc_pt_array, 0.0);
	model.init_fixed_vz_s_bc(vz_bc_pt_array.get_num(), vz_bc_pt_array.get_mem());
	model.init_fixed_vz_f_bc(vz_bc_pt_array.get_num(), vz_bc_pt_array.get_mem());

	//QtApp_Prep_T3D_CHM_mt md_disp(argc, argv);
	//md_disp.set_win_size(900, 900);
	//md_disp.set_view_dir(30.0, -30.0);
	//md_disp.set_light_dir(90.0, -15.0);
	//md_disp.set_model(model);
	//md_disp.set_pts_from_node_id(vx_bc_pt_array.get_mem(), vx_bc_pt_array.get_num(), 0.01);
	////md_disp.set_pts_from_node_id(vy_bc_pt_array.get_mem(), vy_bc_pt_array.get_num(), 0.01);
	////md_disp.set_pts_from_node_id(vz_bc_pt_array.get_mem(), vz_bc_pt_array.get_num(), 0.01);
	////md_disp.set_pts_from_pcl_id(tbc_pcl_array.get_mem(), tbc_pcl_array.get_num(), 0.012);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_chm_tbb_1d_consolidation.h5");

	ModelData_T3D_CHM_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T3D_CHM_TBB_complete out1("consolidation");
	out1.set_interval_num(50);
	out1.set_output_init_state();
	out1.set_res_file(res_file_hdf5);
	TimeHistory_ConsoleProgressBar out_cpb;
	out_cpb.set_interval_num(1000);

	Step_T3D_CHM_TBB step("step1");
	step.set_model(model);
	step.set_step_time(10.0);
	step.set_dtime(1.0e-5);
	step.set_thread_num(5);
	step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();
}

#include "QtApp_Posp_T3D_CHM_mt.h"
#include "test_model_view_omp.h"

void test_t3d_chm_tbb_1d_consolidation_result(int argc, char **argv)
{
	ResultFile_hdf5 rf;
	rf.open("t3d_chm_tbb_1d_consolidation.h5");

	//QtApp_Posp_T3D_CHM_mt app(argc, argv, QtApp_Posp_T3D_ME_s::SingleFrame);
	//app.set_res_file(rf, "compression", 2, Hdf5Field::p);
	QtApp_Posp_T3D_CHM_mt app(argc, argv, QtApp_Posp_T3D_CHM_mt::Animation);
	app.set_ani_time(5.0);
	app.set_win_size(900, 900);
	app.set_view_dir(30.0f, 30.0f);
	app.set_light_dir(90.0f, 30.0f);
	app.set_res_file(rf, "consolidation", Hdf5Field::p);
	app.set_color_map_fld_range(0.0, 10.0);
	app.set_color_map_geometry(0.7f, 0.45f, 0.5f);
	//app.set_png_name("t3d_chm_tbb_1d_consolidation");
	app.set_gif_name("t3d_chm_tbb_1d_consolidation");
	app.start();
}
