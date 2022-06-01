#include "TestsParallel_pcp.h"

#include "ParticleGenerator3D.hpp"
#include "LinearElasticity.h"
#include "Model_T3D_ME_mt.h"
#include "Step_T3D_ME_mt.h"
#include "ModelData_T3D_ME_mt.h"
#include "TimeHistory_T3D_ME_mt_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "QtApp_Prep_T3D_ME_mt.h"
#include "test_parallel_utils.h"
#include "test_simulations_omp.h"

void test_t3d_me_mt_block_sliding(int argc, char** argv)
{
	TetrahedronMesh teh_mesh;
	teh_mesh.load_mesh_from_hdf5("../../Asset/brick_mesh_1.00_20x5x1.h5");
	teh_mesh.init_search_grid(0.5, 0.5, 0.5);

	Model_T3D_ME_mt model;
	model.init_mesh(teh_mesh);
	model.init_search_grid(teh_mesh);

	ParticleGenerator3D<TetrahedronMesh> pcl_generator;
	pcl_generator.generate_pcls_grid(Cube(0.0, 20.0, 0.0, 5.0, 0.75, 1.0), 0.25, 0.25, 0.25);
	model.init_pcls(pcl_generator, 10.0);
	size_t pcl_num = model.get_pcl_num();
	MatModel::MaterialModel** mms = model.get_mat_models();
	MatModel::LinearElasticity* les = model.add_LinearElasticity(pcl_num);
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		les->set_param(1000.0, 0.0);
		mms[pcl_id] = les;
		les = model.following_LinearElasticity(les);
	}

	model.init_rigid_cube(2.5, 2.5, 1.5, 3.0, 3.0, 1.0, 1.0);
	model.set_rigid_cube_force(3.0, 0.0, -3.0);
	model.set_contact_param(200.0, 200.0, 0.2, 0.1);
	//model.set_rough_contact_between_pcl_and_rect();
	//model.set_frictional_contact_between_pcl_and_rect();
	model.set_sticky_contact_between_pcl_and_rect();

	IndexArray all_n_ids(model.get_node_num());
	for (size_t n_id = 0; n_id < model.get_node_num(); ++n_id)
		all_n_ids.add(n_id);
	model.init_fixed_vx_bc(all_n_ids.get_num(), all_n_ids.get_mem());
	model.init_fixed_vy_bc(all_n_ids.get_num(), all_n_ids.get_mem());
	model.init_fixed_vz_bc(all_n_ids.get_num(), all_n_ids.get_mem());

	//QtApp_Prep_T3D_ME_mt md_disp(argc, argv);
	//md_disp.set_win_size(1200, 950);
	//md_disp.set_view_dir(100.0f, 20.0f);
	//md_disp.set_light_dir(120.0f, 20.0f);
	//md_disp.set_model(model);
	//md_disp.set_pts_from_node_id(all_n_ids.get_mem(), all_n_ids.get_num(), 0.1);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_me_mt_block_sliding.h5");

	ModelData_T3D_ME_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T3D_ME_mt_complete out1("sliding");
	out1.set_res_file(res_file_hdf5);
	out1.set_interval_num(100);
	out1.set_output_init_state();
	out1.set_output_final_state();
	TimeHistory_ConsoleProgressBar out_cpb;

	Step_T3D_ME_mt step("step1");
	step.set_model(model);
	step.set_step_time(5.0); // 10.0
	//step.set_step_time(1.0e-5);
	step.set_dtime(1.0e-5);
	//step.set_thread_num(2);
	step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();
}

#include "QtApp_Posp_T3D_ME_mt.h"
#include "test_model_view_omp.h"

void test_t3d_me_mt_block_sliding_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t3d_me_mt_block_sliding.h5");

	QtApp_Posp_T3D_ME_mt app(argc, argv, QtApp_Posp_T3D_ME_mt::Animation);
	app.set_res_file(rf, "sliding", Hdf5Field::s33);
	app.set_ani_time(5.0);
	app.set_win_size(1200, 950);
	app.set_view_dir(90.0f, 20.0f);
	app.set_light_dir(90.0f, 20.0f);
	app.set_color_map_fld_range(-10.0, 0.0);
	//app.set_color_map_geometry(0.7f, 0.45f, 0.5f);
	//app.set_png_name("t3d_me_mt_block_sliding");
	app.set_gif_name("t3d_me_mt_block_sliding");
	app.start();
}
