#include "TestsParallel_pcp.h"

#include "test_parallel_utils.h"
#include "TriangleMesh.h"
#include "ParticleGenerator2D.hpp"
#include "Model_T2D_ME_mt.h"
#include "Model_T2D_ME_mt_hdf5_utilities.h"
#include "ModelData_T2D_ME_mt.h"
#include "Step_T2D_ME_mt.h"
#include "TimeHistory_T2D_ME_mt_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "QtApp_Prep_T2D_ME_mt.h"
#include "test_simulations_omp.h"

void test_t2d_me_mt_block_collision(int argc, char** argv)
{
	TriangleMesh tri_mesh;
	tri_mesh.load_mesh_from_hdf5("../../Asset/rect_mesh_6by5.h5");

	Model_T2D_ME_mt model;
	model.init_mesh(tri_mesh);
	model.init_search_grid(tri_mesh, 0.05, 0.05);

	ParticleGenerator2D<TriangleMesh> pcl_generator;
	pcl_generator.generate_pcls_in_grid_layout(Rect(1.0, 3.0, 1.0, 5.0), 0.2, 0.2);
	model.init_pcls(pcl_generator, 1.0);
	MatModel::MaterialModel** mms = model.get_mat_models();
	MatModel::LinearElasticity* les = model.add_LinearElasticity(model.get_pcl_num());
	for (uint32_t p_id = 0; p_id < model.get_pcl_num(); ++p_id)
	{
		les[p_id].set_param(10000.0, 0.0);
		mms[p_id] = &les[p_id];
	}
	Model_T2D_ME_mt::Velocity *pcl_v = model.get_ini_pcl_v();
	for (uint32_t p_id = 0; p_id < model.get_pcl_num(); ++p_id)
	{
		Model_T2D_ME_mt::Velocity &pv = pcl_v[p_id];
		pv.vx = 0.5;
		pv.vy = 0.0;
	}

	model.init_rigid_rect(5.0, 3.0, 2.0, 4.0, 1.0);
	model.set_rigid_rect_ini_velocity(-0.5, 0.0, 0.0);
	// only smooth contact is considered here
	model.set_contact_param(2000.0, 2000.0, 0.2, 3.0);

	//QtApp_Prep_T2D_ME_mt md_disp(argc, argv);
	//md_disp.set_win_size(900, 900);
	//md_disp.set_model(model);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_me_mt_block_collision.h5");

	ModelData_T2D_ME_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T2D_ME_mt_complete out("collision");
	out.set_interval_num(100);
	out.set_output_init_state();
	out.set_output_final_state();
	out.set_res_file(res_file_hdf5);
	TimeHistory_ConsoleProgressBar out_pb;

	Step_T2D_ME_mt step("step1");
	step.set_model(model);
	step.set_step_time(3.0);
	step.set_dtime(1.0e-5);
	//step.set_thread_num(4);
	step.add_time_history(out);
	step.add_time_history(out_pb);
	step.solve();
}

#include "QtApp_Posp_T2D_ME_mt.h"
#include "test_model_view_omp.h"

void test_t2d_me_mt_block_collision_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t2d_me_mt_block_collision.h5");

	QtApp_Posp_T2D_ME_mt app(argc, argv, QtApp_Posp_T2D_ME_mt::Animation);
	app.set_win_size(900, 900);
	app.set_ani_time(5.0);
	app.set_res_file(rf, "collision", Hdf5Field::s22);
	app.set_color_map_fld_range(-50.0, 0.0);
	//app.set_color_map_geometry(1.0f, 0.45f, 0.5f);
	//app.set_png_name("t2d_me_mt_block_collision");
	app.set_gif_name("t2d_me_mt_block_collision");
	app.start();
}
