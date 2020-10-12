#include "TestsParallel_pcp.h"

#include "TriangleMesh.h"
#include "ParticleGenerator2D.hpp"
#include "Model_T2D_ME_mt.h"
#include "Model_T2D_ME_mt_hdf5_utilities.h"
#include "ModelData_T2D_ME_mt.h"
#include "Step_T2D_ME_mt.h"
#include "TimeHistory_T2D_ME_mt_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"
//#include "QtApp_Prep_T2D_CHM_s.h"

#include "test_simulations_omp.h"

void test_t2d_me_mt_test1(int argc, char** argv)
{
	Model_T2D_ME_mt model;

	TriangleMesh tri_mesh;
	tri_mesh.load_mesh_from_hdf5("../../Asset/square_mesh.h5");
	model.init_mesh(tri_mesh);

	model.init_search_grid(tri_mesh, 0.6, 0.6);

	ParticleGenerator2D<TriangleMesh> pcl_generator;
	//pcl_generator.generate_pcls_in_grid_layout(Rect(0.0, 1.0, 0.0, 1.0), 0.25, 0.25);
	//pcl_generator.generate_pcls_in_grid_layout(Rect(0.0, 1.0, 0.5, 1.5), 0.25, 0.25);
	pcl_generator.generate_pcls_in_grid_layout(Rect(0.0, 1.0, 1.0, 2.0), 0.25, 0.25);
	model.init_pcls(pcl_generator, 10.0);
	MatModel::MaterialModel **mms = model.get_mat_models();
	MatModel::LinearElasticity *les = model.add_LinearElasticity(model.get_pcl_num());
	for (uint32_t p_id = 0; p_id < model.get_pcl_num(); ++p_id)
	{
		les[p_id].set_param(1000.0, 0.0);
		mms[p_id] = &les[p_id];
	}

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_me_mt_test1.h5");

	ModelData_T2D_ME_mt md;
	md.output_model(model, res_file_hdf5);

	//TimeHistory_T2D_ME_mt_complete out("test");
	//out.set_res_file(res_file_hdf5);
	//out.set_output_init_state();
	//out.set_interval_num(10);
	//TimeHistory_ConsoleProgressBar out_pb;

	Step_T2D_ME_mt step("step1");
	step.set_model(model);
	step.set_step_time(1.0e-1);
	step.set_dtime(1.0e-1);
	//step.add_time_history(out);
	//step.add_time_history(out_pb);
	//step.set_thread_num(3);
	//step.init_calculation();
	//step.finalize_calculation();
	step.solve();
}

//#include "QtApp_Posp_T2D_CHM_s.h"
#include "test_model_view_omp.h"

void test_t2d_me_mt_test1_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t2d_me_mt_test1.h5");

	//QtApp_Posp_T2D_CHM_s app(argc, argv, QtApp_Posp_T2D_CHM_s::Animation);
	//app.set_win_size(900, 900);
	//app.set_ani_time(5.0);
	//app.set_res_file(rf, "penetration", Hdf5Field::s22);
	//app.set_display_range(-3.6, 3.6, -5.1, 0.6);
	////app.set_color_map_fld_range(-30000.0, -10000.0); // s22
	////app.set_color_map_fld_range(0, 20000.0); // pore pressure
	//app.set_color_map_fld_range(0, 0.4); // mises strain
	//app.set_color_map_geometry(1.0f, 0.45f, 0.5f);
	////app.set_png_name("t2d_chm_s_pipe_conference_restart1");
	////app.set_gif_name("t2d_chm_s_pipe_conference_restart1");
	//app.start();
}