#include "TestsParallel_pcp.h"

#include "test_parallel_utils.h"
#include "TriangleMesh.h"
#include "ParticleGenerator2D.hpp"
#include "Model_T2D_ME_mt.h"
#include "Model_T2D_ME_mt_hdf5_utilities.h"
#include "ModelData_T2D_ME_mt.h"
#include "Step_T2D_ME_mt_Geo.h"
#include "TimeHistory_T2D_ME_mt_Geo_complete.h"
#include "TimeHistory_T2D_ME_mt_Geo_ratio.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "QtApp_Prep_T2D_ME_mt.h"
#include "test_simulations_omp.h"

void test_t2d_me_mt_geostatic(int argc, char** argv)
{
	TriangleMesh tri_mesh;
	tri_mesh.load_mesh_from_hdf5("../../Asset/rect_mesh.h5");
	
	Model_T2D_ME_mt model;
	model.init_mesh(tri_mesh);
	model.init_search_grid(tri_mesh, 0.05, 0.05);
	
	ParticleGenerator2D<TriangleMesh> pcl_generator;
	pcl_generator.generate_pcls_in_grid_layout(Rect(0.0, 0.2, 0.0, 1.0), 0.02, 0.02);
	model.init_pcls(pcl_generator, 10.0);
	MatModel::MaterialModel** mms = model.get_mat_models();
	MatModel::LinearElasticity *les = model.add_LinearElasticity(model.get_pcl_num());
	for (uint32_t p_id = 0; p_id < model.get_pcl_num(); ++p_id)
	{
		les->set_param(1000.0, 0.0);
		mms[p_id] = les;
		les = model.following_LinearElasticity(les);
	}

	IndexArray tbc_pt_array(50);
	find_2d_pcls(model, tbc_pt_array, Rect(0.0, 0.2, 0.987, 1.0));
	MemoryUtils::ItemArray<double> tys_mem(tbc_pt_array.get_num());
	double ty_mag = 0.02 * -1.0;
	for (size_t t_id = 0; t_id < tbc_pt_array.get_num(); ++t_id)
		tys_mem.add(ty_mag);
	model.init_tys(tbc_pt_array.get_num(), tbc_pt_array.get_mem(), tys_mem.get_mem());

	// vx bc
	IndexArray vx_bc_pt_array(50);
	find_2d_nodes_on_x_line(model, vx_bc_pt_array, 0.0);
	find_2d_nodes_on_x_line(model, vx_bc_pt_array, 0.2, false);
	model.init_fixed_vx_bc(vx_bc_pt_array.get_num(), vx_bc_pt_array.get_mem());

	// vy bc
	IndexArray vy_bc_pt_array(50);
	find_2d_nodes_on_y_line(model, vy_bc_pt_array, 0.0);
	model.init_fixed_vy_bc(vy_bc_pt_array.get_num(), vy_bc_pt_array.get_mem());

	//QtApp_Prep_T2D_ME_mt md_disp(argc, argv);
	//md_disp.set_win_size(900, 900);
	//md_disp.set_model(model);
	////md_disp.set_pts_from_node_id(vx_bc_pt_array.get_mem(), vx_bc_pt_array.get_num(), 0.01);
	////md_disp.set_pts_from_node_id(vy_bc_pt_array.get_mem(), vy_bc_pt_array.get_num(), 0.01);
	//md_disp.set_pts_from_pcl_id(tbc_pt_array.get_mem(), tbc_pt_array.get_num(), 0.01);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_me_mt_geostatic.h5");

	ModelData_T2D_ME_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T2D_ME_mt_Geo_complete out("geostatic");
	out.set_interval_num(100);
	out.set_output_init_state();
	out.set_output_final_state();
	out.set_res_file(res_file_hdf5);
	TimeHistory_T2D_ME_mt_Geo_ratio out2("ratio");
	out2.set_interval_num(200);
	out2.set_res_file("t2d_me_mt_geostatic.csv");
	TimeHistory_ConsoleProgressBar out_pb;
	out_pb.set_interval_num(500);

	Step_T2D_ME_mt_Geo step("step1");
	step.set_model(model);
	step.set_step_time(1.0);
	step.set_dtime(1.0e-5);
	step.set_thread_num(3);
	step.add_time_history(out);
	step.add_time_history(out2);
	step.add_time_history(out_pb);
	step.solve();
}

#include "QtApp_Posp_T2D_ME_mt.h"
#include "test_model_view_omp.h"

void test_t2d_me_mt_geostatic_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t2d_me_mt_geostatic.h5");

	QtApp_Posp_T2D_ME_mt app(argc, argv, QtApp_Posp_T2D_ME_mt::Animation);
	app.set_ani_time(5.0);
	app.set_res_file(rf, "geostatic", Hdf5Field::s22);
	app.set_win_size(900, 900);
	app.set_color_map_fld_range(-2.0, 0.0);
	//app.set_color_map_geometry(1.0f, 0.45f, 0.5f);
	//app.set_png_name("t2d_me_mt_1d_compression");
	//app.set_gif_name("t2d_me_mt_1d_compression");
	app.start();
}
