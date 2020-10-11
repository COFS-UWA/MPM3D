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
	Step_T2D_ME_mt step("step1");

	TriangleMesh tri_mesh;
	tri_mesh.load_mesh_from_hdf5("../../Asset/square_mesh.h5");
	model.init_mesh(tri_mesh);

	ParticleGenerator2D<TriangleMesh> pcl_generator;
	pcl_generator.generate_pcls_in_grid_layout(Rect(0.0, 1.0, 0.0, 1.0), 0.25, 0.25);
	model.init_pcls(pcl_generator, 10.0);
	
	//IndexArray mid_tbc_pt_array, left_right_tbc_pt_array;
	//mid_tbc_pt_array.reserve(100);
	//left_right_tbc_pt_array.reserve(100);
	//for (size_t t_id = 0; t_id < model.get_ty_num(); ++t_id)
	//{
	//	if (model.get_tys()[t_id].t > -500)
	//		mid_tbc_pt_array.add(model.get_tys()[t_id].pcl_id);
	//	else
	//		left_right_tbc_pt_array.add(model.get_tys()[t_id].pcl_id);
	//}
	//
	//IndexArray left_right_bc_pt_array;
	//left_right_bc_pt_array.reserve(model.get_vsx_num());
	//for (size_t v_id = 0; v_id < model.get_vsx_num(); ++v_id)
	//	left_right_bc_pt_array.add(model.get_vsxs()[v_id].node_id);

	//IndexArray bottom_bc_pt_array;
	//bottom_bc_pt_array.reserve(model.get_vsy_num());
	//for (size_t v_id = 0; v_id < model.get_vsy_num(); ++v_id)
	//	bottom_bc_pt_array.add(model.get_vsys()[v_id].node_id);

	//QtApp_Prep_T2D_CHM_s md_disp(argc, argv);
	//md_disp.set_win_size(900, 900);
	//md_disp.set_model(model);
	////md_disp.set_pts_from_node_id(left_right_bc_pt_array.get_mem(), left_right_bc_pt_array.get_num(), 0.05);
	////md_disp.set_pts_from_node_id(bottom_bc_pt_array.get_mem(), bottom_bc_pt_array.get_num(), 0.05);
	////md_disp.set_pts_from_pcl_id(mid_tbc_pt_array.get_mem(), mid_tbc_pt_array.get_num(), 0.01);
	//md_disp.set_pts_from_pcl_id(left_right_tbc_pt_array.get_mem(), left_right_tbc_pt_array.get_num(), 0.015);
	//// all
	//md_disp.set_display_range(-3.6, 3.6, -5.1, 1.1);
	//// left
	////md_disp.set_display_range(-3.8, -2.2, -1.0, 1.0);
	//// middle
	////md_disp.set_display_range(-1.5, 1.5, -0.75, 0.25);
	//// right
	////md_disp.set_display_range(2.2, 3.8, -1.0, 1.0);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_me_mt_test1.h5");

	ModelData_T2D_ME_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T2D_ME_mt_complete out("test");
	out.set_res_file(res_file_hdf5);
	out.set_output_init_state();
	out.set_interval_num(10);
	TimeHistory_ConsoleProgressBar out_pb;

	step.set_model(model);
	step.set_step_time(1.0e-1);
	step.set_dtime(1.0e-1);
	step.add_time_history(out);
	//step.add_time_history(out_pb);
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