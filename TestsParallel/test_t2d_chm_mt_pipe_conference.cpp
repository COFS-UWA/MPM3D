#include "TestsParallel_pcp.h"

#include "test_parallel_utils.h"
#include "Model_T2D_CHM_mt.h"
#include "Step_T2D_CHM_mt.h"
#include "ModelData_T2D_CHM_mt.h"
#include "TimeHistory_T2D_CHM_mt_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "QtApp_Prep_T2D_CHM_mt.h"
#include "test_simulations_omp.h"

void test_t2d_chm_mt_pipe_conference(int argc, char** argv)
{
	TriangleMesh tri_mesh;
	tri_mesh.load_mesh_from_hdf5("../../Asset/rect_pipe_conference_mesh.h5");
	tri_mesh.init_search_grid(0.05, 0.05);

	ParticleGenerator2D<TriangleMesh> pcl_generator;
	pcl_generator.generate_pcls_in_grid_layout(Rect(-3.5, 3.5, -3.5, 0.0), 0.04, 0.04);
	pcl_generator.generate_pcls_in_grid_layout(Rect(-3.5, 3.5, -5.0, -3.5), 0.04, 0.04);
	pcl_generator.replace_with_pcls_in_grid_layout(Rect(-2.5, 2.5, -3.5, 0.0), 0.02, 0.02);
	pcl_generator.adjust_pcl_size_to_fit_elems(tri_mesh);
	
	Model_T2D_CHM_mt model;
	model.init_mesh(tri_mesh);
	model.init_search_grid(tri_mesh, 0.05, 0.05);
	tri_mesh.clear();
	model.init_pcls(pcl_generator, 0.6, 2650.0, 1000.0, 1.0e7, 1.0e-12, 1.0e-3);
	pcl_generator.clear();

	const size_t pcl_num = model.get_pcl_num();
	MatModel::MaterialModel **mms = model.get_mat_models();
	// mcc
	MatModel::ModifiedCamClay* mccs = model.add_ModifiedCamClay(pcl_num);
	double K = 1.0 - sin(23.5 / 180.0 * 3.14159165359);
	double ini_stress[6] = { -12025.0, -20000.0, -12025.0, 0.0, 0.0, 0.0 };
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		//pcl.s11 = ini_stress[0];
		//pcl.s22 = ini_stress[1];
		//pcl.s12 = 0.0;
		MatModel::ModifiedCamClay& mcc = mccs[p_id];
		mcc.set_param_NC(0.3, 0.044, 0.205, 23.5, 3.6677, ini_stress);
		mms[p_id] = &mcc;
	}

	model.init_rigid_circle(1.0e5, 1.0e3, 0.0, 0.5, 0.5);
	model.set_rigid_circle_velocity(0.0, -0.5, 0.0);

	IndexArray left_right_bc_pt_array(100);
	find_2d_nodes_on_x_line<Model_T2D_CHM_mt, Model_T2D_CHM_mt::Position>(model, left_right_bc_pt_array, -3.5);
	find_2d_nodes_on_x_line<Model_T2D_CHM_mt, Model_T2D_CHM_mt::Position>(model, left_right_bc_pt_array, 3.5, false);
	model.init_fixed_vx_s_bc(left_right_bc_pt_array.get_num(), left_right_bc_pt_array.get_mem());
	model.init_fixed_vx_f_bc(left_right_bc_pt_array.get_num(), left_right_bc_pt_array.get_mem());

	IndexArray bottom_bc_pt_array(100);
	find_2d_nodes_on_y_line<Model_T2D_CHM_mt, Model_T2D_CHM_mt::Position>(model, bottom_bc_pt_array, -5.0);
	model.init_fixed_vy_s_bc(bottom_bc_pt_array.get_num(), bottom_bc_pt_array.get_mem());
	model.init_fixed_vy_f_bc(bottom_bc_pt_array.get_num(), bottom_bc_pt_array.get_mem());
	
	//QtApp_Prep_T2D_CHM_mt md_disp(argc, argv);
	//md_disp.set_win_size(1200, 900);
	//md_disp.set_model(model);
	////md_disp.set_pts_from_vx_s_bc(0.03);
	////md_disp.set_pts_from_vx_f_bc(0.03);
	////md_disp.set_pts_from_vy_s_bc(0.03);
	//md_disp.set_pts_from_vy_f_bc(0.03);
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
	res_file_hdf5.create("t2d_chm_mt_pipe_conference1.h5");

	ModelData_T2D_CHM_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T2D_CHM_mt_complete out("penetration");
	out.set_res_file(res_file_hdf5);
	out.set_interval_num(100);
	out.set_output_init_state();
	out.set_output_final_state();
	TimeHistory_ConsoleProgressBar out_pb;

	Step_T2D_CHM_mt step("step1");
	step.set_model(model);
	step.set_step_time(1.0);
	//step.set_step_time(5.0e-4);
	step.set_dtime(2.0e-6);
	step.set_thread_num(7);
	step.add_time_history(out);
	step.add_time_history(out_pb);
	step.solve();
}

#include "QtApp_Posp_T2D_CHM_mt.h"
#include "test_model_view_omp.h"

void test_t2d_chm_mt_pipe_conference_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t2d_chm_mt_pipe_conference1.h5");

	QtApp_Posp_T2D_CHM_mt app(argc, argv, QtApp_Posp_T2D_CHM_mt::Animation);
	app.set_ani_time(5.0);
	app.set_win_size(1200, 950);
	app.set_display_range(-3.6, 3.6, -5.1, 0.6);
	//app.set_res_file(rf, "penetration", Hdf5Field::s22);
	//app.set_color_map_fld_range(-10000.0, 10000.0); // s22
	//app.set_res_file(rf, "penetration", Hdf5Field::p);
	//app.set_color_map_fld_range(0, 20000.0); // pore pressure
	//app.set_res_file(rf, "penetration", Hdf5Field::mises_strain_2d);
	//app.set_color_map_fld_range(0, 0.1); // mises strain
	app.set_res_file(rf, "penetration", Hdf5Field::max_shear_stress);
	app.set_color_map_fld_range(0, 100.0); 
	app.set_color_map_geometry(1.0f, 0.45f, 0.5f);
	//app.set_png_name("t2d_chm_mt_pipe_conference1");
	//app.set_gif_name("t2d_chm_mt_pipe_conference1");
	app.start();
}