#include "TestsParallel_pcp.h"

#include "Model_T2D_CHM_mt.h"
#include "ModelData_T2D_CHM_mt.h"
#include "Step_T2D_CHM_mt_Geo.h"
#include "TimeHistory_T2D_CHM_mt_Geo_complete.h"
#include "Step_T2D_CHM_mt.h"
#include "TimeHistory_T2D_CHM_mt_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "QtApp_Prep_T2D_CHM_mt.h"

#include "test_parallel_utils.h"
#include "test_simulations_omp.h"

void test_t2d_chm_mt_pipe_embedment(int argc, char** argv)
{
	TriangleMesh tri_mesh;
	tri_mesh.load_mesh_from_hdf5("../../Asset/rect_pipe_conference_mesh2_half.h5");
	tri_mesh.init_search_grid(0.02, 0.02);

	ParticleGenerator2D<TriangleMesh> pcl_generator;
	pcl_generator.generate_pcls_in_grid_layout(Rect(0.0, 5.0, -3.5, 0.0), 0.03, 0.03);
	pcl_generator.generate_pcls_in_grid_layout(Rect(0.0, 5.0, -5.0, -3.5), 0.03, 0.03);
	pcl_generator.replace_with_pcls_in_grid_layout(Rect(0.0, 3.5, -3.5, 0.0), 0.01, 0.01);
	pcl_generator.adjust_pcl_size_to_fit_elems(tri_mesh);
	
	Model_T2D_CHM_mt model;
	model.init_mesh(tri_mesh);
	model.init_search_grid(tri_mesh, 0.02, 0.02);
	//model.init_pcls(pcl_generator, 0.6, 2650.0, 1000.0, 1.0e7, 1.0e-12, 1.0e-3); // dt = 1e-6
	model.init_pcls(pcl_generator, 0.6, 2650.0, 1000.0, 1.0e7, 5.0e-13, 1.0e-3); // dt = 1e-6
	//model.init_pcls(pcl_generator, 0.6, 2650.0, 1000.0, 2.0e7, 2.0e-13, 1.0e-3); // dt = 4e-7
	//model.init_pcls(pcl_generator, 0.6, 2650.0, 1000.0, 2.0e7, 1.0e-13, 1.0e-3); // dt = 2e-7
	tri_mesh.clear();
	pcl_generator.clear();

	std::cout << "pcl_num: " << model.get_pcl_num() << ",\n"
		<< "elem_num: " << model.get_elem_num() << ",\n"
		<< "node_num: " << model.get_node_num() << ",\n";

	const size_t pcl_num = model.get_pcl_num();
	MatModel::MaterialModel **mms = model.get_mat_models();
	// mcc
	MatModel::ModifiedCamClay* mccs = model.add_ModifiedCamClay(pcl_num);
	const double K = 1.0 - sin(23.5 / 180.0 * 3.14159165359);
	double ini_stress[6] = { K * -20000.0, -20000.0, K * -20000.0, 0.0, 0.0, 0.0 };
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		//Model_T2D_CHM_mt::Stress &pcl_s = model.get_pcl_stress0()[p_id];
		//pcl_s.s11 = ini_stress[0];
		//pcl_s.s22 = ini_stress[1];
		mccs->set_param_NC(0.3, 0.044, 0.205, 23.5, 3.6677, ini_stress);
		model.add_mat_model(p_id, *mccs, sizeof(MatModel::ModifiedCamClay));
		mccs = model.following_ModifiedCamClay(mccs);
	}

	model.init_rigid_circle(0.0, 0.5, 0.5, 1.0);
	model.set_rigid_circle_velocity(0.0, -0.2, 0.0); // -0.2
	model.set_contact_param(1.0e5 / 0.02, 1.0e5 / 0.02, 0.0, 2.5e3, 2.0e3 / 0.02, 2.0e3 / 0.02);
	//model.set_sticky_contact_between_spcl_and_circle();
	//model.set_rough_contact_between_spcl_and_circle();

	//IndexArray traction_pt_array(100);
	//find_2d_pcls<Model_T2D_CHM_mt>(model, traction_pt_array, Rect(0.0, 3.5, -0.006, 0.0));
	//find_2d_pcls<Model_T2D_CHM_mt>(model, traction_pt_array, Rect(3.5, 5.0, -0.016, 0.0), false);
	//MemoryUtils::ItemArray<double> pt_traction(100);
	//for (size_t p_id = 0; p_id < traction_pt_array.get_num(); p_id++)
	//{
	//	double pcl_t = -20.0e3 * sqrt(model.get_pcl_vol()[p_id]);
	//	pt_traction.add(&pcl_t);
	//}
	//model.init_tys(traction_pt_array.get_num(), traction_pt_array.get_mem(), pt_traction.get_mem());

	IndexArray left_right_bc_pt_array(100);
	find_2d_nodes_on_x_line<Model_T2D_CHM_mt, Model_T2D_CHM_mt::Position>(model, left_right_bc_pt_array, 0.0);
	model.init_fixed_vx_f_bc(left_right_bc_pt_array.get_num(), left_right_bc_pt_array.get_mem());
	find_2d_nodes_on_x_line<Model_T2D_CHM_mt, Model_T2D_CHM_mt::Position>(model, left_right_bc_pt_array, 5.0, false);
	model.init_fixed_vx_s_bc(left_right_bc_pt_array.get_num(), left_right_bc_pt_array.get_mem());

	IndexArray bottom_bc_pt_array(100);
	find_2d_nodes_on_y_line<Model_T2D_CHM_mt, Model_T2D_CHM_mt::Position>(model, bottom_bc_pt_array, -5.0);
	model.init_fixed_vy_s_bc(bottom_bc_pt_array.get_num(), bottom_bc_pt_array.get_mem());
	
	//QtApp_Prep_T2D_CHM_mt md_disp(argc, argv);
	//md_disp.set_win_size(1200, 900);
	//md_disp.set_model(model);
	////md_disp.set_pts_from_vx_s_bc(0.03);
	////md_disp.set_pts_from_vy_s_bc(0.03);
	////md_disp.set_pts_from_vx_f_bc(0.03);
	//md_disp.set_pts_from_vy_f_bc(0.03);
	////md_disp.set_pts_from_pcl_id(traction_pt_array.get_mem(), traction_pt_array.get_num(), 0.01);
	////md_disp.set_display_range(-3.6, 3.6, -5.1, 1.1);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_chm_mt_pipe_embedment.h5");

	ModelData_T2D_CHM_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T2D_CHM_mt_complete out("geostatic");
	//TimeHistory_T2D_CHM_mt_Geo_complete out("geostatic");
	out.set_res_file(res_file_hdf5);
	out.set_interval_num(100);
	out.set_output_init_state();
	out.set_output_final_state();
	TimeHistory_ConsoleProgressBar out_pb;
	out_pb.set_interval_num(2000);

	Step_T2D_CHM_mt step("step1");
	//Step_T2D_CHM_mt_Geo step("step1");
	step.set_model(model);
	step.set_step_time(2.5);
	step.set_dtime(1.0e-6);
	step.set_thread_num(2);
	step.add_time_history(out);
	step.add_time_history(out_pb);
	step.solve();
}

void test_t2d_chm_mt_pipe_embedment_restart(int argc, char** argv)
{
	Model_T2D_CHM_mt model;
	Step_T2D_CHM_mt step("step2");

	using Model_T2D_CHM_mt_hdf5_utilities::load_chm_mt_model_from_hdf5_file;
	load_chm_mt_model_from_hdf5_file(
		model, step,
		"t2d_chm_mt_pipe_embedment.h5",
		"penetration", 5);

	//QtApp_Prep_T2D_CHM_mt md_disp(argc, argv);
	//md_disp.set_win_size(1200, 900);
	//md_disp.set_model(model);
	////md_disp.set_pts_from_vx_s_bc(0.03);
	////md_disp.set_pts_from_vx_f_bc(0.03);
	////md_disp.set_pts_from_vy_s_bc(0.03);
	//md_disp.set_pts_from_vy_f_bc(0.03);
	//// all
	////md_disp.set_display_range(-3.6, 3.6, -5.1, 1.1);
	//// left
	////md_disp.set_display_range(-3.8, -2.2, -1.0, 1.0);
	//// middle
	//md_disp.set_display_range(-1.5, 1.5, -0.75, 0.25);
	//// right
	////md_disp.set_display_range(2.2, 3.8, -1.0, 1.0);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_chm_mt_pipe_embedment2.h5");

	ModelData_T2D_CHM_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T2D_CHM_mt_complete out("penetration");
	out.set_res_file(res_file_hdf5);
	out.set_interval_num(100);
	out.set_output_init_state();
	out.set_output_final_state();
	TimeHistory_ConsoleProgressBar out_pb;

	step.set_model(model);
	//step.set_step_time(5.0);
	step.set_step_time(1.0e-5);
	step.set_dtime(2.0e-6);
	step.set_thread_num(20);
	step.add_time_history(out);
	step.add_time_history(out_pb);
	step.solve();
}

#include "QtApp_Posp_T2D_CHM_mt.h"
#include "test_model_view_omp.h"

void test_t2d_chm_mt_pipe_embedment_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t2d_chm_mt_pipe_embedment.h5");
	//rf.open("t2d_chm_mt_pipe_embedment2.h5");

	QtApp_Posp_T2D_CHM_mt app(argc, argv, QtApp_Posp_T2D_CHM_mt::Animation);
	app.set_ani_time(5.0);
	app.set_win_size(1200, 950);
	app.set_display_range(-3.6, 3.6, -5.1, 0.6);
	app.set_res_file(rf, "penetration", Hdf5Field::s22);
	app.set_color_map_fld_range(-30000.0, -10000.0); // s22
	//app.set_res_file(rf, "penetration", Hdf5Field::p);
	//app.set_color_map_fld_range(0, 20000.0); // pore pressure
	//app.set_res_file(rf, "penetration", Hdf5Field::mises_strain_2d);
	//app.set_color_map_fld_range(0, 0.4); // mises strain
	app.set_color_map_geometry(1.0f, 0.45f, 0.5f);
	//app.set_png_name("t2d_chm_mt_pipe_conference2");
	//app.set_gif_name("t2d_chm_mt_pipe_conference2");
	app.start();
}