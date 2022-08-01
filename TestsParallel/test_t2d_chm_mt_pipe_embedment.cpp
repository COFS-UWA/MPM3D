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
	//tri_mesh.load_mesh_from_hdf5("../../Asset/rect_pipe_conference_mesh2_half.h5");
	tri_mesh.load_mesh_from_hdf5("../../Asset/rect_strip_footing_half.h5");
	tri_mesh.init_search_grid(0.02, 0.02);

	ParticleGenerator2D<TriangleMesh> pcl_generator;
	//pcl_generator.generate_pcls_in_grid_layout(Rect(0.0, 5.0, -3.5, 0.0), 0.03, 0.03);
	//pcl_generator.generate_pcls_in_grid_layout(Rect(0.0, 5.0, -5.0, -3.5), 0.03, 0.03);
	//pcl_generator.replace_with_pcls_in_grid_layout(Rect(0.0, 3.5, -3.5, 0.0), 0.01, 0.01);
	pcl_generator.generate_pcls_in_grid_layout(Rect(0.0, 9.0, -4.5, 0.0), 0.03, 0.03);
	pcl_generator.generate_pcls_in_grid_layout(Rect(0.0, 9.0, -9.0, -4.5), 0.03, 0.03);
	pcl_generator.replace_with_pcls_in_grid_layout(Rect(0.0, 4.5, -4.5, 0.0), 0.01, 0.01);
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
		//pcl_s.s12 = 0.0;
		mccs->set_param_NC(0.3, 0.044, 0.205, 23.5, 3.6677, ini_stress);
		model.add_mat_model(p_id, *mccs, sizeof(MatModel::ModifiedCamClay));
		mccs = model.following_ModifiedCamClay(mccs);
	}

	model.init_rigid_circle(0.0, 1.0, 1.0, 1.0);
	model.set_rigid_circle_velocity(0.0, -0.5, 0.0); // -0.2
	model.set_contact_param(1.0e5 / 0.02, 1.0e5 / 0.02, 0.0, 5.0e3, 2.0e3 / 0.02, 2.0e3 / 0.02);
	//model.set_smooth_contact_between_spcl_and_rb();
	//model.set_sticky_contact_between_spcl_and_circle();
	model.set_rough_contact_between_spcl_and_circle();

	//MemoryUtils::ItemArray<double> pt_traction(100);
	//IndexArray traction_pt_array(100);
	//find_2d_pcls<Model_T2D_CHM_mt>(model, traction_pt_array, Rect(0.0, 4.5, -0.006, 0.0));
	//size_t traction_num_tmp = traction_pt_array.get_num();
	//for (size_t p_id = 0; p_id < traction_num_tmp; p_id++)
	//{
	//	double pcl_t = -20.0e3 * 0.01;
	//	pt_traction.add(&pcl_t);
	//}
	//find_2d_pcls<Model_T2D_CHM_mt>(model, traction_pt_array, Rect(4.5, 9.0, -0.016, 0.0), false);
	//for (size_t p_id = traction_num_tmp; p_id < traction_pt_array.get_num(); p_id++)
	//{
	//	double pcl_t = -20.0e3 * 0.03;
	//	pt_traction.add(&pcl_t);
	//}
	//model.init_tys(traction_pt_array.get_num(), traction_pt_array.get_mem(), pt_traction.get_mem());

	constexpr double boundary_pos = 9.0;
	IndexArray left_right_bc_pt_array(100);
	find_2d_nodes_on_x_line<Model_T2D_CHM_mt, Model_T2D_CHM_mt::Position>(model, left_right_bc_pt_array, 0.0);
	model.init_fixed_vx_f_bc(left_right_bc_pt_array.get_num(), left_right_bc_pt_array.get_mem());
	find_2d_nodes_on_x_line<Model_T2D_CHM_mt, Model_T2D_CHM_mt::Position>(model, left_right_bc_pt_array, boundary_pos, false);
	model.init_fixed_vx_s_bc(left_right_bc_pt_array.get_num(), left_right_bc_pt_array.get_mem());

	IndexArray bottom_bc_pt_array(100);
	find_2d_nodes_on_y_line<Model_T2D_CHM_mt, Model_T2D_CHM_mt::Position>(model, bottom_bc_pt_array, -boundary_pos);
	model.init_fixed_vy_s_bc(bottom_bc_pt_array.get_num(), bottom_bc_pt_array.get_mem());
	
	//QtApp_Prep_T2D_CHM_mt md_disp(argc, argv);
	//md_disp.set_win_size(1200, 900);
	//md_disp.set_model(model);
	////md_disp.set_display_range(4.2, 4.8, -0.25, 0.25);
	////md_disp.set_pts_from_vx_s_bc(0.03);
	////md_disp.set_pts_from_vy_s_bc(0.03);
	////md_disp.set_pts_from_vx_f_bc(0.03);
	////md_disp.set_pts_from_vy_f_bc(0.03);
	//md_disp.set_pts_from_pcl_id(traction_pt_array.get_mem(), traction_pt_array.get_num(), 0.004);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_chm_mt_pipe_embedment_geo.h5");

	ModelData_T2D_CHM_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T2D_CHM_mt_complete out("geostatic");
	//TimeHistory_T2D_CHM_mt_Geo_complete out("geostatic");
	out.set_res_file(res_file_hdf5);
	out.set_interval_num(50);
	out.set_output_init_state();
	out.set_output_final_state();
	TimeHistory_ConsoleProgressBar out_pb;
	out_pb.set_interval_num(2000);

	Step_T2D_CHM_mt step("step1");
	//Step_T2D_CHM_mt_Geo step("step1");
	step.set_model(model);
	step.set_step_time(2.0); //2.0
	step.set_dtime(1.0e-6);
	step.set_thread_num(12);
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
		"t2d_chm_mt_pipe_embedment_geo.h5",
		"geostatic", 50);

	model.set_Kf(1.0e7);

	model.init_rigid_circle(0.0, 1.0, 1.0, 1.0);
	model.set_rigid_circle_velocity(0.0, -0.5, 0.0); // -0.2
	model.set_contact_param(1.0e5 / 0.02, 1.0e5 / 0.02, 0.0, 2.5e3, 1.0e4 / 0.02, 1.0e4 / 0.02);
	//model.set_smooth_contact_between_spcl_and_rb();
	//model.set_sticky_contact_between_spcl_and_circle();
	model.set_rough_contact_between_spcl_and_circle();
	
	//QtApp_Prep_T2D_CHM_mt md_disp(argc, argv);
	//md_disp.set_win_size(1200, 900);
	//md_disp.set_model(model);
	////md_disp.set_display_range(-3.6, 3.6, -5.1, 1.1);
	////md_disp.set_pts_from_vx_s_bc(0.03);
	////md_disp.set_pts_from_vy_s_bc(0.03);
	////md_disp.set_pts_from_vx_f_bc(0.03);
	//md_disp.set_pts_from_vy_f_bc(0.03);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_chm_mt_pipe_embedment.h5");

	ModelData_T2D_CHM_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T2D_CHM_mt_complete out("penetration");
	out.set_res_file(res_file_hdf5);
	out.set_interval_num(100);
	out.set_output_init_state();
	out.set_output_final_state();
	TimeHistory_ConsoleProgressBar out_pb;
	out_pb.set_interval_num(2000);

	step.set_model(model);
	step.set_step_time(1.5);
	step.set_dtime(1.0e-6);
	step.set_thread_num(12);
	step.add_time_history(out);
	step.add_time_history(out_pb);
	step.solve();
}

#include "QtApp_Posp_T2D_CHM_mt.h"
#include "test_model_view_omp.h"

void test_t2d_chm_mt_pipe_embedment_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t2d_chm_mt_pipe_embedment_geo_D2_smh.h5");

	QtApp_Posp_T2D_CHM_mt app(argc, argv, QtApp_Posp_T2D_CHM_mt::Animation);
	app.set_ani_time(10.0);
	app.set_win_size(1600, 950);
	app.set_display_range(-1.0, 4.5, -3.5, 1.0);
	app.set_color_map_geometry(1.5f, 0.45f, 0.5f);
	// s22
	//app.set_res_file(rf, "penetration", Hdf5Field::s22);
	//app.set_color_map_fld_range(-40000.0, 0.0);
	// p
	app.set_res_file(rf, "geostatic", Hdf5Field::p);
	app.set_color_map_fld_range(0, 33000.0);
	// mises_strain
	//app.set_res_file(rf, "penetration", Hdf5Field::mises_strain_2d);
	//app.set_color_map_fld_range(0, 0.1); // mises strain
	//
	//app.set_png_name("t2d_chm_mt_pipe_embedment");
	app.set_gif_name("t2d_chm_mt_pipe_embedment");
	app.start();
}