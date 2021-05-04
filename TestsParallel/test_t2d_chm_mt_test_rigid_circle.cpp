#include "TestsParallel_pcp.h"

#include "test_parallel_utils.h"
#include "TriangleMesh.h"
#include "ParticleGenerator2D.hpp"
#include "Model_T2D_CHM_mt.h"
#include "Step_T2D_CHM_mt.h"
#include "ModelData_T2D_CHM_mt.h"
#include "TimeHistory_T2D_CHM_mt_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "QtApp_Prep_T2D_CHM_mt.h"

#include "test_simulations_omp.h"

void test_t2d_chm_mt_test_rigid_circle(int argc, char** argv)
{
	TriangleMesh tri_mesh;
	tri_mesh.load_mesh_from_hdf5("../../Asset/rect_test_rc_mesh.h5");
	tri_mesh.init_search_grid(0.5, 0.5);

	ParticleGenerator2D<TriangleMesh> pcl_generator;
	pcl_generator.generate_pcls_in_grid_layout(Rect(0.0, 20.0, 0.0, 15.0), 0.2, 0.2);
	pcl_generator.adjust_pcl_size_to_fit_elems(tri_mesh);
	
	Model_T2D_CHM_mt model;
	model.init_mesh(tri_mesh);
	model.init_search_grid(tri_mesh, 0.5, 0.5);
	tri_mesh.clear();
	model.init_pcls(pcl_generator, 0.6, 2650.0, 1000.0, 5.0e6, 5.0e-12, 1.0e-3);
	pcl_generator.clear();

	const size_t pcl_num = model.get_pcl_num();
	MatModel::MaterialModel** mms = model.get_mat_models();
	// elasticity
	//MatModel::LinearElasticity* mms = model.add_LinearElasticity(pcl_num);
	//for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	//{
	//	MatModel::LinearElasticity& mm = mms[pcl_id];
	//	mm.set_param(2.0e6, 0.3);
	//	pcl.set_mat_model(mm);
	//}
	// mcc
	MatModel::ModifiedCamClay *mccs = model.add_ModifiedCamClay(pcl_num);
	double ini_stress[6] = { -40000.0, -24050.0, -24050.0, 0.0, 0.0, 0.0 };
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		MatModel::ModifiedCamClay &mcc = mccs[p_id];
		mcc.set_param_OC(0.3, 0.044, 0.205, 23.5, 3.6677, ini_stress, 39610.0);
		mms[p_id] = &mcc;
	}

	// rigid circle
	model.init_rigid_circle(10.0, 17.5, 2.5, 1.0);
	model.set_rigid_circle_velocity(0.0, -0.25, 0.0);
	model.set_contact_param(1.0e5 / 0.2, 1.0e5 / 0.2, 0.0, 10.0, 1.0e3 / 0.2, 1.0e3 / 0.2);
	model.set_rough_contact_between_spcl_and_circle();
	model.set_rough_contact_between_fpcl_and_circle();

	// vx bc
	IndexArray vx_bc_pt_array(50);
	find_2d_nodes_on_x_line<Model_T2D_CHM_mt, Model_T2D_CHM_mt::Position>(model, vx_bc_pt_array, 0.0);
	find_2d_nodes_on_x_line<Model_T2D_CHM_mt, Model_T2D_CHM_mt::Position>(model, vx_bc_pt_array, 20.0, false);
	model.init_fixed_vx_s_bc(vx_bc_pt_array.get_num(), vx_bc_pt_array.get_mem());

	// vy bc
	IndexArray vy_bc_pt_array(50);
	find_2d_nodes_on_y_line<Model_T2D_CHM_mt, Model_T2D_CHM_mt::Position>(model, vy_bc_pt_array, 0.0);
	model.init_fixed_vy_s_bc(vy_bc_pt_array.get_num(), vy_bc_pt_array.get_mem());

	//QtApp_Prep_T2D_CHM_mt md_disp(argc, argv);
	//md_disp.set_win_size(900, 900);
	//md_disp.set_model(model);
	//md_disp.set_pts_from_vx_s_bc(0.2);
	////md_disp.set_pts_from_vy_s_bc(0.2);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_chm_mt_test_rigid_circle.h5");

	ModelData_T2D_CHM_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T2D_CHM_mt_complete out1("penetration");
	out1.set_res_file(res_file_hdf5);
	out1.set_interval_num(100);
	out1.set_output_init_state();
	TimeHistory_ConsoleProgressBar out_pb;
	out_pb.set_interval_num(1000);

	Step_T2D_CHM_mt step("step1");
	step.set_model(model);
	step.set_step_time(3.0);
	//step.set_step_time(5.0e-4);
	step.set_dtime(5.0e-6);
	step.set_thread_num(10);
	step.add_time_history(out1);
	step.add_time_history(out_pb);
	step.solve();
}

#include "QtApp_Posp_T2D_CHM_mt.h"
#include "test_model_view_omp.h"

void test_t2d_chm_mt_test_rigid_circle_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t2d_chm_mt_test_rigid_circle.h5");

	QtApp_Posp_T2D_CHM_mt app(argc, argv, QtApp_Posp_T2D_CHM_mt::Animation);
	app.set_ani_time(5.0);
	app.set_win_size(900, 900);
	app.set_res_file(rf, "penetration", Hdf5Field::mises_strain_2d);
	app.set_color_map_fld_range(0.0, 0.14);
	//app.set_res_file(rf, "penetration", Hdf5Field::p);
	//app.set_color_map_fld_range(0.0, 20000.0);
	app.set_color_map_geometry(0.82, 0.5, 0.4);
	//app.set_png_name("t2d_chm_mt_test_rigid_circle");
	//app.set_gif_name("t2d_chm_mt_test_rigid_circle");
	app.start();
}