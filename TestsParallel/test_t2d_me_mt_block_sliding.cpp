#include "TestsParallel_pcp.h"

#include "test_parallel_utils.h"
#include "TriangleMesh.h"
#include "ParticleGenerator2D.hpp"
#include "Model_T2D_ME_mt.h"
#include "ModelData_T2D_ME_mt.h"
#include "Step_T2D_ME_mt.h"
#include "TimeHistory_T2D_ME_mt_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "QtApp_Prep_T2D_ME_mt.h"
#include "test_simulations_omp.h"

void test_t2d_me_mt_block_sliding(int argc, char** argv)
{
	TriangleMesh tri_mesh;
	tri_mesh.load_mesh_from_hdf5("../../Asset/rect_mesh_1by10.h5");

	Model_T2D_ME_mt model;
	model.init_mesh(tri_mesh);
	model.init_search_grid(tri_mesh, 0.051, 0.051);

	ParticleGenerator2D<TriangleMesh> pcl_generator;
	ParticleGenerator2D<TriangleMesh>::Particle pcl;
	const double pcl_len = 0.05;
	//const double pcl_len = 0.2;
	pcl.y = pcl_len * 0.5;
	pcl.area = pcl_len * pcl_len;
	const size_t pcl_num = 10.0 / pcl_len;
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		pcl.x = pcl_len * (0.5 + double(pcl_id));
		pcl_generator.add_pcl(pcl);
	}
	model.init_pcls(pcl_generator, 10.0);
	MatModel::MaterialModel** mms = model.get_mat_models();
	MatModel::LinearElasticity* les = model.add_LinearElasticity(model.get_pcl_num());
	for (uint32_t p_id = 0; p_id < model.get_pcl_num(); ++p_id)
	{
		les[p_id].set_param(1000.0, 0.0);
		mms[p_id] = &les[p_id];
	}

	const double Kcn = 10000.0;
	const double Kct = 10000.0;
	//model.init_rigid_rect(1.0, 1.0 + pcl_len - 10.0/2.0/Kcn, 2.0, 2.0, 1.0);
	model.init_rigid_rect(1.0, 1.0 + pcl_len, 2.0, 2.0, 1.0);
	// need a momentum for sticky contact
	model.set_rigid_rect_ext_force(4.0, -10.0);
	//model.set_rigid_rect_ext_force(4.0, -10.0, 3.0); // sticky
	model.set_contact_param(Kcn, Kct, 0.2, 1.5); // pcl_len = 0.05
	//model.set_frictional_contact_between_pcl_and_rect();
	//model.set_sticky_contact_between_pcl_and_rect();
	//model.set_rough_contact_between_pcl_and_rect();

	const size_t node_num = model.get_node_num();
	IndexArray all_n_array(node_num);
	for (size_t n_id = 0; n_id < node_num; ++n_id)
		all_n_array.add(n_id);
	model.init_fixed_vx_bc(node_num, all_n_array.get_mem());
	model.init_fixed_vy_bc(node_num, all_n_array.get_mem());

	QtApp_Prep_T2D_ME_mt md_disp(argc, argv);
	md_disp.set_win_size(1500, 900);
	md_disp.set_model(model);
	md_disp.set_bg_color(1.0, 1.0, 1.0);
	md_disp.set_mesh_color(0.75, 0.75, 0.75);
	md_disp.set_pcl_color(0.26667, 0.44706, 0.76863);
	md_disp.set_rb_color(0.92941, 0.49, 0.19216);
	//md_disp.set_pts_from_node_id(vx_bc_pt_array.get_mem(), vx_bc_pt_array.get_num(), 0.01);
	//md_disp.set_pts_from_node_id(vy_bc_pt_array.get_mem(), vy_bc_pt_array.get_num(), 0.01);
	md_disp.start();
	return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_me_mt_block_sliding.h5");

	ModelData_T2D_ME_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T2D_ME_mt_complete out("slide");
	out.set_interval_num(100);
	out.set_output_init_state();
	out.set_output_final_state();
	out.set_res_file(res_file_hdf5);
	TimeHistory_ConsoleProgressBar out_pb;

	Step_T2D_ME_mt step("step1");
	step.set_model(model);
	step.set_step_time(4.0);
	//step.set_step_time(5.0e-4);
	step.set_dtime(1.0e-5);
	//step.set_thread_num(4);
	step.add_time_history(out);
	step.add_time_history(out_pb);
	step.solve();
}

#include "QtApp_Posp_T2D_ME_mt.h"
#include "test_model_view_omp.h"

void test_t2d_me_mt_block_sliding_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t2d_me_mt_block_sliding.h5");

	// mono color ?
	QtApp_Posp_T2D_ME_mt app(argc, argv, QtApp_Posp_T2D_ME_mt::Animation);
	app.set_win_size(900, 900);
	app.set_ani_time(5.0);
	app.set_res_file(rf, "slide", Hdf5Field::s22);
	app.set_color_map_fld_range(-50.0, 0.0);
	//app.set_color_map_geometry(1.0f, 0.45f, 0.5f);
	//app.set_png_name("t2d_me_mt_block_sliding");
	app.set_gif_name("t2d_me_mt_block_sliding");
	app.start();
}
