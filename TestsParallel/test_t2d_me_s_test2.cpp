#include "TestsParallel_pcp.h"

#include "ItemArray.hpp"
#include "test_parallel_utils2.h"
#include "Model_T2D_ME_s.h"
#include "Step_T2D_ME_s.h"
#include "ModelData_T2D_ME_s.h"
#include "TimeHistory_T2D_ME_s_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "QtApp_Prep_T2D_ME_s.h"

#include "test_simulations_omp.h"

void test_t2d_me_s_test2(int argc, char **argv)
{
	Model_T2D_ME_s model;
	model.load_mesh_from_hdf5("../../Asset/rect_mesh.h5");
	model.init_search_grid(0.05, 0.05);

	Model_T2D_ME_s::Particle pcl;
	pcl.x = 0.0700950846;
	pcl.y = 0.946040322;
	auto *elem = model.find_in_which_element(pcl);
	double a1 = elem->a1;
	double a2 = elem->a2;
	double a3 = elem->a3;
	double b1 = elem->b1;
	double b2 = elem->b2;
	double b3 = elem->b3;
	double c1 = elem->coef1;
	double c2 = elem->coef2;
	double c3 = elem->coef3;

	Point2D pt2d(0.0700950846, 0.946040322);
	elem = model.find_in_which_element(pt2d);

	ParticleGenerator2D<Model_T2D_ME_s> pcl_generator;
	pcl_generator.generate_pcls_in_grid_layout(Rect(0.0, 0.2, 0.0, 1.0), 0.02, 0.02);
	model.init_pcls(pcl_generator, 10.0);

	size_t pcl_num = model.get_pcl_num();
	Model_T2D_ME_s::Particle* pcls = model.get_pcls();
	MatModel::LinearElasticity* mms = model.add_LinearElasticity(pcl_num);
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		Model_T2D_ME_s::Particle& pcl = pcls[pcl_id];
		MatModel::LinearElasticity& mm = mms[pcl_id];
		mm.set_param(1000.0, 0.0);
		pcl.set_mat_model(mm);
	}

	// vx bc
	IndexArray vx_bc_pt_array(50);
	find_2d_nodes_on_x_line(model, vx_bc_pt_array, 0.0);
	find_2d_nodes_on_x_line(model, vx_bc_pt_array, 0.2, false);
	size_t* vx_bc_n_id = vx_bc_pt_array.get_mem();
	model.init_vxs(vx_bc_pt_array.get_num());
	size_t vx_num = model.get_vx_num();
	VelocityBC* vxs = model.get_vxs();
	for (size_t v_id = 0; v_id < vx_num; ++v_id)
	{
		VelocityBC &vbc = vxs[v_id];
		vbc.node_id = vx_bc_n_id[v_id];
		vbc.v = 0.0;
	}

	// vy bc
	IndexArray vy_bc_pt_array(50);
	find_2d_nodes_on_y_line(model, vy_bc_pt_array, 0.0);
	size_t* vy_bc_n_id = vy_bc_pt_array.get_mem();
	model.init_vys(vy_bc_pt_array.get_num());
	size_t vy_num = model.get_vy_num();
	VelocityBC* vys = model.get_vys();
	for (size_t v_id = 0; v_id < vy_num; ++v_id)
	{
		VelocityBC &vbc = vys[v_id];
		vbc.node_id = vy_bc_n_id[v_id];
		vbc.v = 0.0;
	}
	
	// traction bc
	IndexArray tbc_pt_array(50);
	find_2d_pcls(model, tbc_pt_array, Rect(0.0, 0.2, 0.987, 1.0));
	size_t* tbc_pcl_id = tbc_pt_array.get_mem();
	model.init_tys(tbc_pt_array.get_num());
	size_t ty_num = model.get_ty_num();
	TractionBCAtPcl* tys = model.get_tys();
	for (size_t t_id = 0; t_id < ty_num; ++t_id)
	{
		TractionBCAtPcl &tbc = tys[t_id];
		tbc.pcl_id = tbc_pcl_id[t_id];
		//tbc.t = 0.05 * -1.0;
		tbc.t = 0.02 * -10.0;
	}

	//QtApp_Prep_T2D_ME_s md_disp(argc, argv);
	//md_disp.set_win_size(900, 900);
	//md_disp.set_model(model);
	////md_disp.set_pts_from_node_id(vx_bc_pt_array.get_mem(), vx_bc_pt_array.get_num(), 0.01);
	////md_disp.set_pts_from_node_id(vy_bc_pt_array.get_mem(), vy_bc_pt_array.get_num(), 0.01);
	////md_disp.set_pts_from_pcl_id(tbc_pt_array.get_mem(), tbc_pt_array.get_num(), 0.01);
	////size_t disp_n_id = 11;
	////md_disp.set_pts_from_node_id(&disp_n_id, 1, 0.01);
	//Model_T2D_ME_s::Particle disp_pcl;
	//disp_pcl.x = 0.0700950846;
	//disp_pcl.y = 0.946040322;
	//md_disp.set_pts(&disp_pcl, 1, 0.01);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_me_s_1d_compression.h5");

	ModelData_T2D_ME_s md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T2D_ME_s_complete out1("compression");
	out1.set_interval_num(100);
	out1.set_res_file(res_file_hdf5);
	out1.set_output_init_state();

	TimeHistory_ConsoleProgressBar out_pb;

	Step_T2D_ME_s step("step1");
	step.set_model(model);
	step.set_step_time(1.0); // 1.0
	step.set_dtime(1.0e-5);
	step.add_time_history(out1);
	step.add_time_history(out_pb);
	step.solve();
}

#include "QtApp_Posp_T2D_ME_s.h"
#include "test_model_view_omp.h"

void test_t2d_me_s_test2_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t2d_me_s_1d_compression.h5");

	QtApp_Posp_T2D_ME_s app(argc, argv, QtApp_Posp_T2D_ME_s::Animation);
	app.set_ani_time(5.0);
	app.set_win_size(900, 900);
	app.set_res_file(rf, "compression", Hdf5Field::s22);
	app.set_color_map_fld_range(-20.0, 0.0);
	//app.set_png_name("t2d_me_s_1d_compression");
	app.start();
}
