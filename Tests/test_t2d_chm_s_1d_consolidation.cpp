#include "Tests_pcp.h"

#include "ItemArray.hpp"
#include "utils.h"
#include "Model_T2D_CHM_s.h"
#include "Step_T2D_CHM_s.h"
#include "ModelData_T2D_CHM_s.h"
#include "TimeHistory_T2D_CHM_s_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"

#include "test_simulations.h"

void test_t2d_mpm_chm_s_1d_consolidation()
{
	Model_T2D_CHM_s model;
	model.load_mesh_from_hdf5("..\\..\\Asset\\rect_mesh.h5");
	model.init_search_grid(0.05, 0.05);

	ParticleGenerator2D<Model_T2D_CHM_s> pcl_generator;
	pcl_generator.generate_pcls_in_grid_layout(Rect(0.0, 0.2, 0.0, 1.0), 0.02, 0.02);
	model.init_pcls(pcl_generator, 0.4, 20.0, 10.0, 40000.0, 1.0e-4, 1.0);

	size_t pcl_num = model.get_pcl_num();
	Model_T2D_CHM_s::Particle* pcls = model.get_pcls();
	MatModel::LinearElasticity* mms = model.add_LinearElasticity(pcl_num);
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		Model_T2D_CHM_s::Particle& pcl = pcls[pcl_id];
		MatModel::LinearElasticity& mm = mms[pcl_id];
		mm.set_param(1000.0, 0.0);
		pcl.set_mat_model(mm);
	}

	IndexArray pt_array(50);

	find_2d_nodes_on_x_line(model, pt_array, 0.0);
	find_2d_nodes_on_x_line(model, pt_array, 0.2, false);
	//size_t vx_bc_n_id[] = { 0, 3, 15, 16, 17, 18, 19, 20, 21, 22, 23,
	//						1, 2,  5,  6,  7,  8,  9, 10, 11, 12, 13 };
	size_t* vx_bc_n_id = pt_array.get_mem();
	model.init_vsxs(pt_array.get_num());
	for (size_t v_id = 0; v_id < model.vsx_num; ++v_id)
	{
		VelocityBC &vbc = model.vsxs[v_id];
		vbc.node_id = vx_bc_n_id[v_id];
		vbc.v = 0.0;
	}
	model.init_vfxs(pt_array.get_num());
	for (size_t v_id = 0; v_id < model.vfx_num; ++v_id)
	{
		VelocityBC &vbc = model.vfxs[v_id];
		vbc.node_id = vx_bc_n_id[v_id];
		vbc.v = 0.0;
	}

	find_2d_nodes_on_y_line(model, pt_array, 0.0);
	//size_t vy_bc_n_id[] = { 0, 1, 4 };
	size_t* vy_bc_n_id = pt_array.get_mem();
	model.init_vsys(pt_array.get_num());
	for (size_t v_id = 0; v_id < model.vsy_num; ++v_id)
	{
		VelocityBC &vbc = model.vsys[v_id];
		vbc.node_id = vy_bc_n_id[v_id];
		vbc.v = 0.0;
	}
	model.init_vfys(pt_array.get_num());
	for (size_t v_id = 0; v_id < model.vfy_num; ++v_id)
	{
		VelocityBC &vbc = model.vfys[v_id];
		vbc.node_id = vy_bc_n_id[v_id];
		vbc.v = 0.0;
	}
	
	find_2d_pcls(model, pt_array, Rect(0.0, 0.2, 0.987, 1.0));
	size_t *tbc_pcl_id = pt_array.get_mem();
	model.init_tys(pt_array.get_num());
	for (size_t t_id = 0; t_id < model.ty_num; ++t_id)
	{
		TractionBCAtPcl &tbc = model.tys[t_id];
		tbc.pcl_id = tbc_pcl_id[t_id];
		//tbc.t = 0.05 * -1.0;
		tbc.t = 0.02 * -10.0;
	}

	//DisplayModel_T2D disp_model;
	//disp_model.init_win();
	//disp_model.init_model(model);
	//disp_model.init_points(pt_array.get_mem(), pt_array.get_num()/3);
	//disp_model.display(-0.05, 0.25, -0.05, 1.05);
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_1d_consolidation.h5");
	
	ModelData_T2D_CHM_s md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T2D_CHM_s_complete out1("consolidation");
	out1.set_interval_num(100);
	out1.set_res_file(res_file_hdf5);
	out1.set_output_init_state();

	TimeHistory_ConsoleProgressBar out_pb;

	Step_T2D_CHM_s step("step1");
	step.set_model(model);
	step.set_step_time(15.0);
	step.set_dtime(1.0e-5);
	step.add_time_history(out1);
	step.add_time_history(out_pb);
	step.solve();
}

//void test_color_animation_t2d_chm_s_1d_consolidation()
//{
//	double soil_height = 1.0;
//	double soil_width = 0.2;
//	double padding_height = soil_height * 0.05;
//	double padding_width = soil_width * 0.05;
//	ColorGraph::Colori colors[] = {
//		{ 0,   0,   255 },
//		{ 0,   93,  255 },
//		{ 0,   185, 255 },
//		{ 0,   255, 232 },
//		{ 0,   255, 139 },
//		{ 0,   255, 46 },
//		{ 46,  255, 0 },
//		{ 139, 255, 0 },
//		{ 232, 255, 0 },
//		{ 255, 185, 0 },
//		{ 255, 93,  0 },
//		{ 255, 0,   0 }
//	};
//	GA_T2D_CHM_s_hdf5 gen(1000, 1000); // window size
//	gen.init_color_graph(
//		850.0, 500.0, 50.0, 450.0,
//		0.0, 10.0,
//		colors, sizeof(colors) / sizeof(ColorGraph::Colori)
//	);
//	gen.generate(5.0,
//		-padding_width, soil_width + padding_width,
//		-padding_height, soil_height + padding_height,
//		"t2d_mpm_1d_consolidation.h5",
//		"consolidation",
//		"t2d_mpm_1d_consolidation.gif"
//	);
//
//}
