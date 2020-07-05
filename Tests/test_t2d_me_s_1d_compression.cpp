#include "Tests_pcp.h"

#include "ItemArray.hpp"
#include "utils.h"
#include "Model_T2D_ME_s.h"
#include "Step_T2D_ME_s.h"
#include "ModelData_T2D_ME_s.h"
#include "TimeHistory_T2D_ME_s_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "QtApp_Prep_2DMPM.h"

#include "test_simulations.h"

void test_t2d_mpm_me_s_1d_compression(int argc, char **argv)
{
	Model_T2D_ME_s model;
	model.load_mesh_from_hdf5("..\\..\\Asset\\rect_mesh.h5");
	model.init_search_grid(0.05, 0.05);

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
	for (size_t v_id = 0; v_id < model.vx_num; ++v_id)
	{
		VelocityBC &vbc = model.vxs[v_id];
		vbc.node_id = vx_bc_n_id[v_id];
		vbc.v = 0.0;
	}

	// vy bc
	IndexArray vy_bc_pt_array(50);
	find_2d_nodes_on_y_line(model, vy_bc_pt_array, 0.0);
	size_t* vy_bc_n_id = vy_bc_pt_array.get_mem();
	model.init_vys(vy_bc_pt_array.get_num());
	for (size_t v_id = 0; v_id < model.vy_num; ++v_id)
	{
		VelocityBC &vbc = model.vys[v_id];
		vbc.node_id = vy_bc_n_id[v_id];
		vbc.v = 0.0;
	}
	
	// traction bc
	IndexArray tbc_pt_array(50);
	find_2d_pcls(model, tbc_pt_array, Rect(0.0, 0.2, 0.987, 1.0));
	size_t* tbc_pcl_id = tbc_pt_array.get_mem();
	model.init_tys(tbc_pt_array.get_num());
	for (size_t t_id = 0; t_id < model.ty_num; ++t_id)
	{
		TractionBCAtPcl &tbc = model.tys[t_id];
		tbc.pcl_id = tbc_pcl_id[t_id];
		//tbc.t = 0.05 * -1.0;
		tbc.t = 0.02 * -10.0;
	}

	QtApp_Prep_2DMPM md_disp(argc, argv);
	md_disp.set_win_size(1000, 1000);
	md_disp.set_model(model);
	//md_disp.set_pts_from_node_id(vx_bc_pt_array.get_mem(), vx_bc_pt_array.get_num(), 0.01);
	//md_disp.set_pts_from_node_id(vy_bc_pt_array.get_mem(), vy_bc_pt_array.get_num(), 0.01);
	//md_disp.set_pts_from_pcl_id(tbc_pt_array.get_mem(), tbc_pt_array.get_num(), 0.01);
	md_disp.start();
	return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_me_1d_compression.h5");

	// output model
	ModelData_T2D_ME_s md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T2D_ME_s_complete out1("compression");
	out1.set_interval_num(100);
	out1.set_res_file(res_file_hdf5);
	out1.set_output_init_state();

	TimeHistory_ConsoleProgressBar out_pb;

	Step_T2D_ME_s step("step1");
	step.set_model(model);
	step.set_step_time(1.0);
	step.set_dtime(1.0e-5);
	step.add_time_history(out1);
	step.add_time_history(out_pb);
	step.solve();
}

//void test_color_animation_t2d_me_s_1d_compression(void)
//{
//	double soil_height = 1.0;
//	double soil_width = 0.2;
//	double padding_height = soil_height * 0.05;
//	double padding_width = soil_width * 0.05;
//	// Abaqus "rainbow" spectrum scheme
//	ColorGraph::Colori colors[] = {
//		{ 0,   0,   255 },
//		{ 0,   93,  255 },
//		{ 0,   185, 255 },
//		{ 0,   255, 232 },
//		{ 0,   255, 139 },
//		{ 0,   255, 46  },
//		{ 46,  255, 0   },
//		{ 139, 255, 0   },
//		{ 232, 255, 0   },
//		{ 255, 185, 0   },
//		{ 255, 93,  0   },
//		{ 255, 0,   0   }
//	};
//	GA_T2D_ME_s_hdf5 gen;
//	gen.init_color_graph(-1.0, 1.0, colors, sizeof(colors)/sizeof(ColorGraph::Colori));
//	gen.generate(
//		5.0,
//		-padding_width,
//		soil_width + padding_width,
//		-padding_height,
//		soil_height + padding_height,
//		"t2d_mpm_me_1d_compression.hdf5",
//		"compression",
//		"t2d_mpm_me_1d_compression.gif");
//}
