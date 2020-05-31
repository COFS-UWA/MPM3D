#include "Tests_pcp.h"

#include "test_simulations.h"

#include "ItemArray.hpp"

#include "ParticleGenerator3D.hpp"
#include "LinearElasticity.h"
#include "Model_T3D_ME_s.h"
#include "Step_T3D_ME_s.h"

#include "ModelData_T3D_ME_s.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "TimeHistory_T3D_ME_s_complete.h"

#include "PrepMPM3DApp.h"

#include "utils.h"

namespace
{

void display_model(int argc, char** argv,
	float theta, float fai,
	Model_T3D_ME_s &model, Point3DArray &ptlist)
{
	PrepMPM3DApp view_app(argc, argv);
	view_app.set_view_dir(theta, fai);
	view_app.set_model(model);
	if (ptlist.get_num())
		view_app.set_points(ptlist.get_mem(), ptlist.get_num());
	view_app.start();
}

}

void test_t3d_me_s_1d_compression(int argc, char **argv)
{
	Model_T3D_ME_s model;
	model.load_mesh_from_hdf5("..\\..\\Asset\\bar_mesh1.h5");
	std::cout << "node num: " << model.get_node_num() << "\n"
			  << "elem num: " << model.get_elem_num() << "\n";

	model.init_search_grid(0.05, 0.05, 0.05);

	ParticleGenerator3D<Model_T3D_ME_s> pg;
	//Cube bar_box = { 0.0, 0.1, 0.0, 0.1, 0.0, 0.5 };
	//pg.generate_pcls_grid(bar_box, 0.02, 0.02, 0.02);
	pg.generate_pcls_second_order_gauss(model);
	model.init_pcls(pg, 10.0);
	std::cout << "pcl num: " << model.get_pcl_num() << "\n";

	size_t pcl_num = model.get_pcl_num();
	Model_T3D_ME_s::Particle *pcls = model.get_pcls();
	MatModel::LinearElasticity *mms = model.add_LinearElasticity(pcl_num);
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		Model_T3D_ME_s::Particle &pcl = pcls[pcl_id];
		MatModel::LinearElasticity &mm = mms[pcl_id];
		mm.set_param(100.0, 0.0);
		pcl.set_mat_model(mm);
	}

	IndexArray pt_array(100);
	find_nodes_on_x_plane(model, pt_array, 0.0);
	find_nodes_on_x_plane(model, pt_array, 0.1, false);
	size_t *vx_bc_n_id = pt_array.get_mem();
	model.init_vxs(pt_array.get_num());
	for (size_t v_id = 0; v_id < model.vx_num; ++v_id)
	{
		VelocityBC &vbc = model.vxs[v_id];
		vbc.node_id = vx_bc_n_id[v_id];
		vbc.v = 0.0;
	}

	find_nodes_on_y_plane(model, pt_array, 0.0);
	find_nodes_on_y_plane(model, pt_array, 0.1, false);
	size_t *vy_bc_n_id = pt_array.get_mem();
	model.init_vys(pt_array.get_num());
	for (size_t v_id = 0; v_id < model.vy_num; ++v_id)
	{
		VelocityBC &vbc = model.vys[v_id];
		vbc.node_id = vy_bc_n_id[v_id];
		vbc.v = 0.0;
	}
	
	find_nodes_on_z_plane(model, pt_array, 0.0);
	size_t *vz_bc_n_id = pt_array.get_mem();
	model.init_vzs(pt_array.get_num());
	for (size_t v_id = 0; v_id < model.vz_num; ++v_id)
	{
		VelocityBC &vbc = model.vzs[v_id];
		vbc.node_id = vz_bc_n_id[v_id];
		vbc.v = 0.0;
	}

	find_pcls(model, pt_array, Cube(0.0, 0.1, 0.0, 0.1, 0.5-0.02, 0.5));
	size_t *tbc_pcl_id = pt_array.get_mem();
	model.init_tzs(pt_array.get_num());
	for (size_t t_id = 0; t_id < model.tz_num; ++t_id)
	{
		TractionBCAtPcl &tbc = model.tzs[t_id];
		tbc.pcl_id = tbc_pcl_id[t_id];
		tbc.t = 1.6667e-3 * -1.0;
	}
	std::cout << "tz_num: " << model.tz_num << "\n";

	MemoryUtils::ItemArray<Point3D> ptlist(50);
	//init_vx_bcs_display(model, ptlist);
	init_vy_bcs_display(model, ptlist);
	//init_vz_bcs_display(model, ptlist);
	//init_tz_bcs_display(model, ptlist);
	display_model(argc, argv, 0.0f, 45.0f, model, ptlist);
	return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_me_s_1d_compression.h5");

	ModelData_T3D_ME_s md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_ConsoleProgressBar out_cpb;
	TimeHistory_T3D_ME_s_complete out1("compression");
	out1.set_res_file(res_file_hdf5);
	out1.set_output_init_state();
	out1.set_interval_num(10);

	Step_T3D_ME_s step("step1");
	step.set_model(model);
	step.set_step_time(1.0e-4);
	step.set_dtime(1.0e-5);
	step.add_time_history(out_cpb);
	step.add_time_history(out1);
	step.solve();
}


void t3d_me_s_1d_compression_postprocess(int argc, char **argv)
{
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
//	GA_T3D_ME_s_hdf5 gen;
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
}
