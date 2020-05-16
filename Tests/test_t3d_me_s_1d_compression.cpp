#include "Tests_pcp.h"

#include "test_simulations.h"

#include "ItemArray.hpp"

#include "Model_T3D_ME_s.h"
#include "ParticleGenerator3D.hpp"
#include "LinearElasticity.h"
#include "Step_T3D_ME_s.h"

#include "ModelData_T3D_ME_s.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "TimeHistory_T3D_ME_s_complete.h"

#include "ModelViewer3D.h"

namespace
{
typedef MemoryUtils::ItemArray<size_t> IndexArray;

void find_x_bc_nodes(Model_T3D_ME_s &md, IndexArray &pt_array)
{
	pt_array.reset();
	size_t node_num = md.get_node_num();
	Model_T3D_ME_s::Node *nodes = md.get_nodes();
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Model_T3D_ME_s::Node &n = nodes[n_id];
		if (n.x < 1.0e-3 || n.x > 0.2 - 1.0e-3)
			pt_array.add(n.id);
	}
}

void find_y_bc_nodes(Model_T3D_ME_s &md, IndexArray &pt_array)
{
	pt_array.reset();
	size_t node_num = md.get_node_num();
	Model_T3D_ME_s::Node *nodes = md.get_nodes();
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Model_T3D_ME_s::Node &n = nodes[n_id];
		if (n.y < 1.0e-3 || n.y > 0.2 - 1.0e-3)
			pt_array.add(n.id);
	}
}

void find_z_bc_nodes(Model_T3D_ME_s &md, IndexArray &pt_array)
{
	pt_array.reset();
	size_t node_num = md.get_node_num();
	Model_T3D_ME_s::Node *nodes = md.get_nodes();
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Model_T3D_ME_s::Node &n = nodes[n_id];
		if (n.z < 1.0e-3)
			pt_array.add(n.id);
	}
}

void find_top_pcls(Model_T3D_ME_s &md, IndexArray &pt_array)
{
	pt_array.reset();
	size_t pcl_num = md.get_pcl_num();
	Model_T3D_ME_s::Particle *pcls = md.get_pcls();
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		Model_T3D_ME_s::Particle &pcl = pcls[p_id];
		if (pcl.z > 1.0 - 0.011)
			pt_array.add(pcl.id);
	}
}
}

void test_t3d_me_s_1d_compression()
{
	Model_T3D_ME_s model;
	model.load_mesh_from_hdf5("..\\..\\Asset\\bar_mesh.h5");
	std::cout << "node num: " << model.get_node_num() << "\n"
			  << "elem num: " << model.get_elem_num() << "\n";
	
	ParticleGenerator3D<Model_T3D_ME_s> pg;
	Cube bar_box = { 0.0, 0.2, 0.0, 0.2, 0.0, 1.0 };
	pg.generate_pcls_grid(bar_box, 0.02, 0.02, 0.02);
	model.init_pcls(pg, 1.0);

	size_t pcl_num = model.get_pcl_num();
	Model_T3D_ME_s::Particle *pcls = model.get_pcls();
	MatModel::LinearElasticity *mms = model.add_LinearElasticity(pcl_num);
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		Model_T3D_ME_s::Particle &pcl = pcls[pcl_id];
		MatModel::LinearElasticity &mm = mms[pcl_id];
		mm.E = 100.0;
		mm.niu = 0.0;
		pcl.set_mat_model(mm);
	}

	IndexArray pt_array(100);
	find_x_bc_nodes(model, pt_array);
	size_t *vx_bc_n_id = pt_array.get_mem();
	model.init_vxs(pt_array.get_num());
	for (size_t v_id = 0; v_id < model.vx_num; ++v_id)
	{
		VelocityBC &vbc = model.vxs[v_id];
		vbc.node_id = vx_bc_n_id[v_id];
		vbc.v = 0.0;
	}

	find_y_bc_nodes(model, pt_array);
	size_t *vy_bc_n_id = pt_array.get_mem();
	model.init_vys(pt_array.get_num());
	for (size_t v_id = 0; v_id < model.vy_num; ++v_id)
	{
		VelocityBC &vbc = model.vys[v_id];
		vbc.node_id = vy_bc_n_id[v_id];
		vbc.v = 0.0;
	}
	
	find_z_bc_nodes(model, pt_array);
	size_t *vz_bc_n_id = pt_array.get_mem();
	model.init_vzs(pt_array.get_num());
	for (size_t v_id = 0; v_id < model.vz_num; ++v_id)
	{
		VelocityBC &vbc = model.vzs[v_id];
		vbc.node_id = vz_bc_n_id[v_id];
		vbc.v = 0.0;
	}

	find_top_pcls(model, pt_array);
	size_t *tbc_pcl_id = pt_array.get_mem();
	model.init_tzs(pt_array.get_num());
	for (size_t t_id = 0; t_id < model.tz_num; ++t_id)
	{
		TractionBCAtPcl &tbc = model.tzs[t_id];
		tbc.pcl_id = tbc_pcl_id[t_id];
		tbc.t = 0.02 * -10.0;
	}

	MemoryUtils::ItemArray<float> coord_array(300);
	float coord;
	//for (size_t n_id = 0; n_id < sizeof(vx_bc_n_id)/sizeof(vx_bc_n_id[0]); ++n_id)
	//{
	//	Model_T3D_ME_s::Node &n = model.nodes[vx_bc_n_id[n_id]];
	//	pt_coord = double(n.x);
	//	pt_array.add(&pt_coord);
	//	pt_coord = double(n.y);
	//	pt_array.add(&pt_coord);
	//	pt_coord = 0.0f;
	//	pt_array.add(&pt_coord);
	//}
	//for (size_t n_id = 0; n_id < sizeof(vy_bc_n_id) / sizeof(vy_bc_n_id[0]); ++n_id)
	//{
	//	Model_T3D_ME_s::Node &n = model.nodes[vy_bc_n_id[n_id]];
	//	pt_coord = double(n.x);
	//	pt_array.add(&pt_coord);
	//	pt_coord = double(n.y);
	//	pt_array.add(&pt_coord);
	//	pt_coord = 0.0f;
	//	pt_array.add(&pt_coord);
	//}
	//for (size_t p_id = 0; p_id < sizeof(tbc_pcl_id) / sizeof(tbc_pcl_id[0]); ++p_id)
	//{
	//	Model_T3D_ME_s::Particle &pcl = model.pcls[tbc_pcl_id[p_id]];
	//	pt_coord = double(pcl.x);
	//	pt_array.add(&pt_coord);
	//	pt_coord = double(pcl.y);
	//	pt_array.add(&pt_coord);
	//	pt_coord = 0.0f;
	//	pt_array.add(&pt_coord);
	//}

	//// disp only one point
	//MemoryUtilities::ItemArray<GLfloat> pt_array;
	//GLfloat pt_coord;
	//Model_T3D_ME_s::Particle &pt = model.pcls[5];
	////Model_T3D_ME_s::Node &pt = model.nodes[14];
	//pt_array.reserve(3);
	//pt_coord = double(pt.x);
	//pt_array.add(&pt_coord);
	//pt_coord = double(pt.y);
	//pt_array.add(&pt_coord);
	//pt_coord = 0.0f;
	//pt_array.add(&pt_coord);

	view_model_t3d_me(model);
	return;
	//DisplayModel_T2D disp_model;
	//disp_model.init_win();
	//disp_model.init_model(model);
	////disp_model.init_points(pt_array.get_mem(), pt_array.get_num()/3);
	//disp_model.display(-0.05, 0.25, -0.05, 1.05);
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_me_s_1d_compression.h5");

	// output model
	ModelData_T3D_ME_s md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_ConsoleProgressBar out_cpb;
	TimeHistory_T3D_ME_s_complete out1("compression");
	out1.set_res_file(res_file_hdf5);
	out1.set_output_init_state();
	out1.set_interval_num(100);

	Step_T3D_ME_s step("step1");
	step.set_model(model);
	step.set_step_time(1.0);
	step.set_dtime(1.0e-5);
	step.add_time_history(out_cpb);
	step.add_time_history(out1);
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
//}