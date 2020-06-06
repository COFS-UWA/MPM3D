#include "Tests_pcp.h"

#include "ItemArray.hpp"

#include "ParticleGenerator3D.hpp"
#include "LinearElasticity.h"
#include "Model_T3D_CHM_s.h"
#include "Step_T3D_CHM_s.h"
#include "ModelData_T3D_CHM_s.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "TimeHistory_T3D_CHM_s_complete.h"

#include "PrepMPM3DApp.h"

#include "test_simulations.h"

#include "PospMPM3DApp.h"

#include "utils.h"
#include "test_model_view.h"

void test_t3d_chm_s_1d_consolidation(int argc, char **argv)
{
	Model_T3D_CHM_s model;
	model.load_mesh_from_hdf5("..\\..\\Asset\\bar_mesh1.h5");
	std::cout << "node num: " << model.get_node_num() << "\n"
			  << "elem num: " << model.get_elem_num() << "\n";

	model.init_search_grid(0.05, 0.05, 0.05);

	ParticleGenerator3D<Model_T3D_CHM_s> pg;
	//Cube bar_box = { 0.0, 0.1, 0.0, 0.1, 0.0, 0.5 };
	//pg.generate_pcls_grid(bar_box, 0.02, 0.02, 0.02);
	pg.generate_pcls_second_order_gauss(model);
	model.init_pcls(pg, 0.3, 20.0, 10.0, 1.0e4, 1.0, 1.0);
	std::cout << "pcl num: " << model.get_pcl_num() << "\n";

	size_t pcl_num = model.get_pcl_num();
	Model_T3D_CHM_s::Particle *pcls = model.get_pcls();
	MatModel::LinearElasticity *mms = model.add_LinearElasticity(pcl_num);
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		Model_T3D_CHM_s::Particle &pcl = pcls[pcl_id];
		MatModel::LinearElasticity &mm = mms[pcl_id];
		mm.set_param(100.0, 0.0);
		pcl.set_mat_model(mm);
	}

	IndexArray pt_array(100);

	find_nodes_on_x_plane(model, pt_array, 0.0);
	find_nodes_on_x_plane(model, pt_array, 0.1, false);
	size_t *vx_bc_n_id = pt_array.get_mem();
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

	find_nodes_on_y_plane(model, pt_array, 0.0);
	find_nodes_on_y_plane(model, pt_array, 0.1, false);
	size_t *vy_bc_n_id = pt_array.get_mem();
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
		VelocityBC& vbc = model.vfys[v_id];
		vbc.node_id = vy_bc_n_id[v_id];
		vbc.v = 0.0;
	}

	find_nodes_on_z_plane(model, pt_array, 0.0);
	size_t *vz_bc_n_id = pt_array.get_mem();
	model.init_vszs(pt_array.get_num());
	for (size_t v_id = 0; v_id < model.vsz_num; ++v_id)
	{
		VelocityBC &vbc = model.vszs[v_id];
		vbc.node_id = vz_bc_n_id[v_id];
		vbc.v = 0.0;
	}
	model.init_vfzs(pt_array.get_num());
	for (size_t v_id = 0; v_id < model.vfz_num; ++v_id)
	{
		VelocityBC& vbc = model.vfzs[v_id];
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
	//init_vsx_bcs_display(model, ptlist);
	//init_vsy_bcs_display(model, ptlist);
	//init_vsz_bcs_display(model, ptlist);
	//init_vfx_bcs_display(model, ptlist);
	//init_vfy_bcs_display(model, ptlist);
	//init_vfz_bcs_display(model, ptlist);
	//init_tz_bcs_display(model, ptlist);
	//display_model(argc, argv, 60.0, 30.0, 40.0, 20.0, model, ptlist, 1.0e-5);
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_chm_s_1d_consolidation.h5");

	ModelData_T3D_CHM_s md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_ConsoleProgressBar out_cpb;
	TimeHistory_T3D_CHM_s_complete out1("consolidation");
	out1.set_res_file(res_file_hdf5);
	out1.set_output_init_state();
	out1.set_interval_num(50);

	Step_T3D_CHM_s step("step1");
	step.set_model(model);
	step.set_step_time(1.0);
	step.set_dtime(1.0e-5);
	step.add_time_history(out_cpb);
	step.add_time_history(out1);
	step.solve();
}


void test_t3d_chm_s_1d_consolidation_result(int argc, char **argv)
{
	PospMPM3DApp app(argc, argv, PospMPM3DApp::Animation);
	app.set_view_dir(10.0f, 30.0f);
	app.set_light_dir(10.0f, 30.0f);

	app.set_ani_time(5.0);
	app.set_gif_name("1d_consolidation.gif");

	app.init_color_scale(-1.5, 0.0,
		ColorScaleExamples::get_color_scale(),
		ColorScaleExamples::get_color_num());

	ResultFile_hdf5 rf;
	rf.open("t3d_chm_s_1d_consolidation.h5");
	int res = app.set_res_file(
					rf,
					"consolidation",
					"p",
					MPM3DModelView::BallShape
					);

	app.start();
}
