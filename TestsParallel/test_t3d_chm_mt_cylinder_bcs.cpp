#include "TestsParallel_pcp.h"

#include "ParticleGenerator3D.hpp"
#include "LinearElasticity.h"
#include "Model_T3D_CHM_mt.h"
#include "Step_T3D_CHM_mt_Geo.h"
#include "Step_T3D_CHM_mt.h"
#include "ModelData_T3D_CHM_mt.h"
#include "TimeHistory_T3D_CHM_mt_Geo_complete.h"
#include "TimeHistory_T3D_CHM_mt_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "QtApp_Prep_T3D_CHM_mt.h"
#include "test_parallel_utils.h"
#include "test_simulations_omp.h"

void test_t3d_chm_mt_cylinder_bcs(int argc, char **argv)
{
	TetrahedronMesh teh_mesh;
	teh_mesh.load_mesh_from_hdf5("../../Asset/cylinder_2x2_model.h5");
	teh_mesh.rotate_mesh(asin(1.0), 0.0, 0.0);
	teh_mesh.init_search_grid(0.05, 0.05, 0.05);

	Model_T3D_CHM_mt model;
	model.init_mesh(teh_mesh);
	model.init_search_grid(teh_mesh);

	ParticleGenerator3D<TetrahedronMesh> pcl_generator;
	const double hz = pcl_generator.generate_pcls_in_cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 2.0, 0.05, 0.05, 0.05);
	pcl_generator.adjust_pcl_size_to_fit_elems(teh_mesh);
	model.init_pcls(pcl_generator, 0.3, 20.0, 10.0, 30000.0, 1.0e-4, 1.0);
	size_t pcl_num = model.get_pcl_num();
	std::cout << "pcl_num: " << pcl_num << "\n"
		<< "elem_num: " << model.get_elem_num() << "\n"
		<< "node_num: " << model.get_node_num() << "\n";

	// material model
	MatModel::MaterialModel **mms = model.get_mat_models();
	MatModel::LinearElasticity *les = model.add_LinearElasticity(pcl_num);
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		MatModel::LinearElasticity &le = les[pcl_id];
		le.set_param(1000.0, 0.0);
		mms[pcl_id] = &le;
	}

	// body forces
	IndexArray bfz_pcl_array(pcl_num);
	MemoryUtils::ItemArray<double> bfz_array(pcl_num);
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		double bfz = -1.0;
		bfz_pcl_array.add(pcl_id);
		bfz_array.add(bfz);
	}
	model.init_bfz_ss(pcl_num, bfz_pcl_array.get_mem(), bfz_array.get_mem());
	
	// side bcs
	IndexArray vbc_size_pt_array(100);
	find_3d_nodes_on_cylinder(model, vbc_size_pt_array,
		Vector3D(0.0, 0.0, 0.0), Vector3D(0.0, 0.0, 1.0), 1.0, 2.0);
	const auto *node_pos = model.get_node_pos();
	for (size_t n_id = 0; n_id < vbc_size_pt_array.get_num(); n_id++)
	{
		const size_t nid = vbc_size_pt_array[n_id];
		const auto& np = node_pos[nid];
		model.set_vbc_vec_s(nid, np.x, np.y, 0.0);
		model.set_vbc_vec_f(nid, np.x, np.y, 0.0);
	}

	// bottom bcs
	IndexArray vz_bc_pt_array(100);
	find_3d_nodes_on_z_plane(model, vz_bc_pt_array, 0.0);
	model.init_fixed_vz_s_bc(vz_bc_pt_array.get_num(), vz_bc_pt_array.get_mem());
	model.init_fixed_vz_f_bc(vz_bc_pt_array.get_num(), vz_bc_pt_array.get_mem());

	//QtApp_Prep_T3D_ME_mt md_disp(argc, argv);
	//md_disp.set_win_size(1200, 950);
	//md_disp.set_view_dir(60.0f, 30.0f);
	//md_disp.set_light_dir(45.0f, 20.0f);
	//md_disp.set_model(model);
	////md_disp.set_pts_from_node_id(vbc_size_pt_array.get_mem(), vbc_size_pt_array.get_num(), 0.025);
	//md_disp.set_pts_from_node_id(vz_bc_pt_array.get_mem(), vz_bc_pt_array.get_num(), 0.025);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_chm_mt_cylinder_bcs.h5");

	ModelData_T3D_CHM_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T3D_CHM_mt_Geo_complete out1("geostatic");
	out1.set_res_file(res_file_hdf5);
	out1.set_output_init_state();
	out1.set_output_final_state();
	out1.set_interval_num(100);
	TimeHistory_ConsoleProgressBar out_cpb;

	Step_T3D_CHM_mt_Geo step("step1");
	step.set_model(model);
	step.set_step_time(1.0);
	//step.set_step_time(2.0e-4);
	step.set_dtime(1.0e-4);
	step.set_thread_num(5);
	step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();

	//TimeHistory_T3D_CHM_mt_complete out1("geostatic");
	//out1.set_res_file(res_file_hdf5);
	//out1.set_output_init_state();
	//out1.set_output_final_state();
	//out1.set_interval_num(100);
	//TimeHistory_ConsoleProgressBar out_cpb;

	//Step_T3D_CHM_mt step("step1");
	//step.set_model(model);
	//step.set_thread_num(6);
	//step.set_step_time(1.0);
	////step.set_step_time(1.0e-4);
	//step.set_dtime(1.0e-4);
	//step.add_time_history(out1);
	//step.add_time_history(out_cpb);
	//step.solve();
}

void test_t3d_chm_mt_cylinder_bcs2(int argc, char** argv)
{
	Model_T3D_CHM_mt model;
	Step_T3D_CHM_mt step("step1");
	Model_T3D_CHM_mt_hdf5_utilities::load_model_from_hdf5_file(
		model, step, "t3d_chm_mt_cylinder_bcs.h5", "geostatic", 101);

	IndexArray tbc_pcl_array(100);
	find_3d_pcls(model, tbc_pcl_array, Cube(-1.0, 1.0, -1.0, 1.0, 2.0 - 0.03, 2.0));
	MemoryUtils::ItemArray<double> tzs_mem(tbc_pcl_array.get_num());
	double tz_mag = 0.05 * 0.05 * -10.0;
	for (size_t t_id = 0; t_id < tbc_pcl_array.get_num(); ++t_id)
		tzs_mem.add(tz_mag);
	model.init_tzs(tbc_pcl_array.get_num(), tbc_pcl_array.get_mem(), tzs_mem.get_mem());

	//QtApp_Prep_T3D_CHM_mt md_disp(argc, argv);
	//md_disp.set_win_size(1200, 950);
	//md_disp.set_view_dir(30.0f, -30.0f);
	//md_disp.set_light_dir(45.0f, -30.0f);
	//md_disp.set_model(model);
	////md_disp.set_pts_from_vz_bc_s(0.02);
	////md_disp.set_pts_from_vz_bc_f(0.02);
	////md_disp.set_pts_from_vec_bc_s(0.02);
	////md_disp.set_pts_from_vec_bc_f(0.02);
	////md_disp.set_pts_from_pcl_id(tbc_pcl_array.get_mem(), tbc_pcl_array.get_num(), 0.02);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_chm_mt_cylinder_bcs2.h5");

	ModelData_T3D_CHM_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T3D_CHM_mt_complete out1("consolidation");
	out1.set_res_file(res_file_hdf5);
	out1.set_output_init_state();
	out1.set_output_final_state();
	out1.set_interval_num(100);
	TimeHistory_ConsoleProgressBar out_cpb;

	step.set_model(model);
	step.set_thread_num(6);
	step.set_step_time(5.0);
	//step.set_step_time(2.0e-4);
	step.set_dtime(2.0e-5);
	step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();
}

#include "QtApp_Posp_T3D_CHM_mt.h"
#include "test_model_view_omp.h"

void test_t3d_chm_mt_cylinder_bcs_result(int argc, char **argv)
{
	ResultFile_hdf5 rf;
	rf.open("t3d_chm_mt_cylinder_bcs.h5");
	//rf.open("t3d_chm_mt_cylinder_bcs2.h5");

	QtApp_Posp_T3D_CHM_mt app(argc, argv, QtApp_Posp_T3D_CHM_mt::Animation);
	app.set_ani_time(5.0);
	app.set_win_size(900, 900);
	app.set_view_dir(30.0f, 30.0f);
	app.set_light_dir(90.0f, 30.0f);
	app.set_color_map_geometry(0.7f, 0.45f, 0.5f);
	//
	//app.set_res_file(rf, "geostatic", Hdf5Field::s33);
	//app.set_color_map_fld_range(-30.0, 0.0);
	//
	app.set_res_file(rf, "geostatic", Hdf5Field::p);
	app.set_color_map_fld_range(0.0, 2.0);
	//
	//app.set_res_file(rf, "consolidation", Hdf5Field::s33);
	//app.set_color_map_fld_range(-30.0, 0.0);
	//
	//app.set_res_file(rf, "consolidation", Hdf5Field::p);
	//app.set_color_map_fld_range(0.0, 2.0);
	//app.set_png_name("t3d_chm_mt_cylinder_bcs");
	app.set_gif_name("t3d_chm_mt_cylinder_bcs2");
	app.start();
}
