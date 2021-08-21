#include "TestsParallel_pcp.h"

#include "ParticleGenerator3D.hpp"
#include "LinearElasticity.h"
#include "Model_T3D_CHM_mt.h"
#include "Step_T3D_CHM_ud_mt.h"
#include "ModelData_T3D_CHM_mt.h"
#include "TimeHistory_T3D_CHM_mt_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "QtApp_Prep_T3D_CHM_mt_Div.h"

#include "test_parallel_utils.h"
#include "test_simulations_omp.h"

void test_t3d_chm_mt_triaxial_compression(int argc, char **argv)
{
	TetrahedronMesh teh_mesh;
	teh_mesh.load_mesh_from_hdf5("../../Asset/cylinder_2x2_quad_model.h5");
	teh_mesh.rotate_mesh(asin(1.0), 0.0, 0.0);
	teh_mesh.init_search_grid(0.05, 0.05, 0.05);

	ParticleGenerator3D<TetrahedronMesh> pcl_generator;
	pcl_generator.generate_pcls_in_cylinder(
		0.0, 0.0, 0.0,
		0.0, 0.0, 1.0,
		0.0, 1.0,
		asin(1.0), asin(1.0) * 2.0,
		2.0,
		0.05, 0.05, 0.05);
	pcl_generator.adjust_pcl_size_to_fit_elems(teh_mesh);

	constexpr double e0 = 0.55;
	constexpr double den_grain = 2670.0;
	constexpr double den_dry = den_grain / (e0 + 1.0);
	constexpr double den_sat = den_grain / (e0 + 1.0) + 1000 * e0 / (e0 + 1.0);
	constexpr double den_float = den_sat - 1000.0;
	Model_T3D_CHM_mt model;
	model.init_mesh(teh_mesh);
	model.init_search_grid(teh_mesh);
	model.init_pcls(pcl_generator, e0 / (1.0 + e0), den_grain, 1000.0, 2.0e8, 5.0e-9, 1.0);
	MatModel::MaterialModel** mms = model.get_mat_models();
	const size_t pcl_num = model.get_pcl_num();
	//// linear elasticitiy
	//MatModel::LinearElasticity* les = model.add_LinearElasticity(pcl_num);
	//for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	//{
	//	MatModel::LinearElasticity& le = les[pcl_id];
	//	le.set_param(10000.0, 0.3);
	//	mms[pcl_id] = &le;
	//}
	// Tresca
	//MatModel::Tresca* tcs = model.add_Tresca(pcl_num);
	//for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	//{
	//	MatModel::Tresca& tc = tcs[pcl_id];
	//	tc.set_param(10000.0, 0.0, 150.0);
	//	mms[pcl_id] = &tc;
	//}
	// sand hypoplasticity with yield surface
	const double ini_stress[6] = { -10.0e3, -10.0e3, -10.0e3, 0.0, 0.0, 0.0 };
	MatModel::SandHypoplasticityStbWrapper* shps = model.add_SandHypoplasticityStbWrapper(pcl_num);
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		MatModel::SandHypoplasticityStbWrapper& shp = shps[pcl_id];
		shp.set_param(
			ini_stress, 0.55,
			30.0, 1354.0e6, 0.34,
			0.18, 1.27,
			0.49, 0.76, 0.86,
			0.3, 3.6, 200.0,
			200.0, 0.2);
		mms[pcl_id] = &shp;
	}
	
	// rigid cap
	model.init_rigid_cylinder(0.0, 0.0, 2.1, 0.2, 1.2);
	model.set_rigid_cylinder_velocity(0.0, 0.0, -0.1);
	const double Kct = 1.0e5 / (0.025 * 0.025);
	model.set_contact_param(Kct, Kct, 0.1, 0.2, Kct / 100.0, Kct / 100.0);
	
	IndexArray vx_bc_pt_array(100);
	find_3d_nodes_on_x_plane(model, vx_bc_pt_array, 0.0);
	model.init_fixed_vx_s_bc(vx_bc_pt_array.get_num(), vx_bc_pt_array.get_mem());
	model.init_fixed_vx_f_bc(vx_bc_pt_array.get_num(), vx_bc_pt_array.get_mem());

	IndexArray vy_bc_pt_array(100);
	find_3d_nodes_on_y_plane(model, vy_bc_pt_array, 0.0);
	model.init_fixed_vy_s_bc(vy_bc_pt_array.get_num(), vy_bc_pt_array.get_mem());
	model.init_fixed_vy_f_bc(vy_bc_pt_array.get_num(), vy_bc_pt_array.get_mem());

	IndexArray vz_bc_pt_array(100);
	find_3d_nodes_on_z_plane(model, vz_bc_pt_array, 0.0);
	model.init_fixed_vz_s_bc(vz_bc_pt_array.get_num(), vz_bc_pt_array.get_mem());
	model.init_fixed_vz_f_bc(vz_bc_pt_array.get_num(), vz_bc_pt_array.get_mem());

	std::cout << "pcl_num: " << model.get_pcl_num() << "\n"
			  << "elem_num: " << model.get_elem_num() << "\n"
			  << "node_num: " << model.get_node_num() << "\n";
	QtApp_Prep_T3D_CHM_mt_Div<> md_disp(argc, argv);
	md_disp.set_win_size(1200, 950);
	md_disp.set_view_dir(45.0f, -30.0f);
	md_disp.set_light_dir(30.0f, -20.0f);
	//md_disp.set_view_dist_scale(0.5);
	md_disp.set_model(model);
	//md_disp.set_pts_from_node_id(vx_bc_pt_array.get_mem(), vx_bc_pt_array.get_num(), 0.02);
	//md_disp.set_pts_from_node_id(vy_bc_pt_array.get_mem(), vy_bc_pt_array.get_num(), 0.02);
	md_disp.set_pts_from_node_id(vz_bc_pt_array.get_mem(), vz_bc_pt_array.get_num(), 0.02);
	md_disp.start();
	return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_chm_mt_triaxial_compression.h5");

	ModelData_T3D_CHM_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T3D_CHM_mt_complete out1("compression");
	out1.set_interval_num(100);
	out1.set_output_init_state();
	out1.set_output_final_state();
	out1.set_res_file(res_file_hdf5);
	TimeHistory_ConsoleProgressBar out_cpb;
	out_cpb.set_interval_num(500);

	Step_T3D_CHM_ud_mt step("step1");
	step.set_model(model);
	step.set_thread_num(5);
	step.set_step_time(0.8);
	//step.set_step_time(5.0e-5);
	step.set_dtime(5.0e-6);
	step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();
}

#include "QtApp_Posp_T3D_CHM_mt.h"
#include "test_model_view_omp.h"

void test_t3d_chm_mt_triaxial_compression_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t3d_chm_mt_triaxial_compression.h5");

	//QtApp_Posp_T3D_CHM_mt app(argc, argv, QtApp_Posp_T3D_CHM_mt::SingleFrame);
	//app.set_res_file(rf, "compression", 2, Hdf5Field::max_shear_stress);
	
	QtApp_Posp_T3D_CHM_mt app(argc, argv, QtApp_Posp_T3D_CHM_mt::Animation);
	app.set_win_size(1200, 950);
	app.set_ani_time(10.0);
	app.set_view_dir(60.0f, 10.0f);
	app.set_light_dir(45.0f, 10.0f);
	//app.set_view_dist_scale(1.1);
	app.set_color_map_geometry(0.85f, 0.45f, 0.5f);
	// s33
	//app.set_res_file(rf, "compression", Hdf5Field::s33);
	//app.set_color_map_fld_range(-500.0e3, 0.0);
	// shear stress
	//app.set_res_file(rf, "compression", Hdf5Field::max_shear_stress);
	//app.set_color_map_fld_range(0.0, 30.0);
	// p
	app.set_res_file(rf, "compression", Hdf5Field::p);
	app.set_color_map_fld_range(-10.0e3, 10.0e3);
	//
	//app.set_png_name("t3d_chm_mt_triaxial_compression");
	app.set_gif_name("t3d_chm_mt_triaxial_compression");
	app.start();
}
