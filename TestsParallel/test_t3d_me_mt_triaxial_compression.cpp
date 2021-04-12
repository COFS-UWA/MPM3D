#include "TestsParallel_pcp.h"

#include "ParticleGenerator3D.hpp"
#include "LinearElasticity.h"
#include "Model_T3D_ME_mt.h"
#include "Step_T3D_ME_mt.h"
#include "ModelData_T3D_ME_mt.h"
#include "TimeHistory_T3D_ME_mt_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "QtApp_Prep_T3D_ME_mt_Div.h"

#include "test_parallel_utils.h"
#include "test_simulations_omp.h"

void test_t3d_me_mt_triaxial_compression(int argc, char **argv)
{
	TetrahedronMesh teh_mesh;
	teh_mesh.load_mesh_from_hdf5("../../Asset/brick_mesh_1.00_2x2x10.h5");
	teh_mesh.init_search_grid(0.05, 0.05, 0.05);

	ParticleGenerator3D<TetrahedronMesh> pcl_generator;
	pcl_generator.generate_pcls_grid(Cube(0.0, 0.2, 0.0, 0.2, 0.0, 1.0), 0.025, 0.025, 0.025);
	pcl_generator.adjust_pcl_size_to_fit_elems(teh_mesh);

	Model_T3D_ME_mt model;
	model.init_mesh(teh_mesh);
	model.init_search_grid(teh_mesh);
	model.init_pcls(pcl_generator, 10.0);
	MatModel::MaterialModel** mms = model.get_mat_models();
	const size_t pcl_num = model.get_pcl_num();
	//// linear elasticitiy
	//MatModel::LinearElasticity* les = model.add_LinearElasticity(pcl_num);
	//for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	//{
	//	MatModel::LinearElasticity& le = les[pcl_id];
	//	le.set_param(1000.0, 0.3);
	//	mms[pcl_id] = &le;
	//}
	// sand hypoplasticity
	MatModel::SandHypoplasticityByUmat* shps = model.add_SandHypoplasticityByUmat(pcl_num);
	const double ini_stress[6] = { -100.0, -100.0, -100.0, 0.0, 0.0, 0.0 };
	constexpr double R = 1.0e-4;
	const double ig_strain[6] = { -R / sqrt(3.0), -R / sqrt(3.0), -R / sqrt(3.0), 0.0, 0.0, 0.0 };
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		MatModel::SandHypoplasticityByUmat &shp = shps[pcl_id];
		// intergranular strain
		//shp.set_param(ini_stress, 0.817,
		//	33.1, 4.0e9, 0.27, 0.677, 1.054, 1.212, 0.14, 2.5,
		//	2.2, 1.1, R, 0.1, 5.5, ig_strain);
		// no intergranular strain
		shp.set_param(ini_stress, 0.817,
			33.1, 4.0e9, 0.27, 0.677, 1.054, 1.212, 0.14, 2.5,
			0.0, 1.1, R, 0.1, 5.5, ig_strain);
		shp.set_integration_step_ratio(0.5);
		mms[pcl_id] = &shp;
	}
	// sand hypoplasticity
	//const double ini_stress[6] = { -100.0, -100.0, -100.0, 0.0, 0.0, 0.0 };
	//MatModel::SandHypoplasticityWrapper* shps = model.add_SandHypoplasticityWrapper(pcl_num);
	//for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	//{
	//	MatModel::SandHypoplasticityWrapper& shp = shps[pcl_id];
	//	shp.set_param(ini_stress, 0.817,
	//		33.1, 4.0e9, 0.27,
	//		0.677, 1.054, 1.212,
	//		0.14, 2.5);
	//	mms[pcl_id] = &shp;
	//}
	
	// rigid cap
	model.init_rigid_cylinder(0.1, 0.1, 1.025, 0.05, 0.2);
	model.set_rigid_cylinder_velocity(0.0, 0.0, -0.1);
	model.set_contact_param(200000.0, 200000.0, 0.1);

	// top pressure
	//IndexArray tbc_pcl_array(100);
	//find_3d_pcls(model, tbc_pcl_array, Cube(0.0, 0.2, 0.0, 0.2, 1.0 - 0.013, 1.0));
	//MemoryUtils::ItemArray<double> tzs_mem(tbc_pcl_array.get_num());
	//double tz_mag = 0.025 * 0.025 * -10.0;
	//for (size_t t_id = 0; t_id < tbc_pcl_array.get_num(); ++t_id)
	//	tzs_mem.add(tz_mag);
	//model.init_tzs(tbc_pcl_array.get_num(), tbc_pcl_array.get_mem(), tzs_mem.get_mem());

	IndexArray vx_bc_pt_array(100);
	find_3d_nodes_on_x_plane(model, vx_bc_pt_array, 0.0);
	//find_3d_nodes_on_x_plane(model, vx_bc_pt_array, 0.2, false);
	model.init_fixed_vx_bc(vx_bc_pt_array.get_num(), vx_bc_pt_array.get_mem());

	IndexArray vy_bc_pt_array(100);
	find_3d_nodes_on_y_plane(model, vy_bc_pt_array, 0.0);
	//find_3d_nodes_on_y_plane(model, vy_bc_pt_array, 0.2, false);
	model.init_fixed_vy_bc(vy_bc_pt_array.get_num(), vy_bc_pt_array.get_mem());

	IndexArray vz_bc_pt_array(100);
	find_3d_nodes_on_z_plane(model, vz_bc_pt_array, 0.0);
	model.init_fixed_vz_bc(vz_bc_pt_array.get_num(), vz_bc_pt_array.get_mem());

	//QtApp_Prep_T3D_ME_mt_Div<> md_disp(argc, argv);
	//md_disp.set_win_size(1200, 950);
	//md_disp.set_view_dir(200.0f, 30.0f);
	//md_disp.set_light_dir(200.0f, 30.0f);
	////md_disp.set_view_dist_scale(0.5);
	//md_disp.set_model(model);
	////md_disp.set_pts_from_node_id(vx_bc_pt_array.get_mem(), vx_bc_pt_array.get_num(), 0.01);
	//md_disp.set_pts_from_node_id(vy_bc_pt_array.get_mem(), vy_bc_pt_array.get_num(), 0.01);
	////md_disp.set_pts_from_node_id(vz_bc_pt_array.get_mem(), vz_bc_pt_array.get_num(), 0.01);
	//md_disp.start();
	//return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_me_mt_triaxial_compression.h5");

	ModelData_T3D_ME_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T3D_ME_mt_complete out1("compression");
	out1.set_interval_num(100);
	out1.set_output_init_state();
	out1.set_output_final_state();
	out1.set_res_file(res_file_hdf5);
	TimeHistory_ConsoleProgressBar out_cpb;
	out_cpb.set_interval_num(500);

	Step_T3D_ME_mt step("step1");
	step.set_model(model);
	step.set_step_time(0.5);
	//step.set_step_time(1.0e-4);
	step.set_dtime(1.0e-5);
	step.set_thread_num(5);
	step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();
}

#include "QtApp_Posp_T3D_ME_mt.h"
#include "test_model_view_omp.h"

void test_t3d_me_mt_triaxial_compression_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t3d_me_mt_triaxial_compression.h5");

	//QtApp_Posp_T3D_ME_mt app(argc, argv, QtApp_Posp_T3D_ME_mt::SingleFrame);
	//app.set_res_file(rf, "compression", 2, Hdf5Field::max_shear_stress);
	QtApp_Posp_T3D_ME_mt app(argc, argv, QtApp_Posp_T3D_ME_mt::Animation);
	app.set_res_file(rf, "compression", Hdf5Field::s33);
	//app.set_res_file(rf, "compression", Hdf5Field::max_shear_stress);
	app.set_win_size(1200, 950);
	app.set_ani_time(10.0);
	app.set_view_dir(30.0f, 0.0f);
	app.set_light_dir(30.0f, 20.0f);
	//app.set_view_dist_scale(1.1);
	app.set_color_map_fld_range(-500.0, 0.0);
	//app.set_color_map_fld_range(0.0, 30.0);
	app.set_color_map_geometry(0.85f, 0.45f, 0.5f);
	//app.set_png_name("t3d_me_mt_triaxial_compression");
	app.set_gif_name("t3d_me_mt_triaxial_compression");
	app.start();
}
