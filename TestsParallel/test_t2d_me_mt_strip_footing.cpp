#include "TestsParallel_pcp.h"

#include <cmath>
#include <iomanip>

#include "TriangleMesh.h"
#include "ParticleGenerator2D.hpp"
#include "Model_T2D_ME_mt.h"
#include "Model_T2D_ME_mt_hdf5_utilities.h"
#include "ModelData_T2D_ME_mt.h"
#include "Step_T2D_ME_mt.h"
#include "TimeHistory_T2D_ME_mt_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "QtApp_Prep_T2D_ME_mt.h"

#include "test_parallel_utils.h"
#include "test_simulations_omp.h"

void test_t2d_me_mt_strip_footing(int argc, char** argv)
{
	TriangleMesh tri_mesh;
	//tri_mesh.load_mesh_from_hdf5("../../Asset/rect_pipe_conference_mesh2_half.h5");
	tri_mesh.load_mesh_from_hdf5("../../Asset/rect_strip_footing_half.h5");
	tri_mesh.init_search_grid(0.02, 0.02);

	ParticleGenerator2D<TriangleMesh> pcl_generator;
	//pcl_generator.generate_pcls_in_grid_layout(Rect(0.0, 5.0, -3.5, 0.0), 0.03, 0.03);
	//pcl_generator.generate_pcls_in_grid_layout(Rect(0.0, 5.0, -5.0, -3.5), 0.03, 0.03);
	//pcl_generator.replace_with_pcls_in_grid_layout(Rect(0.0, 3.5, -3.5, 0.0), 0.01, 0.01);
	pcl_generator.generate_pcls_in_grid_layout(Rect(0.0, 9.0, -4.5, 0.0), 0.03, 0.03);
	pcl_generator.generate_pcls_in_grid_layout(Rect(0.0, 9.0, -9.0, -4.5), 0.03, 0.03);
	pcl_generator.replace_with_pcls_in_grid_layout(Rect(0.0, 4.5, -4.5, 0.0), 0.01, 0.01);
	pcl_generator.adjust_pcl_size_to_fit_elems(tri_mesh);

	Model_T2D_ME_mt model;
	model.init_mesh(tri_mesh);
	model.init_search_grid(tri_mesh, 0.02, 0.02);
	constexpr double density = 1800.0;
	model.init_pcls(pcl_generator, density);
	tri_mesh.clear();
	pcl_generator.clear();

	const size_t pcl_num = model.get_pcl_num();
	std::cout << "pcl_num: " << pcl_num << ",\n"
		<< "elem_num: " << model.get_elem_num() << ",\n"
		<< "node_num: " << model.get_node_num() << ",\n";

	MatModel::MaterialModel** mms = model.get_mat_models();
	// Tresca
	//MatModel::Tresca* tes = model.add_Tresca(model.get_pcl_num());
	//for (size_t p_id = 0; p_id < model.get_pcl_num(); ++p_id)
	//{
	//	tes->set_param(4000.0e3, 0.45, 5.0e3);
	//	model.add_mat_model(p_id, *tes, sizeof(MatModel::Tresca));
	//	tes = model.following_Tresca(tes);
	//}
	// Mohr-Coulomb
	double mc_stress[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	auto* pcl_stresses = model.get_pcl_stress0();
	MatModel::MohrCoulombWrapper* mcs = model.add_MohrCoulombWrapper(model.get_pcl_num());
	for (size_t p_id = 0; p_id < model.get_pcl_num(); ++p_id)
	{
		double depth = model.get_pcl_pos()[p_id].y;
		//
		auto& pcl_s = pcl_stresses[p_id];
		pcl_s.s22 = depth * 9.81 * (density - 1000.0);
		pcl_s.s11 = pcl_s.s22 * 0.5;
		//
		if (depth > -0.05)
			depth = -0.05;
		mc_stress[1] = depth * 9.81 * (density - 1000.0);
		mc_stress[0] = mc_stress[1] * 0.5; // 30.0
		mc_stress[2] = mc_stress[0];
		mcs->set_param(mc_stress, 30.0, 0.0, 0.1, 50.0e6, 0.3);
		model.add_mat_model(p_id, *mcs, sizeof(MatModel::MohrCoulombWrapper));
		mcs = model.following_MohrCoulombWrapper(mcs);
	}

	const double contact_fric_ang = 20.0;
	model.init_rigid_rect(0.0, 0.1, 1.0, 0.2, 1.0);
	model.set_rigid_rect_velocity(0.0, -0.1, 0.0);
	model.set_contact_param(1.0e5 / 0.01, 1.0e5 / 0.01, tan(contact_fric_ang/180.0*3.14159265359), 1.5);
	//model.set_frictional_contact_between_pcl_and_rect();
	//model.set_rough_contact_between_pcl_and_rect();

	// gravity force, float unit weight
	IndexArray bfy_pcl_array(pcl_num);
	MemoryUtils::ItemArray<double> bfy_array(pcl_num);
	double bfy = -9.81 * (density - 1000.0) / density;
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		bfy_pcl_array.add(pcl_id);
		bfy_array.add(bfy);
	}
	model.init_bfys(pcl_num, bfy_pcl_array.get_mem(), bfy_array.get_mem());

	// vx bc
	IndexArray vx_bc_pt_array(50);
	find_2d_nodes_on_x_line(model, vx_bc_pt_array, 0.0);
	find_2d_nodes_on_x_line(model, vx_bc_pt_array, 5.0, false);
	model.init_fixed_vx_bc(vx_bc_pt_array.get_num(), vx_bc_pt_array.get_mem());

	// vy bc
	IndexArray vy_bc_pt_array(50);
	find_2d_nodes_on_y_line(model, vy_bc_pt_array, -5.0);
	model.init_fixed_vy_bc(vy_bc_pt_array.get_num(), vy_bc_pt_array.get_mem());

	QtApp_Prep_T2D_ME_mt md_disp(argc, argv);
	md_disp.set_win_size(1500, 950);
	md_disp.set_model(model);
	md_disp.set_display_range(-0.6, 1.0, -1.5, 0.6);
	//md_disp.set_pts_from_node_id(vx_bc_pt_array.get_mem(), vx_bc_pt_array.get_num(), 0.02);
	md_disp.set_pts_from_node_id(vy_bc_pt_array.get_mem(), vy_bc_pt_array.get_num(), 0.02);
	md_disp.start();
	return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t2d_me_mt_strip_footing.h5");

	ModelData_T2D_ME_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T2D_ME_mt_complete out1("loading");
	out1.set_output_init_state();
	out1.set_interval_num(100);
	out1.set_res_file(res_file_hdf5);
	TimeHistory_ConsoleProgressBar out_pb;
	out_pb.set_interval_num(2000);

	Step_T2D_ME_mt step("step1");
	step.set_model(model);
	step.set_step_time(0.5); // 1.0
	step.set_dtime(5.0e-6);
	step.set_thread_num(5);
	step.add_time_history(out1);
	step.add_time_history(out_pb);
	step.solve();
}

#include "QtApp_Posp_T2D_ME_mt.h"
#include "test_model_view_omp.h"

void test_t2d_me_mt_strip_footing_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t2d_me_mt_strip_footing.h5");

	QtApp_Posp_T2D_ME_mt app(argc, argv, QtApp_Posp_T2D_ME_mt::Animation);
	app.set_win_size(1400, 1000);
	app.set_ani_time(10.0);
	//app.set_display_range(0.0, 2.0, -3.5, 0.3);
	app.set_color_map_geometry(1.05, 0.3f, 0.5f);
	// s22
	//app.set_res_file(rf, "loading", Hdf5Field::s22);
	//app.set_color_map_fld_range(-1.0e5, 0.0);
	// mat s22
	app.set_res_file(rf, "loading", Hdf5Field::mat_s22);
	app.set_color_map_fld_range(-1.0e5, 0.0);
	// pe
	//app.set_res_file(rf, "loading", Hdf5Field::plastic_mises_strain_2d);
	//app.set_color_map_fld_range(0.0, 0.05);
	//
	//app.set_gif_name("t2d_me_mt_strip_footing");

	//QtApp_Posp_T2D_ME_mt app(argc, argv, QtApp_Posp_T2D_ME_mt::SingleFrame);
	//app.set_win_size(2000, 900);
	////app.set_display_range(0.0, 4.0, -3.5, 0.3);
	////app.set_color_map_geometry(1.6, 0.3f, 0.5f);
	////app.set_res_file(rf, "loading", 86, Hdf5Field::s22);
	////app.set_color_map_fld_range(-25.7e3, 0.0);
	////
	//app.set_display_range(0.0, 2.2, -1.5, 0.3);
	//app.set_color_map_geometry(1.75f, 0.3f, 0.5f);
	//app.set_res_file(rf, "loading", 86, Hdf5Field::plastic_mises_strain_2d);
	//app.set_color_map_fld_range(0.0, 0.05);
	////
	////app.set_res_file(rf, "loading", Hdf5Field::vx);
	////app.set_color_map_fld_range(-0.003, 0.003);

	////app.set_png_name("t2d_me_mt_strip_footing");

	app.start();
}
