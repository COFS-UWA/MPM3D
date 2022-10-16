#include "TestsParallel_pcp.h"

#include "ParticleGenerator3D.hpp"
#include "LinearElasticity.h"
#include "Model_T3D_ME_mt.h"
#include "Step_T3D_ME_mt.h"
#include "ModelData_T3D_ME_mt.h"
#include "TimeHistory_T3D_ME_mt_complete.h"
#include "TimeHistory_ConsoleProgressBar.h"
#include "QtApp_Prep_T3D_ME_mt.h"
#include "test_parallel_utils.h"
#include "test_simulations_omp.h"

void test_contact_model_3d(int argc, char** argv)
{
	TetrahedronMesh teh_mesh;
	teh_mesh.load_mesh_from_hdf5("../../Asset/brick_mesh_0.10_5x5x1.h5");
	teh_mesh.init_search_grid(0.05, 0.05, 0.05);

	Model_T3D_ME_mt model;
	model.init_mesh(teh_mesh);
	model.init_search_grid(teh_mesh);

	ParticleGenerator3D<TetrahedronMesh> pcl_generator;
	pcl_generator.generate_pcls_grid(Cube(0.0, 0.5, 0.0, 0.5, 0.0, 0.1), 0.025, 0.025, 0.025);
	model.init_pcls(pcl_generator, 10.0);
	size_t pcl_num = model.get_pcl_num();
	MatModel::MaterialModel** mms = model.get_mat_models();
	MatModel::LinearElasticity* les = model.add_LinearElasticity(pcl_num);
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		MatModel::LinearElasticity& le = les[pcl_id];
		le.set_param(1000.0, 0.0);
		mms[pcl_id] = &le;
	}

	// penetrate 0.05
	model.init_rigid_cylinder(0.25, 0.25, 0.12, 0.05, 0.2);
	model.set_rigid_cylinder_velocity(0.02, 0.0, 0.0);
	model.set_contact_param(20000.0, 200000.0, 0.01, 0.1);

	IndexArray all_n_ids(model.get_node_num());
	for (size_t n_id = 0; n_id < model.get_node_num(); ++n_id)
		all_n_ids.add(n_id);
	model.init_fixed_vx_bc(all_n_ids.get_num(), all_n_ids.get_mem());
	model.init_fixed_vy_bc(all_n_ids.get_num(), all_n_ids.get_mem());
	model.init_fixed_vz_bc(all_n_ids.get_num(), all_n_ids.get_mem());

	QtApp_Prep_T3D_ME_mt md_disp(argc, argv);
	md_disp.set_win_size(1200, 950);
	md_disp.set_view_dir(30.0f, 30.0f);
	md_disp.set_light_dir(90.0f, 30.0f);
	md_disp.set_model(model);
	//md_disp.set_pts_from_node_id(all_n_ids.get_mem(), all_n_ids.get_num(), 0.01);
	md_disp.start();
	return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_me_mt_cylinder_sliding.h5");

	ModelData_T3D_ME_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T3D_ME_mt_complete out1("slide");
	out1.set_res_file(res_file_hdf5);
	out1.set_output_init_state();
	out1.set_interval_num(50);
	TimeHistory_ConsoleProgressBar out_cpb;

	Step_T3D_ME_mt step("step1");
	step.set_model(model);
	//step.set_step_time(1.0);
	step.set_step_time(1.0e-3);
	step.set_dtime(1.0e-5);
	step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();
}

#include <chrono>

void test_contact_3d_rigid_ball_pap2(int argc, char** argv)
{
	RigidObjectByT3DMesh rb_mh;
	rb_mh.init(1.0, "../../Asset/ba_pap2.h5", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.1);
	rb_mh.init_max_dist(0.5);

	Point3D pt(0.0, 0.0, -1.2);
	double dist;
	Vector3D lnorm;
	Point3D lcontpos;

	std::chrono::high_resolution_clock::time_point t0, t1;
	const size_t repeat_num = 1000000;
	size_t pcl_se_time = 0;
	for (size_t i = 0; i < repeat_num; ++i)
	{
		t0 = std::chrono::high_resolution_clock::now();
		bool res = rb_mh.detect_collision_with_point(pt.x, pt.y, pt.z, 0.3, dist, lnorm, lcontpos);
		t1 = std::chrono::high_resolution_clock::now();
		pcl_se_time += (t1 - t0).count();
	}

	std::cout << "dist: " << dist << ", lnorm: "
		<< lnorm.x << ", " << lnorm.y << ", " << lnorm.z << ",\n";
	const double avg_find_time = double(pcl_se_time) / double(repeat_num);
	std::cout << "search time: " << avg_find_time << ",\n";
}

void test_contact_3d_rigid_cylinder_pap2(int argc, char** argv)
{
	RigidObjectByT3DMesh rb_mh;
	//rb_mh.init(1.0, "../../Asset/cy_pap2.h5", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.1);
	//rb_mh.init(1.0, "C:\\MyData\\Work\\Contact algorithm paper2\\pap2_md\\cy_pap2.h5", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 0.05, 0.05);
	rb_mh.init(1.0, "C:\\MyData\\Work\\Contact algorithm paper2\\pap2_md\\ba_pap2.h5", 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 0.05, 0.05);
	std::cout << "Complete initing model.\n";

	rb_mh.init_max_dist(0.15);
	std::cout << "Complete loading model.\n";

	//Point3D pt(0.0, 0.0, 0.0);
	Point3D pt(0.0, 0.0, -1.0);
	double dist;
	Vector3D lnorm;
	Point3D lcontpos;

	std::chrono::high_resolution_clock::time_point t0, t1;
	constexpr size_t repeat_num = 1000000;
	size_t pcl_se_time;
	
	for (size_t j = 0; j < 10; ++j)
	{
		pcl_se_time = 0;
		for (size_t i = 0; i < repeat_num; ++i)
		{
			t0 = std::chrono::high_resolution_clock::now();
			bool res = rb_mh.detect_collision_with_point(pt.x, pt.y, pt.z, 0.01, dist, lnorm, lcontpos);
			t1 = std::chrono::high_resolution_clock::now();
			pcl_se_time += (t1 - t0).count();
		}
		std::cout << "dist: " << dist << ", lnorm: "
			<< lnorm.x << ", " << lnorm.y << ", " << lnorm.z << ",\n";
		const double avg_find_time = double(pcl_se_time) / double(repeat_num);
		std::cout << "search time: " << avg_find_time << ",\n";
	}
}
