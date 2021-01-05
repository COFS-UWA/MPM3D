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

namespace
{
	class RigidObjectByT3DMeshTestVer :
		public RigidObjectByT3DMesh
	{
	public:
		int init(double _density, const char* filename,
			double dx, double dy, double dz,
			double dx_ang, double dy_ang, double dz_ang,
			double ghx,	double ghy, double ghz)
		{
			RigidObjectByT3DMesh::init(_density, filename,
				dx, dy, dz, dx_ang, dy_ang, dz_ang, ghx, ghy, ghz);
			size_t z_id = 2;
			for (size_t y_id = 0; y_id < grid.y_num; ++y_id)
			{
				for (size_t x_id = 0; x_id < grid.x_num; ++x_id)
				{
					GridPosType& g = grid_pos_type[offset_from_xyz_id(x_id, y_id, z_id)];
					std::cout << size_t(g);
				}
				std::cout << "\n";
			}
			std::cout << "\n";
			return 0;
		}

		void init_max_dist(double dist)
		{
			RigidObjectByT3DMesh::init_max_dist(dist);
			size_t z_id = 2;
			for (size_t y_id = 0; y_id < grid.y_num; ++y_id)
			{
				for (size_t x_id = 0; x_id < grid.x_num; ++x_id)
				{
					GridPosType &g = grid_pos_type[offset_from_xyz_id(x_id, y_id, z_id)];
					std::cout << size_t(g);
				}
				std::cout << "\n";
			}
		}
	};
}

void test_t3d_rigid_mesh(int argc, char** argv)
{
	RigidObjectByT3DMeshTestVer rb;
	//rb.init(1.0, "../../Asset/brick_mesh_1.00_1x1x1.h5",
	//	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.2);
	rb.init(1.0, "../../Asset/cylinder_model.h5",
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0125, 0.0125, 0.0125);
	rb.init_max_dist(0.01);
	
	const double g_xl = rb.get_grid_xl();
	const double g_yl = rb.get_grid_yl();
	const double g_zl = rb.get_grid_zl();
	const double g_xu = rb.get_grid_xu();
	const double g_yu = rb.get_grid_yu();
	const double g_zu = rb.get_grid_zu();

	double dist;
	Vector3D norm;
	Point3D contpos;
	//bool res = rb.detect_collision_with_point(
	//	0.1, 0.1, 0.015, 0.0, dist, norm, contpos);
	bool res = rb.detect_collision_with_point(
		0.1515, -0.042275, -0.00505 + rb.get_pos().z, 0.00505, dist, norm, contpos);

	int efe = 0;
}

void test_t3d_me_mt_test_rigid_mesh(int argc, char** argv)
{
	TetrahedronMesh teh_mesh;
	teh_mesh.load_mesh_from_hdf5("../../Asset/brick_mesh_1.00_2x2x10.h5");
	teh_mesh.init_search_grid(0.05, 0.05, 0.05);

	ParticleGenerator3D<TetrahedronMesh> pcl_generator;
	pcl_generator.generate_pcls_grid(Cube(0.0, 0.2, 0.0, 0.2, 0.0, 1.0), 0.025, 0.025, 0.025);
	
	Model_T3D_ME_mt model;
	model.init_mesh(teh_mesh);
	model.init_search_grid(teh_mesh);
	model.init_pcls(pcl_generator, 10.0);
	const size_t pcl_num = model.get_pcl_num();
	MatModel::MaterialModel** mms = model.get_mat_models();
	MatModel::LinearElasticity* les = model.add_LinearElasticity(pcl_num);
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		MatModel::LinearElasticity& le = les[pcl_id];
		le.set_param(1000.0, 0.0);
		mms[pcl_id] = &le;
	}

	model.init_t3d_rigid_mesh(1.0, "../../Asset/cylinder_model.h5",
		0.1, 0.1, 1.0, 0.0, 0.0, 0.0, 0.0125, 0.0125, 0.0125);
	model.set_t3d_rigid_mesh_velocity(0.0, 0.0, -0.1);
	model.set_contact_param(20000.0, 20000.0, 0.1);

	IndexArray vx_bc_pt_array(100);
	find_3d_nodes_on_x_plane(model, vx_bc_pt_array, 0.0);
	find_3d_nodes_on_x_plane(model, vx_bc_pt_array, 0.2, false);
	model.init_fixed_vx_bc(vx_bc_pt_array.get_num(), vx_bc_pt_array.get_mem());

	IndexArray vy_bc_pt_array(100);
	find_3d_nodes_on_y_plane(model, vy_bc_pt_array, 0.0);
	find_3d_nodes_on_y_plane(model, vy_bc_pt_array, 0.2, false);
	model.init_fixed_vy_bc(vy_bc_pt_array.get_num(), vy_bc_pt_array.get_mem());

	IndexArray vz_bc_pt_array(100);
	find_3d_nodes_on_z_plane(model, vz_bc_pt_array, 0.0);
	model.init_fixed_vz_bc(vz_bc_pt_array.get_num(), vz_bc_pt_array.get_mem());

	QtApp_Prep_T3D_ME_mt md_disp(argc, argv);
	md_disp.set_win_size(1200, 950);
	md_disp.set_view_dir(30.0f, 30.0f);
	md_disp.set_light_dir(90.0f, 30.0f);
	md_disp.set_model(model);
	//md_disp.set_pts_from_node_id(vx_bc_pt_array.get_mem(), vx_bc_pt_array.get_num(), 0.01);
	//md_disp.set_pts_from_node_id(vy_bc_pt_array.get_mem(), vy_bc_pt_array.get_num(), 0.01);
	//md_disp.set_pts_from_node_id(vz_bc_pt_array.get_mem(), vz_bc_pt_array.get_num(), 0.01);
	//md_disp.set_pts_from_pcl_id(tbc_pcl_array.get_mem(), tbc_pcl_array.get_num(), 0.012);
	md_disp.start();
	return;

	ResultFile_hdf5 res_file_hdf5;
	res_file_hdf5.create("t3d_me_mt_rigid_mesh_compression.h5");

	ModelData_T3D_ME_mt md;
	md.output_model(model, res_file_hdf5);

	TimeHistory_T3D_ME_mt_complete out1("compression");
	out1.set_res_file(res_file_hdf5);
	out1.set_output_init_state();
	out1.set_interval_num(50);
	TimeHistory_ConsoleProgressBar out_cpb;

	Step_T3D_ME_mt step("step1");
	step.set_model(model);
	step.set_step_time(0.5);
	//step.set_step_time(1.0e-5);
	step.set_dtime(1.0e-5);
	step.set_thread_num(4);
	step.add_time_history(out1);
	step.add_time_history(out_cpb);
	step.solve();
}

#include "QtApp_Posp_T3D_ME_mt.h"
#include "test_model_view_omp.h"

void test_t3d_me_mt_test_rigid_mesh_result(int argc, char** argv)
{
	ResultFile_hdf5 rf;
	rf.open("t3d_me_mt_rigid_mesh_compression.h5");

	//QtApp_Posp_T3D_ME_mt app(argc, argv, QtApp_Posp_T3D_ME_mt::SingleFrame);
	//app.set_res_file(rf, "compression", 0, Hdf5Field::s33);
	QtApp_Posp_T3D_ME_mt app(argc, argv, QtApp_Posp_T3D_ME_mt::Animation);
	app.set_ani_time(5.0);
	app.set_res_file(rf, "compression", Hdf5Field::s33);
	app.set_win_size(1200, 800);
	app.set_view_dir(0.0f, 10.0f);
	app.set_light_dir(0.0f, 10.0f);
	app.set_color_map_fld_range(-20.0, 0.0);
	app.set_color_map_geometry(1.2f, 0.4f, 0.45f);
	//app.set_png_name("t3d_me_mt_rigid_mesh_compression");
	//app.set_gif_name("t3d_me_mt_rigid_mesh_compression");
	app.start();
}
