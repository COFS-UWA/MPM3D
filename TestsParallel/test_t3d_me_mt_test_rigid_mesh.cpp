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
			//size_t z_id = 2;
			//for (size_t y_id = 0; y_id < grid.y_num; ++y_id)
			//{
			//	for (size_t x_id = 0; x_id < grid.x_num; ++x_id)
			//	{
			//		GridPosType& g = grid_pos_type[offset_from_xyz_id(x_id, y_id, z_id)];
			//		std::cout << size_t(g);
			//	}
			//	std::cout << "\n";
			//}
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
	rb.init(1.0, "../../Asset/brick_mesh_1.00_1x1x1.h5",
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.2, 0.2);
	rb.init_max_dist(0.1);
	
	//bool res = rb.detect_collision_with_point(0.5, 0.3, 0.05);
}

void test_t3d_me_mt_test_rigid_mesh(int argc, char** argv)
{
	TetrahedronMesh teh_mesh;
	teh_mesh.load_mesh_from_hdf5("../../Asset/brick_mesh_1.00_2x2x10.h5");
	teh_mesh.init_search_grid(0.05, 0.05, 0.05);

	Model_T3D_ME_mt model;
	model.init_mesh(teh_mesh);
	model.init_search_grid(teh_mesh);

	ParticleGenerator3D<TetrahedronMesh> pcl_generator;
	pcl_generator.generate_pcls_grid(Cube(0.0, 0.2, 0.0, 0.2, 0.0, 1.0), 0.025, 0.025, 0.025);
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

	model.init_t3d_rigid_mesh(1.0, "../../Asset/spudcan_model.h5",
		0.0, 0.0, 0.0, 90.0, 0.0, 0.0, 0.05, 0.05, 0.05);
	model.set_t3d_rigid_mesh_velocity(0.0, 0.0, 0.0);

	IndexArray tbc_pcl_array(100);
	find_3d_pcls(model, tbc_pcl_array, Cube(0.0, 0.2, 0.0, 0.2, 1.0 - 0.013, 1.0));
	MemoryUtils::ItemArray<double> tzs_mem(tbc_pcl_array.get_num());
	double tz_mag = 0.025 * 0.025 * -10.0;
	for (size_t t_id = 0; t_id < tbc_pcl_array.get_num(); ++t_id)
		tzs_mem.add(tz_mag);
	model.init_tzs(tbc_pcl_array.get_num(), tbc_pcl_array.get_mem(), tzs_mem.get_mem());

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
}
