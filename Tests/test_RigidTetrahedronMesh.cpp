#include "Tests_pcp.h"

#include <iostream>
#include "RigidBody/RigidTetrahedronMesh.h"

#include "test_simulations.h"

namespace
{
	class TestRigidTetrahedronMesh : public RigidTetrahedronMesh
	{
	public:
		TestRigidTetrahedronMesh() {}
		~TestRigidTetrahedronMesh() {}

		inline void init_teh_aabb_collision(Element& e)
		{ return RigidTetrahedronMesh::init_teh_aabb_collision(e); }
		inline bool detect_teh_aabb_collision(Cube& box)
		{ return RigidTetrahedronMesh::detect_teh_aabb_collision(box); }
	
		inline void init_tri_aabb_collision(Face& f)
		{ return RigidTetrahedronMesh::init_tri_aabb_collision(f); }
		inline bool detect_tri_aabb_collision(Cube& box)
		{ return RigidTetrahedronMesh::detect_tri_aabb_collision(box); }
		
		inline Grid& grid_by_id(size_t x_id, size_t y_id, size_t z_id)
		{ return RigidTetrahedronMesh::grid_by_id(x_id, y_id, z_id); }
	
		inline const Cube& get_grid_bbox() { return RigidTetrahedronMesh::get_grid_bbox(); }
	};
}

void test_RigidTetrahedronMesh_intersection(int argc, char** argv)
{
	TestRigidTetrahedronMesh rb;
	rb.init_mesh("../../Asset/teh_mesh.h5", 0.0, 0.0, 0.0);
	std::cout << "node num: " << rb.get_node_num()
			<< ", elem num: " << rb.get_elem_num()
			<< ", face num: " << rb.get_bface_num() << "\n";
	
	bool res;
	Cube test_box;
	rb.init_teh_aabb_collision(rb.get_elems()[0]);
	
	//test_box.xl = -3.0;
	//test_box.xu = -2.0;
	//test_box.yl = -1.0;
	//test_box.yu = 1.0;
	//test_box.zl = 0.0;
	//test_box.zu = 1.0;
	//res = rb.detect_teh_aabb_collision(test_box);

	//test_box.xl = -2.0;
	//test_box.xu = -1.0;
	//test_box.yl = -1.5;
	//test_box.yu = -0.5;
	//test_box.zl = -0.5;
	//test_box.zu = 0.5;
	//res = rb.detect_teh_aabb_collision(test_box);

	//test_box.xl = 1.5;
	//test_box.xu = 2.0;
	//test_box.yl = -1.5;
	//test_box.yu = -1.0;
	//test_box.zl = -0.5;
	//test_box.zu = 0.5;
	//res = rb.detect_teh_aabb_collision(test_box);

	//test_box.xl = 0.25;
	//test_box.xu = 0.75;
	//test_box.yl = -0.8;
	//test_box.yu = 0.2;
	//test_box.zl = 0.1;
	//test_box.zu = 0.5;
	//res = rb.detect_teh_aabb_collision(test_box);
	
	//test_box.xl = -0.5;
	//test_box.xu = 0.5;
	//test_box.yl = -1.1;
	//test_box.yu = 1.5;
	//test_box.zl = -0.1;
	//test_box.zu = 2.0;
	//res = rb.detect_teh_aabb_collision(test_box);
	
	rb.init_tri_aabb_collision(rb.get_bfaces()[3]);

	//test_box.xl = -3.0;
	//test_box.xu = -2.0;
	//test_box.yl = -2.0;
	//test_box.yu = -1.0;
	//test_box.zl = -0.5;
	//test_box.zu = 0.5;
	//res = rb.detect_tri_aabb_collision(test_box);

	//test_box.xl = -2.5;
	//test_box.xu = -1.5;
	//test_box.yl = -1.5;
	//test_box.yu = -0.5;
	//test_box.zl = -0.5;
	//test_box.zu = 0.5;
	//res = rb.detect_tri_aabb_collision(test_box);

	//test_box.xl = 1.5;
	//test_box.xu = 2.5;
	//test_box.yl = -1.5;
	//test_box.yu = -0.5;
	//test_box.zl = -0.5;
	//test_box.zu = 0.5;
	//res = rb.detect_tri_aabb_collision(test_box);
	
	test_box.xl = -0.5;
	test_box.xu = 0.5;
	test_box.yl = -1.5;
	test_box.yu = -0.5;
	test_box.zl = -0.5;
	test_box.zu = 0.5;
	res = rb.detect_tri_aabb_collision(test_box);

	int efe = 0;
}

void test_RigidTetrahedronMesh_bg_grid(int argc, char** argv)
{
	TestRigidTetrahedronMesh rb;
	rb.init_mesh("../../Asset/square_cap_mesh.h5", 0.0, 0.0, 0.0);
	std::cout << "node num: " << rb.get_node_num()
			  << ", elem num: " << rb.get_elem_num()
			  << ", face num: " << rb.get_bface_num() << "\n";

	int res = rb.init_bg_grids(0.05, 0.07); // 0.06, 0.075
	const Cube& g_box = rb.get_grid_bbox();
	size_t g_x_num = rb.get_grid_x_num();
	size_t g_y_num = rb.get_grid_y_num();
	size_t g_z_num = rb.get_grid_z_num();
	std::cout << "bg grid: " << g_x_num
					 << ", " << g_y_num
					 << ", " << g_z_num << "\n";

	size_t disp_z_id = 2;
	for (size_t y_id = 0; y_id < g_y_num; ++y_id)
	{
		for (size_t x_id = 0; x_id < g_x_num; ++x_id)
		{
			auto& g = rb.grid_by_id(x_id, y_id, disp_z_id);
			switch (g.pos_type)
			{
			case TestRigidTetrahedronMesh::PosType::Inside:
				std::cout << "I";
				break;
			case TestRigidTetrahedronMesh::PosType::AtBoundary:
				std::cout << "A";
				break;
			case TestRigidTetrahedronMesh::PosType::Outside:
				std::cout << "O";
				break;
			default:
				break;
			}
		}
		std::cout << "\n";
	}

	int efe = 2;
}

// close to boundary
void test_RigidTetrahedronMesh_close_to_boundary(int argc, char** argv)
{
	TestRigidTetrahedronMesh rb;
	rb.init_mesh("../../Asset/brick_mesh_0.50_5x5x5.h5", 0.0, 0.0, 0.0);
	std::cout << "node num: " << rb.get_node_num()
			<< ", elem num: " << rb.get_elem_num()
			<< ", face num: " << rb.get_bface_num() << "\n";

	int res = rb.init_bg_grids(0.05, 0.220); // 0.06, 0.075
	const Cube& g_box = rb.get_grid_bbox();
	size_t g_x_num = rb.get_grid_x_num();
	size_t g_y_num = rb.get_grid_y_num();
	size_t g_z_num = rb.get_grid_z_num();
	std::cout << "bg grid: " << g_x_num
			<< ", " << g_y_num
			<< ", " << g_z_num << "\n";

	rb.set_dist_max(0.05);

	size_t disp_z_id = 7;
	for (size_t y_id = 0; y_id < g_y_num; ++y_id)
	{
		for (size_t x_id = 0; x_id < g_x_num; ++x_id)
		{
			auto& g = rb.grid_by_id(x_id, y_id, disp_z_id);
			if (g.close_to_boundary)
				std::cout << "O";
			else
				std::cout << "X";
		}
		std::cout << "\n";
	}

	int efe = 2;
}

void test_RigidTetrahedronMesh_search_dist(int argc, char** argv)
{
	TestRigidTetrahedronMesh rb;
	rb.init_mesh("../../Asset/brick_mesh_0.50_5x5x5.h5", 0.0, 0.0, 1.0);
	std::cout << "node num: " << rb.get_node_num()
			<< ", elem num: " << rb.get_elem_num()
			<< ", face num: " << rb.get_bface_num() << "\n";

	int res = rb.init_bg_grids(0.05, 0.07); // 0.06, 0.075

	const Cube& g_box = rb.get_grid_bbox();
	size_t g_x_num = rb.get_grid_x_num();
	size_t g_y_num = rb.get_grid_y_num();
	size_t g_z_num = rb.get_grid_z_num();
	std::cout << "bg grid: "
		<< g_box.xl << ", " << g_box.xu << ", "
		<< g_box.yl << ", " << g_box.yu << ", "
		<< g_box.zl << ", " << g_box.zu << ", "
		<< g_x_num  << ", " << g_y_num << ", " << g_z_num << "\n";
	
	//auto& g = rb.grid_by_id(1, 1, 1);
	//for (auto iter = g.bfaces; iter; iter = iter->next)
	//{
	//	auto& f = *(iter->pface);
	//	std::cout << f.id << ", ";
	//}
	//std::cout << "\n";

	//size_t disp_z_id = 1;
	//for (size_t y_id = 0; y_id < g_y_num; ++y_id)
	//{
	//	for (size_t x_id = 0; x_id < g_x_num; ++x_id)
	//	{
	//		auto& g = rb.grid_by_id(x_id, y_id, disp_z_id);
	//		switch (g.pos_type)
	//		{
	//		case TestRigidTetrahedronMesh::PosType::Inside:
	//			std::cout << "I";
	//			break;
	//		case TestRigidTetrahedronMesh::PosType::AtBoundary:
	//			std::cout << "A";
	//			break;
	//		case TestRigidTetrahedronMesh::PosType::Outside:
	//			std::cout << "O";
	//			break;
	//		default:
	//			break;
	//		}
	//	}
	//	std::cout << "\n";
	//}

	double dist, nx, ny, nz;
	bool res2;

	//rb.set_dist_max(0.0); // 0 grid
	//Point3D pt(0.01, 0.01, -0.01);
	//res2 = rb.cal_distance_to_boundary(pt, dist, nx, ny, nz);
	//std::cout << "pt: (" << pt.x << ", " << pt.y << ", " << pt.z
	//		  << ")\ndist: " << dist << "\nnormal: (" << nx << ", " << ny << ", " << nz << ")\n";

	//rb.set_dist_max(0.045); // 1 grid
	//Point3D pt(0.01, 0.01, -0.01);
	//res2 = rb.cal_distance_to_boundary(pt, dist, nx, ny, nz);
	//std::cout << "pt: (" << pt.x << ", " << pt.y << ", " << pt.z
	//	<< ")\ndist: " << dist << "\nnormal: (" << nx << ", " << ny << ", " << nz << ")\n";

	//rb.set_dist_max(0.075); // 2 grid
	//Point3D pt(0.01, 0.01, -0.01);
	//res2 = rb.cal_distance_to_boundary(pt, dist, nx, ny, nz);
	//std::cout << "pt: (" << pt.x << ", " << pt.y << ", " << pt.z
	//	<< ")\ndist: " << dist << "\nnormal: (" << nx << ", " << ny << ", " << nz << ")\n";

	rb.set_dist_max(0.075); // 2 grid
	Point3D pt(0.01, 0.01, 1.0 - 0.01);
	res2 = rb.cal_dist_and_dir_to_pt(pt, dist, nx, ny, nz);
	std::cout << "pt: (" << pt.x << ", " << pt.y << ", " << pt.z
		<< ")\ndist: " << dist << "\nnormal: (" << nx << ", " << ny << ", " << nz << ")\n";

	int efe = 2;
}