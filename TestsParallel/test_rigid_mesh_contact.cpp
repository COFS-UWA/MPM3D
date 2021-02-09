#include "TestsParallel_pcp.h"

#include "RigidObject/RigidMeshT3D.h"
#include "Model_T3D_ME_mt.h"

#include "test_simulations_omp.h"

void test_rigid_mesh_contact(int argc, char** argv)
{
	TetrahedronMesh mh;
	mh.load_mesh_from_hdf5("../../Asset/brick_mesh_1.00_2x2x2.h5");

	RigidMeshT3D rm;
	rm.init_from_mesh(mh, 0.5, 0.5, 0.5);
	std::cout << rm.get_grid_x_num() << ", "
		<< rm.get_grid_y_num() << ", "
		<< rm.get_grid_z_num() << "\n";
	const size_t fig_num = rm.get_face_in_grid_list_len();
	const size_t *fig_list = rm.get_face_in_grid_list();
	for (size_t i = 0; i < fig_num; i++)
		std::cout << fig_list[i] << ", ";

	//rm.init_max_dist(1.0);
	Point3D pt;
	double pt_r, dist;
	Vector3D norm;
	pt.x = 0.5;
	pt.y = 0.5;
	pt.z = -0.25;
	pt_r = 0.3;
	rm.detect_collision_with_point(pt, pt_r, dist, norm);
	std::cout << "dist: " << dist
		<< ", norm: (" << norm.x
		<< ", " << norm.y
		<< ", " << norm.z << ")\n";
}

void test_rigid_mesh_contact2(int argc, char** argv)
{
	Model_T3D_ME_mt md;
	//md.init_t3d_rigid_mesh(1.0, "../../Asset/brick_mesh_1.00_2x2x2.h5",
	//	0.0, 0.0, 0.0, 0.0, -90.0, 0.0, 0.3, 0.3, 0.3);
	md.init_t3d_rigid_mesh(1.0, "../../Asset/brick_mesh_1.00_2x2x2.h5",
		0.0, 0.0, 0.0, 90.0, 0.0, 0.0, 0.3, 0.3, 0.3);

	auto& rb = md.get_t3d_rigid_mesh();
	
	Point3D pt;
	double pt_r, dist;
	Vector3D norm;
	Point3D ct_pos;

	//pt.x = -0.5;
	//pt.y = 0.5;
	//pt.z = -0.25;
	//pt_r = 0.3;
	//rb.detect_collision_with_point(
	//	pt.x, pt.y, pt.z, pt_r, dist, norm, ct_pos);
	//std::cout << "dist: " << dist
	//	<< ", norm: (" << norm.x
	//	<< ", " << norm.y
	//	<< ", " << norm.z << ")\n";

	pt.x = 0.5;
	pt.y = -0.5;
	pt.z = -0.25;
	pt_r = 0.3;
	rb.detect_collision_with_point(
		pt.x, pt.y, pt.z, pt_r, dist, norm, ct_pos);
	std::cout << "dist: " << dist
		<< ", norm: (" << norm.x
		<< ", " << norm.y
		<< ", " << norm.z << ")\n";

}
