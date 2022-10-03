#include "TestsParallel_pcp.h"

#include "RigidObject/RigidMeshT2D.h"
#include "RigidObject/RigidMeshT3D.h"
#include "RigidObject/RigidObjectByT2DMesh.h"
#include "DetectCollisionSAT.hpp"
#include "Model_T3D_ME_mt.h"
#include "SymMatEigen.h"

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
	//std::cout << "dist: " << dist
	//	<< ", norm: (" << norm.x
	//	<< ", " << norm.y
	//	<< ", " << norm.z << ")\n";

	const double* rb_moi = rb.get_moi();
	std::cout << rb_moi[0] << ", " << rb_moi[3] << ", " << rb_moi[5] << ",\n"
		<< rb_moi[3] << ", " << rb_moi[1] << ", " << rb_moi[4] << ",\n"
		<< rb_moi[5] << ", " << rb_moi[4] << ", " << rb_moi[2] << ",\n";
}

void test_rigid_mesh_contact_2d(int argc, char** argv)
{
	Point2D pt, lpt1, lpt2, lpt3;
	double dist_to_pt;
	Vector2D norm_to_pt;

	lpt1.x = 0.0;
	lpt1.y = 0.0;
	lpt2.x = 1.0;
	lpt2.y = 0.0;

	pt.x = 2.0;
	pt.y = 1.0;

	PointToLineDistance ptl;
	ptl.init_line(lpt1, lpt2);
	unsigned char dist_type = ptl.cal_distance_to_point(pt, dist_to_pt);
	ptl.cal_normal_to_point(pt, dist_type, norm_to_pt);
	std::cout << "type: " << int(dist_type) << "\ndist: " << dist_to_pt
			  << "\nnorm: " << norm_to_pt.x << ", " << norm_to_pt.y << ",\n";

	Rect bbox;
	lpt1.x = 0.0;
	lpt1.y = -1.0;
	lpt2.x = 2.0;
	lpt2.y = 0.0;
	DetectLineAABBCollisionSAT laabb;
	laabb.init_line(lpt1, lpt2);
	laabb.get_ln_bbox(bbox);
	std::cout << "bbox: " << bbox.xl << ", " << bbox.xu << ", "
		<< bbox.yl << ", " << bbox.yu << ",\n";
	bbox.xl = -0.5;
	bbox.xu = 1.0001;
	bbox.yl = -0.5;
	bbox.yu = 0.5;
	bool res = laabb.detect(bbox);
	std::cout << "ln lp: " << int(res) << "\n";

	lpt3.x = -1.0;
	lpt3.y = -1.5;
	DetectTriangleAABBCollisionSAT taabb;
	taabb.init_triangle(lpt1, lpt2, lpt3);
	taabb.get_tri_bbox(bbox);
	std::cout << "bbox: " << bbox.xl << ", " << bbox.xu << ", "
		<< bbox.yl << ", " << bbox.yu << ",\n";
	bbox.xl = -0.5;
	bbox.xu = 0.99;
	bbox.yl = -0.5;
	bbox.yu = 0.5;
	res = laabb.detect(bbox);
	std::cout << "tri lp: " << int(res) << "\n";

	// test rigid mesh
	TriangleMesh tri_mh;
	tri_mh.load_mesh_from_hdf5("../../Asset/rect_mesh_1by4.h5");

	RigidMeshT2D rb;
	rb.init_from_mesh(tri_mh, 0.3, 0.3);

	//rb.init_max_dist();

	double pt_r, dist;
	Vector2D norm;

	pt.x = 1.3;
	pt.y = -0.2;
	pt_r = 0.3;
	res = rb.detect_collision_with_point(pt, pt_r, dist, norm);

	std::cout << "is op: " << res << "\ndist: " << dist << "\nnorm: " << norm.x << ", " << norm.y << std::endl;

	RigidObjectByT2DMesh rmh;
	rmh.init(1.0, "../../Asset/rect_mesh_1by4.h5", 0.0, 0.0, 0.0, 0.3, 0.3);
	std::cout << "moi: " << rmh.get_moi() << std::endl;
}

void test_stress_rotate()
{
	double stress[6] = {
		10.0,
		66.5,
		-711.0,
		12.0,
		-6561.54,
		426.0
	};
	std::cout << stress[0] << ", " << stress[1] << ", " << stress[2] << ", "
		<< stress[3] << ", " << stress[4] << ", " << stress[5] << ",\n";

	double prin_s[3], prin_vecs[3][3];

	cal_sym_mat_eigen(stress, prin_s, prin_vecs);

	prin_s[0] += 1.0;
	prin_s[1] += 1.0;

	double stress2[6];
	rotate_eigen_mat_to_sym_mat(prin_s, prin_vecs, stress2);

	std::cout << stress2[0] << ", " << stress2[1] << ", " << stress2[2] << ", "
			  << stress2[3] << ", " << stress2[4] << ", " << stress2[5] << ",\n";
}
