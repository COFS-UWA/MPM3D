#include "TestsParallel_pcp.h"

#include "test_simulations_omp.h"

#include "RigidObject/RigidCylinder.h"

void test_rigid_cylinder(int argc, char** argv)
{
	RigidCylinder rc;
	rc.init(0.0, 0.0, 0.0, 0.5, 1.0);
	
	double h = rc.get_h();
	double r = rc.get_r();
	const Point3D &cen = rc.get_centre();

	rc.set_vbc(0.0, 1.0, 2.0);
	const Vector3D &rc_v = rc.get_velocity();

	rc.reset_cont_force();
	const ContactForce3D &cf = rc.get_cont_force();
	
	bool res; double dist; Vector3D norm; Point3D cont_pos;

	// inside
	// r
	res = rc.detect_collision_with_point(
		0.7, 0.7, 0.0, 0.0,
		dist, norm, cont_pos);
	// up
	res = rc.detect_collision_with_point(
		0.6, 0.6, 0.23, 0.0,
		dist, norm, cont_pos);
	// down
	res = rc.detect_collision_with_point(
		0.6, 0.6, -0.23, 0.0,
		dist, norm, cont_pos);

	// up
	res = rc.detect_collision_with_point(
		0.6, 0.6, 0.26, 0.05,
		dist, norm, cont_pos);
	// up, r
	res = rc.detect_collision_with_point(
		0.72, -0.72, 0.26, 0.05,
		dist, norm, cont_pos);

	// down
	res = rc.detect_collision_with_point(
		0.6, 0.6, -0.26, 0.05,
		dist, norm, cont_pos);
	// down, r
	res = rc.detect_collision_with_point(
		-0.72, 0.72, -0.26, 0.05,
		dist, norm, cont_pos);

	res = true;
}

#include "RigidObject/RigidCone.h"

void test_rigid_cone(int argc, char** argv)
{
	RigidCone rc;
	rc.init(0.0, 0.0, 0.0, 0.5, 1.0, 2.0);

	double r = rc.get_r();
	double h_tip = rc.get_h_tip();
	double h_shaft = rc.get_h_shaft();
	const Point3D& cen = rc.get_centre();

	rc.set_vbc(0.0, 1.0, 2.0);
	const Vector3D& rc_v = rc.get_velocity();

	rc.reset_cont_force();
	const ContactForce3D& cf = rc.get_cont_force();

	bool res; double dist; Vector3D norm; Point3D cont_pos;

	// outside
	// 1
	res = rc.detect_collision_with_point(
		0.355, 0.355, 2.01, 0.05,
		dist, norm, cont_pos);
	// 2
	res = rc.detect_collision_with_point(
		0.3, 0.0, 2.03, 0.05,
		dist, norm, cont_pos);
	// 3
	res = rc.detect_collision_with_point(
		0.52, 0.0, 1.0, 0.05,
		dist, norm, cont_pos);
	// 4
	res = rc.detect_collision_with_point(
		0.36, 0.36, -0.004, 0.05,
		dist, norm, cont_pos);
	// 5
	res = rc.detect_collision_with_point(
		0.19, 0.19, -0.5, 0.05,
		dist, norm, cont_pos);
	// 6
	res = rc.detect_collision_with_point(
		0.01, 0.01, -1.02, 0.05,
		dist, norm, cont_pos);
	// 7 at tip
	res = rc.detect_collision_with_point(
		0.0, 0.0, -1.0, 0.05,
		dist, norm, cont_pos);

	// inside
	// 1
	res = rc.detect_collision_with_point(
		0.3, 0.3, 1.98, 0.0,
		dist, norm, cont_pos);
	// 2
	res = rc.detect_collision_with_point(
		0.35, 0.35, 1.0, 0.0,
		dist, norm, cont_pos);
	// 3
	res = rc.detect_collision_with_point(
		0.15, 0.15, -0.38, 0.0,
		dist, norm, cont_pos);

	res = true;
}
