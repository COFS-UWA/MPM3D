#ifndef __Tetrahedron_Utils_h__
#define __Tetrahedron_Utils_h__

#include "Geometry.h"

template <typename Point3D>
Cube get_tetrahedron_bounding_box(
	Point3D &n1, Point3D &n2, Point3D &n3, Point3D &n4)
{
	Cube res;
	res.xl = n1.x;
	if (res.xl > n2.x)
		res.xl = n2.x;
	if (res.xl > n3.x)
		res.xl = n3.x;
	if (res.xl > n4.x)
		res.xl = n4.x;
	res.xu = n1.x;
	if (res.xu < n2.x)
		res.xu = n2.x;
	if (res.xu < n3.x)
		res.xu = n3.x;
	if (res.xu < n4.x)
		res.xu = n4.x;
	res.yl = n1.y;
	if (res.yl > n2.y)
		res.yl = n2.y;
	if (res.yl > n3.y)
		res.yl = n3.y;
	if (res.yl > n4.y)
		res.yl = n4.y;
	res.yu = n1.y;
	if (res.yu < n2.y)
		res.yu = n2.y;
	if (res.yu < n3.y)
		res.yu = n3.y;
	if (res.yu < n4.y)
		res.yu = n4.y;
	res.zl = n1.z;
	if (res.zl > n2.z)
		res.zl = n2.z;
	if (res.zl > n3.z)
		res.zl = n3.z;
	if (res.zl > n4.z)
		res.zl = n4.z;
	res.zu = n1.z;
	if (res.zu < n2.z)
		res.zu = n2.z;
	if (res.zu < n3.z)
		res.zu = n3.z;
	if (res.zu < n4.z)
		res.zu = n4.z;
	return res;
}

template <typename Point3D>
Cube get_3Dtriangle_bounding_box(
	Point3D& n1, Point3D& n2, Point3D& n3)
{
	Cube res;
	res.xl = n1.x;
	if (res.xl > n2.x)
		res.xl = n2.x;
	if (res.xl > n3.x)
		res.xl = n3.x;
	res.xu = n1.x;
	if (res.xu < n2.x)
		res.xu = n2.x;
	if (res.xu < n3.x)
		res.xu = n3.x;
	res.yl = n1.y;
	if (res.yl > n2.y)
		res.yl = n2.y;
	if (res.yl > n3.y)
		res.yl = n3.y;
	res.yu = n1.y;
	if (res.yu < n2.y)
		res.yu = n2.y;
	if (res.yu < n3.y)
		res.yu = n3.y;
	res.zl = n1.z;
	if (res.zl > n2.z)
		res.zl = n2.z;
	if (res.zl > n3.z)
		res.zl = n3.z;
	res.zu = n1.z;
	if (res.zu < n2.z)
		res.zu = n2.z;
	if (res.zu < n3.z)
		res.zu = n3.z;
	return res;
}

template <typename Point3D>
bool test_3Dtriangle_aabb_intersection(
	Point3D& n1, Point3D& n2, Point3D& n3,
	Cube &aabb)
{

	return true;
}

template <typename Point3D>
bool test_tetrahedron_aabb_intersection(
	Point3D& n1, Point3D& n2, Point3D& n3, Point3D& n4,
	Cube &aabb)
{

	return true;
}

template <typename Point3D>
bool test_3Dtriangle_aabb_intersection(
	Point3D& n1, Point3D& n2, Point3D& n3,
	Cube&& aabb)
{

	return true;
}

#endif