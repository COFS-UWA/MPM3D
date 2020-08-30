#ifndef __Tetrahedron_Utils_h__
#define __Tetrahedron_Utils_h__

#include <assert.h>
#include "Geometry.h"

template <typename Point3D>
Cube get_tetrahedron_bounding_box(
	Point3D& n1, Point3D& n2, Point3D& n3, Point3D& n4)
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

template <typename Node3D, typename Point3D>
inline double cal_tetrahedron_vol(Node3D& n1, Node3D& n2, Node3D& n3, Point3D& p4)
{
	double v21_x, v21_y, v21_z;
	double v31_x, v31_y, v31_z;
	double v41_x, v41_y, v41_z;
	double cp_x, cp_y, cp_z;
	v21_x = n2.x - n1.x;
	v21_y = n2.y - n1.y;
	v21_z = n2.z - n1.z;
	v31_x = n3.x - n1.x;
	v31_y = n3.y - n1.y;
	v31_z = n3.z - n1.z;
	v41_x = p4.x - n1.x;
	v41_y = p4.y - n1.y;
	v41_z = p4.z - n1.z;
	cp_x = v21_y * v31_z - v31_y * v21_z;
	cp_y = v21_z * v31_x - v31_z * v21_x;
	cp_z = v21_x * v31_y - v31_x * v21_y;
	return (v41_x * cp_x + v41_y * cp_y + v41_z * cp_z) / 6.0;
}

template <typename Node3D>
struct PointInTetrahedron
{
protected:
	double a1, b1, c1, coef1;
	double a2, b2, c2, coef2;
	double a3, b3, c3, coef3;
	double a4, b4, c4, coef4;

public:
	inline void init_tetrahedron(Node3D& n1, Node3D& n2, Node3D& n3, Node3D& n4, double vol)
	{
		double inv_vol_6 = 1.0 / (6.0 * vol);
		// N1
		a1 = ((n4.y - n2.y) * (n3.z - n2.z) - (n3.y - n2.y) * (n4.z - n2.z)) * inv_vol_6;
		b1 = ((n4.z - n2.z) * (n3.x - n2.x) - (n3.z - n2.z) * (n4.x - n2.x)) * inv_vol_6;
		c1 = ((n4.x - n2.x) * (n3.y - n2.y) - (n3.x - n2.x) * (n4.y - n2.y)) * inv_vol_6;
		coef1 = a1 * n2.x + b1 * n2.y + c1 * n2.z;
		// N2
		a2 = ((n4.y - n3.y) * (n1.z - n3.z) - (n1.y - n3.y) * (n4.z - n3.z)) * inv_vol_6;
		b2 = ((n4.z - n3.z) * (n1.x - n3.x) - (n1.z - n3.z) * (n4.x - n3.x)) * inv_vol_6;
		c2 = ((n4.x - n3.x) * (n1.y - n3.y) - (n1.x - n3.x) * (n4.y - n3.y)) * inv_vol_6;
		coef2 = a2 * n3.x + b2 * n3.y + c2 * n3.z;
		// N3
		a3 = ((n2.y - n4.y) * (n1.z - n4.z) - (n1.y - n4.y) * (n2.z - n4.z)) * inv_vol_6;
		b3 = ((n2.z - n4.z) * (n1.x - n4.x) - (n1.z - n4.z) * (n2.x - n4.x)) * inv_vol_6;
		c3 = ((n2.x - n4.x) * (n1.y - n4.y) - (n1.x - n4.x) * (n2.y - n4.y)) * inv_vol_6;
		coef3 = a3 * n4.x + b3 * n4.y + c3 * n4.z;
		// N4
		a4 = ((n2.y - n1.y) * (n3.z - n1.z) - (n3.y - n1.y) * (n2.z - n1.z)) * inv_vol_6;
		b4 = ((n2.z - n1.z) * (n3.x - n1.x) - (n3.z - n1.z) * (n2.x - n1.x)) * inv_vol_6;
		c4 = ((n2.x - n1.x) * (n3.y - n1.y) - (n3.x - n1.x) * (n2.y - n1.y)) * inv_vol_6;
		coef4 = a4 * n1.x + b4 * n1.y + c4 * n1.z;
	}
	inline void init_tetrahedron(Node3D& n1, Node3D& n2, Node3D& n3, Node3D& n4)
	{ init_tetrahedron(n1, n2, n3, n4, cal_tetrahedron_vol<Node3D, Node3D>(n1, n2, n3, n4)); }

	inline bool is_in_tetrahedron(double x, double y, double z)
	{
		double N1v = N1(x, y, z);
		double N2v = N2(x, y, z);
		double N3v = N3(x, y, z);
		double N4v = 1.0 - N1v - N2v - N3v;
		return !(N1v < 0.0 || N1v > 1.0 || N2v < 0.0 || N2v > 1.0 ||
				 N3v < 0.0 || N3v > 1.0 || N4v < 0.0 || N4v > 1.0);
	}
	template <typename Point3D>
	inline bool is_in_tetrahedron(Point3D& p)
	{ return is_in_tetrahedron(p.x, p.y, p.z); }

	// shape functions
	inline double N1(double x, double y, double z) const noexcept
	{ return a1 * x + b1 * y + c1 * z - coef1; }
	inline double N2(double x, double y, double z) const noexcept
	{ return a2 * x + b2 * y + c2 * z - coef2; }
	inline double N3(double x, double y, double z) const noexcept
	{ return a3 * x + b3 * y + c3 * z - coef3; }
	inline double N4(double x, double y, double z) const noexcept
	{ return a4 * x + b4 * y + c4 * z - coef4; }
	inline void cal_N(double x, double y, double z,
		double &N1v, double &N2v, double &N3v, double &N4v) const noexcept
	{
		N1v = N1(x, y, z); N2v = N2(x, y, z); N3v = N3(x, y, z);
		N4v = 1.0 - N1v - N2v - N3v;
	}

	// shape function derivatives
	inline double dN1_dx() const noexcept { return a1; }
	inline double dN1_dy() const noexcept { return b1; }
	inline double dN1_dz() const noexcept { return c1; }
	inline double dN2_dx() const noexcept { return a2; }
	inline double dN2_dy() const noexcept { return b2; }
	inline double dN2_dz() const noexcept { return c2; }
	inline double dN3_dx() const noexcept { return a3; }
	inline double dN3_dy() const noexcept { return b3; }
	inline double dN3_dz() const noexcept { return c3; }
	inline double dN4_dx() const noexcept { return a4; }
	inline double dN4_dy() const noexcept { return b4; }
	inline double dN4_dz() const noexcept { return c4; }
};

template <typename Node3D>
struct TetrahedronAABBCollisionSAT
{
protected:
	double hx, hy, hz;
	Point3D n1, n2, n3, n4;
	// 4 face normal + 3 * 6 edge cross product
	Vector3D axes[22];

	inline bool is_seperating_axis(Vector3D& axis,
		Point3D& p1, Point3D& p2, Point3D& p3, Point3D& p4)
	{
#define Norm_Tol 1.0e-10
		if (axis.norm() < Norm_Tol)
			return false;
		double box_range = 0.5 * (hx * abs(axis.x) + hy * abs(axis.y) + hz * abs(axis.z));
		double p1_proj = p1.x * axis.x + p1.y * axis.y + p1.z * axis.z;
		double p2_proj = p2.x * axis.x + p2.y * axis.y + p2.z * axis.z;
		double p3_proj = p3.x * axis.x + p3.y * axis.y + p3.z * axis.z;
		double p4_proj = p4.x * axis.x + p4.y * axis.y + p4.z * axis.z;
		return ((p1_proj >  box_range && p2_proj >  box_range &&
				 p3_proj >  box_range && p4_proj >  box_range) ||
				(p1_proj < -box_range && p2_proj < -box_range &&
				 p3_proj < -box_range && p4_proj < -box_range));
#undef Norm_Tol
	}

public:
	void init_tetrahedron(Node3D& _n1, Node3D& _n2, Node3D& _n3, Node3D& _n4)
	{
		n1.x = _n1.x;
		n1.y = _n1.y;
		n1.z = _n1.z;
		n2.x = _n2.x;
		n2.y = _n2.y;
		n2.z = _n2.z;
		n3.x = _n3.x;
		n3.y = _n3.y;
		n3.z = _n3.z;
		n4.x = _n4.x;
		n4.y = _n4.y;
		n4.z = _n4.z;
		double e12_x = n1.x - n2.x;
		double e12_y = n1.y - n2.y;
		double e12_z = n1.z - n2.z;
		double e13_x = n1.x - n3.x;
		double e13_y = n1.y - n3.y;
		double e13_z = n1.z - n3.z;
		double e14_x = n1.x - n4.x;
		double e14_y = n1.y - n4.y;
		double e14_z = n1.z - n4.z;
		double e23_x = n2.x - n3.x;
		double e23_y = n2.y - n3.y;
		double e23_z = n2.z - n3.z;
		double e24_x = n2.x - n4.x;
		double e24_y = n2.y - n4.y;
		double e24_z = n2.z - n4.z;
		double e34_x = n3.x - n4.x;
		double e34_y = n3.y - n4.y;
		double e34_z = n3.z - n4.z;
		// 4 face normal
		axes[0].cross(e13_x, e13_y, e13_z, e12_x, e12_y, e12_z);
		axes[1].cross(e12_x, e12_y, e12_z, e14_x, e14_y, e14_z);
		axes[2].cross(e14_x, e14_y, e14_z, e13_x, e13_y, e13_z);
		axes[3].cross(e23_x, e23_y, e23_z, e24_x, e24_y, e24_z);
		// 3 * 6 edge cross product
		axes[4].cross(1.0, 0.0, 0.0, e12_x, e12_y, e12_z);
		axes[5].cross(0.0, 1.0, 0.0, e12_x, e12_y, e12_z);
		axes[6].cross(0.0, 0.0, 1.0, e12_x, e12_y, e12_z);
		axes[7].cross(1.0, 0.0, 0.0, e13_x, e13_y, e13_z);
		axes[8].cross(0.0, 1.0, 0.0, e13_x, e13_y, e13_z);
		axes[9].cross(0.0, 0.0, 1.0, e13_x, e13_y, e13_z);
		axes[10].cross(1.0, 0.0, 0.0, e14_x, e14_y, e14_z);
		axes[11].cross(0.0, 1.0, 0.0, e14_x, e14_y, e14_z);
		axes[12].cross(0.0, 0.0, 1.0, e14_x, e14_y, e14_z);
		axes[13].cross(1.0, 0.0, 0.0, e23_x, e23_y, e23_z);
		axes[14].cross(0.0, 1.0, 0.0, e23_x, e23_y, e23_z);
		axes[15].cross(0.0, 0.0, 1.0, e23_x, e23_y, e23_z);
		axes[16].cross(1.0, 0.0, 0.0, e24_x, e24_y, e24_z);
		axes[17].cross(0.0, 1.0, 0.0, e24_x, e24_y, e24_z);
		axes[18].cross(0.0, 0.0, 1.0, e24_x, e24_y, e24_z);
		axes[19].cross(1.0, 0.0, 0.0, e34_x, e34_y, e34_z);
		axes[20].cross(0.0, 1.0, 0.0, e34_x, e34_y, e34_z);
		axes[21].cross(0.0, 0.0, 1.0, e34_x, e34_y, e34_z);
	}

	bool detect_collision_with_aabb(Cube& aabb)
	{
		hx = aabb.xu - aabb.xl;
		hy = aabb.yu - aabb.yl;
		hz = aabb.zu - aabb.zl;
		double box_xc = (aabb.xl + aabb.xu) * 0.5;
		double box_yc = (aabb.yl + aabb.yu) * 0.5;
		double box_zc = (aabb.zl + aabb.zu) * 0.5;
		Point3D n1_m, n2_m, n3_m, n4_m; // moved tetrahedron nodes
		n1_m.x = n1.x - box_xc;
		n1_m.y = n1.y - box_yc;
		n1_m.z = n1.z - box_zc;
		n2_m.x = n2.x - box_xc;
		n2_m.y = n2.y - box_yc;
		n2_m.z = n2.z - box_zc;
		n3_m.x = n3.x - box_xc;
		n3_m.y = n3.y - box_yc;
		n3_m.z = n3.z - box_zc;
		n4_m.x = n4.x - box_xc;
		n4_m.y = n4.y - box_yc;
		n4_m.z = n4.z - box_zc;
		// if there is one seperating axis, there is no collision
		if (is_seperating_axis(axes[0], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(axes[1], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(axes[2], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(axes[3], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(axes[4], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(axes[5], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(axes[6], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(axes[7], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(axes[8], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(axes[9], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(axes[10], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(axes[11], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(axes[12], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(axes[13], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(axes[14], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(axes[15], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(axes[16], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(axes[17], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(axes[18], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(axes[19], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(axes[20], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(axes[21], n1_m, n2_m, n3_m, n4_m))
			return false;
		return true;
	}

	inline const Point3D& get_n1() const noexcept { return n1; }
	inline const Point3D& get_n2() const noexcept { return n2; }
	inline const Point3D& get_n3() const noexcept { return n3; }
	inline const Point3D& get_n4() const noexcept { return n4; }
};


template <typename Node3D>
struct TriangleAABBCollisionSAT
{
protected:
	double hx, hy, hz;
	Point3D n1, n2, n3;
	// 1 face normal + 3 * 3 edge cross product
	Vector3D axes[10];

	inline bool is_seperating_axis(Vector3D& axis,
		Point3D& p1, Point3D& p2, Point3D& p3)
	{
#define Norm_Tol 1.0e-10
		if (axis.norm() < Norm_Tol)
			return false;
		double box_range = 0.5 * (hx * abs(axis.x) + hy * abs(axis.y) + hz * abs(axis.z));
		double p1_proj = p1.x * axis.x + p1.y * axis.y + p1.z * axis.z;
		double p2_proj = p2.x * axis.x + p2.y * axis.y + p2.z * axis.z;
		double p3_proj = p3.x * axis.x + p3.y * axis.y + p3.z * axis.z;
		return ((p1_proj >  box_range && p2_proj >  box_range && p3_proj >  box_range) ||
				(p1_proj < -box_range && p2_proj < -box_range && p3_proj < -box_range));
#undef Norm_Tol
	}

public:
	void init_triangle(Node3D& _n1, Node3D& _n2, Node3D& _n3)
	{
		n1.x = _n1.x;
		n1.y = _n1.y;
		n1.z = _n1.z;
		n2.x = _n2.x;
		n2.y = _n2.y;
		n2.z = _n2.z;
		n3.x = _n3.x;
		n3.y = _n3.y;
		n3.z = _n3.z;
		double e12_x = n1.x - n2.x;
		double e12_y = n1.y - n2.y;
		double e12_z = n1.z - n2.z;
		double e13_x = n1.x - n3.x;
		double e13_y = n1.y - n3.y;
		double e13_z = n1.z - n3.z;
		double e23_x = n2.x - n3.x;
		double e23_y = n2.y - n3.y;
		double e23_z = n2.z - n3.z;
		// 1 face normal
		axes[0].cross(e13_x, e13_y, e13_z, e12_x, e12_y, e12_z);
		// 3 * 3 edge cross product
		axes[1].cross(1.0, 0.0, 0.0, e12_x, e12_y, e12_z);
		axes[2].cross(0.0, 1.0, 0.0, e12_x, e12_y, e12_z);
		axes[3].cross(0.0, 0.0, 1.0, e12_x, e12_y, e12_z);
		axes[4].cross(1.0, 0.0, 0.0, e13_x, e13_y, e13_z);
		axes[5].cross(0.0, 1.0, 0.0, e13_x, e13_y, e13_z);
		axes[6].cross(0.0, 0.0, 1.0, e13_x, e13_y, e13_z);
		axes[7].cross(1.0, 0.0, 0.0, e23_x, e23_y, e23_z);
		axes[8].cross(0.0, 1.0, 0.0, e23_x, e23_y, e23_z);
		axes[9].cross(0.0, 0.0, 1.0, e23_x, e23_y, e23_z);
	}

	bool detect_collision_with_aabb(Cube& aabb)
	{
		hx = aabb.xu - aabb.xl;
		hy = aabb.yu - aabb.yl;
		hz = aabb.zu - aabb.zl;
		double box_xc = (aabb.xl + aabb.xu) * 0.5;
		double box_yc = (aabb.yl + aabb.yu) * 0.5;
		double box_zc = (aabb.zl + aabb.zu) * 0.5;
		Point3D n1_m, n2_m, n3_m;
		n1_m.x = n1.x - box_xc;
		n1_m.y = n1.y - box_yc;
		n1_m.z = n1.z - box_zc;
		n2_m.x = n2.x - box_xc;
		n2_m.y = n2.y - box_yc;
		n2_m.z = n2.z - box_zc;
		n3_m.x = n3.x - box_xc;
		n3_m.y = n3.y - box_yc;
		n3_m.z = n3.z - box_zc;
		// if there is one seperating axis, there is no collision
		if (is_seperating_axis(axes[0], n1_m, n2_m, n3_m) ||
			is_seperating_axis(axes[1], n1_m, n2_m, n3_m) ||
			is_seperating_axis(axes[2], n1_m, n2_m, n3_m) ||
			is_seperating_axis(axes[3], n1_m, n2_m, n3_m) ||
			is_seperating_axis(axes[4], n1_m, n2_m, n3_m) ||
			is_seperating_axis(axes[5], n1_m, n2_m, n3_m) ||
			is_seperating_axis(axes[6], n1_m, n2_m, n3_m) ||
			is_seperating_axis(axes[7], n1_m, n2_m, n3_m) ||
			is_seperating_axis(axes[8], n1_m, n2_m, n3_m) ||
			is_seperating_axis(axes[9], n1_m, n2_m, n3_m))
			return false;
		return true;
	}

	inline const Point3D& get_n1() const noexcept { return n1; }
	inline const Point3D& get_n2() const noexcept { return n2; }
	inline const Point3D& get_n3() const noexcept { return n3; }
};

template <typename Node3D>
struct PointToTriangleDistance
{
protected:
	Point3D n1, n2, n3;
	union
	{
		struct { Vector3D ix1, iy1, iz1; };
		double T1[3][3];
	};
	union
	{
		struct { Vector3D ix2, iy2, iz2; };
		double T2[3][3];
	};
	union
	{
		struct { Vector3D ix3, iy3, iz3; };
		double T3[3][3];
	};
	double a1, a2, a3;

public:
	void init_triangle(Node3D& _n1, Node3D& _n2, Node3D& _n3)
	{
		n1.x = _n1.x;
		n1.y = _n1.y;
		n1.z = _n1.z;
		n2.x = _n2.x;
		n2.y = _n2.y;
		n2.z = _n2.z;
		n3.x = _n3.x;
		n3.y = _n3.y;
		n3.z = _n3.z;
		Vector3D tmp1, tmp2;
		double coef_tmp;
		// line 1 at x axis
		ix1.substract<Point3D>(n2, n1);
		a1 = ix1.norm();
		ix1.scale(1.0 / a1);
		tmp1.substract<Point3D>(n3, n1);
		coef_tmp = tmp1.dot(ix1);
		tmp2.x = coef_tmp * ix1.x;
		tmp2.y = coef_tmp * ix1.y;
		tmp2.z = coef_tmp * ix1.z;
		iy1.substract<Vector3D>(tmp1, tmp2);
		iy1.normalize();
		iz1.cross<Vector3D>(ix1, iy1);
		// line 2 at x axis
		ix2.substract<Point3D>(n3, n2);
		a2 = ix2.norm();
		ix2.scale(1.0 / a2);
		tmp1.substract<Point3D>(n1, n2);
		coef_tmp = tmp1.dot(ix2);
		tmp2.x = coef_tmp * ix2.x;
		tmp2.y = coef_tmp * ix2.y;
		tmp2.z = coef_tmp * ix2.z;
		iy2.substract<Vector3D>(tmp1, tmp2);
		iy2.normalize();
		iz2.cross<Vector3D>(ix2, iy2);
		// line 3 at x axis
		ix3.substract<Point3D>(n1, n3);
		a3 = ix3.norm();
		ix3.scale(1.0 / a3);
		tmp1.substract<Point3D>(n2, n3);
		coef_tmp = tmp1.dot(ix3);
		tmp2.x = coef_tmp * ix3.x;
		tmp2.y = coef_tmp * ix3.y;
		tmp2.z = coef_tmp * ix3.z;
		iy3.substract<Vector3D>(tmp1, tmp2);
		iy3.normalize();
		iz3.cross<Vector3D>(ix3, iy3);
	}

	unsigned char cal_distance_to_point(Point3D& p, double& dist)
	{
		Vector3D pj1, pj2, pj3, tmp1;
		tmp1.substract<Point3D>(p, n1);
		pj1.x = T1[0][0] * tmp1.x + T1[0][1] * tmp1.y + T1[0][2] * tmp1.z;
		pj1.y = T1[1][0] * tmp1.x + T1[1][1] * tmp1.y + T1[1][2] * tmp1.z;
		pj1.z = T1[2][0] * tmp1.x + T1[2][1] * tmp1.y + T1[2][2] * tmp1.z;
		tmp1.substract<Point3D>(p, n2);
		pj2.x = T2[0][0] * tmp1.x + T2[0][1] * tmp1.y + T2[0][2] * tmp1.z;
		pj2.y = T2[1][0] * tmp1.x + T2[1][1] * tmp1.y + T2[1][2] * tmp1.z;
		tmp1.substract<Point3D>(p, n3);
		pj3.x = T3[0][0] * tmp1.x + T3[0][1] * tmp1.y + T3[0][2] * tmp1.z;
		pj3.y = T3[1][0] * tmp1.x + T3[1][1] * tmp1.y + T3[1][2] * tmp1.z;

		if (pj1.y >= 0.0 && pj2.y >= 0.0 && pj3.y >= 0.0)
		{
			dist = pj1.z;
			return 0;
		}

		if (pj1.y < 0.0 && pj1.x >= 0.0 && pj1.x <= a1)
		{
			dist = sqrt(pj1.y * pj1.y + pj1.z * pj1.z);
			if (pj1.z < 0)
				dist = -dist;
			return 1;
		}

		if (pj2.y < 0.0 && pj2.x >= 0.0 && pj2.x <= a2)
		{
			dist = sqrt(pj2.y * pj2.y + pj1.z * pj1.z);
			if (pj1.z < 0)
				dist = -dist;
			return 2;
		}

		if (pj3.y < 0.0 && pj3.x >= 0.0 && pj3.x <= a3)
		{
			dist = sqrt(pj3.y * pj3.y + pj1.z * pj1.z);
			if (pj1.z < 0)
				dist = -dist;
			return 3;
		}

		//double dist1;
		if (pj1.y < 0.0 && pj1.x > a1 || pj2.y < 0.0 && pj2.x < 0.0)
		{
			//dist1 = sqrt((pj1.x - a1) * (pj1.x - a1) + pj1.y * pj1.y + pj1.z * pj1.z);
			dist = sqrt(pj2.x * pj2.x + pj2.y * pj2.y + pj1.z * pj1.z);
			if (pj1.z < 0)
				dist = -dist;
			return 4;
		}

		if (pj2.y < 0.0 && pj2.x > a2 || pj3.y < 0.0 && pj3.x < 0.0)
		{
			//dist1 = sqrt((pj2.x - a2) * (pj2.x - a2) + pj2.y * pj2.y + pj1.z * pj1.z);
			dist = sqrt(pj3.x * pj3.x + pj3.y * pj3.y + pj1.z * pj1.z);
			if (pj1.z < 0)
				dist = -dist;
			return 5;
		}

		if (pj3.y < 0.0 && pj3.x > a3 || pj1.y < 0.0 && pj1.x < 0.0)
		{
			//dist1 = sqrt((pj3.x - a3) * (pj3.x - a3) + pj3.y * pj3.y + pj1.z * pj1.z);
			dist = sqrt(pj1.x * pj1.x + pj1.y * pj1.y + pj1.z * pj1.z);
			if (pj1.z < 0)
				dist = -dist;
			return 6;
		}

		assert(0);
		return 7;
	}

	void cal_normal_to_point(Point3D& pt, unsigned char norm_type, Vector3D& normal)
	{
		Vector3D v1, v2, tmp1;
		double coef;
		switch (norm_type)
		{
		case 0:
			v1.substract(n2, n1);
			v2.substract(n3, n1);
			normal.cross(v1, v2);
			normal.normalize();
			tmp1.substract<Point3D>(pt, n1);
			coef = T1[2][0] * tmp1.x + T1[2][1] * tmp1.y + T1[2][2] * tmp1.z;
			if (coef < 0.0)
				normal.reverse();
			return;
		case 1:
			v1.substract(n2, n1);
			v1.normalize();
			v2.substract(pt, n1);
			coef = v1.dot(v2);
			normal.substract(v2, v1.scale(coef));
			break;
		case 2:
			v1.substract(n3, n2);
			v1.normalize();
			v2.substract(pt, n2);
			coef = v1.dot(v2);
			normal.substract(v2, v1.scale(coef));
			break;
		case 3:
			v1.substract(n1, n3);
			v1.normalize();
			v2.substract(pt, n3);
			coef = v1.dot(v2);
			v1.scale(coef);
			normal.substract(v2, v1);
			break;
		case 4:
			normal.substract(pt, n2);
			break;
		case 5:
			normal.substract(pt, n3);
			break;
		case 6:
			normal.substract(pt, n1);
			break;
		default:
			assert(0);
			return;
		}
		double norm = normal.norm();
		if (norm != 0.0)
			normal.scale(1.0 / norm);
		else
		{
			v1.substract(n2, n1);
			v2.substract(n3, n1);
			normal.cross(v1, v2);
			normal.normalize();
		}
		return;
	}
};

#endif