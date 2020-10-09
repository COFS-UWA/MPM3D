#ifndef __Triangle_Utils_h__
#define __Triangle_Utils_h__

#include <assert.h>
#include "Geometry2D.h"

template <typename Point2D>
inline void get_triangle_bounding_box(
	const Point2D& n1,
	const Point2D& n2,
	const Point2D& n3,
	Rect &res
	) noexcept
{
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
}

template <typename Node2D, typename Point2D>
inline double cal_triangle_area(
	const Node2D& p1,
	const Node2D& p2,
	const Point2D& p3
	) noexcept
{ return 0.5 * ((p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y)); }

struct PointInTriangle
{
protected:
	double a1, b1, coef1;
	double a2, b2, coef2;
	double a3, b3, coef3;

public:
	template <typename Node2D>
	inline void init_triangle(
		const Node2D& n1,
		const Node2D& n2,
		const Node2D& n3,
		double area
		) noexcept
	{
		double area2 = area * 2.0;
		a1 = (n2.y - n3.y) / area2;
		b1 = (n3.x - n2.x) / area2;
		coef1 = (n2.x * n3.y - n3.x * n2.y) / area2;
		a2 = (n3.y - n1.y) / area2;
		b2 = (n1.x - n3.x) / area2;
		coef2 = (n3.x * n1.y - n1.x * n3.y) / area2;
		a3 = (n1.y - n2.y) / area2;
		b3 = (n2.x - n1.x) / area2;
		coef3 = (n1.x * n2.y - n2.x * n1.y) / area2;
	}
	template <typename Node2D>
	inline void init_triangle(
		const Node2D& n1,
		const Node2D& n2,
		const Node2D& n3
		) noexcept
	{ init_triangle(n1, n2, n3, cal_triangle_area<Node2D, Node2D>(n1, n2, n3)); }

	inline bool is_in_triangle(double x, double y) const noexcept
	{
		double N1v = N1(x, y);
		double N2v = N2(x, y);
		double N3v = 1.0 - N1v - N2v;
		return !(N1v < 0.0 || N1v > 1.0 ||
				 N2v < 0.0 || N2v > 1.0 ||
				 N3v < 0.0 || N3v > 1.0);
	}
	template <typename Point2D>
	inline bool is_in_triangle(const Point2D& p) const noexcept
	{ return is_in_triangle(p.x, p.y); }

	// shape functions
	inline double N1(double x, double y) const noexcept
	{ return a1 * x + b1 * y - coef1; }
	inline double N2(double x, double y) const noexcept
	{ return a2 * x + b2 * y - coef2; }
	inline double N3(double x, double y) const noexcept
	{ return a3 * x + b3 * y - coef3; }
	inline void cal_N(double x, double y,
		double &N1v, double &N2v, double &N3v) const noexcept
	{ N1v = N1(x, y); N2v = N2(x, y); N3v = 1.0 - N1v - N2v; }

	// shape function derivatives
	inline double dN1_dx() const noexcept { return a1; }
	inline double dN1_dy() const noexcept { return b1; }
	inline double dN2_dx() const noexcept { return a2; }
	inline double dN2_dy() const noexcept { return b2; }
	inline double dN3_dx() const noexcept { return a3; }
	inline double dN3_dy() const noexcept { return b3; }

	inline double get_a1() const noexcept { return a1; }
	inline double get_b1() const noexcept { return b1; }
	inline double get_coef1() const noexcept { return coef1; }
	inline double get_a2() const noexcept { return a2; }
	inline double get_b2() const noexcept { return b2; }
	inline double get_coef2() const noexcept { return coef2; }
	inline double get_a3() const noexcept { return a3; }
	inline double get_b3() const noexcept { return b3; }
	inline double get_coef3() const noexcept { return coef3; }
};


//template <typename Node3D>
//struct TriangleAABBCollisionSAT
//{
//protected:
//	double hx, hy, hz;
//	Point3D n1, n2, n3;
//	// 1 face normal + 3 * 3 edge cross product
//	Vector3D axes[10];
//
//	inline bool is_seperating_axis(Vector3D& axis,
//		Point3D& p1, Point3D& p2, Point3D& p3)
//	{
//#define Norm_Tol 1.0e-6
//		if (axis.norm() < Norm_Tol)
//			return false;
//		double box_range = 0.5 * (hx * abs(axis.x) + hy * abs(axis.y) + hz * abs(axis.z)) * (1.0 + Norm_Tol);
//		double p1_proj = p1.x * axis.x + p1.y * axis.y + p1.z * axis.z;
//		double p2_proj = p2.x * axis.x + p2.y * axis.y + p2.z * axis.z;
//		double p3_proj = p3.x * axis.x + p3.y * axis.y + p3.z * axis.z;
//		return ((p1_proj >  box_range && p2_proj >  box_range && p3_proj >  box_range) ||
//				(p1_proj < -box_range && p2_proj < -box_range && p3_proj < -box_range));
//#undef Norm_Tol
//	}
//
//public:
//	void init_triangle(Node3D& _n1, Node3D& _n2, Node3D& _n3)
//	{
//		n1.x = _n1.x;
//		n1.y = _n1.y;
//		n1.z = _n1.z;
//		n2.x = _n2.x;
//		n2.y = _n2.y;
//		n2.z = _n2.z;
//		n3.x = _n3.x;
//		n3.y = _n3.y;
//		n3.z = _n3.z;
//		double e12_x = n1.x - n2.x;
//		double e12_y = n1.y - n2.y;
//		double e12_z = n1.z - n2.z;
//		double e13_x = n1.x - n3.x;
//		double e13_y = n1.y - n3.y;
//		double e13_z = n1.z - n3.z;
//		double e23_x = n2.x - n3.x;
//		double e23_y = n2.y - n3.y;
//		double e23_z = n2.z - n3.z;
//		// 1 face normal
//		axes[0].cross(e13_x, e13_y, e13_z, e12_x, e12_y, e12_z);
//		// 3 * 3 edge cross product
//		axes[1].cross(1.0, 0.0, 0.0, e12_x, e12_y, e12_z);
//		axes[2].cross(0.0, 1.0, 0.0, e12_x, e12_y, e12_z);
//		axes[3].cross(0.0, 0.0, 1.0, e12_x, e12_y, e12_z);
//		axes[4].cross(1.0, 0.0, 0.0, e13_x, e13_y, e13_z);
//		axes[5].cross(0.0, 1.0, 0.0, e13_x, e13_y, e13_z);
//		axes[6].cross(0.0, 0.0, 1.0, e13_x, e13_y, e13_z);
//		axes[7].cross(1.0, 0.0, 0.0, e23_x, e23_y, e23_z);
//		axes[8].cross(0.0, 1.0, 0.0, e23_x, e23_y, e23_z);
//		axes[9].cross(0.0, 0.0, 1.0, e23_x, e23_y, e23_z);
//	}
//
//	bool detect_collision_with_cube(Cube& cube)
//	{
//		hx = cube.xu - cube.xl;
//		hy = cube.yu - cube.yl;
//		hz = cube.zu - cube.zl;
//		double box_xc = (cube.xl + cube.xu) * 0.5;
//		double box_yc = (cube.yl + cube.yu) * 0.5;
//		double box_zc = (cube.zl + cube.zu) * 0.5;
//		Point3D n1_m, n2_m, n3_m;
//		n1_m.x = n1.x - box_xc;
//		n1_m.y = n1.y - box_yc;
//		n1_m.z = n1.z - box_zc;
//		n2_m.x = n2.x - box_xc;
//		n2_m.y = n2.y - box_yc;
//		n2_m.z = n2.z - box_zc;
//		n3_m.x = n3.x - box_xc;
//		n3_m.y = n3.y - box_yc;
//		n3_m.z = n3.z - box_zc;
//		// if there is one seperating axis, there is no collision
//		if (is_seperating_axis(axes[0], n1_m, n2_m, n3_m) ||
//			is_seperating_axis(axes[1], n1_m, n2_m, n3_m) ||
//			is_seperating_axis(axes[2], n1_m, n2_m, n3_m) ||
//			is_seperating_axis(axes[3], n1_m, n2_m, n3_m) ||
//			is_seperating_axis(axes[4], n1_m, n2_m, n3_m) ||
//			is_seperating_axis(axes[5], n1_m, n2_m, n3_m) ||
//			is_seperating_axis(axes[6], n1_m, n2_m, n3_m) ||
//			is_seperating_axis(axes[7], n1_m, n2_m, n3_m) ||
//			is_seperating_axis(axes[8], n1_m, n2_m, n3_m) ||
//			is_seperating_axis(axes[9], n1_m, n2_m, n3_m))
//			return false;
//		return true;
//	}
//
//	inline const Point3D& get_n1() const noexcept { return n1; }
//	inline const Point3D& get_n2() const noexcept { return n2; }
//	inline const Point3D& get_n3() const noexcept { return n3; }
//};

//template <typename Node3D>
//struct PointToLineDistance
//{
//protected:
//	Point3D n1, n2, n3;
//	union
//	{
//		struct { Vector3D ix1, iy1, iz1; };
//		double T1[3][3];
//	};
//	union
//	{
//		struct { Vector3D ix2, iy2, iz2; };
//		double T2[3][3];
//	};
//	union
//	{
//		struct { Vector3D ix3, iy3, iz3; };
//		double T3[3][3];
//	};
//	double a1, a2, a3;
//
//public:
//	PointToTriangleDistance() {}
//	void init_triangle(Node3D& _n1, Node3D& _n2, Node3D& _n3)
//	{
//		n1.x = _n1.x;
//		n1.y = _n1.y;
//		n1.z = _n1.z;
//		n2.x = _n2.x;
//		n2.y = _n2.y;
//		n2.z = _n2.z;
//		n3.x = _n3.x;
//		n3.y = _n3.y;
//		n3.z = _n3.z;
//		Vector3D tmp1, tmp2;
//		double coef_tmp;
//		// line 1 at x axis
//		ix1.substract<Point3D>(n2, n1);
//		a1 = ix1.norm();
//		ix1.scale(1.0 / a1);
//		tmp1.substract<Point3D>(n3, n1);
//		coef_tmp = tmp1.dot(ix1);
//		tmp2.x = coef_tmp * ix1.x;
//		tmp2.y = coef_tmp * ix1.y;
//		tmp2.z = coef_tmp * ix1.z;
//		iy1.substract<Vector3D>(tmp1, tmp2);
//		iy1.normalize();
//		iz1.cross<Vector3D>(ix1, iy1);
//		// line 2 at x axis
//		ix2.substract<Point3D>(n3, n2);
//		a2 = ix2.norm();
//		ix2.scale(1.0 / a2);
//		tmp1.substract<Point3D>(n1, n2);
//		coef_tmp = tmp1.dot(ix2);
//		tmp2.x = coef_tmp * ix2.x;
//		tmp2.y = coef_tmp * ix2.y;
//		tmp2.z = coef_tmp * ix2.z;
//		iy2.substract<Vector3D>(tmp1, tmp2);
//		iy2.normalize();
//		iz2.cross<Vector3D>(ix2, iy2);
//		// line 3 at x axis
//		ix3.substract<Point3D>(n1, n3);
//		a3 = ix3.norm();
//		ix3.scale(1.0 / a3);
//		tmp1.substract<Point3D>(n2, n3);
//		coef_tmp = tmp1.dot(ix3);
//		tmp2.x = coef_tmp * ix3.x;
//		tmp2.y = coef_tmp * ix3.y;
//		tmp2.z = coef_tmp * ix3.z;
//		iy3.substract<Vector3D>(tmp1, tmp2);
//		iy3.normalize();
//		iz3.cross<Vector3D>(ix3, iy3);
//	}
//
//	unsigned char cal_distance_to_point(Point3D& p, double& dist)
//	{
//		Vector3D pj1, pj2, pj3, tmp1;
//		tmp1.substract<Point3D>(p, n1);
//		pj1.x = T1[0][0] * tmp1.x + T1[0][1] * tmp1.y + T1[0][2] * tmp1.z;
//		pj1.y = T1[1][0] * tmp1.x + T1[1][1] * tmp1.y + T1[1][2] * tmp1.z;
//		pj1.z = T1[2][0] * tmp1.x + T1[2][1] * tmp1.y + T1[2][2] * tmp1.z;
//		tmp1.substract<Point3D>(p, n2);
//		pj2.x = T2[0][0] * tmp1.x + T2[0][1] * tmp1.y + T2[0][2] * tmp1.z;
//		pj2.y = T2[1][0] * tmp1.x + T2[1][1] * tmp1.y + T2[1][2] * tmp1.z;
//		tmp1.substract<Point3D>(p, n3);
//		pj3.x = T3[0][0] * tmp1.x + T3[0][1] * tmp1.y + T3[0][2] * tmp1.z;
//		pj3.y = T3[1][0] * tmp1.x + T3[1][1] * tmp1.y + T3[1][2] * tmp1.z;
//
//		if (pj1.y >= 0.0 && pj2.y >= 0.0 && pj3.y >= 0.0)
//		{
//			dist = pj1.z;
//			return 0;
//		}
//
//		if (pj1.y < 0.0 && pj1.x >= 0.0 && pj1.x <= a1)
//		{
//			dist = sqrt(pj1.y * pj1.y + pj1.z * pj1.z);
//			if (pj1.z < 0)
//				dist = -dist;
//			return 1;
//		}
//
//		if (pj2.y < 0.0 && pj2.x >= 0.0 && pj2.x <= a2)
//		{
//			dist = sqrt(pj2.y * pj2.y + pj1.z * pj1.z);
//			if (pj1.z < 0)
//				dist = -dist;
//			return 2;
//		}
//
//		if (pj3.y < 0.0 && pj3.x >= 0.0 && pj3.x <= a3)
//		{
//			dist = sqrt(pj3.y * pj3.y + pj1.z * pj1.z);
//			if (pj1.z < 0)
//				dist = -dist;
//			return 3;
//		}
//
//		//double dist1;
//		if (pj1.y < 0.0 && pj1.x > a1 || pj2.y < 0.0 && pj2.x < 0.0)
//		{
//			//dist1 = sqrt((pj1.x - a1) * (pj1.x - a1) + pj1.y * pj1.y + pj1.z * pj1.z);
//			dist = sqrt(pj2.x * pj2.x + pj2.y * pj2.y + pj1.z * pj1.z);
//			if (pj1.z < 0)
//				dist = -dist;
//			return 4;
//		}
//
//		if (pj2.y < 0.0 && pj2.x > a2 || pj3.y < 0.0 && pj3.x < 0.0)
//		{
//			//dist1 = sqrt((pj2.x - a2) * (pj2.x - a2) + pj2.y * pj2.y + pj1.z * pj1.z);
//			dist = sqrt(pj3.x * pj3.x + pj3.y * pj3.y + pj1.z * pj1.z);
//			if (pj1.z < 0)
//				dist = -dist;
//			return 5;
//		}
//
//		if (pj3.y < 0.0 && pj3.x > a3 || pj1.y < 0.0 && pj1.x < 0.0)
//		{
//			//dist1 = sqrt((pj3.x - a3) * (pj3.x - a3) + pj3.y * pj3.y + pj1.z * pj1.z);
//			dist = sqrt(pj1.x * pj1.x + pj1.y * pj1.y + pj1.z * pj1.z);
//			if (pj1.z < 0)
//				dist = -dist;
//			return 6;
//		}
//
//		assert(0);
//		return 7;
//	}
//
//	void cal_normal_to_point(Point3D& pt, unsigned char norm_type, Vector3D& normal)
//	{
//		Vector3D v1, v2, tmp1;
//		double coef;
//		switch (norm_type)
//		{
//		case 0:
//			v1.substract(n2, n1);
//			v2.substract(n3, n1);
//			normal.cross(v1, v2);
//			normal.normalize();
//			tmp1.substract<Point3D>(pt, n1);
//			coef = T1[2][0] * tmp1.x + T1[2][1] * tmp1.y + T1[2][2] * tmp1.z;
//			if (coef < 0.0)
//				normal.reverse();
//			return;
//		case 1:
//			v1.substract(n2, n1);
//			v1.normalize();
//			v2.substract(pt, n1);
//			coef = v1.dot(v2);
//			normal.substract(v2, v1.scale(coef));
//			break;
//		case 2:
//			v1.substract(n3, n2);
//			v1.normalize();
//			v2.substract(pt, n2);
//			coef = v1.dot(v2);
//			normal.substract(v2, v1.scale(coef));
//			break;
//		case 3:
//			v1.substract(n1, n3);
//			v1.normalize();
//			v2.substract(pt, n3);
//			coef = v1.dot(v2);
//			v1.scale(coef);
//			normal.substract(v2, v1);
//			break;
//		case 4:
//			normal.substract(pt, n2);
//			break;
//		case 5:
//			normal.substract(pt, n3);
//			break;
//		case 6:
//			normal.substract(pt, n1);
//			break;
//		default:
//			assert(0);
//			return;
//		}
//		double norm = normal.norm();
//		if (norm != 0.0)
//			normal.scale(1.0 / norm);
//		else
//		{
//			v1.substract(n2, n1);
//			v2.substract(n3, n1);
//			normal.cross(v1, v2);
//			normal.normalize();
//		}
//		return;
//	}
//};
//
//template <typename Point3D>
//void cal_triangle_moi(
//	double xc,
//	double yc,
//	Point3D &p1,
//	Point3D &p2,
//	Point3D &p3, 
//	double vol,
//	double &moi
//	)
//{
//	moi_mat[0] = vol / 10.0 *
//		 (p1.y * p1.y + p1.y * p2.y + p2.y * p2.y + p1.y * p3.y + p2.y * p3.y
//		+ p3.y * p3.y + p1.y * p4.y + p2.y * p4.y + p3.y * p4.y + p4.y * p4.y
//		+ p1.z * p1.z + p1.z * p2.z + p2.z * p2.z + p1.z * p3.z + p2.z * p3.z
//		+ p3.z * p3.z + p1.z * p4.z + p2.z * p4.z + p3.z * p4.z + p4.z * p4.z);
//	moi_mat[1] = vol / 10.0 *
//		 (p1.x * p1.x + p1.x * p2.x + p2.x * p2.x + p1.x * p3.x + p2.x * p3.x
//		+ p3.x * p3.x + p1.x * p4.x + p2.x * p4.x + p3.x * p4.x + p4.x * p4.x
//		+ p1.z * p1.z + p1.z * p2.z + p2.z * p2.z + p1.z * p3.z + p2.z * p3.z
//		+ p3.z * p3.z + p1.z * p4.z + p2.z * p4.z + p3.z * p4.z + p4.z * p4.z);
//	moi_mat[2] = vol / 10.0 *
//		 (p1.x * p1.x + p1.x * p2.x + p2.x * p2.x + p1.x * p3.x + p2.x * p3.x
//		+ p3.x * p3.x + p1.x * p4.x + p2.x * p4.x + p3.x * p4.x + p4.x * p4.x
//		+ p1.y * p1.y + p1.y * p2.y + p2.y * p2.y + p1.y * p3.y + p2.y * p3.y
//		+ p3.y * p3.y + p1.y * p4.y + p2.y * p4.y + p3.y * p4.y + p4.y * p4.y);
//	moi_mat[3] = -vol / 20.0 *
//		 (2.0 * p1.y * p1.z + p2.y * p1.z + p3.y * p1.z + p4.y * p1.z + p1.y * p2.z
//		+ 2.0 * p2.y * p2.z + p3.y * p2.z + p4.y * p2.z + p1.y * p3.z + p2.y * p3.z
//		+ 2.0 * p3.y * p3.z + p4.y * p3.z + p1.y * p4.z + p2.y * p4.z + p3.y * p4.z
//		+ 2.0 * p4.y * p4.z);
//	moi_mat[4] = -vol / 20.0 *
//		 (2.0 * p1.x * p1.z + p2.x * p1.z + p3.x * p1.z + p4.x * p1.z + p1.x * p2.z
//		+ 2.0 * p2.x * p2.z + p3.x * p2.z + p4.x * p2.z + p1.x * p3.z + p2.x * p3.z
//		+ 2.0 * p3.x * p3.z + p4.x * p3.z + p1.x * p4.z + p2.x * p4.z + p3.x * p4.z
//		+ 2.0 * p4.x * p4.z);
//	moi_mat[5] = -vol / 20.0 *
//		 (2.0 * p1.x * p1.y + p2.x * p1.y + p3.x * p1.y + p4.x * p1.y + p1.x * p2.y
//		+ 2.0 * p2.x * p2.y + p3.x * p2.y + p4.x * p2.y + p1.x * p3.y + p2.x * p3.y
//		+ 2.0 * p3.x * p3.y + p4.x * p3.y + p1.x * p4.y + p2.x * p4.y + p3.x * p4.y
//		+ 2.0 * p4.x * p4.y);
//}

#endif