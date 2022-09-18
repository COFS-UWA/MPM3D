#ifndef __Triangle_Utils_h__
#define __Triangle_Utils_h__

#include <assert.h>
#include "Geometry2D.h"

template <typename Point2D>
inline void get_triangle_bounding_box(
	const Point2D& n1,
	const Point2D& n2,
	const Point2D& n3,
	Rect &res) noexcept
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
	const Point2D& p3) noexcept
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
		double area) noexcept
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
		const Node2D& n3) noexcept
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

struct PointToLineDistance
{
	Point2D n1, n2;
	union
	{
		struct { Vector2D ix, iy; };
		double T[2][2];
	};
	double len;

	PointToLineDistance() {}
	~PointToLineDistance() {}

	template <typename Node2D>
	void init_line(const Node2D& _n1, const Node2D& _n2)
	{
		n1.x = _n1.x;
		n1.y = _n1.y;
		n2.x = _n2.x;
		n2.y = _n2.y;
		ix.substract<Point2D>(n2, n1);
		len = ix.norm();
		ix.normalize();
		iy.x = -ix.y;
		iy.y =  ix.x;
	}

	unsigned char cal_distance_to_point(const Point2D &p, double& dist) const noexcept
	{
		Vector2D lp, pj;
		lp.x = p.x - n1.x;
		lp.y = p.y - n1.y;
		pj.x = T[0][0] * lp.x + T[0][1] * lp.y;
		pj.y = T[1][0] * lp.x + T[1][1] * lp.y;

		if (pj.x < 0.0)
		{
			dist = lp.norm();
			if (pj.y < 0.0)
				dist = -dist;
			return 0;
		}

		if (pj.x > len)
		{
			Vector2D lp1;
			lp1.x = p.x - n2.x;
			lp1.y = p.y - n2.y;
			dist = lp1.norm();
			if (pj.y < 0.0)
				dist = -dist;
			return 2;
		}

		dist = pj.y;
		return 1;
	}

	template <typename Point2D>
	void cal_normal_to_point(
		const Point2D& pt,
		unsigned char norm_type,
		Vector2D& normal) const
	{
		double side;
		switch (norm_type)
		{
		case 0:
			normal.substract<Point2D>(pt, n1).normalize();
			side = normal.x * iy.x + normal.y * iy.y;
			if (side < 0.0)
				normal.reverse();
			return;
		case 1:
			normal.x = iy.x;
			normal.y = iy.y;
			return;
		case 2:
			normal.substract<Point2D>(pt, n1);
			side = normal.x * iy.x + normal.y * iy.y;
			normal.substract<Point2D>(pt, n2).normalize();
			if (side < 0.0)
				normal.reverse();
			return;
		}
		assert(0);
		return;
	}
};

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