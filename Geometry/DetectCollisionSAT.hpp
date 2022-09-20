#ifndef __Detect_Collision_SAT_hpp__
#define __Detect_Collision_SAT_hpp__

#include "Geometry2D.h"
#include "TriangleUtils.h"
#include "Geometry3D.h"
#include "TetrahedronUtils.h"

struct DetectLineAABBCollisionSAT
{
protected:
	Rect ln_bbox;
	double hx, hy;
	Point2D n1, n2;
	// normal of the line
	Vector2D axes;

	inline bool is_seperating_axis(Vector2D& axis, Point2D& p1, Point2D& p2)
	{
#define Norm_Tol 1.0e-6
		if (axis.norm() < Norm_Tol)
			return false;
		double box_range = 0.5 * (hx * abs(axis.x) + hy * abs(axis.y)) * (1.0 + Norm_Tol);
		double p1_proj = p1.x * axis.x + p1.y * axis.y;
		double p2_proj = p2.x * axis.x + p2.y * axis.y;
		return ((p1_proj > box_range && p2_proj > box_range) ||
			(p1_proj < -box_range && p2_proj < -box_range));
#undef Norm_Tol
	}

public:
	inline void get_ln_bbox(Rect& box) const noexcept
	{
		box.xl = ln_bbox.xl; box.xu = ln_bbox.xu;
		box.yl = ln_bbox.yl; box.yu = ln_bbox.yu;
	}

	template <typename Node2D>
	void init_line(Node2D& _n1, Node2D& _n2)
	{
		n1.x = _n1.x;
		n1.y = _n1.y;
		n2.x = _n2.x;
		n2.y = _n2.y;
		const double e12_x = n1.x - n2.x;
		const double e12_y = n1.y - n2.y;
		axes.x = e12_y;
		axes.y = -e12_x;
		//
		if (n1.x < n2.x)
		{
			ln_bbox.xl = n1.x;
			ln_bbox.xu = n2.x;
		}
		else
		{
			ln_bbox.xl = n2.x;
			ln_bbox.xu = n1.x;
		}
		if (n1.y < n2.y)
		{
			ln_bbox.yl = n1.y;
			ln_bbox.yu = n2.y;
		}
		else
		{
			ln_bbox.yl = n2.y;
			ln_bbox.yu = n1.y;
		}
	}

	bool detect(const Rect& rect)
	{
		hx = rect.xu - rect.xl;
		hy = rect.yu - rect.yl;
		double box_xc = (rect.xl + rect.xu) * 0.5;
		double box_yc = (rect.yl + rect.yu) * 0.5;
		Point2D n1_m, n2_m;
		n1_m.x = n1.x - box_xc;
		n1_m.y = n1.y - box_yc;
		n2_m.x = n2.x - box_xc;
		n2_m.y = n2.y - box_yc;
		// if there is one seperating axis, there is no collision
		if (is_seperating_axis(axes, n1_m, n2_m))
			return false;
		return true;
	}

	inline const Point2D& get_n1() const noexcept { return n1; }
	inline const Point2D& get_n2() const noexcept { return n2; }
};

struct DetectTriangleAABBCollisionSAT
{
protected:
	double hx, hy;
	Point2D n1, n2, n3;
	Rect tri_bbox;
	// 3 face normal
	Vector2D axes[3];

	inline bool is_seperating_axis(
		Vector2D& axis,
		Point2D &p1,
		Point2D &p2,
		Point2D &p3)
	{
#define Norm_Tol 1.0e-6
		if (axis.norm() < Norm_Tol)
			return false;
		double box_range = 0.5 * (hx * abs(axis.x) + hy * abs(axis.y)) * (1.0 + Norm_Tol);
		double p1_proj = p1.x * axis.x + p1.y * axis.y;
		double p2_proj = p2.x * axis.x + p2.y * axis.y;
		double p3_proj = p3.x * axis.x + p3.y * axis.y;
		return ((p1_proj >  box_range && p2_proj >  box_range && p3_proj >  box_range) ||
				(p1_proj < -box_range && p2_proj < -box_range && p3_proj < -box_range));
#undef Norm_Tol
	}

public:
	inline void get_tri_bbox(Rect& box) const noexcept
	{
		box.xl = tri_bbox.xl; box.xu = tri_bbox.xu;
		box.yl = tri_bbox.yl; box.yu = tri_bbox.yu;
	}
	
	template <typename Node2D>
	void init_triangle(Node2D& _n1, Node2D& _n2, Node2D& _n3)
	{
		n1.x = _n1.x;
		n1.y = _n1.y;
		n2.x = _n2.x;
		n2.y = _n2.y;
		n3.x = _n3.x;
		n3.y = _n3.y;
		const double e12_x = n1.x - n2.x;
		const double e12_y = n1.y - n2.y;
		axes[0].x =  e12_y;
		axes[0].y = -e12_x;
		const double e13_x = n1.x - n3.x;
		const double e13_y = n1.y - n3.y;
		axes[1].x =  e13_y;
		axes[1].y = -e13_x;
		const double e23_x = n2.x - n3.x;
		const double e23_y = n2.y - n3.y;
		axes[2].x =  e23_y;
		axes[2].y = -e23_x;
		get_triangle_bounding_box(_n1, _n2, _n3, tri_bbox);
	}

	bool detect(const Rect &rect)
	{
		hx = rect.xu - rect.xl;
		hy = rect.yu - rect.yl;
		double box_xc = (rect.xl + rect.xu) * 0.5;
		double box_yc = (rect.yl + rect.yu) * 0.5;
		Point2D n1_m, n2_m, n3_m;
		n1_m.x = n1.x - box_xc;
		n1_m.y = n1.y - box_yc;
		n2_m.x = n2.x - box_xc;
		n2_m.y = n2.y - box_yc;
		n3_m.x = n3.x - box_xc;
		n3_m.y = n3.y - box_yc;
		// if there is one seperating axis, there is no collision
		if (is_seperating_axis(axes[0], n1_m, n2_m, n3_m) ||
			is_seperating_axis(axes[1], n1_m, n2_m, n3_m) ||
			is_seperating_axis(axes[2], n1_m, n2_m, n3_m))
			return false;
		return true;
	}

	inline const Point2D& get_n1() const noexcept { return n1; }
	inline const Point2D& get_n2() const noexcept { return n2; }
	inline const Point2D& get_n3() const noexcept { return n3; }
};

class DetectTetrahedronAABBCollisionSAT
{
public:
	inline DetectTetrahedronAABBCollisionSAT() {}
	~DetectTetrahedronAABBCollisionSAT() {}

	inline void get_teh_bbox(Cube& box) const noexcept
	{
		box.xl = teh_bbox.xl; box.xu = teh_bbox.xu;
		box.yl = teh_bbox.yl; box.yu = teh_bbox.yu;
		box.zl = teh_bbox.zl; box.zu = teh_bbox.zu;
	}
	inline void get_teh_n1(Point3D& n) const noexcept
	{ n.x = teh_n1.x; n.y = teh_n1.y; n.z = teh_n1.z; }
	inline void get_teh_n2(Point3D& n) const noexcept
	{ n.x = teh_n2.x; n.y = teh_n2.y; n.z = teh_n2.z; }
	inline void get_teh_n3(Point3D &n) const noexcept
	{ n.x = teh_n3.x; n.y = teh_n3.y; n.z = teh_n3.z; }
	inline void get_teh_n4(Point3D &n) const noexcept
	{ n.x = teh_n4.x; n.y = teh_n4.y; n.z = teh_n4.z; }

	template <typename Node3D>
	void init_tetrahedron(
		const Node3D &n1,
		const Node3D &n2,
		const Node3D &n3,
		const Node3D &n4
		) noexcept
	{
		get_tetrahedron_bounding_box(n1, n2, n3, n4, teh_bbox);
		teh_n1.x = n1.x; teh_n1.y = n1.y; teh_n1.z = n1.z;
		teh_n2.x = n2.x; teh_n2.y = n2.y; teh_n2.z = n2.z;
		teh_n3.x = n3.x; teh_n3.y = n3.y; teh_n3.z = n3.z;
		teh_n4.x = n4.x; teh_n4.y = n4.y; teh_n4.z = n4.z;
		pt_in_teh.init_tetrahedron(n1, n2, n3, n4);
		cal_seperating_axes(n1, n2, n3, n4);
	}

	bool detect(const Cube& aabb) const noexcept
	{
		if (aabb.xu < teh_bbox.xl || aabb.xl > teh_bbox.xu ||
			aabb.yu < teh_bbox.yl || aabb.yl > teh_bbox.yu ||
			aabb.zu < teh_bbox.zl || aabb.zl > teh_bbox.zu)
			return false;

		// whether tetrahedron nodes locate in box
		// efficient when tetrahedron is much smaller than grid
		if (aabb.is_in_box(teh_n1) || aabb.is_in_box(teh_n2) ||
			aabb.is_in_box(teh_n3) || aabb.is_in_box(teh_n4))
			return true;

		// whether box corners locate in tetrahedron
		// efficient when grid is much smaller than tetrahedron
		if (pt_in_teh.is_in_tetrahedron(aabb.xl, aabb.yl, aabb.zl) ||
			pt_in_teh.is_in_tetrahedron(aabb.xl, aabb.yl, aabb.zu) ||
			pt_in_teh.is_in_tetrahedron(aabb.xl, aabb.yu, aabb.zl) ||
			pt_in_teh.is_in_tetrahedron(aabb.xl, aabb.yu, aabb.zu) ||
			pt_in_teh.is_in_tetrahedron(aabb.xu, aabb.yl, aabb.zl) ||
			pt_in_teh.is_in_tetrahedron(aabb.xu, aabb.yl, aabb.zu) ||
			pt_in_teh.is_in_tetrahedron(aabb.xu, aabb.yu, aabb.zl) ||
			pt_in_teh.is_in_tetrahedron(aabb.xu, aabb.yu, aabb.zu))
			return true;

		// applied separating axis theory
		// for case when grid and tetrahedron is of comparable size
		// move origin to centre of the Cube
		const double box_xc = (aabb.xl + aabb.xu) * 0.5;
		const double box_yc = (aabb.yl + aabb.yu) * 0.5;
		const double box_zc = (aabb.zl + aabb.zu) * 0.5;
		const double box_h[3] = { aabb.xu - aabb.xl, aabb.yu - aabb.yl, aabb.zu - aabb.zl };
		// moved tetrahedron nodes
		const Point3D n1_m(teh_n1.x - box_xc, teh_n1.y - box_yc, teh_n1.z - box_zc);
		const Point3D n2_m(teh_n2.x - box_xc, teh_n2.y - box_yc, teh_n2.z - box_zc);
		const Point3D n3_m(teh_n3.x - box_xc, teh_n3.y - box_yc, teh_n3.z - box_zc);
		const Point3D n4_m(teh_n4.x - box_xc, teh_n4.y - box_yc, teh_n4.z - box_zc);
		// if there is one seperating axis, there is no collision
		return !(is_seperating_axis(seperating_axes[0], box_h, n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[1], box_h, n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[2], box_h, n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[3], box_h, n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[4], box_h, n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[5], box_h, n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[6], box_h, n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[7], box_h, n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[8], box_h, n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[9], box_h, n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[10], box_h, n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[11], box_h, n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[12], box_h, n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[13], box_h, n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[14], box_h, n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[15], box_h, n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[16], box_h, n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[17], box_h, n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[18], box_h, n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[19], box_h, n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[20], box_h, n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[21], box_h, n1_m, n2_m, n3_m, n4_m));
	}

protected:
	Cube teh_bbox;
	Point3D teh_n1, teh_n2, teh_n3, teh_n4;
	PointInTetrahedron pt_in_teh;
	// 22 seperating axises:
	// 4 face normals and 3 * 6 edge cross-products
	Vector3D seperating_axes[22];
	
	inline bool is_seperating_axis(
		const Vector3D& sa,
		const double box_len[3] /* hx, hy, hz */,
		const Point3D &p1,
		const Point3D &p2,
		const Point3D &p3,
		const Point3D &p4
		) const noexcept
	{
		constexpr double norm_tol = 1.0e-8;
		if (sa.norm() < norm_tol)
			return false;
		const double box_range = 0.5 * (box_len[0] * abs(sa.x)
			+ box_len[1] * abs(sa.y) + box_len[2] * abs(sa.z));
		const double p1_proj = p1.x * sa.x + p1.y * sa.y + p1.z * sa.z;
		const double p2_proj = p2.x * sa.x + p2.y * sa.y + p2.z * sa.z;
		const double p3_proj = p3.x * sa.x + p3.y * sa.y + p3.z * sa.z;
		const double p4_proj = p4.x * sa.x + p4.y * sa.y + p4.z * sa.z;
		return (p1_proj >  box_range && p2_proj >  box_range && p3_proj >  box_range && p4_proj >  box_range) ||
			   (p1_proj < -box_range && p2_proj < -box_range && p3_proj < -box_range && p4_proj < -box_range);
	}

	template <typename Node3D>
	void cal_seperating_axes(
		const Node3D &n1,
		const Node3D &n2,
		const Node3D &n3,
		const Node3D &n4
		) noexcept
	{
		const double e12_x = n1.x - n2.x;
		const double e12_y = n1.y - n2.y;
		const double e12_z = n1.z - n2.z;
		const double e13_x = n1.x - n3.x;
		const double e13_y = n1.y - n3.y;
		const double e13_z = n1.z - n3.z;
		const double e14_x = n1.x - n4.x;
		const double e14_y = n1.y - n4.y;
		const double e14_z = n1.z - n4.z;
		const double e23_x = n2.x - n3.x;
		const double e23_y = n2.y - n3.y;
		const double e23_z = n2.z - n3.z;
		const double e24_x = n2.x - n4.x;
		const double e24_y = n2.y - n4.y;
		const double e24_z = n2.z - n4.z;
		const double e34_x = n3.x - n4.x;
		const double e34_y = n3.y - n4.y;
		const double e34_z = n3.z - n4.z;
		// 4 face normal
		seperating_axes[0].cross(e13_x, e13_y, e13_z, e12_x, e12_y, e12_z);
		seperating_axes[1].cross(e12_x, e12_y, e12_z, e14_x, e14_y, e14_z);
		seperating_axes[2].cross(e14_x, e14_y, e14_z, e13_x, e13_y, e13_z);
		seperating_axes[3].cross(e23_x, e23_y, e23_z, e24_x, e24_y, e24_z);
		// 3 * 6 edge cross product
		seperating_axes[4].cross(1.0, 0.0, 0.0, e12_x, e12_y, e12_z);
		seperating_axes[5].cross(0.0, 1.0, 0.0, e12_x, e12_y, e12_z);
		seperating_axes[6].cross(0.0, 0.0, 1.0, e12_x, e12_y, e12_z);
		seperating_axes[7].cross(1.0, 0.0, 0.0, e13_x, e13_y, e13_z);
		seperating_axes[8].cross(0.0, 1.0, 0.0, e13_x, e13_y, e13_z);
		seperating_axes[9].cross(0.0, 0.0, 1.0, e13_x, e13_y, e13_z);
		seperating_axes[10].cross(1.0, 0.0, 0.0, e14_x, e14_y, e14_z);
		seperating_axes[11].cross(0.0, 1.0, 0.0, e14_x, e14_y, e14_z);
		seperating_axes[12].cross(0.0, 0.0, 1.0, e14_x, e14_y, e14_z);
		seperating_axes[13].cross(1.0, 0.0, 0.0, e23_x, e23_y, e23_z);
		seperating_axes[14].cross(0.0, 1.0, 0.0, e23_x, e23_y, e23_z);
		seperating_axes[15].cross(0.0, 0.0, 1.0, e23_x, e23_y, e23_z);
		seperating_axes[16].cross(1.0, 0.0, 0.0, e24_x, e24_y, e24_z);
		seperating_axes[17].cross(0.0, 1.0, 0.0, e24_x, e24_y, e24_z);
		seperating_axes[18].cross(0.0, 0.0, 1.0, e24_x, e24_y, e24_z);
		seperating_axes[19].cross(1.0, 0.0, 0.0, e34_x, e34_y, e34_z);
		seperating_axes[20].cross(0.0, 1.0, 0.0, e34_x, e34_y, e34_z);
		seperating_axes[21].cross(0.0, 0.0, 1.0, e34_x, e34_y, e34_z);
	}
};

class Detect3DTriangleAABBCollisionSAT
{
public:
	inline Detect3DTriangleAABBCollisionSAT() {}
	~Detect3DTriangleAABBCollisionSAT() {}
	
	inline void get_tri_bbox(Cube& box) const noexcept
	{
		box.xl = tri_bbox.xl; box.xu = tri_bbox.xu;
		box.yl = tri_bbox.yl; box.yu = tri_bbox.yu;
		box.zl = tri_bbox.zl; box.zu = tri_bbox.zu;
	}
	inline void get_tri_n1(Point3D& n) const noexcept
	{ n.x = tri_n1.x; n.y = tri_n1.y; n.z = tri_n1.z; }
	inline void get_tri_n2(Point3D& n) const noexcept
	{ n.x = tri_n2.x; n.y = tri_n2.y; n.z = tri_n2.z; }
	inline void get_tri_n3(Point3D &n) const noexcept
	{ n.x = tri_n3.x; n.y = tri_n3.y; n.z = tri_n3.z; }
	inline const Point3D& get_tri_n1() const noexcept { return tri_n1; }
	inline const Point3D& get_tri_n2() const noexcept { return tri_n2; }
	inline const Point3D& get_tri_n3() const noexcept { return tri_n3; }

	template <typename Node3D>
	void init_triangle(
		const Node3D& n1,
		const Node3D& n2,
		const Node3D& n3
		) noexcept
	{
		get_3Dtriangle_bounding_box(n1, n2, n3, tri_bbox);
		tri_n1.x = n1.x; tri_n1.y = n1.y; tri_n1.z = n1.z;
		tri_n2.x = n2.x; tri_n2.y = n2.y; tri_n2.z = n2.z;
		tri_n3.x = n3.x; tri_n3.y = n3.y; tri_n3.z = n3.z;
		cal_seperating_axes(n1, n2, n3);
	}

	bool detect(const Cube& aabb) const noexcept
	{
		if (aabb.xu < tri_bbox.xl || aabb.xl > tri_bbox.xu ||
			aabb.yu < tri_bbox.yl || aabb.yl > tri_bbox.yu ||
			aabb.zu < tri_bbox.zl || aabb.zl > tri_bbox.zu)
			return false;

		// whether tetrahedron nodes locate in box
		// efficient when tetrahedron is much smaller than grid
		if (aabb.is_in_box(tri_n1) || aabb.is_in_box(tri_n2) ||
			aabb.is_in_box(tri_n3))
			return true;

		const double box_h[3] = { aabb.xu - aabb.xl,
			   aabb.yu - aabb.yl, aabb.zu - aabb.zl };
		const double box_xc = (aabb.xl + aabb.xu) * 0.5;
		const double box_yc = (aabb.yl + aabb.yu) * 0.5;
		const double box_zc = (aabb.zl + aabb.zu) * 0.5;
		Point3D n1_m(tri_n1.x - box_xc, tri_n1.y - box_yc, tri_n1.z - box_zc);
		Point3D n2_m(tri_n2.x - box_xc, tri_n2.y - box_yc, tri_n2.z - box_zc);
		Point3D n3_m(tri_n3.x - box_xc, tri_n3.y - box_yc, tri_n3.z - box_zc);
		// if there is one seperating axis, there is no collision
		return !(is_seperating_axis(seperating_axes[0], box_h, n1_m, n2_m, n3_m) ||
			is_seperating_axis(seperating_axes[1], box_h, n1_m, n2_m, n3_m) ||
			is_seperating_axis(seperating_axes[2], box_h, n1_m, n2_m, n3_m) ||
			is_seperating_axis(seperating_axes[3], box_h, n1_m, n2_m, n3_m) ||
			is_seperating_axis(seperating_axes[4], box_h, n1_m, n2_m, n3_m) ||
			is_seperating_axis(seperating_axes[5], box_h, n1_m, n2_m, n3_m) ||
			is_seperating_axis(seperating_axes[6], box_h, n1_m, n2_m, n3_m) ||
			is_seperating_axis(seperating_axes[7], box_h, n1_m, n2_m, n3_m) ||
			is_seperating_axis(seperating_axes[8], box_h, n1_m, n2_m, n3_m) ||
			is_seperating_axis(seperating_axes[9], box_h, n1_m, n2_m, n3_m));
	}

protected:
	Cube tri_bbox;
	Point3D tri_n1, tri_n2, tri_n3;
	// 10 seperating axises:
	// 1 face normals and 3 * 3 edge cross-products
	Vector3D seperating_axes[10];

	inline bool is_seperating_axis(
		const Vector3D &sa,
		const double box_len[3] /* hx, hy, hz */,
		const Point3D &p1,
		const Point3D &p2,
		const Point3D &p3
		) const noexcept
	{
		constexpr double norm_tol = 1.0e-8;
		if (sa.norm() < norm_tol)
			return false;
		const double box_range = 0.5 * (box_len[0] * abs(sa.x)
			+ box_len[1] * abs(sa.y) + box_len[2] * abs(sa.z));
		const double p1_proj = p1.x * sa.x + p1.y * sa.y + p1.z * sa.z;
		const double p2_proj = p2.x * sa.x + p2.y * sa.y + p2.z * sa.z;
		const double p3_proj = p3.x * sa.x + p3.y * sa.y + p3.z * sa.z;
		return (p1_proj >  box_range && p2_proj >  box_range && p3_proj >  box_range) ||
			   (p1_proj < -box_range && p2_proj < -box_range && p3_proj < -box_range);
	}

	template <typename Node3D>
	void cal_seperating_axes(
		const Node3D& n1,
		const Node3D& n2,
		const Node3D& n3
		) noexcept
	{
		const double e12_x = n1.x - n2.x;
		const double e12_y = n1.y - n2.y;
		const double e12_z = n1.z - n2.z;
		const double e13_x = n1.x - n3.x;
		const double e13_y = n1.y - n3.y;
		const double e13_z = n1.z - n3.z;
		const double e23_x = n2.x - n3.x;
		const double e23_y = n2.y - n3.y;
		const double e23_z = n2.z - n3.z;
		// 1 face normal
		seperating_axes[0].cross(e13_x, e13_y, e13_z, e12_x, e12_y, e12_z);
		// 3 * 3 edge cross product
		seperating_axes[1].cross(1.0, 0.0, 0.0, e12_x, e12_y, e12_z);
		seperating_axes[2].cross(0.0, 1.0, 0.0, e12_x, e12_y, e12_z);
		seperating_axes[3].cross(0.0, 0.0, 1.0, e12_x, e12_y, e12_z);
		seperating_axes[4].cross(1.0, 0.0, 0.0, e13_x, e13_y, e13_z);
		seperating_axes[5].cross(0.0, 1.0, 0.0, e13_x, e13_y, e13_z);
		seperating_axes[6].cross(0.0, 0.0, 1.0, e13_x, e13_y, e13_z);
		seperating_axes[7].cross(1.0, 0.0, 0.0, e23_x, e23_y, e23_z);
		seperating_axes[8].cross(0.0, 1.0, 0.0, e23_x, e23_y, e23_z);
		seperating_axes[9].cross(0.0, 0.0, 1.0, e23_x, e23_y, e23_z);
	}
};

struct OBB3DAABBCollisionSAT
{
	double xhlen, yhlen, zhlen;
	Vector3D ix, iy, iz;
	double ix_min, ix_max;
	double iy_min, iy_max;
	double iz_min, iz_max;
	double it1x_min, it1x_max;
	double it1y_min, it1y_max;
	double it1z_min, it1z_max;
	double it2x_min, it2x_max;
	double it2y_min, it2y_max;
	double it2z_min, it2z_max;
	double it3x_min, it3x_max;
	double it3y_min, it3y_max;
	double it3z_min, it3z_max;
	bool it1x_not_zero, it1y_not_zero, it1z_not_zero;
	bool it2x_not_zero, it2y_not_zero, it2z_not_zero;
	bool it3x_not_zero, it3y_not_zero, it3z_not_zero;

	inline void init_obb3d(OBB3D& obb)
	{
		constexpr double norm_tol = 1.0e-8;
		xhlen = obb.hx * 0.5;
		yhlen = obb.hy * 0.5;
		zhlen = obb.hz * 0.5;
		ix = obb.ix;
		iy = obb.iy;
		iz = obb.iz;
		double cen, radius;
		// ix
		cen = ix.x * obb.xo + ix.y * obb.yo + ix.z * obb.zo;
		ix_min = cen - xhlen;
		ix_max = cen + xhlen;
		// iy
		cen = iy.x * obb.xo + iy.y * obb.yo + iy.z * obb.zo;
		iy_min = cen - yhlen;
		iy_max = cen + yhlen;
		// iz
		cen = iz.x * obb.xo + iz.y * obb.yo + iz.z * obb.zo;
		iz_min = cen - zhlen;
		iz_max = cen + zhlen;
		// it1x
		if ((ix.z * ix.z + ix.y * ix.y) > (norm_tol * norm_tol))
		{
			it1x_not_zero = true;
			cen = -ix.z * obb.yo + ix.y * obb.zo;
			radius = xhlen * abs(-ix.z * ix.y + ix.y * ix.z)
				+ yhlen * abs(-ix.z * iy.y + ix.y * iy.z)
				+ zhlen * abs(-ix.z * iz.y + ix.y * iz.z);
			it1x_min = cen - radius;
			it1x_max = cen + radius;
		}
		else
		{
			it1x_not_zero = false;
		}
		// it1y
		if ((iy.z * iy.z + iy.y * iy.y) > (norm_tol * norm_tol))
		{
			it1y_not_zero = true;
			cen = -iy.z * obb.yo + iy.y * obb.zo;
			radius = xhlen * abs(-iy.z * ix.y + iy.y * ix.z)
				+ yhlen * abs(-iy.z * iy.y + iy.y * iy.z)
				+ zhlen * abs(-iy.z * iz.y + iy.y * iz.z);
			it1y_min = cen - radius;
			it1y_max = cen + radius;
		}
		else
		{
			it1y_not_zero = false;
		}
		// it1z
		if ((iz.z * iz.z + iz.y * iz.y) > (norm_tol * norm_tol))
		{
			it1z_not_zero = true;
			cen = -iz.z * obb.yo + iz.y * obb.zo;
			radius = xhlen * abs(-iz.z * ix.y + iz.y * ix.z)
				+ yhlen * abs(-iz.z * iy.y + iz.y * iy.z)
				+ zhlen * abs(-iz.z * iz.y + iz.y * iz.z);
			it1z_min = cen - radius;
			it1z_max = cen + radius;
		}
		else
		{
			it1z_not_zero = false;
		}
		// it2x
		if ((ix.z * ix.z + ix.x * ix.x) > (norm_tol * norm_tol))
		{
			it2x_not_zero = true;
			cen = ix.z * obb.xo - ix.x * obb.zo;
			radius = xhlen * abs(ix.z * ix.x - ix.x * ix.z)
				+ yhlen * abs(ix.z * iy.x - ix.x * iy.z)
				+ zhlen * abs(ix.z * iz.x - ix.x * iz.z);
			it2x_min = cen - radius;
			it2x_max = cen + radius;
		}
		else
		{
			it2x_not_zero = false;
		}
		// it2y
		if ((iy.z * iy.z + iy.x * iy.x) > (norm_tol * norm_tol))
		{
			it2y_not_zero = true;
			cen = iy.z * obb.xo - iy.x * obb.zo;
			radius = xhlen * abs(iy.z * ix.x - iy.x * ix.z)
				+ yhlen * abs(iy.z * iy.x - iy.x * iy.z)
				+ zhlen * abs(iy.z * iz.x - iy.x * iz.z);
			it2y_min = cen - radius;
			it2y_max = cen + radius;
		}
		else
		{
			it2y_not_zero = false;
		}
		// it2z
		if ((iz.z * iz.z + iz.x * iz.x) > (norm_tol * norm_tol))
		{
			it2z_not_zero = true;
			cen = iz.z * obb.xo - iz.x * obb.zo;
			radius = xhlen * abs(iz.z * ix.x - iz.x * ix.z)
				+ yhlen * abs(iz.z * iy.x - iz.x * iy.z)
				+ zhlen * abs(iz.z * iz.x - iz.x * iz.z);
			it2z_min = cen - radius;
			it2z_max = cen + radius;
		}
		else
		{
			it2z_not_zero = false;
		}
		// it3x
		if ((ix.y * ix.y + ix.x * ix.x) > (norm_tol * norm_tol))
		{
			it3x_not_zero = true;
			cen = -ix.y * obb.xo + ix.x * obb.yo;
			radius = xhlen * abs(-ix.y * ix.x + ix.x * ix.y)
				+ yhlen * abs(-ix.y * iy.x + ix.x * iy.y)
				+ zhlen * abs(-ix.y * iz.x + ix.x * iz.y);
			it3x_min = cen - radius;
			it3x_max = cen + radius;
		}
		else
		{
			it3x_not_zero = false;
		}
		// it3y
		if ((iy.y * iy.y + iy.x * iy.x) > (norm_tol * norm_tol))
		{
			it3y_not_zero = true;
			cen = -iy.y * obb.xo + iy.x * obb.yo;
			radius = xhlen * abs(-iy.y * ix.x + iy.x * ix.y)
				+ yhlen * abs(-iy.y * iy.x + iy.x * iy.y)
				+ zhlen * abs(-iy.y * iz.x + iy.x * iz.y);
			it3y_min = cen - radius;
			it3y_max = cen + radius;
		}
		else
		{
			it3y_not_zero = false;
		}
		// it3z
		if ((iz.y * iz.y + iz.x * iz.x) > (norm_tol * norm_tol))
		{
			it3z_not_zero = true;
			cen = -iz.y * obb.xo + iz.x * obb.yo;
			radius = xhlen * abs(-iz.y * ix.x + iz.x * ix.y)
				+ yhlen * abs(-iz.y * iy.x + iz.x * iy.y)
				+ zhlen * abs(-iz.y * iz.x + iz.x * iz.y);
			it3z_min = cen - radius;
			it3z_max = cen + radius;
		}
		else
		{
			it3z_not_zero = false;
		}
	}
	inline bool detect(Cube& cube)
	{
		double cube_xm = (cube.xl + cube.xu) * 0.5;
		double cube_ym = (cube.yl + cube.yu) * 0.5;
		double cube_zm = (cube.zl + cube.zu) * 0.5;
		double cube_xhlen = (cube.xu - cube.xl) * 0.5;
		double cube_yhlen = (cube.yu - cube.yl) * 0.5;
		double cube_zhlen = (cube.zu - cube.zl) * 0.5;
		double ix_cube_cen = ix.x * cube_xm + ix.y * cube_ym + ix.z * cube_zm;
		double ix_cube_range = abs(ix.x * cube_xhlen) + abs(ix.y * cube_yhlen) + abs(ix.z * cube_zhlen);
		double iy_cube_cen = iy.x * cube_xm + iy.y * cube_ym + iy.z * cube_zm;
		double iy_cube_range = abs(iy.x * cube_xhlen) + abs(iy.y * cube_yhlen) + abs(iy.z * cube_zhlen);
		double iz_cube_cen = iz.x * cube_xm + iz.y * cube_ym + iz.z * cube_zm;
		double iz_cube_range = abs(iz.x * cube_xhlen) + abs(iz.y * cube_yhlen) + abs(iz.z * cube_zhlen);
		double it1x_cube_cen = -ix.z * cube_ym + ix.y * cube_zm;
		double it1x_cube_range = abs(-ix.z * cube_yhlen) + abs(ix.y * cube_zhlen);
		double it1y_cube_cen = -iy.z * cube_ym + iy.y * cube_zm;
		double it1y_cube_range = abs(-iy.z * cube_yhlen) + abs(iy.y * cube_zhlen);
		double it1z_cube_cen = -iz.z * cube_ym + iz.y * cube_zm;
		double it1z_cube_range = abs(-iz.z * cube_yhlen) + abs(iz.y * cube_zhlen);
		double it2x_cube_cen = ix.z * cube_xm - ix.x * cube_zm;
		double it2x_cube_range = abs(ix.z * cube_xhlen) + abs(-ix.x * cube_zhlen);
		double it2y_cube_cen = iy.z * cube_xm - iy.x * cube_zm;
		double it2y_cube_range = abs(iy.z * cube_xhlen) + abs(-iy.x * cube_zhlen);
		double it2z_cube_cen = iz.z * cube_xm - iz.x * cube_zm;
		double it2z_cube_range = abs(iz.z * cube_xhlen) + abs(-iz.x * cube_zhlen);
		double it3x_cube_cen = -ix.y * cube_xm + ix.x * cube_ym;
		double it3x_cube_range = abs(-ix.y * cube_xhlen) + abs(ix.x * cube_yhlen);
		double it3y_cube_cen = -iy.y * cube_xm + iy.x * cube_ym;
		double it3y_cube_range = abs(-iy.y * cube_xhlen) + abs(iy.x * cube_yhlen);
		double it3z_cube_cen = -iz.y * cube_xm + iz.x * cube_ym;
		double it3z_cube_range = abs(-iz.y * cube_xhlen) + abs(iz.x * cube_yhlen);

		return !(is_seperating_axis(ix_cube_range, ix_cube_cen, ix_min, ix_max) ||
			is_seperating_axis(iy_cube_range, iy_cube_cen, iy_min, iy_max) ||
			is_seperating_axis(iz_cube_range, iz_cube_cen, iz_min, iz_max) ||
			(it1x_not_zero && is_seperating_axis(it1x_cube_range, it1x_cube_cen, it1x_min, it1x_max)) ||
			(it1y_not_zero && is_seperating_axis(it1y_cube_range, it1y_cube_cen, it1y_min, it1y_max)) ||
			(it1z_not_zero && is_seperating_axis(it1z_cube_range, it1z_cube_cen, it1z_min, it1z_max)) ||
			(it2x_not_zero && is_seperating_axis(it2x_cube_range, it2x_cube_cen, it2x_min, it2x_max)) ||
			(it2y_not_zero && is_seperating_axis(it2y_cube_range, it2y_cube_cen, it2y_min, it2y_max)) ||
			(it2z_not_zero && is_seperating_axis(it2z_cube_range, it2z_cube_cen, it2z_min, it2z_max)) ||
			(it3x_not_zero && is_seperating_axis(it3x_cube_range, it3x_cube_cen, it3x_min, it3x_max)) ||
			(it3y_not_zero && is_seperating_axis(it3y_cube_range, it3y_cube_cen, it3y_min, it3y_max)) ||
			(it3z_not_zero && is_seperating_axis(it3z_cube_range, it3z_cube_cen, it3z_min, it3z_max)));
	}

protected:
	inline bool is_seperating_axis(double cube_range, double cube_cen,
		double obb_min, double obb_max)
	{
		return ((obb_min - cube_cen) > cube_range ||
			(obb_max + cube_cen) < -cube_range);
	}
};

#endif