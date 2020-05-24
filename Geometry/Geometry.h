#ifndef __Geometry_h__
#define __Geometry_h__

#include <cmath>

struct Point2D { double x, y; };
struct Rect { double xl, xu, yl, yu; };

struct Point3D { double x, y, z; };
struct Cube
{
	double xl, xu, yl, yu, zl, zu;
	Cube() {}
	Cube(double _xl, double _xu,
		 double _yl, double _yu,
		 double _zl, double _zu) :
		xl(_xl), xu(_xu), yl(_yl),
		yu(_yu), zl(_zl), zu(_zu) {}
	inline bool is_in_box(double x, double y, double z)
	{
		if (x < xl || x > xu || y < yl || y > yu || z < zl || z > zu)
			return false;
		return true;
	}
	template <typename Point3D>
	inline bool is_in_box(Point3D &point)
	{
		return is_in_box(point.x, point.y, point.z);
	}
};

struct Vector3D
{
	double x, y, z;
	inline double norm() { return sqrt(x*x + y*y + z*z); }
	inline void normalize()
	{
		double len = norm();
		if (len != 0.0)
		{
			x /= len;
			y /= len;
			z /= len;
		}
	}
	inline void cross(double e1_x, double e1_y, double e1_z,
					  double e2_x, double e2_y, double e2_z)
	{
		x = e1_y * e2_z - e1_z * e2_y;
		y = e1_z * e2_x - e1_x * e2_z;
		z = e1_x * e2_y - e1_y * e2_x;
	}
};

// distance between 2D points
template <typename Node2D, typename Point2D>
inline double cal_distance_2D(Node2D &p1, Point2D &p2) noexcept
{
	double dx = p1.x - p2.x;
	double dy = p1.y - p2.y;
	return sqrt(dx * dx + dy * dy);
}

// distance between 3D points
template <typename Node3D, typename Point3D>
inline double cal_distance_3D(Node3D &n, Point3D &p)
{
	double dx = n.x - p.x;
	double dy = n.y - p.y;
	double dz = n.z - p.z;
	return sqrt(dx*dx + dy*dy + dz*dz);
}

// distance between rectangle and 2D point
inline double cal_distance_2D(Rect &rec, Point2D &p) noexcept
{
	double cx, cy, x_diff, y_diff;
	cx = p.x < rec.xu ? p.x : rec.xu;
	cx = rec.xl > cx ? rec.xl : cx;
	cy = p.y < rec.yu ? p.y : rec.yu;
	cy = rec.yl > cy ? rec.yl : cy;
	x_diff = p.x - cx;
	y_diff = p.y - cy;
	return sqrt(x_diff * x_diff + y_diff * y_diff);
}

template <typename Node3D, typename Point3D>
inline double cal_tetrahedron_vol(Node3D &n1, Node3D &n2, Node3D &n3, Point3D &p4)
{
	double v21_x, v21_y, v21_z;
	double v31_x, v31_y, v31_z;
	double v41_x, v41_y, v41_z;
	v21_x = n2.x - n1.x;
	v21_y = n2.y - n1.y;
	v21_z = n2.z - n1.z;
	v31_x = n3.x - n1.x;
	v31_y = n3.y - n1.y;
	v31_z = n3.z - n1.z;
	v41_x = p4.x - n1.x;
	v41_y = p4.y - n1.y;
	v41_z = p4.z - n1.z;
	double cp_x, cp_y, cp_z;
	cp_x = v21_y * v31_z - v31_y * v21_z;
	cp_y = v21_z * v31_x - v31_z * v21_x;
	cp_z = v21_x * v31_y - v31_x * v21_y;
	return (v41_x * cp_x + v41_y * cp_y + v41_z * cp_z) / 6.0;
}

#endif