#ifndef __Geometry_h__
#define __Geometry_h__

#include <cmath>

struct Point2D { double x, y; };

struct Rect
{
	double xl, xu, yl, yu;
	Rect() {}
	Rect(double _xl, double _xu, double _yl, double _yu) :
		xl(_xl), xu(_xu), yl(_yl), yu(_yu) {}
	inline bool is_in_box(double x, double y)
	{
		if (x < xl || x > xu ||	y < yl || y > yu)
			return false;
		return true;
	}
	template <typename Point2D>
	inline bool is_in_box(Point2D& point)
	{ return is_in_box(point.x, point.y); }
};

struct Vector2D
{
	double x, y;

	inline double norm() { return sqrt(x*x + y*y); }
	inline Vector2D &normalize()
	{
		double len = norm();
		if (len != 0.0)
		{
			x /= len;
			y /= len;
		}
		return *this;
	}

	inline Vector2D& scale(double fac) noexcept
	{ x *= fac; y *= fac; return *this; }
	inline Vector2D& reverse() noexcept
	{ x = -x; y = -y; return *this; }

	inline Vector2D& substract(double e1_x, double e1_y,
							   double e2_x, double e2_y)
	{ x = e1_x - e2_x; y = e1_y - e2_y; return *this; }
	inline double dot(double e2_x, double e2_y)
	{ return x * e2_x + y * e2_y; }

	template <typename Point2D>
	inline Vector2D& substract(Point2D& p1, Point2D& p2)
	{ x = p1.x - p2.x; y = p1.y - p2.y; return *this; }
	template <typename Point2D>
	inline double dot(Point2D& p2)
	{ return x * p2.x + y * p2.y; }
};

struct Point3D { double x, y, z; };

struct Cube
{
	double xl, xu, yl, yu, zl, zu;

	Cube() {}
	Cube(double _xl, double _xu, double _yl, double _yu, double _zl, double _zu) :
		xl(_xl), xu(_xu), yl(_yl), yu(_yu), zl(_zl), zu(_zu) {}

	inline bool is_in_box(double x, double y, double z) const noexcept
	{ return !(x < xl || x > xu || y < yl || y > yu || z < zl || z > zu); }
	template <typename Point3D>
	inline bool is_in_box(Point3D &point) const noexcept
	{ return is_in_box(point.x, point.y, point.z); }
};

struct Vector3D
{
	double x, y, z;

	inline double norm() { return sqrt(x*x + y*y + z*z); }
	inline Vector3D &normalize()
	{
		double len = norm();
		if (len != 0.0)
		{
			x /= len;
			y /= len;
			z /= len;
		}
		return *this;
	}

	inline Vector3D &scale(double fac) noexcept
	{ x *= fac; y *= fac; z *= fac; return *this; }
	inline Vector3D &reverse() noexcept
	{ x = -x; y = -y; z = -z; return *this; }
	
	inline Vector3D& substract(double e1_x, double e1_y, double e1_z,
							   double e2_x, double e2_y, double e2_z)
	{ x = e1_x - e2_x; y = e1_y - e2_y; z = e1_z - e2_z; return *this; }
	inline double dot(double e2_x, double e2_y, double e2_z)
	{ return x * e2_x + y * e2_y + z * e2_z; }
	inline Vector3D& cross(double e1_x, double e1_y, double e1_z,
						   double e2_x, double e2_y, double e2_z)
	{
		x = e1_y * e2_z - e1_z * e2_y;
		y = e1_z * e2_x - e1_x * e2_z;
		z = e1_x * e2_y - e1_y * e2_x;
		return *this;
	}

	template <typename Point3D>
	inline Vector3D &substract(Point3D &p1, Point3D &p2)
	{ x = p1.x - p2.x; y = p1.y - p2.y; z = p1.z - p2.z; return *this; }
	template <typename Point3D>
	inline double dot(Point3D& p2)
	{ return x * p2.x + y * p2.y + z * p2.z; }
	template <typename Point3D>
	inline Vector3D &cross(Point3D& p1, Point3D& p2)
	{
		x = p1.y * p2.z - p1.z * p2.y;
		y = p1.z * p2.x - p1.x * p2.z;
		z = p1.x * p2.y - p1.y * p2.x;
		return *this;
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
inline double cal_distance_3D(Node3D &n, Point3D &p) noexcept
{
	double dx = n.x - p.x;
	double dy = n.y - p.y;
	double dz = n.z - p.z;
	return sqrt(dx*dx + dy*dy + dz*dz);
}

// distance between rectangle and point
inline double cal_rect_point_distance(Rect& rec, Point2D& p) noexcept
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

inline bool detect_rect_collision(Rect &c1, Rect &c2) noexcept
{
	return !(c1.xu < c2.xl || c1.xl > c2.xu ||
			 c1.yu < c2.yl || c1.yl > c2.yu);
}

template <typename Node2D, typename Point2D>
inline double cal_triangle_area(Node2D &p1, Node2D &p2, Point2D &p3)
{ return 0.5 * ((p1.x-p3.x)*(p2.y-p3.y) - (p2.x-p3.x)*(p1.y-p3.y)); }

// limit theta to [-pi, pi]
#define PI 3.14159265359
inline void trim_to_pi(double& theta)
{
	if (theta > PI)
		theta -= (2.0 * PI) * long long((theta + PI) / (2.0 * PI));
	else if (theta < -PI)
		theta -= (2.0 * PI) * long long((theta - PI) / (2.0 * PI));
}

template <typename Item>
inline void sort_array_3_acc(Item ids[3])
{
	Item tmp, min_id;
	min_id = 0;
	if (ids[1] < ids[0])
		min_id = 1;
	if (ids[2] < ids[min_id])
		min_id = 2;
	if (min_id != 0)
	{
		tmp = ids[0];
		ids[0] = ids[min_id];
		ids[min_id] = tmp;
	}
	if (ids[2] < ids[1])
	{
		tmp = ids[1];
		ids[1] = ids[2];
		ids[2] = tmp;
	}
}

template <typename Item>
inline void sort_array_4_acc(Item ids[4])
{
	Item tmp;
	size_t min_id;
	for (size_t i = 0; i < 3; ++i)
	{
		min_id = i;
		for (size_t j = i + 1; j < 4; ++j)
		{
			if (ids[j] < ids[min_id])
				min_id = j;
		}
		if (min_id != i)
		{
			tmp = ids[min_id];
			ids[min_id] = ids[i];
			ids[i] = tmp;
		}
	}
}

#undef PI
#endif