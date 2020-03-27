#ifndef __Geometry_h__
#define __Geometry_h__

#include <cmath>

struct Point2D { double x, y; };
struct Point3D { double x, y, z; };
struct Rect { double xl, xu, yl, yu; };
struct Cube { double xl, xu, yl, yu, zl, zu; };

inline double distance(Point2D &p1, Point2D &p2) noexcept
{
	double x_diff = p1.x - p2.x;
	double y_diff = p1.y - p2.y;
	return sqrt(x_diff * x_diff + y_diff * y_diff);
}

inline double distance(Rect &rec, Point2D &p) noexcept
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

#endif