#ifndef __Geometry_3D_h__
#define __Geometry_3D_h__

#include "Geometry.h"

inline bool detect_cube_collision(Cube& c1, Cube& c2) noexcept
{
	return !(c1.xu < c2.xl || c1.xl > c2.xu ||
			 c1.yu < c2.yl || c1.yl > c2.yu ||
			 c1.zu < c2.zl || c1.zl > c2.zu);
}

inline double cal_cube_point_distance(Cube &box, Point3D &p) noexcept
{
	double cx, cy, cz, x_diff, y_diff, z_diff;
	cx = p.x < box.xu ? p.x : box.xu;
	cx = box.xl > cx ? box.xl : cx;
	cy = p.y < box.yu ? p.y : box.yu;
	cy = box.yl > cy ? box.yl : cy;
	cz = p.z < box.zu ? p.z : box.zu;
	cz = box.zl > cz ? box.zl : cz;
	x_diff = p.x - cx;
	y_diff = p.y - cy;
	z_diff = p.z - cz;
	return sqrt(x_diff * x_diff + y_diff * y_diff + z_diff * z_diff);
}

#endif