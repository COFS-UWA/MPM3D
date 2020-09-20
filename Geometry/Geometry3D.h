#ifndef __Geometry_3D_h__
#define __Geometry_3D_h__

#include <cmath>

struct Point3D
{
	union
	{
		struct { double x, y, z; };
		double data[3];
	};

	inline Point3D() {}
	inline Point3D(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
	inline Point3D(const Point3D& other) : x(other.x), y(other.y), z(other.z) {}
};

struct Vector3D
{
	union
	{
		struct { double x, y, z; };
		double data[3];
	};

	inline Vector3D() {}
	inline Vector3D(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
	inline Vector3D(const Vector3D& other) : x(other.x), y(other.y), z(other.z) {}

	inline double norm() const noexcept { return sqrt(x * x + y * y + z * z); }
	inline Vector3D& normalize()
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

	inline Vector3D& reverse() noexcept
	{
		x = -x; y = -y; z = -z; return *this;
	}
	inline Vector3D& scale(double fac) noexcept
	{
		x *= fac; y *= fac; z *= fac; return *this;
	}
	inline Vector3D& add(double e1_x, double e1_y, double e1_z)
	{
		x += e1_x; y += e1_y; z += e1_z; return *this;
	}
	inline Vector3D& add(double e1_x, double e1_y, double e1_z,
		double e2_x, double e2_y, double e2_z)
	{
		x = e1_x + e2_x; y = e1_y + e2_x; z = e1_z + e2_z; return *this;
	}
	inline Vector3D& substract(double e1_x, double e1_y, double e1_z)
	{
		x -= e1_x; y -= e1_y; z -= e1_z; return *this;
	}
	inline Vector3D& substract(double e1_x, double e1_y, double e1_z,
		double e2_x, double e2_y, double e2_z)
	{
		x = e1_x - e2_x; y = e1_y - e2_y; z = e1_z - e2_z; return *this;
	}
	inline double dot(double e2_x, double e2_y, double e2_z)
	{
		return x * e2_x + y * e2_y + z * e2_z;
	}
	inline Vector3D& cross(double e1_x, double e1_y, double e1_z,
		double e2_x, double e2_y, double e2_z)
	{
		double _x = e1_y * e2_z - e1_z * e2_y;
		double _y = e1_z * e2_x - e1_x * e2_z;
		double _z = e1_x * e2_y - e1_y * e2_x;
		x = _x;
		y = _y;
		z = _z;
		return *this;
	}

	template <typename Point3D>
	inline Vector3D& add(Point3D& p)
	{
		x += p.x; y += p.y; z += p.z; return *this;
	}
	template <typename Point3D>
	inline Vector3D& add(Point3D& p1, Point3D& p2)
	{
		x = p1.x + p2.x; y = p1.y + p2.y; z = p1.z + p2.z; return *this;
	}
	template <typename Point3D>
	inline Vector3D& substract(Point3D& p)
	{
		x -= p.x; y -= p.y; z -= p.z; return *this;
	}
	template <typename Point3D>
	inline Vector3D& substract(Point3D& p1, Point3D& p2)
	{
		x = p1.x - p2.x; y = p1.y - p2.y; z = p1.z - p2.z; return *this;
	}
	template <typename Point3D>
	inline double dot(Point3D& p2)
	{
		return x * p2.x + y * p2.y + z * p2.z;
	}
	template <typename Point3D>
	inline Vector3D& cross(Point3D& p1, Point3D& p2)
	{
		double _x = p1.y * p2.z - p1.z * p2.y;
		double _y = p1.z * p2.x - p1.x * p2.z;
		double _z = p1.x * p2.y - p1.y * p2.x;
		x = _x;
		y = _y;
		z = _z;
		return *this;
	}
};


struct Cube
{
	double xl, xu, yl, yu, zl, zu;

	inline Cube() {}
	inline Cube(double _xl, double _xu, double _yl, double _yu, double _zl, double _zu) :
		xl(_xl), xu(_xu), yl(_yl), yu(_yu), zl(_zl), zu(_zu) {}
	inline Cube(const Cube& other) :
		xl(other.xl), xu(other.xu), yl(other.yl), yu(other.yu), zl(other.zl), zu(other.zu) {}
	
	inline bool is_in_box(double x, double y, double z) const noexcept
	{
		return !(x < xl || x > xu || y < yl || y > yu || z < zl || z > zu);
	}
	template <typename Point3D>
	inline bool is_in_box(Point3D& point) const noexcept
	{
		return is_in_box(point.x, point.y, point.z);
	}
};

// distance between 3D points
template <typename Node3D, typename Point3D>
inline double cal_distance_3D(Node3D& n, Point3D& p) noexcept
{
	double dx = n.x - p.x;
	double dy = n.y - p.y;
	double dz = n.z - p.z;
	return sqrt(dx * dx + dy * dy + dz * dz);
}


struct IdCube
{
	long long xl_id, xu_id;
	long long yl_id, yu_id;
	long long zl_id, zu_id;

	inline IdCube() {}
	inline IdCube(long long _xl_id, long long _xu_id,
		long long _yl_id, long long _yu_id,
		long long _zl_id, long long _zu_id) :
		xl_id(_xl_id), xu_id(_xu_id), yl_id(_yl_id), yu_id(_yu_id),
		zl_id(_zl_id), zu_id(_zu_id) {}

	inline bool does_not_overlap(const IdCube& ic) const noexcept
	{
		return xu_id <= ic.xl_id || xl_id >= ic.xu_id ||
			yu_id <= ic.yl_id || yl_id >= ic.yu_id ||
			zu_id <= ic.zl_id || zl_id >= ic.zu_id;
	}
	inline void trim_by(const IdCube& id)
	{
		if (xl_id < id.xl_id)
			xl_id = id.xl_id;
		if (xu_id > id.xu_id)
			xu_id = id.xu_id;
		if (yl_id < id.yl_id)
			yl_id = id.yl_id;
		if (yu_id > id.yu_id)
			yu_id = id.yu_id;
		if (zl_id < id.zl_id)
			zl_id = id.zl_id;
		if (zu_id > id.zu_id)
			zu_id = id.zu_id;
	}
	inline void from_cube(const Cube& c,
		double xl, double yl, double zl,
		double hx, double hy, double hz, double tol = 0.0)
	{
		xl_id = long long(floor((c.xl - tol - xl) / hx));
		xu_id = long long(ceil((c.xu + tol - xl) / hx));
		yl_id = long long(floor((c.yl - tol - yl) / hy));
		yu_id = long long(ceil((c.yu + tol - yl) / hy));
		zl_id = long long(floor((c.zl - tol - zl) / hz));
		zu_id = long long(ceil((c.zu + tol - zl) / hz));
	}
};

inline bool detect_cube_collision(Cube& c1, Cube& c2) noexcept
{
	return !(c1.xu < c2.xl || c1.xl > c2.xu ||
		c1.yu < c2.yl || c1.yl > c2.yu ||
		c1.zu < c2.zl || c1.zl > c2.zu);
}

inline double cal_cube_point_distance(Cube& box, Point3D& p) noexcept
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

struct OBB3D
{
	double xo, yo, zo;
	double xlen, ylen, zlen;
	Vector3D ix, iy, iz;
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
#define Norm_Tol (1.0e-10)
		xhlen = obb.xlen * 0.5;
		yhlen = obb.ylen * 0.5;
		zhlen = obb.zlen * 0.5;
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
		if ((ix.z * ix.z + ix.y * ix.y) > (Norm_Tol * Norm_Tol))
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
		if ((iy.z * iy.z + iy.y * iy.y) > (Norm_Tol * Norm_Tol))
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
		if ((iz.z * iz.z + iz.y * iz.y) > (Norm_Tol * Norm_Tol))
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
		if ((ix.z * ix.z + ix.x * ix.x) > (Norm_Tol * Norm_Tol))
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
		if ((iy.z * iy.z + iy.x * iy.x) > (Norm_Tol * Norm_Tol))
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
		if ((iz.z * iz.z + iz.x * iz.x) > (Norm_Tol * Norm_Tol))
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
		if ((ix.y * ix.y + ix.x * ix.x) > (Norm_Tol * Norm_Tol))
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
		if ((iy.y * iy.y + iy.x * iy.x) > (Norm_Tol * Norm_Tol))
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
		if ((iz.y * iz.y + iz.x * iz.x) > (Norm_Tol * Norm_Tol))
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
#undef Norm_Tol
	}
	inline bool detect_collision_with_cube(Cube& cube)
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

template <typename Point3DType1, typename Point3DType2>
inline void point_from_global_to_local_coordinate(
	const Point3D& loc_cen,
	const Vector3D& loc_ix,
	const Vector3D& loc_iy,
	const Vector3D& loc_iz,
	const Point3DType1& gp,
	Point3DType2& lp
) noexcept
{
	double dx = gp.x - loc_cen.x;
	double dy = gp.y - loc_cen.y;
	double dz = gp.z - loc_cen.z;
	lp.x = loc_ix.x * dx + loc_ix.y * dy + loc_ix.z * dz;
	lp.y = loc_iy.x * dx + loc_iy.y * dy + loc_iy.z * dz;
	lp.z = loc_iz.x * dx + loc_iz.y * dy + loc_iz.z * dz;
}

template <typename Point3DType1, typename Point3DType2>
inline void point_from_local_to_global_coordinate(
	const Point3D& loc_cen,
	const Vector3D& loc_ix,
	const Vector3D& loc_iy,
	const Vector3D& loc_iz,
	const Point3DType1& lp,
	Point3DType2& gp
) noexcept
{
	gp.x = loc_ix.x * lp.x + loc_iy.x * lp.y + loc_iz.x * lp.z + loc_cen.x;
	gp.y = loc_ix.y * lp.x + loc_iy.y * lp.y + loc_iz.y * lp.z + loc_cen.y;
	gp.z = loc_ix.z * lp.x + loc_iy.z * lp.y + loc_iz.z * lp.z + loc_cen.z;
}

template <typename Vector3DType1, typename Vector3DType2>
inline void vector_from_global_to_local_coordinate(
	const Vector3D& loc_ix,
	const Vector3D& loc_iy,
	const Vector3D& loc_iz,
	const Vector3DType1& gp,
	Vector3DType2& lp
	) noexcept
{
	lp.x = loc_ix.x * gp.x + loc_ix.y * gp.y + loc_ix.z * gp.z;
	lp.y = loc_iy.x * gp.x + loc_iy.y * gp.y + loc_iy.z * gp.z;
	lp.z = loc_iz.x * gp.x + loc_iz.y * gp.y + loc_iz.z * gp.z;
}

template <typename Vector3DType1, typename Vector3DType2>
inline void vector_from_local_to_global_coordinate(
	const Vector3D& loc_ix,
	const Vector3D& loc_iy,
	const Vector3D& loc_iz,
	const Vector3DType1& lp,
	Vector3DType2& gp
	) noexcept
{
	gp.x = loc_ix.x * lp.x + loc_iy.x * lp.y + loc_iz.x * lp.z;
	gp.y = loc_ix.y * lp.x + loc_iy.y * lp.y + loc_iz.y * lp.z;
	gp.z = loc_ix.z * lp.x + loc_iy.z * lp.y + loc_iz.z * lp.z;
}

// rotate coordinates ix, iy, iz by ang
// pure geometric operations
//namespace Geometry3D_Internal
//{
//	inline void rotate_axis(
//		Vector3D& rix,
//		Vector3D& riy,
//		Vector3D& riz,
//		double sin_ang,
//		double cos_ang,
//		Vector3D& axis
//	)
//	{
//		riy.cross(riz, axis);
//		double riy_norm = riy.norm();
//		if (riy_norm != 0.0)
//		{
//			riy.scale(1.0 / riy_norm);
//			rix.cross(riy, riz);
//			double pj_len = axis.dot(rix);
//			double pj_ht = axis.dot(riz);
//			axis.x = (rix.x * cos_ang + riy.x * sin_ang) * pj_len + riz.x * pj_ht;
//			axis.y = (rix.y * cos_ang + riy.y * sin_ang) * pj_len + riz.y * pj_ht;
//			axis.z = (rix.z * cos_ang + riy.z * sin_ang) * pj_len + riz.z * pj_ht;
//			//axis.add(riz.x * pj_ht, riz.y * pj_ht, riz.z * pj_ht);
//		}
//	}
//}
//
//inline void rotate_axses_by_angle(
//	const Vector3D& ang,
//	Vector3D& ix,
//	Vector3D& iy,
//	Vector3D& iz
//) noexcept
//{
//	Vector3D riz(ang.x, ang.y, ang.z);
//	double riz_norm = riz.norm();
//	if (riz_norm != 0.0)
//	{
//		riz.scale(1.0 / riz_norm);
//		double sin_ang = sin(riz_norm);
//		double cos_ang = cos(riz_norm);
//		Vector3D rix, riy;
//		Geometry3D_Internal::rotate_axis(rix, riy, riz, sin_ang, cos_ang, ix);
//		Geometry3D_Internal::rotate_axis(rix, riy, riz, sin_ang, cos_ang, iy);
//		Geometry3D_Internal::rotate_axis(rix, riy, riz, sin_ang, cos_ang, iz);
//	}
//}

// use quaternion
inline void rotate_axses_by_angle(
	const Vector3D& ang,
	Vector3D& ix,
	Vector3D& iy,
	Vector3D& iz
) noexcept
{
	double theta, Kx, Ky, Kz;
	theta = ang.norm();
	if (theta != 0.0)
	{
		Kx = ang.x / theta;
		Ky = ang.y / theta;
		Kz = ang.z / theta;
	}
	else
	{
		Kx = 0.0;
		Ky = 0.0;
		Kz = 0.0;
	}

	// quaternion
	double q0, q1, q2, q3;
	q0 = cos(0.5 * theta);
	q1 = Kx * sin(0.5 * theta);
	q2 = Ky * sin(0.5 * theta);
	q3 = Kz * sin(0.5 * theta);

	ix.x = 1.0 - 2.0 * (q2 * q2 + q3 * q3);
	ix.y = 2.0 * (q1 * q2 + q0 * q3);
	ix.z = 2.0 * (q1 * q3 - q0 * q2);
	iy.x = 2.0 * (q1 * q2 - q0 * q3);
	iy.y = 1.0 - 2.0 * (q3 * q3 + q1 * q1);
	iy.z = 2.0 * (q2 * q3 + q0 * q1);
	iz.x = 2.0 * (q1 * q3 + q0 * q2);
	iz.y = 2.0 * (q2 * q3 - q0 * q1);
	iz.z = 1.0 - 2.0 * (q1 * q1 + q2 * q2);
}

#endif