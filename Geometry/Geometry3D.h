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
	inline Vector3D& add(const Point3D& p)
	{
		x += p.x; y += p.y; z += p.z; return *this;
	}
	template <typename Point3D>
	inline Vector3D& add(const Point3D& p1, const Point3D& p2)
	{
		x = p1.x + p2.x; y = p1.y + p2.y; z = p1.z + p2.z; return *this;
	}
	template <typename Point3D>
	inline Vector3D& substract(const Point3D& p)
	{
		x -= p.x; y -= p.y; z -= p.z; return *this;
	}
	template <typename Point3D>
	inline Vector3D& substract(const Point3D& p1, const Point3D& p2)
	{
		x = p1.x - p2.x; y = p1.y - p2.y; z = p1.z - p2.z; return *this;
	}
	template <typename Point3D>
	inline double dot(const Point3D& p2)
	{
		return x * p2.x + y * p2.y + z * p2.z;
	}
	template <typename Point3D>
	inline Vector3D& cross(const Point3D& p1, const Point3D& p2)
	{
		double _x = p1.y * p2.z - p1.z * p2.y;
		double _y = p1.z * p2.x - p1.x * p2.z;
		double _z = p1.x * p2.y - p1.y * p2.x;
		x = _x;	y = _y; z = _z;
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

	inline void envelop(const Cube &other)
	{
		if (xl > other.xl)
			xl = other.xl;
		if (xu < other.xu)
			xu = other.xu;
		if (yl > other.yl)
			yl = other.yl;
		if (yu < other.yu)
			yu = other.yu;
		if (zl > other.zl)
			zl = other.zl;
		if (zu < other.zu)
			zu = other.zu;
	}
	inline void envelop(const Cube &cube1, const Cube &cube2)
	{
		xl = cube1.xl < cube2.xl ? cube1.xl : cube2.xl;
		xu = cube1.xu > cube2.xu ? cube1.xu : cube2.xu;
		yl = cube1.yl < cube2.yl ? cube1.yl : cube2.yl;
		yu = cube1.yu > cube2.yu ? cube1.yu : cube2.yu;
		zl = cube1.zl < cube2.zl ? cube1.zl : cube2.zl;
		zu = cube1.zu > cube2.zu ? cube1.zu : cube2.zu;
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

inline double cal_cube_point_distance(const Cube& box, const Point3D& p) noexcept
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
	double hx, hy, hz;
	Vector3D ix, iy, iz;
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

// use quaternion
inline void rotate_axses_by_angle(
	const Vector3D& ang,
	Vector3D& ix,
	Vector3D& iy,
	Vector3D& iz
) noexcept
{
	double Kx, Ky, Kz;
	const double theta = ang.norm();
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
	const double q0 = cos(0.5 * theta);
	const double q1 = Kx * sin(0.5 * theta);
	const double q2 = Ky * sin(0.5 * theta);
	const double q3 = Kz * sin(0.5 * theta);

	ix.x = 1.0 - 2.0 * (q2 * q2 + q3 * q3);
	ix.y = 2.0 * (q1 * q2 - q0 * q3);
	ix.z = 2.0 * (q1 * q3 + q0 * q2);
	iy.x = 2.0 * (q1 * q2 + q0 * q3);
	iy.y = 1.0 - 2.0 * (q3 * q3 + q1 * q1);
	iy.z = 2.0 * (q2 * q3 - q0 * q1);
	iz.x = 2.0 * (q1 * q3 - q0 * q2);
	iz.y = 2.0 * (q2 * q3 + q0 * q1);
	iz.z = 1.0 - 2.0 * (q1 * q1 + q2 * q2);
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

#endif