#ifndef __Geometry_2D_h__
#define __Geometry_2D_h__

#include <cmath>

struct Point2D
{
	union
	{
		struct { double x, y; };
		double data[2];
	};

	inline Point2D() {}
	inline Point2D(double _x, double _y) : x(_x), y(_y) {}
	inline Point2D(const Point2D &other) : x(other.x), y(other.y) {}
};

struct Rect
{
	double xl, xu, yl, yu;
	
	inline Rect() {}
	inline Rect(double _xl, double _xu, double _yl, double _yu) :
		xl(_xl), xu(_xu), yl(_yl), yu(_yu) {}
	inline Rect(const Rect& other) :
		xl(other.xl), xu(other.xu), yl(other.yl), yu(other.yu) {}

	inline bool is_in_box(double x, double y) const noexcept
	{ return !(x < xl || x > xu || y < yl || y > yu); }
	template <typename Point2D>
	inline bool is_in_box(Point2D& point) const noexcept
	{ return is_in_box(point.x, point.y); }

	inline void envelop(const Rect& other) noexcept
	{
		if (xl > other.xl)
			xl = other.xl;
		if (xu < other.xu)
			xu = other.xu;
		if (yl > other.yl)
			yl = other.yl;
		if (yu < other.yu)
			yu = other.yu;
	}
	inline void envelop(
		const Rect &rect1,
		const Rect &rect2
		) noexcept
	{
		xl = rect1.xl < rect2.xl ? rect1.xl : rect2.xl;
		xu = rect1.xu > rect2.xu ? rect1.xu : rect2.xu;
		yl = rect1.yl < rect2.yl ? rect1.yl : rect2.yl;
		yu = rect1.yu > rect2.yu ? rect1.yu : rect2.yu;
	}
};

struct Vector2D
{
	double x, y;

	inline Vector2D() {}
	inline Vector2D(double _x, double _y) : x(_x), y(_y) {}
	inline Vector2D(const Vector2D &other) : x(other.x), y(other.y) {}

	inline double norm() const noexcept { return sqrt(x*x + y*y); }
	inline Vector2D &normalize() noexcept
	{
		double len = norm();
		if (len != 0.0)
		{
			x /= len;
			y /= len;
		}
		return *this;
	}

	inline Vector2D &reverse() noexcept
	{ x = -x; y = -y; return *this; }
	inline Vector2D &scale(double fac) noexcept
	{ x *= fac; y *= fac; return *this; }
	inline Vector2D& add(double e_x, double e_y) noexcept
	{ x += e_x; y += e_y; return *this; }
	inline Vector2D& add(double e1_x, double e1_y,
						 double e2_x, double e2_y) noexcept
	{ x = e1_x + e2_x; y = e1_y + e2_y; return *this; }
	inline Vector2D& substract(double e_x, double e_y) noexcept
	{ x -= e_x; y -= e_y; return *this; }
	inline Vector2D& substract(double e1_x, double e1_y,
							   double e2_x, double e2_y) noexcept
	{ x = e1_x - e2_x; y = e1_y - e2_y; return *this; }
	inline double dot(double e2_x, double e2_y) const noexcept
	{ return x * e2_x + y * e2_y; }

	template <typename Point2D>
	inline Vector2D& add(const Point2D& p) noexcept
	{ x += p.x; y += p.y; return *this; }
	template <typename Point2D>
	inline Vector2D& add(const Point2D& p1, const Point2D& p2) noexcept
	{ x = p1.x + p2.x; y = p1.y + p2.y; return *this; }
	template <typename Point2D>
	inline Vector2D& substract(const Point2D& p) noexcept
	{ x -= p.x; y -= p.y; return *this; }
	template <typename Point2D>
	inline Vector2D& substract(const Point2D& p1, const Point2D& p2) noexcept
	{ x = p1.x - p2.x; y = p1.y - p2.y; return *this; }
	template <typename Point2D>
	inline double dot(const Point2D& p2) const noexcept
	{ return x * p2.x + y * p2.y; }
};

// distance between 2D points
template <typename Node2D, typename Point2D>
inline double cal_point2d_distance(const Node2D &p1, const Point2D &p2) noexcept
{
	double dx = p1.x - p2.x;
	double dy = p1.y - p2.y;
	return sqrt(dx * dx + dy * dy);
}

// distance between rectangle and point
inline double cal_rect_point_distance(const Rect& rec, const Point2D& p) noexcept
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

inline bool detect_rect_collision(const Rect &c1, const Rect &c2) noexcept
{
	return !(c1.xu < c2.xl || c1.xl > c2.xu ||
			 c1.yu < c2.yl || c1.yl > c2.yu);
}

struct OBB2D
{
	double xo, yo;
	double xlen, ylen;
	Vector2D ix, iy;
};

struct OBB2DAABBCollisionSAT
{
	double xhlen, yhlen;
	Vector2D ix, iy;
	double ix_min, ix_max;
	double iy_min, iy_max;
	double it1x_min, it1x_max;
	double it1y_min, it1y_max;
	double it2x_min, it2x_max;
	double it2y_min, it2y_max;

	inline void init_obb2d(const OBB2D& obb) noexcept
	{
		xhlen = obb.xlen * 0.5;
		yhlen = obb.ylen * 0.5;
		ix = obb.ix;
		iy = obb.iy;
		double cen;
		// ix
		cen = ix.x * obb.xo + ix.y * obb.yo;
		ix_min = cen - xhlen;
		ix_max = cen + xhlen;
		// iy
		cen = iy.x * obb.xo + iy.y * obb.yo;
		iy_min = cen - yhlen;
		iy_max = cen + yhlen;
	}
	inline bool detect_collision_with_rect(const Rect& rect) const noexcept
	{
		double rect_xm = (rect.xl + rect.xu) * 0.5;
		double rect_ym = (rect.yl + rect.yu) * 0.5;
		double rect_xhlen = (rect.xu - rect.xl) * 0.5;
		double rect_yhlen = (rect.yu - rect.yl) * 0.5;
		double ix_rect_cen = ix.x * rect_xm + ix.y * rect_ym;
		double ix_rect_range = abs(ix.x * rect_xhlen) + abs(ix.y * rect_yhlen);
		double iy_rect_cen = iy.x * rect_xm + iy.y * rect_ym;
		double iy_rect_range = abs(iy.x * rect_xhlen) + abs(iy.y * rect_yhlen);
		return !(is_seperating_axis(ix_rect_range, ix_rect_cen, ix_min, ix_max) ||
				 is_seperating_axis(iy_rect_range, iy_rect_cen, iy_min, iy_max));
	}

protected:
	inline bool is_seperating_axis(
		double rect_range,
		double rect_cen,
		double obb_min,
		double obb_max
		) const noexcept
	{
		return ((obb_min - rect_cen) >  rect_range ||
				(obb_max + rect_cen) < -rect_range);
	}
};

inline bool detect_obb2d_cube_collision(
	const OBB2D &obb,
	const Rect &rect
	) noexcept
{
	double rect_xl = rect.xl - obb.xo;
	double rect_xu = rect.xu - obb.xo;
	double rect_yl = rect.yl - obb.yo;
	double rect_yu = rect.yu - obb.yo;
	double obb_hxlen = obb.xlen * 0.5;
	double obb_hylen = obb.ylen * 0.5;

	double tmp1, tmp2, proj_min, proj_max;
	// ix as seperating axis
	const Vector2D &obb_ix = obb.ix;
	tmp1 = rect_xl * obb_ix.x;
	tmp2 = rect_xu * obb_ix.x;
	if (tmp1 > tmp2)
	{
		proj_max = tmp1;
		proj_min = tmp2;
	}
	else
	{
		proj_max = tmp2;
		proj_min = tmp1;
	}
	tmp1 = rect_yl * obb_ix.y;
	tmp2 = rect_yu * obb_ix.y;
	if (tmp1 > tmp2)
	{
		proj_max += tmp1;
		proj_min += tmp2;
	}
	else
	{
		proj_max += tmp2;
		proj_min += tmp1;
	}
	// whether ix is the seperating axis
	if (proj_max < -obb_hxlen || proj_min > obb_hxlen)
		return false;

	// iy as seperating axis
	const Vector2D& obb_iy = obb.iy;
	tmp1 = rect_xl * obb_iy.x;
	tmp2 = rect_xu * obb_iy.x;
	if (tmp1 > tmp2)
	{
		proj_max = tmp1;
		proj_min = tmp2;
	}
	else
	{
		proj_max = tmp2;
		proj_min = tmp1;
	}
	tmp1 = rect_yl * obb_iy.y;
	tmp2 = rect_yu * obb_iy.y;
	if (tmp1 > tmp2)
	{
		proj_max += tmp1;
		proj_min += tmp2;
	}
	else
	{
		proj_max += tmp2;
		proj_min += tmp1;
	}
	// whether iy is the seperating axis
	if (proj_max < -obb_hylen || proj_min > obb_hylen)
		return false;

	return true;
}

template <typename Point2DType1, typename Point2DType2>
inline void point_from_global_to_local_coordinate(
	const Point2D &loc_cen,
	const Vector2D &loc_ix,
	const Vector2D &loc_iy,
	const Point2DType1& gp,
	Point2DType2& lp
	) noexcept
{
	double dx = gp.x - loc_cen.x;
	double dy = gp.y - loc_cen.y;
	lp.x = loc_ix.x * dx + loc_ix.y * dy;
	lp.y = loc_iy.x * dx + loc_iy.y * dy;
}

template <typename Point2DType1, typename Point2DType2>
inline void point_from_local_to_global_coordinate(
	const Point2D& loc_cen,
	const Vector2D& loc_ix,
	const Vector2D& loc_iy,
	const Point2DType1& lp,
	Point2DType2& gp
	) noexcept
{
	gp.x = loc_ix.x * lp.x + loc_iy.x * lp.y + loc_cen.x;
	gp.y = loc_ix.y * lp.x + loc_iy.y * lp.y + loc_cen.y;
}

template <typename Vector2DType1, typename Vector2DType2>
inline void vector_from_global_to_local_coordinate(
	const Vector2D& loc_ix,
	const Vector2D& loc_iy,
	const Vector2DType1& gp,
	Vector2DType2& lp
	) noexcept
{
	lp.x = loc_ix.x * gp.x + loc_ix.y * gp.y;
	lp.y = loc_iy.x * gp.x + loc_iy.y * gp.y;
}

template <typename Vector2DType1, typename Vector2DType2>
inline void vector_from_local_to_global_coordinate(
	const Vector2D& loc_ix,
	const Vector2D& loc_iy,
	const Vector2DType1& lp,
	Vector2DType2& gp
	) noexcept
{
	gp.x = loc_ix.x * lp.x + loc_iy.x * lp.y;
	gp.y = loc_ix.y * lp.x + loc_iy.y * lp.y;
}

template <typename Point2DType1, typename Point2DType2>
inline void point_from_global_to_local_coordinate(
	const Point2D& loc_cen,
	const double angle,
	const Point2DType1& gp,
	Point2DType2& lp
	) noexcept
{
	double dx = gp.x - loc_cen.x;
	double dy = gp.y - loc_cen.y;
	double sin_ang = sin(angle);
	double cos_ang = cos(angle);
	lp.x =  cos_ang * dx + sin_ang * dy;
	lp.y = -sin_ang * dx + cos_ang * dy;
}

template <typename Point2DType1, typename Point2DType2>
inline void point_from_local_to_global_coordinate(
	const Point2D& loc_cen,
	const double angle,
	const Point2DType1& lp,
	Point2DType2& gp
	) noexcept
{
	double sin_ang = sin(angle);
	double cos_ang = cos(angle);
	gp.x = cos_ang * lp.x - sin_ang * lp.y + loc_cen.x;
	gp.y = sin_ang * lp.x + cos_ang * lp.y + loc_cen.y;
}

template <typename Vector2DType1, typename Vector2DType2>
inline void vector_from_global_to_local_coordinate(
	const double angle,
	const Vector2DType1& gp,
	Vector2DType2& lp
	) noexcept
{
	double sin_ang = sin(angle);
	double cos_ang = cos(angle);
	lp.x =  cos_ang * gp.x + sin_ang * gp.y;
	lp.y = -sin_ang * gp.x + cos_ang * gp.y;
}

template <typename Vector2DType1, typename Vector2DType2>
inline void vector_from_local_to_global_coordinate(
	const double angle,
	const Vector2DType1& lp,
	Vector2DType2& gp
	) noexcept
{
	double sin_ang = sin(angle);
	double cos_ang = cos(angle);
	gp.x = cos_ang * lp.x - sin_ang * lp.y;
	gp.y = sin_ang * lp.x + cos_ang * lp.y;
}

#endif