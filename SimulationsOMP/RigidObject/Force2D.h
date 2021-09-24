#ifndef __Force_2D_h__
#define __Force_2D_h__

struct Force2D
{
	double fx, fy, m;

	inline void reset() noexcept { fx = 0.0; fy = 0.0; m = 0.0; }

	inline void combine(const Force2D &other) noexcept
	{ fx += other.fx; fy += other.fy; m += other.m; }

	inline void add_force(
		double p_x,
		double p_y,
		double _fx,
		double _fy,
		double cen_x,
		double cen_y
		) noexcept
	{
		fx += _fx; fy += _fy;
		double dx = p_x - cen_x;
		double dy = p_y - cen_y;
		m += dx * _fy - dy * _fx;
	}

	inline Force2D& operator=(const Force2D& other)
	{
		fx = other.fx; fy = other.fy; m = other.m;
		return *this;
	}

	inline Force2D& operator+= (const Force2D& other) noexcept
	{ fx += other.fx; fy += other.fy; m += other.m; return *this; }
};

#endif