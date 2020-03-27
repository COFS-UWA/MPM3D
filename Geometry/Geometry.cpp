#include "Geometry_pcp.h"

#include "Geometry.h"

inline bool clip_test(double p, double q, double &t_min, double &t_max)
{
	bool is_in_win = true;
	double t;
	if (p < 0.0) //line entry point 
	{
		t = q / p;
		// line is outside window 
		if (t > t_max)
			return false;
		else if (t > t_min)
			t_min = t;
	}
	else if (p > 0.0) //line leaving point
	{
		t = q / p;
		// line is outside window   
		if (t < t_min)
			return false;
		else if (t < t_max)
			t_max = t;
	}
	else //line is parallel to window edge
	{
		// line is outside window  
		if (q < 0.0)
			return false;
	}
	return true;
}

bool clip_line(double xl, double xu, double yl, double yu,
			   double &x1, double &y1, double &x2, double &y2)
{
	double t_min = 0.0, t_max = 1.0;
	double dx = x2 - x1;
	if (clip_test(-dx, x1 - x1, t_min, t_max) && // left edge
		clip_test( dx, xu - x1, t_min, t_max))   // right edge
	{
		double dy = y2 - y1;
		if (clip_test(-dy, y1 - yl, t_min, t_max) && // bottom edge
			clip_test( dy, yu - y1, t_min, t_max))   // top edge
			{
				if (t_max != 1.0)
				{
					x2 = x1 + t_max * dx;
					y2 = y1 + t_max * dy;
				}
				if (t_min != 0.0) 
				{
					x1 += t_min*dx;
					y1 += t_min*dy;
				}
				return true;
			}
	}
	return false;
}

inline void swap(double &a, double &b) { double c = a; a = b; b = c; }

bool test_AABB_triangle_intersection(double xl, double xu, double yl, double yu,
	double x0, double y0, double x1, double y1, double x2, double y2)
{
	double hx, hy, xc, yc;
	hx = xu - xl;
	hy = yu - yl;
	xc = (xl + xu) * 0.5;
	yc = (yl + yu) * 0.5;
	// translate triangle
	// take grid centre as origin
	x0 -= xc;
	y0 -= yc;
	x1 -= xc;
	y1 -= yc;
	x2 -= xc;
	y2 -= yc;

	double r, v_min, v_max;
	// a31
	r = (hx * abs(y1 - y0) + hy * abs(x1 - x0)) * 0.5;
	v_min = x1 * y0 - x0 * y1;
	v_max = (y0 - y1) * x2 + (x1 - x0) * y2;
	if (v_min > v_max)
		swap(v_min, v_max);
	if (v_min >= r || v_max <= -r)
		return false;
	// a32
	r = (hx * abs(y2 - y1) + hy * abs(x2 - x1)) * 0.5;
	v_min = x2 * y1 - x1 * y2;
	v_max = (y1 - y2) * x0 + (x2 - x1) * y0;
	if (v_min > v_max)
		swap(v_min, v_max);
	if (v_min >= r || v_max <= -r)
		return false;
	// a33
	r = (hx * abs(y0 - y2) + hy * abs(x0 - x2)) * 0.5;
	v_min = x0 * y2 - x2 * y0;
	v_max = (y2 - y0) * x1 + (x0 - x2) * y1;
	if (v_min > v_max)
		swap(v_min, v_max);
	if (v_min >= r || v_max <= -r)
		return false;

	return true;
}
