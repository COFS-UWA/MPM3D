#include "Geometry_pcp.h"

#include <math.h>

#include "GeometryUtils.h"

namespace
{

	bool clip_test(double p, double q, double& t_min, double& t_max)
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
		double& x1, double& y1, double& x2, double& y2)
	{
		double t_min = 0.0, t_max = 1.0;
		double dx = x2 - x1;
		if (clip_test(-dx, x1 - x1, t_min, t_max) && // left edge
			clip_test(dx, xu - x1, t_min, t_max))   // right edge
		{
			double dy = y2 - y1;
			if (clip_test(-dy, y1 - yl, t_min, t_max) && // bottom edge
				clip_test(dy, yu - y1, t_min, t_max))   // top edge
			{
				if (t_max != 1.0)
				{
					x2 = x1 + t_max * dx;
					y2 = y1 + t_max * dy;
				}
				if (t_min != 0.0)
				{
					x1 += t_min * dx;
					y1 += t_min * dy;
				}
				return true;
			}
		}
		return false;
	}
}
