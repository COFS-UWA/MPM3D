#include "Simulations_pcp.h"

#include "RigidCircle.h"

int RigidCircle::init(double _r, double _x, double _y, double max_pcl_size)
{
	r = _r;
	r2 = r * r;
	cen_x = _x;
	cen_y = _y;
	theta = 0.0;
	rfx = 0.0;
	rfy = 0.0;
	rm  = 0.0;
	return 0;
}
