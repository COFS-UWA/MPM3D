#include "SimulationsOMP_pcp.h"

#include "RigidCone.h"

void RigidCone::init(
	double _x,
	double _y,
	double _z,
	double _r,
	double _tip_h,
	double _shaft_h
	)
{
	x = _x;
	y = _y;
	z = _z;
	r = _r;
	tip_h = _tip_h;
	shaft_h = _shaft_h;
}
