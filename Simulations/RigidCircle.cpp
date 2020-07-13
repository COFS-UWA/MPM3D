#include "Simulations_pcp.h"

#include "RigidCircle.h"

void RigidCircle::init(
	double _x, double _y,
	double _r, double _density
	)
{
	x = _x;
	y = _y;
	ang = 0.0;
	r = _r;
	density = _density;
	init_cal_var();
}

void RigidCircle::set_init_state(
	double _r, double _density,
	double _ax, double _ay, double _a_ang,
	double _vx, double _vy, double _v_ang,
	double _x, double _y, double _ang
	)
{
	r = _r;
	density = _density;
	ax = _ax;
	ay = _ay;
	a_ang = _a_ang;
	vx = _vx;
	vy = _vy;
	v_ang = _v_ang;
	x = _x;
	y = _y;
	ang = _ang;
	init_cal_var();
}
