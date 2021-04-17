#include "Simulations_pcp.h"

#include "RigidRect.h"

void RigidRect::init(
	double _x, double _y,
	double _hx, double _hy,
	double _density)
{
	hx = _hx;
	hy = _hy;
	density = _density;
	x = _x;
	y = _y;
	ang = 0.0;
	init_cal_var();
}

void RigidRect::set_init_state(
	double _hx, double _hy, double _density,
	double _fx_cont, double _fy_cont, double _m_cont,
	double _ax, double _ay, double _a_ang,
	double _vx, double _vy, double _v_ang,
	double _x, double _y, double _ang)
{
	hx = _hx;
	hy = _hy;
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
	fx_con = _fx_cont;
	fy_con = _fy_cont;
	m_con = _m_cont;
	init_cal_var();
}

const Vector2D RigidRect::edge_normals[4] = {
	{ -1.0, 0.0 },
	{  1.0, 0.0 },
	{ 0.0, -1.0 },
	{ 0.0,  1.0 }
};
