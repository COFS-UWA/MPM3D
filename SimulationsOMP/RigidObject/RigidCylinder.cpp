#include "SimulationsOMP_pcp.h"

#include "RigidCylinder.h"

RigidCylinder::RigidCylinder() : 
	vx(0.0), vy(0.0), vz(0.0),
	fx_cont(0.0), fy_cont(0.0), fz_cont(0.0),
	mx_cont(0.0), my_cont(0.0), mz_cont(0.0)
{
	norm[0][0] = 1.0;
	norm[0][1] = 0.0;
	norm[0][2] = 0.0;
	norm[1][0] = 0.0;
	norm[1][1] = 1.0;
	norm[1][0] = 0.0;
	norm[2][2] = 0.0;
}

RigidCylinder::~RigidCylinder()
{

}

void RigidCylinder::init(
	double _x,
	double _y,
	double _z,
	double _h,
	double _r
	) noexcept
{
	x = _x;
	y = _y;
	z = _z;
	h = _h;
	r = _r;
	r2 = r * r;
	h_div_2 = 0.5 * h;
	bbox.xl = x - r;
	bbox.xu = x + r;
	bbox.yl = y - r;
	bbox.yu = y + r;
	bbox.zl = z - h_div_2;
	bbox.zu = z + h_div_2;
}
