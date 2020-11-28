#include "SimulationsOMP_pcp.h"

#include "RigidCube.h"

const Vector3D RigidCube::norms[6] = {
	{ 1.0, 0.0, 0.0 },
	{ -1.0, 0.0, 0.0 },
	{ 0.0, 1.0, 0.0 },
	{ 0.0, -1.0, 0.0 },
	{ 0.0, 0.0, 1.0 },
	{ 0.0, 0.0, -1.0 },
};

RigidCube::RigidCube() :
	ax(0.0), ay(0.0), az(0.0),
	vx(0.0), vy(0.0), vz(0.0),
	fx_cont(0.0), fy_cont(0.0), fz_cont(0.0),
	mx_cont(0.0), my_cont(0.0), mz_cont(0.0),
	vx_bc(0.0), vy_bc(0.0), vz_bc(0.0),
	has_vx_bc(false), has_vy_bc(false), has_vz_bc(false),
	fx_ext(0.0), fy_ext(0.0), fz_ext(0.0),
	mx_ext(0.0), my_ext(0.0), mz_ext(0.0) {}

RigidCube::~RigidCube() {}

void RigidCube::init(
	double _x,
	double _y,
	double _z,
	double _hx,
	double _hy,
	double _hz,
	double _density
	) noexcept
{
	x = _x;
	y = _y;
	z = _z;
	hx = _hx;
	hy = _hy;
	hz = _hz;
	density = _density;
	m = hx * hy * hz * density;
	hhx = 0.5 * hx;
	hhy = 0.5 * hy;
	hhz = 0.5 * hz;
}

void RigidCube::set_cont_force(
	double fx,
	double fy,
	double fz,
	double mx,
	double my,
	double mz
) noexcept
{
	fx_cont = fx; fy_cont = fy; fz_cont = fz;
	mx_cont = mx; my_cont = my; mz_cont = mz;
}

void RigidCube::set_cont_force(const Force3D &cf) noexcept
{
	fx_cont = cf.fx; fy_cont = cf.fy; fz_cont = cf.fz;
	mx_cont = cf.mx; my_cont = cf.my; mz_cont = cf.mz;
}

void RigidCube::set_ext_force(
	double fx,
	double fy,
	double fz,
	double mx,
	double my,
	double mz
	) noexcept
{
	fx_ext = fx; fy_ext = fy; fz_ext = fz;
	mx_ext = mx; my_ext = my; mz_ext = mz;
}

bool RigidCube::detect_collision_with_point(
	double p_x,
	double p_y,
	double p_z,
	double p_r,
	double& dist,
	Vector3D& lnorm,
	Point3D& lcontpos
	) noexcept
{
	double lp_x = p_x - x;
	double lp_y = p_y - y;
	double lp_z = p_z - z;
	double p_hhx = hhx + p_r;
	double p_hhy = hhy + p_r;
	double p_hhz = hhz + p_r;
	if (lp_x < -p_hhx || lp_x > p_hhx ||
		lp_y < -p_hhy || lp_y > p_hhy ||
		lp_z < -p_hhz || lp_z > p_hhz)
		return false;

	unsigned char type;
	double tmp;
	if (lp_x >= 0.0)
	{
		type = 0;
		dist = p_hhx - lp_x;
	}
	else
	{
		type = 1;
		dist = lp_x + p_hhx;
	}

	if (lp_y > 0.0)
	{
		tmp = p_hhy - lp_y;
		if (tmp < dist)
		{
			type = 2;
			dist = tmp;
		}
	}
	else
	{
		tmp = lp_y + p_hhy;
		if (tmp < dist)
		{
			type = 3;
			dist = tmp;
		}
	}

	if (lp_z > 0.0)
	{
		tmp = p_hhz - lp_z;
		if (tmp < dist)
		{
			type = 4;
			dist = tmp;
		}
	}
	else
	{
		tmp = lp_z + p_hhz;
		if (tmp < dist)
		{
			type = 5;
			dist = tmp;
		}
	}

	const Vector3D &norm = norms[type];
	lnorm.x = norm.x;
	lnorm.y = norm.y;
	lnorm.z = norm.z;
	lcontpos.x = lp_x - p_r * lnorm.x;
	lcontpos.y = lp_y - p_r * lnorm.y;
	lcontpos.z = lp_z - p_r * lnorm.z;
	return true;
}
