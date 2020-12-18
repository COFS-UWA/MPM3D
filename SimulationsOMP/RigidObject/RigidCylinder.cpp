#include "SimulationsOMP_pcp.h"

#include "RigidCylinder.h"

RigidCylinder::RigidCylinder() : 
	vx(0.0), vy(0.0), vz(0.0),
	fx_cont(0.0), fy_cont(0.0), fz_cont(0.0),
	mx_cont(0.0), my_cont(0.0), mz_cont(0.0) {}

RigidCylinder::~RigidCylinder() {}

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
	lbbox.xl = -r;
	lbbox.xu = r;
	lbbox.yl = -r;
	lbbox.yu = r;
	lbbox.zl = -h_div_2;
	lbbox.zu = h_div_2;
}

void RigidCylinder::set_vbc(
	double _vx,
	double _vy,
	double _vz
	) noexcept
{
	vx = _vx;
	vy = _vy;
	vz = _vz;
}

void RigidCylinder::set_cont_force(
	double fx,
	double fy,
	double fz,
	double mx,
	double my,
	double mz
	) noexcept
{
	fx_cont = fx;
	fy_cont = fy;
	fz_cont = fz;
	mx_cont = mx;
	my_cont = my;
	mz_cont = mz;
}

void RigidCylinder::set_cont_force(
	const Force3D& cf
	) noexcept
{
	fx_cont = cf.fx;
	fy_cont = cf.fy;
	fz_cont = cf.fz;
	mx_cont = cf.mx;
	my_cont = cf.my;
	mz_cont = cf.mz;
}

bool RigidCylinder::detect_collision_with_point(
	double p_x,
	double p_y,
	double p_z,
	double p_r,
	double& dist,
	Vector3D& lnorm,
	Point3D& lcontpos
	) const noexcept
{
	double lp_x = p_x - x;
	double lp_y = p_y - y;
	double lp_z = p_z - z;
	if (lp_x < lbbox.xl - p_r || lp_x > lbbox.xu + p_r ||
		lp_y < lbbox.yl - p_r || lp_y > lbbox.yu + p_r ||
		lp_z < lbbox.zl - p_r || lp_z > lbbox.zu + p_r)
		return false;

	double rxy2 = lp_x * lp_x + lp_y * lp_y;
	double rxy, tmp, z_diff;
	if (lp_z > lbbox.zu)
	{
		if (rxy2 > r2)
		{
			rxy = sqrt(rxy2);
			tmp = rxy - r;
			z_diff = lp_z - lbbox.zu;
			dist = -sqrt(tmp * tmp + z_diff * z_diff) + p_r;
			tmp /= rxy;
			lnorm.x = tmp * lp_x;
			lnorm.y = tmp * lp_y;
			lnorm.z = z_diff;
			tmp = sqrt(lnorm.x * lnorm.x
					 + lnorm.y * lnorm.y
					 + lnorm.z * lnorm.z);
			lnorm.x /= tmp;
			lnorm.y /= tmp;
			lnorm.z /= tmp;
		}
		else
		{
			dist = lbbox.zu - lp_z + p_r;
			lnorm.x = 0.0;
			lnorm.y = 0.0;
			lnorm.z = 1.0;
		}
	}
	else if (lp_z < lbbox.zl)
	{
		if (rxy2 > r2)
		{
			rxy = sqrt(rxy2);
			tmp = rxy - r;
			z_diff = lp_z - lbbox.zl;
			dist = -sqrt(tmp * tmp + z_diff * z_diff) + p_r;
			tmp /= rxy;
			lnorm.x = tmp * lp_x;
			lnorm.y = tmp * lp_y;
			lnorm.z = z_diff;
			tmp = sqrt(lnorm.x * lnorm.x
					 + lnorm.y * lnorm.y
					 + lnorm.z * lnorm.z);
			lnorm.x /= tmp;
			lnorm.y /= tmp;
			lnorm.z /= tmp;
		}
		else
		{
			dist = lp_z + p_r - lbbox.zl;
			lnorm.x = 0.0;
			lnorm.y = 0.0;
			lnorm.z = -1.0;
		}
	}
	else
	{
		if (rxy2 > r2)
		{
			rxy = sqrt(rxy2);
			dist = r - rxy + p_r;
			lnorm.x = lp_x / rxy;
			lnorm.y = lp_y / rxy;
			lnorm.z = 0.0;
		}
		else // inside cylinder
		{
			dist = lbbox.zu - p_z + p_r;
			lnorm.x = 0.0;
			lnorm.y = 0.0;
			lnorm.z = 1.0;
			tmp = p_z + p_r - lbbox.zl;
			if (tmp < dist)
			{
				dist = tmp;
				lnorm.x = 0.0;
				lnorm.y = 0.0;
				lnorm.z = -1.0;
			}
			rxy = sqrt(rxy2);
			tmp = r - rxy + p_r;
			if (tmp < dist)
			{
				dist = tmp;
				if (rxy2 != 0.0)
				{
					lnorm.x = lp_x / rxy;
					lnorm.y = lp_y / rxy;
					lnorm.z = 0.0;
				}
				else
				{
					lnorm.x = 0.0;
					lnorm.y = 0.0;
					lnorm.z = 0.0;
				}
			}
		}
	}

	if (dist < 0.0)
		return false;

	lcontpos.x = lp_x - p_r * lnorm.x;
	lcontpos.y = lp_y - p_r * lnorm.y;
	lcontpos.z = lp_z - p_r * lnorm.z;
	return true;
}
