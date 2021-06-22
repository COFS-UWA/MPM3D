#include "SimulationsOMP_pcp.h"

#include "RigidCylinder.h"

RigidCylinder::RigidCylinder() :
	r(0.0), h(0.0), density(1.0),
	x(0.0), y(0.0), z(0.0),
	ax(0.0), ay(0.0), az(0.0),
	vx(0.0), vy(0.0), vz(0.0),
	fx_cont(0.0), fy_cont(0.0), fz_cont(0.0),
	fx_ext(0.0), fy_ext(0.0), fz_ext(0.0),
	ax_bc_mask(0), ay_bc_mask(0), az_bc_mask(0),
	vx_bc_mask(0), vy_bc_mask(0), vz_bc_mask(0),
	ax_bc(0.0), ay_bc(0.0), az_bc(0.0),
	vx_bc(0.0), vy_bc(0.0), vz_bc(0.0),
	r2(0.0), m(0.0), inv_m(0.0) {}

RigidCylinder::~RigidCylinder() {}

void RigidCylinder::init(
	double _x,
	double _y,
	double _z,
	double _h,
	double _r,
	double _den
	) noexcept
{
	x = _x;
	y = _y;
	z = _z;
	h = _h;
	r = _r;
	density = _den;
	r2 = r * r;
	h_div_2 = 0.5 * h;
	m = 3.14159265359 * r2 * h * density;
	inv_m = 1.0 / m;
	lbbox.xl = -r;
	lbbox.xu = r;
	lbbox.yl = -r;
	lbbox.yu = r;
	lbbox.zl = -h_div_2;
	lbbox.zu = h_div_2;
}

void RigidCylinder::update_motion(double dt) noexcept
{
	// only translation, no rotation
	ax = (fx_cont + fx_ext) * inv_m;
	ay = (fy_cont + fy_ext) * inv_m;
	az = (fz_cont + fz_ext) * inv_m;
	iax = (iax & ~ax_bc_mask) | (iax_bc & ax_bc_mask);
	iay = (iay & ~ay_bc_mask) | (iay_bc & ay_bc_mask);
	iaz = (iaz & ~az_bc_mask) | (iaz_bc & az_bc_mask);
	// update velocity
	vx += ax * dt;
	vy += ay * dt;
	vz += az * dt;
	ivx = (ivx & ~vx_bc_mask) | (ivx_bc & vx_bc_mask);
	ivy = (ivy & ~vy_bc_mask) | (ivy_bc & vy_bc_mask);
	ivz = (ivz & ~vz_bc_mask) | (ivz_bc & vz_bc_mask);
	// update position
	x += vx * dt;
	y += vy * dt;
	z += vz * dt;
}

void RigidCylinder::set_init_v(
	double _vx,
	double _vy,
	double _vz
	) noexcept
{
	vx = _vx;
	vy = _vy;
	vz = _vz;
}

void RigidCylinder::set_vx_bc(double _vx) noexcept
{
	vx_bc_mask = SIZE_MAX;
	vx_bc = _vx;
}

void RigidCylinder::set_vy_bc(double _vy) noexcept
{
	vy_bc_mask = SIZE_MAX;
	vy_bc = _vy;
}

void RigidCylinder::set_vz_bc(double _vz) noexcept
{
	vz_bc_mask = SIZE_MAX;
	vz_bc = _vz;
}

void RigidCylinder::set_vbc(
	double _vx,
	double _vy,
	double _vz
	) noexcept
{
	vx_bc_mask = SIZE_MAX;
	vy_bc_mask = SIZE_MAX;
	vz_bc_mask = SIZE_MAX;
	vx_bc = _vx;
	vy_bc = _vy;
	vz_bc = _vz;
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

void RigidCylinder::set_ext_force(
	double fx, double fy, double fz,
	double mx, double my, double mz
	) noexcept
{
	fx_ext = fx; fy_ext = fy; fz_ext = fz;
	mx_ext = mx; my_ext = my; mz_ext = mz;
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
			dist = lbbox.zu - lp_z + p_r;
			lnorm.x = 0.0;
			lnorm.y = 0.0;
			lnorm.z = 1.0;
			tmp = lp_z + p_r - lbbox.zl;
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

	lcontpos.x = lp_x;
	lcontpos.y = lp_y;
	lcontpos.z = lp_z;
	return true;
}
