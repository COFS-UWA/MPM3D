#include "SimulationsOMP_pcp.h"

#include "RigidCone.h"

RigidCone::RigidCone() :
	vx(0.0), vy(0.0), vz(0.0),
	fx_cont(0.0), fy_cont(0.0), fz_cont(0.0),
	mx_cont(0.0), my_cont(0.0), mz_cont(0.0)
{
	Vector3D& v1 = res_norms[1];
	v1.z = 0.0;
	Vector3D& v2 = res_norms[2];
	v2.x = 0.0;
	v2.y = 0.0;
	v2.z = 1.0;
}

RigidCone::~RigidCone()
{

}

void RigidCone::init(
	double _x,
	double _y,
	double _z,
	double _r,
	double _tip_h,
	double _shaft_h
	) noexcept
{
	x = _x;
	y = _y;
	z = _z;
	r = _r;
	h_tip = _tip_h;
	ht_div_r = h_tip / r;
	r_div_ht = r / h_tip;
	r2_div_ht = r * r / h_tip;
	sqrt_one_ht2_div_r2 = sqrt(1.0 + h_tip * h_tip / (r * r));
	h_shaft = _shaft_h;
	lbbox.xl = x - r;
	lbbox.xu = x + r;
	lbbox.yl = y - r;
	lbbox.yu = y + r;
	lbbox.zl = z - h_tip;
	lbbox.zu = z + h_shaft;
}

void RigidCone::set_vbc(
	double _vx,
	double _vy,
	double _vz
	) noexcept
{
	vx = _vx; vy = _vy; vz = _vz;
}

void RigidCone::set_cont_force(
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

bool RigidCone::detect_collision_with_point(
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
	if (lp_x < lbbox.xl - p_r || lp_x > lbbox.xu + p_r ||
		lp_y < lbbox.yl - p_r || lp_y > lbbox.yu + p_r ||
		lp_z < lbbox.zl - p_r || lp_z > lbbox.zu + p_r)
		return false;

	double lp_r = sqrt(lp_x * lp_x + lp_y * lp_y);
	double tip_line = lp_z - ht_div_r * lp_r + h_tip;
	unsigned char norm_type;
	double tmp, norm_tmp;
	// isnide cone
	if (lp_z < h_shaft && lp_r < r && tip_line > 0.0)
	{
		// inside cone
		norm_type = 0;
		dist = abs(tip_line) / sqrt_one_ht2_div_r2;
		Vector3D& v0 = res_norms[0];
		v0.x = lp_x;
		v0.y = lp_y;
		v0.z = -r_div_ht * lp_r;
		norm_tmp = sqrt(v0.x * v0.x + v0.y * v0.y + v0.z * v0.z);
		if (norm_tmp != 0.0)
		{
			v0.x /= norm_tmp;
			v0.y /= norm_tmp;
			v0.z /= norm_tmp;
		}
		else
		{
			v0.x = 0.0;
			v0.y = 0.0;
			v0.z = -1.0;
		}
		tmp = r - lp_r;
		if (tmp < dist)
		{
			norm_type = 1;
			dist = tmp;
			Vector3D& v1 = res_norms[1];
			v1.x = lp_x;
			v1.y = lp_y;
			norm_tmp = sqrt(v1.x * v1.x + v1.y * v1.y);
			v1.x /= norm_tmp;
			v1.y /= norm_tmp;
		}
		tmp = h_shaft - lp_z;
		if (tmp < dist)
		{
			norm_type = 2;
			dist = tmp;
		}
		dist += p_r;
		Vector3D& ln = res_norms[norm_type];
		lnorm.x = ln.x;
		lnorm.y = ln.y;
		lnorm.z = ln.z;
		lcontpos.x = lp_x - p_r * lnorm.x;
		lcontpos.y = lp_y - p_r * lnorm.y;
		lcontpos.z = lp_z - p_r * lnorm.z;
		return true;
	}

	// outside cone
	double cone_line_up = lp_z + r_div_ht * lp_r - r2_div_ht;
	double cone_line_down = lp_z + r_div_ht * lp_r + h_tip;
	if (lp_z < 0.0)
	{
		if (cone_line_down <= 0.0)
		{
			tmp = lp_z + h_tip;
			dist = -sqrt(lp_r * lp_r + tmp * tmp);
			lnorm.x = lp_x;
			lnorm.y = lp_y;
			lnorm.z = tmp;
			norm_tmp = sqrt(lnorm.x * lnorm.x
				+ lnorm.y * lnorm.y
				+ lnorm.z * lnorm.z);
			if (norm_tmp != 0.0)
			{
				lnorm.x /= norm_tmp;
				lnorm.y /= norm_tmp;
				lnorm.z /= norm_tmp;
			}
			else // at the tip
			{
				lnorm.x = 0.0;
				lnorm.y = 0.0;
				lnorm.z = -1.0;
			}
		}
		else if (cone_line_up < 0.0)
		{
			dist = -abs(tip_line) / sqrt_one_ht2_div_r2;
			lnorm.x = lp_x;
			lnorm.y = lp_y;
			lnorm.z = -r_div_ht * lp_r;
			norm_tmp = sqrt(lnorm.x * lnorm.x
				+ lnorm.y * lnorm.y
				+ lnorm.z * lnorm.z);
			lnorm.x /= norm_tmp;
			lnorm.y /= norm_tmp;
			lnorm.z /= norm_tmp;
		}
		else
		{
			tmp = lp_r - r;
			dist = -sqrt(tmp * tmp + lp_z * lp_z);
			tmp /= lp_r;
			lnorm.x = tmp * lp_x;
			lnorm.y = tmp * lp_y;
			lnorm.z = lp_z;
			norm_tmp = sqrt(lnorm.x * lnorm.x
						  + lnorm.y * lnorm.y
						  + lnorm.z * lnorm.z);
			lnorm.x /= norm_tmp;
			lnorm.y /= norm_tmp;
			lnorm.z /= norm_tmp;
		}
	}
	else // lp_z >= 0.0
	{
		if (lp_z <= h_shaft)
		{
			dist = r - lp_r;
			lnorm.x = lp_x;
			lnorm.y = lp_y;
			norm_tmp = sqrt(lnorm.x * lnorm.x + lnorm.y * lnorm.y);
			lnorm.x /= norm_tmp;
			lnorm.y /= norm_tmp;
			lnorm.z = 0.0;
		}
		else if (lp_r < r)
		{
			dist = h_shaft - lp_z;
			lnorm.x = 0.0;
			lnorm.y = 0.0;
			lnorm.z = 1.0;
		}
		else
		{
			tmp = lp_r - r;
			lnorm.z = lp_z - h_shaft;
			dist = -sqrt(tmp  * tmp + lnorm.z * lnorm.z);
			tmp /= lp_r;
			lnorm.x = tmp * lp_x;
			lnorm.y = tmp * lp_y;
			norm_tmp = sqrt(lnorm.x * lnorm.x
						  + lnorm.y * lnorm.y
						  + lnorm.z * lnorm.z);
			lnorm.x /= norm_tmp;
			lnorm.y /= norm_tmp;
			lnorm.z /= norm_tmp;
		}
	}

	dist += p_r;
	if (dist < 0.0)
		return false;
	lcontpos.x = lp_x - p_r * lnorm.x;
	lcontpos.y = lp_y - p_r * lnorm.y;
	lcontpos.z = lp_z - p_r * lnorm.z;
	return true;
}
