#ifndef __Rigid_Cylinder_h__
#define __Rigid_Cylinder_h__

#include "Geometry3D.h"
#include "ContactForce3D.h"

class RigidCylinder
{
protected:
	double h, r;
	
	double x, y, z;

	double vx, vy, vz;

	double fx_cont, fy_cont, fz_cont;
	double mx_cont, my_cont, mz_cont;
	
	Cube bbox;
	double h_div_2, r2;
	double norm[3][3];

public:
	explicit RigidCylinder();
	~RigidCylinder();

	inline void reset_f_cont() noexcept
	{
		fx_cont = 0.0;
		fy_cont = 0.0;
		fz_cont = 0.0;
		mx_cont = 0.0;
		my_cont = 0.0;
		mz_cont = 0.0;
	}

	inline void combine_f_cont(const ContactForce3D& other) noexcept
	{
		fx_cont += other.fx;
		fy_cont += other.fy;
		fz_cont += other.fz;
		mx_cont += other.mx;
		my_cont += other.my;
		mz_cont += other.mz;
	}

	void init(double _x, double _y, double _z,
			  double _h, double _r) noexcept;

	inline bool detect_collision_with_point(
		double p_x,
		double p_y,
		double p_z,
		double p_vol,
		double &dist,
		double &norm_x,
		double &norm_y,
		double &norm_z
		) noexcept
	{
		double p_r = 0.5 * pow(p_vol, 0.33333333);
		double xl = bbox.xl - p_r;
		double xu = bbox.xu + p_r;
		double yl = bbox.yl - p_r;
		double yu = bbox.yu + p_r;
		double zl = bbox.zl - p_r;
		double zu = bbox.zu + p_r;
		if (p_x < xl || p_x > xu ||
			p_y < yl || p_y > yu ||
			p_z < zl || p_z > zu)
			return false;

		double x_diff = p_x - x;
		double y_diff = p_y - y;
		double rxy2 = x_diff * x_diff + y_diff * y_diff;
		double z_diff, rxy, tmp;
		if (p_z > bbox.zu)
		{
			if (rxy2 > r2)
			{
				rxy = sqrt(rxy2);
				tmp = rxy - r;
				z_diff = p_z - bbox.zu;
				dist = -sqrt(tmp * tmp + z_diff * z_diff) + p_r;
				norm_x = x_diff;
				norm_y = y_diff;
				norm_z = z_diff;
				tmp = sqrt(norm_x * norm_x + norm_y * norm_y + norm_z * norm_z);
				norm_x /= tmp;
				norm_y /= tmp;
				norm_z /= tmp;
			}
			else
			{
				dist = bbox.zu - p_z + p_r;
				norm_x = 0.0;
				norm_y = 0.0;
				norm_z = 1.0;
			}
		}
		else if (p_z < bbox.zl)
		{
			if (rxy2 > r2)
			{
				rxy = sqrt(rxy2);
				tmp = rxy - r;
				z_diff = bbox.zl - p_z;
				dist = -sqrt(tmp * tmp + z_diff * z_diff) + p_r;
				norm_x = x_diff;
				norm_y = y_diff;
				norm_z = z_diff;
				tmp = sqrt(norm_x * norm_x + norm_y * norm_y + norm_z * norm_z);
				norm_x /= tmp;
				norm_y /= tmp;
				norm_z /= tmp;
			}
			else
			{
				dist = p_z + p_r - bbox.zl;
				norm_x = 0.0;
				norm_y = 0.0;
				norm_z = -1.0;
			}
		}
		else
		{
			if (rxy2 > r2)
			{
				rxy = sqrt(rxy2);
				dist = r - rxy + p_r;
				norm_x = x_diff / rxy;
				norm_y = y_diff / rxy;
				norm_z = 0.0;
			}
			else // inside cylinder
			{
				unsigned char type = 0;
				dist = bbox.zu - p_z + p_r;
				tmp = p_z + p_r - bbox.zl;
				if (tmp < dist)
				{
					type = 1;
					dist = tmp;
				}
				rxy = sqrt(rxy2);
				tmp = r - rxy + p_r;
				if (tmp < dist)
				{
					type = 2;
					dist = tmp;
					if (rxy2 != 0.0)
					{
						norm[2][0] = x_diff / rxy;
						norm[2][1] = y_diff / rxy;
					}
					else
					{
						norm[2][0] = 0.0;
						norm[2][1] = 0.0;
					}
				}
				norm_x = norm[type][0];
				norm_y = norm[type][1];
				norm_z = norm[type][2];
			}
		}

		return dist < 0.0 ? false : true;
	}

	inline void update_motion(double dt) noexcept
	{
		x += vx * dt;
		y += vy * dt;
		z += vz * dt;
		bbox.xl = x - r;
		bbox.xu = x + r;
		bbox.yl = y - r;
		bbox.yu = y + r;
		bbox.zl = z - h_div_2;
		bbox.zu = z + h_div_2;
	}
};

#endif