#ifndef __Rigid_Cylinder_h__
#define __Rigid_Cylinder_h__

#include "Geometry3D.h"
#include "ContactForce3D.h"

class RigidCylinder
{
protected:
	double h, r;
	
	union
	{
		struct { double x, y, z; };
		Point3D centre;
	};

	union
	{
		struct { double vx, vy, vz; };
		Vector3D velocity;
	};

	union
	{
		ContactForce3D cont_force;
		struct
		{
			double fx_cont, fy_cont, fz_cont;
			double mx_cont, my_cont, mz_cont;
		};
	};
	
	Cube lbbox;
	double h_div_2, r2;
	Vector3D res_norms[3];

public:
	explicit RigidCylinder();
	~RigidCylinder();

	double get_h() const noexcept { return h; }
	double get_r() const noexcept { return r; }
	const Point3D& get_centre() const noexcept { return centre; }
	const Vector3D& get_velocity() const noexcept { return velocity; }
	const ContactForce3D& get_cont_force() const noexcept { return cont_force; }
	
	void init(double _x, double _y, double _z,
		double _h, double _r) noexcept;
	void set_vbc(double _vx, double _vy, double _vz) noexcept;
	void set_cont_force(double fx, double fy, double fz,
						double mx, double my, double mz) noexcept;

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

	inline void get_global_point(
		const Point3D &lp,
		Point3D &gp
		) const noexcept
	{
		gp.x = lp.x + x;
		gp.y = lp.y + y;
		gp.z = lp.z + z;
	}

	inline void get_local_point(
		const Point3D& gp,
		Point3D& lp
		) const noexcept
	{
		lp.x = gp.x + x;
		lp.y = gp.y + y;
		lp.z = gp.z + z;
	}

	inline void get_global_vector(
		const Vector3D &lv,
		Vector3D &gv
		) const noexcept
	{
		gv.x = lv.x;
		gv.y = lv.y;
		gv.z = lv.z;
	}

	inline void get_local_vector(
		const Vector3D& gv,
		Vector3D& lv
		) const noexcept
	{
		lv.x = gv.x;
		lv.y = gv.y;
		lv.z = gv.z;
	}

	bool detect_collision_with_point(
		double p_x,	double p_y,	double p_z, double p_r,
		double &dist, Vector3D& lnorm, Point3D& lcontpos
		) noexcept;

	inline void update_motion(double dt) noexcept
	{
		x += vx * dt;
		y += vy * dt;
		z += vz * dt;
	}
};

#endif