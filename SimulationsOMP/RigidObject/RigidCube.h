#ifndef __Rigid_Cube_h__
#define __Rigid_Cube_h__

#include "Geometry3D.h"
#include "ContactForce3D.h"

class RigidCube
{
protected:
	double hx, hy, hz;
	double density;

	union
	{
		struct { double ax, ay, az; };
		struct { size_t iax, iay, iaz; };
		Vector3D acceleration;
	};

	union
	{
		struct { double vx, vy, vz; };
		struct { size_t ivx, ivy, ivz; };
		Vector3D velocity;
	};

	union
	{
		struct { double x, y, z; };
		Point3D centre;
	};

	union
	{
		struct
		{
			double fx_cont, fy_cont, fz_cont;
			double mx_cont, my_cont, mz_cont;
		};
		ContactForce3D cont_force;
	};
	
	double m;

	bool has_vx_bc, has_vy_bc, has_vz_bc;
	union
	{
		struct { double vx_bc, vy_bc, vz_bc; };
		struct { size_t ivx_bc, ivy_bc, ivz_bc; };
	};

	double fx_ext, fy_ext, fz_ext;
	double mx_ext, my_ext, mz_ext;

	double hhx, hhy, hhz;

	static const Vector3D norms[6];

public:
	explicit RigidCube();
	~RigidCube();

	inline double get_hx() const noexcept { return hx; }
	inline double get_hy() const noexcept { return hy; }
	inline double get_hz() const noexcept { return hz; }
	inline double get_density() const noexcept { return density; }
	const Point3D& get_centre() const noexcept { return centre; }
	const Vector3D& get_acceleration() const noexcept { return acceleration; }
	const Vector3D& get_velocity() const noexcept { return velocity; }
	const ContactForce3D& get_cont_force() const noexcept { return cont_force; }
	Cube get_bbox() const noexcept
	{ return Cube(-hhx + x, hhx + x, -hhy + y, hhy + y, -hhz + z, hhz + z); }

	inline void set_fx_ext(double f) noexcept { fx_ext = f; }
	inline void set_fy_ext(double f) noexcept { fy_ext = f; }
	inline void set_fz_ext(double f) noexcept { fz_ext = f; }
	inline void set_vx_bc(double vbc) noexcept { has_vx_bc = true; vx_bc = vbc; }
	inline void set_vy_bc(double vbc) noexcept { has_vy_bc = true; vy_bc = vbc; }
	inline void set_vz_bc(double vbc) noexcept { has_vz_bc = true; vz_bc = vbc; }
	inline void set_acceleration(double _ax, double _ay, double _az) noexcept
	{ ax = _ax; ay = _ay; az = _az; }
	inline void set_velocity(double _vx, double _vy, double _vz) noexcept
	{ vx = _vx; vy = _vy; vz = _vz; }

	void init(double _x, double _y, double _z,
			  double _hx,	double _hy, double _hz,
			  double _density) noexcept;
	void set_cont_force(double fx, double fy, double fz,
						double mx, double my, double mz) noexcept;
	void set_cont_force(const ContactForce3D &cf) noexcept;

	inline void reset_cont_force() noexcept
	{
		fx_cont = 0.0; fy_cont = 0.0; fz_cont = 0.0;
		mx_cont = 0.0; my_cont = 0.0; mz_cont = 0.0;
	}

	inline void get_global_point(
		const Point3D& lp,
		Point3D& gp
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
		lp.x = gp.x - x;
		lp.y = gp.y - y;
		lp.z = gp.z - z;
	}

	inline void get_global_vector(
		const Vector3D& lv,
		Vector3D& gv
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
		double p_x, double p_y, double p_z, double p_r,
		double& dist, Vector3D& lnorm, Point3D& lcontpos
	) noexcept;

	inline void update_motion(double dt) noexcept
	{
		ax = (fx_cont + fx_ext) / m;
		ay = (fy_cont + fy_ext) / m;
		az = (fz_cont + fz_ext) / m;
		vx += ax * dt;
		vy += ay * dt;
		vz += az * dt;
		size_t bc_mask;
		bc_mask = size_t(has_vx_bc) + SIZE_MAX;
		iax &= bc_mask;
		ivx = (ivx & bc_mask) | (ivx_bc & ~bc_mask);
		bc_mask = size_t(has_vy_bc) + SIZE_MAX;
		iay &= bc_mask;
		ivy = (ivy & bc_mask) | (ivy_bc & ~bc_mask);
		bc_mask = size_t(has_vz_bc) + SIZE_MAX;
		iaz &= bc_mask;
		ivz = (ivz & bc_mask) | (ivz_bc & ~bc_mask);
		x += vx * dt;
		y += vy * dt;
		z += vz * dt;
	}
};

#endif