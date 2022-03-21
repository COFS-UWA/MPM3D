#ifndef __Rigid_Cylinder_h__
#define __Rigid_Cylinder_h__

#include "Geometry3D.h"
#include "Force3D.h"
#include "Ratio.h"

class RigidCylinder
{
protected:
	double h, r;
	double density;

	union
	{
		struct { double x, y, z; };
		Point3D centre;
	};

	union // acceleration
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
		struct
		{
			double fx_cont, fy_cont, fz_cont;
			double mx_cont, my_cont, mz_cont;
		};
		Force3D cont_force;
	};

	union
	{
		struct
		{
			double fx_ext, fy_ext, fz_ext;
			double mx_ext, my_ext, mz_ext;
		};
		Force3D ext_force;
	};

	size_t ax_bc_mask, ay_bc_mask, az_bc_mask;
	size_t vx_bc_mask, vy_bc_mask, vz_bc_mask;

	union
	{
		struct { double ax_bc, ay_bc, az_bc; };
		struct { size_t iax_bc, iay_bc, iaz_bc; };
		Vector3D acceleration_bc;
	};

	OneRatio one_ratio;
	QuadraticRampUpRatio qru_x_ratio, qru_y_ratio, qru_z_ratio;
	Ratio* pvx_bc_ratio, * pvy_bc_ratio, * pvz_bc_ratio;
	double cur_time;
	union
	{
		struct { double vx_bc, vy_bc, vz_bc; };
		Vector3D velocity_bc;
	};

	double h_div_2, r2;
	double m, inv_m;
	Cube lbbox;
	
public:
	explicit RigidCylinder();
	~RigidCylinder();

	inline double get_h() const noexcept { return h; }
	inline double get_r() const noexcept { return r; }
	inline double get_density() const noexcept { return density; }
	inline const Point3D &get_centre() const noexcept { return centre; }
	inline const Vector3D& get_acceleration() const noexcept { return acceleration; }
	inline const Vector3D &get_velocity() const noexcept { return velocity; }
	inline const size_t get_vx_bc_mask() const noexcept { return vx_bc_mask; }
	inline const size_t get_vy_bc_mask() const noexcept { return vy_bc_mask; }
	inline const size_t get_vz_bc_mask() const noexcept { return vz_bc_mask; }
	inline const Vector3D& get_acceleration_bc() const noexcept { return acceleration_bc; }
	inline const Vector3D& get_velocity_bc() const noexcept { return velocity_bc; }
	inline const Force3D &get_cont_force() const noexcept { return cont_force; }
	inline const Force3D& get_ext_force() const noexcept { return ext_force; }
	inline const double get_m() const noexcept { return m; }
	Cube get_bbox() const noexcept
	{
		return Cube(lbbox.xl + x, lbbox.xu + x,
					lbbox.yl + y, lbbox.yu + y,
					lbbox.zl + z, lbbox.zu + z);
	}

	void init(double _x, double _y, double _z,
			  double _h, double _r, double _den = 1.0) noexcept;
	void set_init_v(double _vx, double _vy, double _vz) noexcept;
	void set_vx_bc(double _vx) noexcept;
	void set_vy_bc(double _vy) noexcept;
	void set_vz_bc(double _vz) noexcept;
	void set_vbc(double _vx, double _vy, double _vz) noexcept;
	void set_cont_force(double fx, double fy, double fz,
						double mx, double my, double mz) noexcept;
	void set_cont_force(const Force3D &cf) noexcept;
	void set_ext_force(double fx, double fy, double fz,
			double mx, double my, double mz) noexcept;
	inline void set_const_vx_bc() noexcept { pvx_bc_ratio = &one_ratio; }
	inline void set_const_vy_bc() noexcept { pvy_bc_ratio = &one_ratio; }
	inline void set_const_vz_bc() noexcept { pvz_bc_ratio = &one_ratio; }
	inline void set_ramp_up_vx_bc(double ru_time) noexcept { pvx_bc_ratio = &qru_x_ratio; qru_x_ratio.set_ramp_up_range(ru_time); }
	inline void set_ramp_up_vy_bc(double ru_time) noexcept { pvy_bc_ratio = &qru_y_ratio; qru_y_ratio.set_ramp_up_range(ru_time); }
	inline void set_ramp_up_vz_bc(double ru_time) noexcept { pvz_bc_ratio = &qru_z_ratio; qru_z_ratio.set_ramp_up_range(ru_time); }

	inline void reset_cont_force() noexcept
	{
		fx_cont = 0.0; fy_cont = 0.0; fz_cont = 0.0;
		mx_cont = 0.0; my_cont = 0.0; mz_cont = 0.0;
	}

	void update_motion(double dt) noexcept;
	
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
		lp.x = gp.x - x;
		lp.y = gp.y - y;
		lp.z = gp.z - z;
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
		double &dist, Vector3D &lnorm, Point3D &lcontpos
		) const noexcept;
};

#endif