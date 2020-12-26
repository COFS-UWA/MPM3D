#ifndef __Rigid_Object_Motion_3D_h__
#define __Rigid_Object_Motion_3D_h__

#include "GeometryUtils.h"
#include "Geometry3D.h"
#include "Force3D.h"

class RigidObjectMotion3D
{
protected:
	union
	{
		struct { double ax, ay, az; };
		struct { size_t iax, iay, iaz; };
		Vector3D acceleration;
	};

	union
	{
		struct { double ax_ang, ay_ang, az_ang; };
		struct { size_t iax_ang, iay_ang, iaz_ang; };
		Vector3D acceleration_ang;
	};

	union
	{
		struct { double vx, vy, vz; };
		struct { size_t ivx, ivy, ivz; };
		Vector3D velocity;
	};
	
	union
	{
		struct { double vx_ang, vy_ang, vz_ang; };
		struct { size_t ivx_ang, ivy_ang, ivz_ang; };
		Vector3D velocity_ang;
	};
	
	double x_ori, y_ori, z_ori;
	double ux, uy, uz;

	union
	{
		struct { double x, y, z; };
		Point3D pos;
	};

	union
	{
		struct { double x_ang, y_ang, z_ang; };
		Vector3D pos_ang;
	};

	size_t ax_bc_mask, ay_bc_mask, az_bc_mask;
	union
	{
		struct { double ax_bc, ay_bc, az_bc; };
		struct { size_t iax_bc, iay_bc, iaz_bc; };
	};

	size_t ax_ang_bc_mask, ay_ang_bc_mask, az_ang_bc_mask;
	union
	{
		struct { double ax_ang_bc, ay_ang_bc, az_ang_bc; };
		struct { size_t iax_ang_bc, iay_ang_bc, iaz_ang_bc; };
	};
	
	size_t vx_bc_mask, vy_bc_mask, vz_bc_mask;
	union
	{
		struct { double vx_bc, vy_bc, vz_bc; };
		struct { size_t ivx_bc, ivy_bc, ivz_bc; };
	};

	size_t vx_ang_bc_mask, vy_ang_bc_mask, vz_ang_bc_mask;
	union
	{
		struct { double vx_ang_bc, vy_ang_bc, vz_ang_bc; };
		struct { size_t ivx_ang_bc, ivy_ang_bc, ivz_ang_bc; };
	};	
	
	union
	{
		struct
		{
			double fx_ext, fy_ext, fz_ext;
			double mx_ext, my_ext, mz_ext;
		};
		Force3D force_ext;
	};
	
	union
	{
		struct
		{
			double fx_cont, fy_cont, fz_cont;
			double mx_cont, my_cont, mz_cont;
		};
		Force3D force_contact;
	};

	double m, moi[6], inv_moi[6];

	union
	{
		struct { Vector3D ix, iy, iz; };
		double T_mat[3][3];
		double T_mat_data[9];
	};

public:
	RigidObjectMotion3D();
	~RigidObjectMotion3D();

	inline double get_m() const noexcept { return m; }
	inline const double *get_moi() const noexcept { return moi; }
	inline const double *get_inv_moi() const noexcept { return inv_moi; }
	inline const double *get_T_mat() const noexcept { return T_mat_data; }
	inline const Vector3D& get_ix() const noexcept { return ix; }
	inline const Vector3D& get_iy() const noexcept { return iy; }
	inline const Vector3D& get_iz() const noexcept { return iz; }
	inline double get_ax() const noexcept { return ax; }
	inline double get_ay() const noexcept { return ay; }
	inline double get_az() const noexcept { return az; }
	inline const Vector3D& get_a() const noexcept { return acceleration; }
	inline double get_ax_ang() const noexcept { return ax_ang; }
	inline double get_ay_ang() const noexcept { return ay_ang; }
	inline double get_az_ang() const noexcept { return az_ang; }
	inline const Vector3D& get_a_ang() const noexcept { return acceleration_ang; }
	inline double get_vx() const noexcept { return vx; }
	inline double get_vy() const noexcept { return vy; }
	inline double get_vz() const noexcept { return vz; }
	inline const Vector3D& get_v() const noexcept { return velocity; }
	inline double get_vx_ang() const noexcept { return vx_ang; }
	inline double get_vy_ang() const noexcept { return vy_ang; }
	inline double get_vz_ang() const noexcept { return vz_ang; }
	inline const Vector3D& get_v_ang() const noexcept { return velocity_ang; }
	inline double get_x() const noexcept { return x; }
	inline double get_y() const noexcept { return y; }
	inline double get_z() const noexcept { return z; }
	inline Point3D get_pos() const noexcept
	{ return Point3D(x_ori + ux, y_ori + uy, z_ori + uz); }
	inline double get_x_ang() const noexcept { return x_ang; }
	inline double get_y_ang() const noexcept { return y_ang; }
	inline double get_z_ang() const noexcept { return z_ang; }
	inline const Vector3D& get_pos_ang() const noexcept { return pos_ang; }
	inline bool has_ax_bc() const noexcept { return ax_bc_mask != 0; }
	inline bool has_ay_bc() const noexcept { return ay_bc_mask != 0; }
	inline bool has_az_bc() const noexcept { return az_bc_mask != 0; }
	inline double get_ax_bc() const noexcept { return ax_bc; }
	inline double get_ay_bc() const noexcept { return ay_bc; }
	inline double get_az_bc() const noexcept { return az_bc; }
	inline bool has_ax_ang_bc() const noexcept { return ax_ang_bc_mask != 0; }
	inline bool has_ay_ang_bc() const noexcept { return ay_ang_bc_mask != 0; }
	inline bool has_az_ang_bc() const noexcept { return az_ang_bc_mask != 0; }
	inline double get_ax_ang_bc() const noexcept { return ax_ang_bc; }
	inline double get_ay_ang_bc() const noexcept { return ay_ang_bc; }
	inline double get_az_ang_bc() const noexcept { return az_ang_bc; }
	inline bool has_vx_bc() const noexcept { return vx_bc_mask != 0; }
	inline bool has_vy_bc() const noexcept { return vy_bc_mask != 0; }
	inline bool has_vz_bc() const noexcept { return vz_bc_mask != 0; }
	inline double get_vx_bc() const noexcept { return vx_bc; }
	inline double get_vy_bc() const noexcept { return vy_bc; }
	inline double get_vz_bc() const noexcept { return vz_bc; }
	inline bool has_vx_ang_bc() const noexcept { return vx_ang_bc_mask != 0; }
	inline bool has_vy_ang_bc() const noexcept { return vy_ang_bc_mask != 0; }
	inline bool has_vz_ang_bc() const noexcept { return vz_ang_bc_mask != 0; }
	inline double get_vx_ang_bc() const noexcept { return vx_ang_bc; }
	inline double get_vy_ang_bc() const noexcept { return vy_ang_bc; }
	inline double get_vz_ang_bc() const noexcept { return vz_ang_bc; }
	inline double get_fx_ext() const noexcept { return fx_ext; }
	inline double get_fy_ext() const noexcept { return fy_ext; }
	inline double get_fz_ext() const noexcept { return fz_ext; }
	inline double get_mx_ext() const noexcept { return mx_ext; }
	inline double get_my_ext() const noexcept { return my_ext; }
	inline double get_mz_ext() const noexcept { return mz_ext; }
	inline const Force3D& get_force_ext() const noexcept { return force_ext; }
	inline double get_fx_contact() const noexcept { return fx_cont; }
	inline double get_fy_contact() const noexcept { return fy_cont; }
	inline double get_fz_contact() const noexcept { return fz_cont; }
	inline double get_mx_contact() const noexcept { return mx_cont; }
	inline double get_my_contact() const noexcept { return my_cont; }
	inline double get_mz_contact() const noexcept { return mz_cont; }
	inline const Force3D& get_force_contact() const noexcept { return force_contact; }

	inline void set_ax_bc(double _ax) noexcept { ax_bc = _ax; ax_bc_mask = SIZE_MAX; }
	inline void set_ay_bc(double _ay) noexcept { ay_bc = _ay; ay_bc_mask = SIZE_MAX; }
	inline void set_az_bc(double _az) noexcept { az_bc = _az; az_bc_mask = SIZE_MAX; }
	inline void set_ax_ang_bc(double _ax_ang) noexcept { ax_ang_bc = _ax_ang; ax_ang_bc_mask = SIZE_MAX; }
	inline void set_ay_ang_bc(double _ay_ang) noexcept { ay_ang_bc = _ay_ang; ay_ang_bc_mask = SIZE_MAX; }
	inline void set_az_ang_bc(double _az_ang) noexcept { az_ang_bc = _az_ang; az_ang_bc_mask = SIZE_MAX; }
	inline void set_vx_bc(double _vx) noexcept { vx_bc = _vx; vx_bc_mask = SIZE_MAX; }
	inline void set_vy_bc(double _vy) noexcept { vy_bc = _vy; vy_bc_mask = SIZE_MAX; }
	inline void set_vz_bc(double _vz) noexcept { vz_bc = _vz; vz_bc_mask = SIZE_MAX; }
	inline void set_vx_ang_bc(double _vx_ang) noexcept { vx_ang_bc = _vx_ang; vx_ang_bc_mask = SIZE_MAX; }
	inline void set_vy_ang_bc(double _vy_ang) noexcept { vy_ang_bc = _vy_ang; vy_ang_bc_mask = SIZE_MAX; }
	inline void set_vz_ang_bc(double _vz_ang) noexcept { vz_ang_bc = _vz_ang; vz_ang_bc_mask = SIZE_MAX; }

	inline void update_motion(double dt) noexcept
	{
		// update a
		ax = (fx_ext + fx_cont) / m;
		ay = (fy_ext + fy_cont) / m;
		az = (fz_ext + fz_cont) / m;
		// apply abc
		iax = (iax & ~ax_bc_mask) | (iax_bc & ax_bc_mask);
		iay = (iay & ~ay_bc_mask) | (iay_bc & ay_bc_mask);
		iaz = (iaz & ~az_bc_mask) | (iaz_bc & az_bc_mask);
		// update velocity
		vx += ax * dt;
		vy += ay * dt;
		vz += az * dt;
		// apply vbc
		ivx = (ivx & ~vx_bc_mask) | (ivx_bc & vx_bc_mask);
		ivy = (ivy & ~vy_bc_mask) | (ivy_bc & vy_bc_mask);
		ivz = (ivz & ~vz_bc_mask) | (ivz_bc & vz_bc_mask);
		// update position
		ux += vx * dt;
		uy += vy * dt;
		uz += vz * dt;
		x = x_ori + ux;
		y = y_ori + uy;
		z = z_ori + uz;

		// 3D rotation (Newton-Euler equation)
		//Matrix3x3 cur_moi = T_mat.transpose() * moi_mat * T_mat;
		//Vector3 a_ang = cur_moi.partialPivLu().solve(m_vec - v_ang_vec.cross(cur_moi * v_ang_vec));
		double vec_tmp1[3], vec_tmp2[3];
		// Tt*moi*T
		vec_tmp1[0] = T_mat[0][0] * vx_ang + T_mat[0][1] * vy_ang + T_mat[0][2] * vz_ang;
		vec_tmp1[1] = T_mat[1][0] * vx_ang + T_mat[1][1] * vy_ang + T_mat[1][2] * vz_ang;
		vec_tmp1[2] = T_mat[2][0] * vx_ang + T_mat[2][1] * vy_ang + T_mat[2][2] * vz_ang;
		vec_tmp2[0] = moi[0] * vec_tmp1[0] + moi[3] * vec_tmp1[1] + moi[5] * vec_tmp1[2];
		vec_tmp2[1] = moi[3] * vec_tmp1[0] + moi[1] * vec_tmp1[1] + moi[4] * vec_tmp1[2];
		vec_tmp2[2] = moi[5] * vec_tmp1[0] + moi[4] * vec_tmp1[1] + moi[2] * vec_tmp1[2];
		vec_tmp1[0] = T_mat[0][0] * vec_tmp2[0] + T_mat[1][0] * vec_tmp2[1] + T_mat[2][0] * vec_tmp2[2];
		vec_tmp1[1] = T_mat[0][1] * vec_tmp2[0] + T_mat[1][1] * vec_tmp2[1] + T_mat[2][1] * vec_tmp2[2];
		vec_tmp1[2] = T_mat[0][2] * vec_tmp2[0] + T_mat[1][2] * vec_tmp2[1] + T_mat[2][2] * vec_tmp2[2];
		// w cross I*w
		vec_tmp2[0] = vec_tmp1[1] * vz_ang - vec_tmp1[2] * vy_ang;
		vec_tmp2[1] = vec_tmp1[2] * vx_ang - vec_tmp1[0] * vz_ang;
		vec_tmp2[2] = vec_tmp1[0] * vy_ang - vec_tmp1[1] * vx_ang;
		// m - w cross I*w
		vec_tmp2[0] = mx_ext + mx_cont - vec_tmp2[0];
		vec_tmp2[1] = my_ext + my_cont - vec_tmp2[1];
		vec_tmp2[2] = mz_ext + mz_cont - vec_tmp2[2];
		// Tt*inv_moi*T*vec2
		vec_tmp1[0] = T_mat[0][0] * vec_tmp2[0] + T_mat[0][1] * vec_tmp2[1] + T_mat[0][2] * vec_tmp2[2];
		vec_tmp1[1] = T_mat[1][0] * vec_tmp2[0] + T_mat[1][1] * vec_tmp2[1] + T_mat[1][2] * vec_tmp2[2];
		vec_tmp1[2] = T_mat[2][0] * vec_tmp2[0] + T_mat[2][1] * vec_tmp2[1] + T_mat[2][2] * vec_tmp2[2];
		vec_tmp2[0] = inv_moi[0] * vec_tmp1[0] + inv_moi[3] * vec_tmp1[1] + inv_moi[5] * vec_tmp1[2];
		vec_tmp2[1] = inv_moi[3] * vec_tmp1[0] + inv_moi[1] * vec_tmp1[1] + inv_moi[4] * vec_tmp1[2];
		vec_tmp2[2] = inv_moi[5] * vec_tmp1[0] + inv_moi[4] * vec_tmp1[1] + inv_moi[2] * vec_tmp1[2];
		ax_ang = T_mat[0][0] * vec_tmp2[0] + T_mat[1][0] * vec_tmp2[1] + T_mat[2][0] * vec_tmp2[2];
		ay_ang = T_mat[0][1] * vec_tmp2[1] + T_mat[1][1] * vec_tmp2[1] + T_mat[2][1] * vec_tmp2[2];
		az_ang = T_mat[0][2] * vec_tmp2[2] + T_mat[1][2] * vec_tmp2[1] + T_mat[2][2] * vec_tmp2[2];
		iax_ang = (iax_ang & ~ax_ang_bc_mask) | (iax_ang_bc & ax_ang_bc_mask);
		iay_ang = (iay_ang & ~ay_ang_bc_mask) | (iay_ang_bc & ay_ang_bc_mask);
		iaz_ang = (iaz_ang & ~az_ang_bc_mask) | (iaz_ang_bc & az_ang_bc_mask);
		vx_ang += ax_ang * dt;
		vy_ang += ay_ang * dt;
		vz_ang += az_ang * dt;
		ivx_ang = (ivx_ang & ~vx_ang_bc_mask) | (ivx_ang_bc & vx_ang_bc_mask);
		ivy_ang = (ivy_ang & ~vy_ang_bc_mask) | (ivy_ang_bc & vy_ang_bc_mask);
		ivz_ang = (ivz_ang & ~vz_ang_bc_mask) | (ivz_ang_bc & vz_ang_bc_mask);
		// adjust local axises
		x_ang += vx_ang * dt;
		trim_to_pi(x_ang);
		y_ang += vy_ang * dt;
		trim_to_pi(y_ang);
		z_ang += vz_ang * dt;
		trim_to_pi(z_ang);
		// update ix, iy, iz
		ix.x = 1.0, ix.y = 0.0, ix.z = 0.0;
		iy.x = 0.0, iy.y = 1.0, iy.z = 0.0;
		iz.x = 0.0, iz.y = 0.0, iz.z = 1.0;
		rotate_axses_by_angle(pos_ang, ix, iy, iz);
	}

	inline void get_local_point(const Point3D& gp, Point3D& lp) const noexcept
	{ point_from_global_to_local_coordinate<Point3D, Point3D>(pos, ix, iy, iz, gp, lp); }
	
	inline void get_global_point(const Point3D& lp, Point3D& gp) const noexcept
	{ point_from_local_to_global_coordinate<Point3D, Point3D>(pos, ix, iy, iz, lp, gp); }

	inline void get_local_vector(const Vector3D& gv, Vector3D& lv) const noexcept
	{ vector_from_global_to_local_coordinate<Vector3D, Vector3D>(ix, iy, iz, gv, lv); }
	
	inline void get_global_vector(const Vector3D& lv, Vector3D& gv) const noexcept
	{ vector_from_local_to_global_coordinate<Vector3D, Vector3D>(ix, iy, iz, lv, gv); }

	void init(double _x, double _y, double _z, double _m, double _moi_data[6]);
	void set_angle(double _x_ang, double _y_ang, double _z_ang); // in radius
	
	void set_translation_velocity_bc(double _vx, double _vy, double _vz);
};

#endif