#ifndef __Rigid_Object_Motion_2D_h__
#define __Rigid_Object_Motion_2D_h__

#include "GeometryUtils.h"
#include "Geometry2D.h"
#include "Force2D.h"
#include "Ratio.h"

class RigidObjectMotion2D
{
protected:
	union
	{
		struct { double ax, ay; };
		struct { size_t iax, iay; };
		Vector2D acceleration;
	};

	double a_ang;

	union
	{
		struct { double vx, vy; };
		struct { size_t ivx, ivy; };
		Vector2D velocity;
	};

	double v_ang;

	double x_ori, y_ori;
	double ux, uy;

	union
	{
		struct { double x, y; };
		Point2D pos;
	};

	double ang;

	size_t ax_bc_mask, ay_bc_mask;
	union
	{
		struct { double ax_bc, ay_bc; };
		struct { size_t iax_bc, iay_bc; };
	};

	size_t a_ang_bc_mask;
	union
	{
		double a_ang_bc;
		size_t ia_ang_bc;
	};

	size_t vx_bc_mask, vy_bc_mask;
	union
	{
		struct { double vx_bc, vy_bc; };
		struct { size_t ivx_bc, ivy_bc; };
	};

	size_t v_ang_bc_mask;
	union
	{
		double v_ang_bc;
		size_t iv_ang_bc;
	};

	union
	{
		struct { double fx_ext, fy_ext, m_ext; };
		Force2D force_ext;
	};
	
	union
	{
		struct { double fx_cont, fy_cont, m_cont; };
		Force2D force_cont;
	};

	double m, moi, inv_moi;

	union
	{
		struct { Vector2D ix, iy; };
		double T_mat[2][2];
		double T_mat_data[4];
	};

public:
	RigidObjectMotion2D();
	~RigidObjectMotion2D();

	inline double get_m() const noexcept { return m; }
	inline const double get_moi() const noexcept { return moi; }
	inline const double get_inv_moi() const noexcept { return inv_moi; }
	inline const double *get_T_mat() const noexcept { return T_mat_data; }
	inline const Vector2D& get_ix() const noexcept { return ix; }
	inline const Vector2D& get_iy() const noexcept { return iy; }
	inline double get_ax() const noexcept { return ax; }
	inline double get_ay() const noexcept { return ay; }
	inline const Vector2D& get_a() const noexcept { return acceleration; }
	inline double get_a_ang() const noexcept { return a_ang; }
	inline double get_vx() const noexcept { return vx; }
	inline double get_vy() const noexcept { return vy; }
	inline const Vector2D& get_v() const noexcept { return velocity; }
	inline double get_v_ang() const noexcept { return v_ang; }
	inline double get_x() const noexcept { return x; }
	inline double get_y() const noexcept { return y; }
	inline Point2D get_pos() const noexcept
	{ return Point2D(x_ori + ux, y_ori + uy); }
	inline double get_pos_ang() const noexcept { return ang; }
	inline bool has_ax_bc() const noexcept { return ax_bc_mask != 0; }
	inline bool has_ay_bc() const noexcept { return ay_bc_mask != 0; }
	inline double get_ax_bc() const noexcept { return ax_bc; }
	inline double get_ay_bc() const noexcept { return ay_bc; }
	inline bool has_a_ang_bc() const noexcept { return a_ang_bc_mask != 0; }
	inline double get_a_ang_bc() const noexcept { return a_ang_bc; }
	inline bool has_vx_bc() const noexcept { return vx_bc_mask != 0; }
	inline bool has_vy_bc() const noexcept { return vy_bc_mask != 0; }
	inline double get_vx_bc() const noexcept { return vx_bc; }
	inline double get_vy_bc() const noexcept { return vy_bc; }
	inline bool has_v_ang_bc() const noexcept { return v_ang_bc_mask != 0; }
	inline double get_v_ang_bc() const noexcept { return v_ang_bc; }
	inline double get_fx_ext() const noexcept { return fx_ext; }
	inline double get_fy_ext() const noexcept { return fy_ext; }
	inline double get_m_ext() const noexcept { return m_ext; }
	inline const Force2D& get_force_ext() const noexcept { return force_ext; }
	inline double get_fx_contact() const noexcept { return fx_cont; }
	inline double get_fy_contact() const noexcept { return fy_cont; }
	inline double get_m_contact() const noexcept { return m_cont; }
	inline const Force2D& get_force_contact() const noexcept { return force_cont; }

	inline void set_ax_bc(double _ax) noexcept { ax_bc = _ax; ax_bc_mask = SIZE_MAX; }
	inline void set_ay_bc(double _ay) noexcept { ay_bc = _ay; ay_bc_mask = SIZE_MAX; }
	inline void set_a_ang_bc(double _ax_ang) noexcept { a_ang_bc = _ax_ang; a_ang_bc_mask = SIZE_MAX; }
	inline void set_vx_bc(double _vx) noexcept { vx_bc = _vx; vx_bc_mask = SIZE_MAX; }
	inline void set_vy_bc(double _vy) noexcept { vy_bc = _vy; vy_bc_mask = SIZE_MAX; }
	inline void set_v_ang_bc(double _vx_ang) noexcept { v_ang_bc = _vx_ang; v_ang_bc_mask = SIZE_MAX; }

	inline void get_local_point(const Point2D& gp, Point2D& lp) const noexcept
	{ point_from_global_to_local_coordinate<Point2D, Point2D>(pos, ix, iy, gp, lp); }
	
	inline void get_global_point(const Point2D& lp, Point2D& gp) const noexcept
	{ point_from_local_to_global_coordinate<Point2D, Point2D>(pos, ix, iy, lp, gp); }

	inline void get_local_vector(const Vector2D& gv, Vector2D& lv) const noexcept
	{ vector_from_global_to_local_coordinate<Vector2D, Vector2D>(ix, iy, gv, lv); }
	
	inline void get_global_vector(const Vector2D& lv, Vector2D& gv) const noexcept
	{ vector_from_local_to_global_coordinate<Vector2D, Vector2D>(ix, iy, lv, gv); }

	void init(double _x, double _y, double _m, double _moi);
	void set_angle(double _ang); // in radius
	
	void set_translation_velocity_bc(double _vx, double _vy);

	void update_motion(double dt) noexcept;
};

#endif