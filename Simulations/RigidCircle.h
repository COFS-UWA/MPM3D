#ifndef __Rigid_Circle_h__
#define __Rigid_Circle_h__

#include <cmath>

#include "Geometry.h"

#define PI 3.14159265359

class RigidCircle;

struct RigidCircleForce
{
protected:
	friend class RigidCircle;
	double fx, fy, m;

public:
	inline void reset_rf() noexcept { fx = 0.0; fy = 0.0; m = 0.0; }
	inline void add_rf(
		double _x,   double _y,
		double _fx,  double _fy,
		double rc_x, double rc_y
		) noexcept
	{
		fx += _fx;
		fy += _fy;
		m += (_x - rc_x) * _fy - (_y - rc_y) * _fx;
	}
	inline void combine(const RigidCircleForce &other) noexcept
	{
		fx += other.fx;
		fy += other.fy;
		m += other.m;
	}
};

// for t-bar penetration and pipe embedment
class RigidCircle
{
protected:
	double r; // radius
	double density;
	// motion
	double rfx, rfy, rm; // reaction force
	double ax, ay, a_ang; // acceleration
	double vx, vy, v_ang; // velocity
	double x, y, ang; // position and angle

	double *pax, *pay, *pa_ang;
	double *pvx, *pvy, *pv_ang;

	// boundary condition
	double rfx_bc;
	double rfy_bc;
	double rm_bc;
	double ax_bc;
	double ay_bc;
	double a_ang_bc;
	double vx_bc;
	double vy_bc;
	double v_ang_bc;

	// calculation variables
	double r2, m, moi; // moment of inertia
	inline void init_cal_var()
	{
		r2 = r * r;
		m = PI * r2 * density;
		moi = 0.5 * m * r2;
	}

public:
	RigidCircle() : r(0.0), density(1.0),
		x(0.0), y(0.0), ang(0.0),
		ax(0.0), ay(0.0), a_ang(0.0),
		vx(0.0), vy(0.0), v_ang(0.0),
		rfx(0.0), rfy(0.0), rm(0.0),
		r2(0.0), m(0.0), moi(0.0),
		pax(&ax), pay(&ay), pa_ang(&a_ang),
		pvx(&vx), pvy(&vy), pv_ang(&v_ang),
		rfx_bc(0.0), rfy_bc(0.0), rm_bc(0.0)
	{ init_cal_var(); }
	~RigidCircle() {}

	inline double get_radius() { return r; }
	inline double get_density() { return density; }
	inline double get_x() { return x; }
	inline double get_y() { return y; }
	inline double get_ang() { return ang; }
	inline double get_ax() { return ax; }
	inline double get_ay() { return ay; }
	inline double get_a_ang() { return a_ang; }
	inline double get_vx() { return vx; }
	inline double get_vy() { return vy; }
	inline double get_v_ang() { return v_ang; }
	inline double get_rfx() { return rfx; }
	inline double get_rfy() { return rfy; }
	inline double get_rm() { return rm; }
	inline double get_m() { return m; }
	inline double get_moi() { return moi; }

	inline bool has_ax_bc() { return pax == &ax_bc; }
	inline bool has_ay_bc() { return pay == &ay_bc; }
	inline bool has_a_ang_bc() { return pa_ang == &a_ang_bc; }
	
	inline bool has_vx_bc() { return pvx == &vx_bc; }
	inline bool has_vy_bc() { return pvy == &vy_bc; }
	inline bool has_v_ang_bc() { return pv_ang == &v_ang_bc; }

	inline double get_ax_bc() { return ax_bc; }
	inline double get_ay_bc() { return ay_bc; }
	inline double get_a_ang_bc() { return a_ang_bc; }

	inline double get_vx_bc() { return vx_bc; }
	inline double get_vy_bc() { return vy_bc; }
	inline double get_v_ang_bc() { return v_ang_bc; }

	inline double get_rfx_bc() { return rfx_bc; }
	inline double get_rfy_bc() { return rfy_bc; }
	inline double get_rm_bc() { return rm_bc; }

	inline void set_ax_bc(double _a) { pax = &ax_bc; ax_bc = _a; }
	inline void set_ay_bc(double _a) { pay = &ay_bc; ay_bc = _a; }
	inline void set_a_ang_bc(double _a) { pa_ang = &a_ang_bc; a_ang_bc = _a; }

	inline void set_vx_bc(double _v) { pvx = &vx_bc; vx_bc = _v; }
	inline void set_vy_bc(double _v) { pvy = &vy_bc; vy_bc = _v; }
	inline void set_v_ang_bc(double _v) { pv_ang = &v_ang_bc; v_ang_bc = _v; }
	inline void set_v_bc(double _vx, double _vy, double _v_ang)
	{ set_vx_bc(_vx); set_vy_bc(_vy); set_v_ang_bc(_v_ang); }

	inline void add_rfx_bc(double _rf) { rfx_bc += _rf; }
	inline void add_rfy_bc(double _rf) { rfy_bc += _rf; }
	inline void add_rm_bc(double _rf) { rm_bc += _rf; }
	inline void add_rf_bc(double _x, double _y, double _fx, double _fy)
	{
		rfx_bc += _fx; rfy_bc += _fy;
		rm_bc += (_x - x) * _fy - (_y - y) * _fx;
	}

	void init(double _x, double _y, double _r, double _density = 1.0);

	void set_init_state(double _r, double _density,
		double _rfx, double _rfy, double _rm,
		double _ax, double _ay, double _a_ang,
		double _vx, double _vy, double _v_ang,
		double _x, double _y, double _ang);

public: // helper function for calculation
	inline Rect get_bbox(double exp_size = 0.0)
	{
		return Rect(x - r - exp_size, x + r + exp_size,
					y - r - exp_size, y + r + exp_size);
	}

	inline void reset_rf() { rfx = 0.0; rfy = 0.0; rm = 0.0; }

	inline void add_rf(double _x, double _y, double fx, double fy)
	{
		rfx += fx; rfy += fy;
		rm += (_x - x) * fy - (_y - y) * fx;
	}

	// for parallelsim
	inline void add_rcf(RigidCircleForce& rcf)
	{
		rfx += rcf.fx; rfy += rcf.fy; rm += rcf.m;
	}

	inline void set_rcf(RigidCircleForce& rcf)
	{
		rfx = rcf.fx; rfy = rcf.fy; rm = rcf.m;
	}

	inline bool is_in_circle(double _x, double _y)
	{
		double x_diff = _x - x;
		double y_diff = _y - y;
		return (x_diff * x_diff + y_diff * y_diff) <= r2;
	}

	// return true if point(x, y) collides with circle
	// also calculate overlapping distance and norm direction at that point
	// return false if point is outside circle
	inline bool detect_collision_with_point(double _x, double _y, double pcl_vol,
		double& overlap_dist, double& norm_x, double& norm_y)
	{
		double x_diff = _x - x;
		double y_diff = _y - y;
		double pcl_radius = sqrt(pcl_vol / PI);
		double dist = sqrt(x_diff * x_diff + y_diff * y_diff) - pcl_radius;
		if (dist > r) // not overlapping
			return false;

		norm_x = x_diff / dist;
		norm_y = y_diff / dist;
		overlap_dist = r - dist;
		if (overlap_dist < 0.0)
			overlap_dist = 0.0;
		return true;
	}

	// update rigid circle motion and position
	void update_motion(double dt)
	{
		// update a
		ax = (rfx + rfx_bc) / m;
		ay = (rfy + rfy_bc) / m;
		a_ang = (rm + rm_bc) / moi;
		// apply abc
		ax = *pax;
		ay = *pay;
		a_ang = *pa_ang;
		// update velocity
		vx += ax * dt;
		vy += ay * dt;
		v_ang += a_ang * dt;
		// apply vbc
		vx = *pvx;
		vy = *pvy;
		v_ang = *pv_ang;
		// update position
		x += vx * dt;
		y += vy * dt;
		ang += v_ang * dt;
	}
};

#undef PI

#endif