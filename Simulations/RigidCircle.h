#ifndef __Rigid_Circle_h__
#define __Rigid_Circle_h__

#include <cmath>

// For t-bar penetration and pipe embedment
class RigidCircle
{
public:
	struct State
	{
		double r, r2; // radius
		double cen_x, cen_y, theta; // position
		double vx, vy, w; // velocity
		double rfx, rfy, rm; // reaction force
	};
protected:
	union
	{
		struct
		{
			double r, r2;
			double cen_x, cen_y, theta;
			double vx, vy, w;
			double rfx, rfy, rm;
		};
		State state;
	};

public:
	RigidCircle() : r(0.0), r2(0.0),
		cen_x(0.0), cen_y(0.0), theta(0.0),
		vx(0.0), vy(0.0), w(0.0),
		rfx(0.0), rfy(0.0), rm(0.0) {}
	~RigidCircle() {}

	inline const State& get_state(void) noexcept { return state; }

	int init(double _r, double _x, double _y, double max_pcl_size);

	inline void reset_rf() { rfx = 0.0; rfy = 0.0; rm = 0.0; }

	inline void set_velocity(double _vx, double _vy, double _w)
	{
		vx = _vx; vy = _vy; w = _w;
	}

	inline bool is_in_circle(double x, double y)
	{
		double x_diff = x - cen_x;
		double y_diff = y - cen_y;
		return x_diff * x_diff + y_diff * y_diff <= r2;
	}

	// return true if point(x, y) collides with circle
	// also calculate overlapping distance and norm force at that point
	// return false if point is outside circle
	inline bool detect_collision_with_point(double x, double y,
		double& overlap_dist, double& norm_x, double& norm_y)
	{
		double x_diff = x - cen_x;
		double y_diff = y - cen_y;
		double dist = x_diff * x_diff + y_diff * y_diff;
		if (dist > r2) // not overlapping
			return false;

		dist = sqrt(dist);
		norm_x = x_diff / dist;
		norm_y = y_diff / dist;
		overlap_dist = r - dist;
		return true;
	}

	// add reaction force
	inline void add_rf(double x, double y, double fx, double fy)
	{
		rfx += fx;
		rfy += fy;
		rm += (x - cen_x) * fy - (y - cen_y) * fx;
	}

	// update rigid circle position
	// currently is displacement control
	void update(double dt)
	{
		double dx, dy, dtheta;
		dx = vx * dt;
		dy = vy * dt;
		dtheta = w * dt;
		// update position
		cen_x += dx;
		cen_y += dy;
		theta += dtheta;
	}
};

#endif