#ifndef __Rigid_Cone_h__
#define __Rigid_Cone_h__

#include "Geometry3D.h"
#include "ContactForce3D.h"

class RigidCone
{
protected:
	double x, y, z;

	double r, tip_h, shaft_h;

	double fx_cont, fy_cont, fz_cont;
	double mx_cont, my_cont, mz_cont;
	
public:
	explicit RigidCone() {}
	~RigidCone() {}

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
		double _r, double _tip_h, double _shaft_h);

	inline bool detect_collision_with_point(
		double p_x,
		double p_y,
		double p_z,
		double p_vol,
		double& dist,
		double& norm_x,
		double& norm_y,
		double& norm_z
		)
	{

	}

	inline void update_motion(double dt) noexcept
	{

	}
};

#endif