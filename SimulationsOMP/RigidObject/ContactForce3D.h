#ifndef __Contact_Force_3D_h__
#define __Contact_Force_3D_h__

struct ContactForce3D
{
	double fx, fy, fz;
	double mx, my, mz;

	inline void reset() noexcept
	{
		fx = 0.0;
		fy = 0.0;
		fz = 0.0;
		mx = 0.0;
		my = 0.0;
		mz = 0.0;
	}

	inline void combine(const ContactForce3D &other) noexcept
	{
		fx += other.fx;
		fy += other.fy;
		fz += other.fz;
		mx += other.mx;
		my += other.my;
		mz += other.mz;
	}

	inline void add_force(
		double p_x,
		double p_y,
		double p_z,
		double _fx,
		double _fy,
		double _fz,
		double cen_x,
		double cen_y,
		double cen_z
		) noexcept
	{
		fx += _fx;
		fy += _fy;
		fz += _fz;
		double dx = p_x - cen_x;
		double dy = p_y - cen_y;
		double dz = p_z - cen_z;
		mx += dy * _fz - dz * _fy;
		my += dz * _fx - dx * _fz;
		mz += dx * _fy - dy * _fx;
	}
};

#endif