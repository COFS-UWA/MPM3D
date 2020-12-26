#include "SimulationsOMP_pcp.h"

#include "RigidCircle.h"

namespace RigidObject
{

	static constexpr double PI = 3.14159265359;

	void RigidCircle::init(
		double _x,
		double _y,
		double _r,
		double _density
	) noexcept
	{
		x = _x;
		y = _y;
		ang = 0.0;
		r = _r;
		density = _density;
		// cal var
		lbbox.xl = -r;
		lbbox.xu = r;
		lbbox.yl = -r;
		lbbox.yu = r;
		r2 = r * r;
		m = PI * r2 * density;
		moi = 0.5 * m * r2;
		inv_m = 1.0 / m;
		inv_moi = 1.0 / moi;
	}

	void RigidCircle::set_vbc(
		double _vx,
		double _vy,
		double _vang
	) noexcept
	{
		vx_bc_mask = SIZE_MAX;
		vy_bc_mask = SIZE_MAX;
		v_ang_bc_mask = SIZE_MAX;
		vx_bc = _vx;
		vy_bc = _vy;
		v_ang_bc = _vang;
	}

	void RigidCircle::set_cont_force(
		Force2D& cont_force
	) noexcept
	{
		fx_cont = cont_force.fx;
		fy_cont = cont_force.fy;
		m_cont = cont_force.m;
	}

}
