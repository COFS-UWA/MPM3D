#include "SimulationsOMP_pcp.h"

#include <Eigen/Eigen>

#include "RigidObjectMotion2D.h"

RigidObjectMotion2D::RigidObjectMotion2D() :
	ax(0.0), ay(0.0), a_ang(0.0),
	vx(0.0), vy(0.0), v_ang(0.0),
	ux(0.0), uy(0.0), ang(0.0),
	ix(1.0, 0.0), iy(0.0, 1.0),
	ax_bc_mask(0), ay_bc_mask(0),
	ax_bc(0.0), ay_bc(0.0),
	a_ang_bc_mask(0), a_ang_bc(0.0),
	vx_bc_mask(0), vy_bc_mask(0),
	vx_bc(0.0), vy_bc(0.0),
	v_ang_bc_mask(0), v_ang_bc(0.0),
	fx_ext(0.0), fy_ext(0.0), m_ext(0.0),
	fx_cont(0.0), fy_cont(0.0), m_cont(0.0) {}

RigidObjectMotion2D::~RigidObjectMotion2D() {}

void RigidObjectMotion2D::init(double _x, double _y, double _m, double _moi)
{
	x_ori = _x;	y_ori = _y;
	x = _x; y = _y;
	m = _m;
	moi = _moi;
	inv_moi = 1.0 / moi;
}

void RigidObjectMotion2D::set_angle(double _ang)
{
	ang = _ang;
	trim_to_pi(_ang);
	rotate_axses_by_angle(_ang, ix, iy);
}

void RigidObjectMotion2D::set_translation_velocity_bc(double _vx, double _vy)
{
	set_vx_bc(_vx);
	set_vy_bc(_vy);
	set_v_ang_bc(0.0);
}

void RigidObjectMotion2D::update_motion(double dt) noexcept
{
	// update a
	ax = (fx_ext + fx_cont) / m;
	ay = (fy_ext + fy_cont) / m;
	// apply abc
	iax = (iax & ~ax_bc_mask) | (iax_bc & ax_bc_mask);
	iay = (iay & ~ay_bc_mask) | (iay_bc & ay_bc_mask);
	// update velocity
	vx += ax * dt;
	vy += ay * dt;
	ivx = (ivx & ~vx_bc_mask) | (ivx_bc & vx_bc_mask);
	ivy = (ivy & ~vy_bc_mask) | (ivy_bc & vy_bc_mask);
	// update position
	ux += vx * dt;
	uy += vy * dt;
	x = x_ori + ux;
	y = y_ori + uy;

	// 2D rotation
	ang = 0.0;
	//ang += v_ang * dt;
	trim_to_pi(ang);
	rotate_axses_by_angle(ang, ix, iy);
}
