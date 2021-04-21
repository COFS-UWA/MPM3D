#include "Simulations_pcp.h"

#include "RigidRect.h"

void RigidRect::init(
	double _x, double _y,
	double _hx, double _hy,
	double _density)
{
	hx = _hx;
	hy = _hy;
	density = _density;
	x = _x;
	y = _y;
	ang = 0.0;
	init_cal_var();
}

void RigidRect::set_init_state(
	double _hx, double _hy, double _density,
	double _fx_cont, double _fy_cont, double _m_cont,
	double _ax, double _ay, double _a_ang,
	double _vx, double _vy, double _v_ang,
	double _x, double _y, double _ang)
{
	hx = _hx;
	hy = _hy;
	density = _density;
	ax = _ax;
	ay = _ay;
	a_ang = _a_ang;
	vx = _vx;
	vy = _vy;
	v_ang = _v_ang;
	x = _x;
	y = _y;
	ang = _ang;
	fx_con = _fx_cont;
	fy_con = _fy_cont;
	m_con = _m_cont;
	init_cal_var();
}

const Vector2D RigidRect::edge_normals[4] = {
	{ -1.0, 0.0 },
	{  1.0, 0.0 },
	{ 0.0, -1.0 },
	{ 0.0,  1.0 }
};

void RigidRect::update_motion(double dt)
{
	// update a
	ax = (fx_ext + fx_con) / m;
	ay = (fy_ext + fy_con) / m;
	a_ang = (m_ext + m_con) / moi;
	// apply abc
	ax_ull = (ax_ull & ~ax_bc_mask) | (ax_bc_ull & ax_bc_mask);
	ay_ull = (ay_ull & ~ay_bc_mask) | (ay_bc_ull & ay_bc_mask);
	a_ang_ull = (a_ang_ull & ~a_ang_bc_mask) | (a_ang_bc_ull & a_ang_bc_mask);
	// update velocity
	vx += ax * dt;
	vy += ay * dt;
	v_ang += a_ang * dt;
	// apply vbc
	vx_ull = (vx_ull & ~vx_bc_mask) | (vx_bc_ull & vx_bc_mask);
	vy_ull = (vy_ull & ~vy_bc_mask) | (vy_bc_ull & vy_bc_mask);
	v_ang_ull = (v_ang_ull & ~v_ang_bc_mask) | (v_ang_bc_ull & v_ang_bc_mask);
	// update position
	x += vx * dt;
	y += vy * dt;
	//ang += v_ang * dt;
	ang = 0.0;
	sin_ang = sin(ang);
	cos_ang = cos(ang);
}

bool RigidRect::detect_collision_with_point(
	double p_x, double p_y, double p_r,
	double& overlap_dist,
	Vector2D& lnorm, Point2D& lcont_pos)
{
	const Vector2D ix(cos_ang, sin_ang), iy(-sin_ang, cos_ang);
	point_from_global_to_local_coordinate(centre,
		ix, iy, Point2D(p_x, p_y), lcont_pos);

	const double bbox_hhx = hhx + p_r;
	const double bbox_hhy = hhy + p_r;
	if (lcont_pos.x < -bbox_hhx || lcont_pos.x > bbox_hhx ||
		lcont_pos.y < -bbox_hhy || lcont_pos.y > bbox_hhy)
		return false;

	// four corners
	double dx, dy, dist_tmp;
	// emperical
#define p_r_adjust_ratio 1.0
	if (lcont_pos.x < -hhx && lcont_pos.y < -hhy)
	{
		dx = lcont_pos.x + hhx;
		dy = lcont_pos.y + hhy;
		dist_tmp = sqrt(dx * dx + dy * dy);
		overlap_dist = p_r_adjust_ratio * p_r - dist_tmp;
		if (overlap_dist < 0.0)
			return false;
		lnorm.x = dx / dist_tmp;
		lnorm.y = dy / dist_tmp;
		return true;
	}
	if (lcont_pos.x > hhx && lcont_pos.y < -hhy)
	{
		dx = lcont_pos.x - hhx;
		dy = lcont_pos.y + hhy;
		dist_tmp = sqrt(dx * dx + dy * dy);
		overlap_dist = p_r_adjust_ratio * p_r - dist_tmp;
		if (overlap_dist < 0.0)
			return false;
		lnorm.x = dx / dist_tmp;
		lnorm.y = dy / dist_tmp;
		return true;
	}
	if (lcont_pos.x > hhx && lcont_pos.y > hhy)
	{
		dx = lcont_pos.x - hhx;
		dy = lcont_pos.y - hhy;
		dist_tmp = sqrt(dx * dx + dy * dy);
		overlap_dist = p_r_adjust_ratio * p_r - dist_tmp;
		if (overlap_dist < 0.0)
			return false;
		lnorm.x = dx / dist_tmp;
		lnorm.y = dy / dist_tmp;
		return true;
	}
	if (lcont_pos.x < -hhx && lcont_pos.y > hhy)
	{
		dx = lcont_pos.x + hhx;
		dy = lcont_pos.y - hhy;
		dist_tmp = sqrt(dx * dx + dy * dy);
		overlap_dist = p_r_adjust_ratio * p_r - dist_tmp;
		if (overlap_dist < 0.0)
			return false;
		lnorm.x = dx / dist_tmp;
		lnorm.y = dy / dist_tmp;
		return true;
	}

	overlap_dist = lcont_pos.x + bbox_hhx;
	size_t norm_id = 0;
	dist_tmp = bbox_hhx - lcont_pos.x;
	if (overlap_dist > dist_tmp)
	{
		overlap_dist = dist_tmp;
		norm_id = 1;
	}
	dist_tmp = lcont_pos.y + bbox_hhy;
	if (overlap_dist > dist_tmp)
	{
		overlap_dist = dist_tmp;
		norm_id = 2;
	}
	dist_tmp = bbox_hhy - lcont_pos.y;
	if (overlap_dist > dist_tmp)
	{
		overlap_dist = dist_tmp;
		norm_id = 3;
	}

	const Vector2D& en = edge_normals[norm_id];
	lnorm.x = en.x;
	lnorm.y = en.y;
	return true;
}

//bool RigidRect::detect_collision_with_point(
//	double p_x, double p_y, double p_r,
//	double& overlap_dist,
//	Vector2D& lnorm, Point2D& lcont_pos)
//{
//	const Vector2D ix(cos_ang, sin_ang), iy(-sin_ang, cos_ang);
//	point_from_global_to_local_coordinate(centre,
//		ix, iy, Point2D(p_x, p_y), lcont_pos);
//
//	// need rotation when change coordinate
//	const double bbox_hhx = hhx + p_r;
//	const double bbox_hhy = hhy + p_r;
//	if (lcont_pos.x < -bbox_hhx || lcont_pos.x > bbox_hhx ||
//		lcont_pos.y < -bbox_hhy || lcont_pos.y > bbox_hhy)
//		return false;
//
//	double dist_tmp;
//	overlap_dist = lcont_pos.x + bbox_hhx;
//	size_t norm_id = 0;
//	dist_tmp = bbox_hhx - lcont_pos.x;
//	if (overlap_dist > dist_tmp)
//	{
//		overlap_dist = dist_tmp;
//		norm_id = 1;
//	}
//	dist_tmp = lcont_pos.y + bbox_hhy;
//	if (overlap_dist > dist_tmp)
//	{
//		overlap_dist = dist_tmp;
//		norm_id = 2;
//	}
//	dist_tmp = bbox_hhy - lcont_pos.y;
//	if (overlap_dist > dist_tmp)
//	{
//		overlap_dist = dist_tmp;
//		norm_id = 3;
//	}
//
//	const Vector2D& en = edge_normals[norm_id];
//	lnorm.x = en.x;
//	lnorm.y = en.y;
//	return true;
//}
