#include "SimulationsOMP_pcp.h"

#include <Eigen/Eigen>

#include "RigidObjectMotion3D.h"

RigidObjectMotion3D::RigidObjectMotion3D() :
	ax(0.0), ay(0.0), az(0.0), ax_ang(0.0), ay_ang(0.0), az_ang(0.0),
	vx(0.0), vy(0.0), vz(0.0), vx_ang(0.0), vy_ang(0.0), vz_ang(0.0),
	ux(0.0), uy(0.0), uz(0.0), x_ang(0.0), y_ang(0.0), z_ang(0.0),
	ix(1.0, 0.0, 0.0), iy(0.0, 1.0, 0.0), iz(0.0, 0.0, 1.0),
	ax_bc_mask(0), ay_bc_mask(0), az_bc_mask(0),
	ax_bc(0.0), ay_bc(0.0), az_bc(0.0),
	ax_ang_bc_mask(0), ay_ang_bc_mask(0), az_ang_bc_mask(0),
	ax_ang_bc(0.0), ay_ang_bc(0.0), az_ang_bc(0.0),
	vx_bc_mask(0), vy_bc_mask(0), vz_bc_mask(0),
	vx_bc(0.0), vy_bc(0.0), vz_bc(0.0),
	vx_ang_bc_mask(0), vy_ang_bc_mask(0), vz_ang_bc_mask(0),
	vx_ang_bc(0.0), vy_ang_bc(0.0), vz_ang_bc(0.0),
	pvx_bc_ratio(&one_ratio), pvy_bc_ratio(&one_ratio), pvz_bc_ratio(&one_ratio),
	fx_ext(0.0), fy_ext(0.0), fz_ext(0.0),
	mx_ext(0.0), my_ext(0.0), mz_ext(0.0),
	fx_cont(0.0), fy_cont(0.0), fz_cont(0.0),
	mx_cont(0.0), my_cont(0.0), mz_cont(0.0) {}

RigidObjectMotion3D::~RigidObjectMotion3D() {}

void RigidObjectMotion3D::init(
	double _x,
	double _y,
	double _z,
	double _m,
	double _moi_data[6])
{
	x_ori = _x;	y_ori = _y;	z_ori = _z;
	x = _x; y = _y;	z = _z;
	m = _m;
	moi[0] = _moi_data[0];	moi[1] = _moi_data[1];
	moi[2] = _moi_data[2];	moi[3] = _moi_data[3];
	moi[4] = _moi_data[4];	moi[5] = _moi_data[5];
	
	cur_time = 0.0;

	Eigen::Matrix3d moi_mat;
	moi_mat << moi[0], moi[3], moi[5],
			   moi[3], moi[1], moi[4],
			   moi[5], moi[4], moi[2];

	Eigen::Matrix3d inv_moi_mat = moi_mat.inverse();
	inv_moi[0] = inv_moi_mat(0, 0);
	inv_moi[1] = inv_moi_mat(1, 1);
	inv_moi[2] = inv_moi_mat(2, 2);
	inv_moi[3] = inv_moi_mat(0, 1);
	inv_moi[4] = inv_moi_mat(1, 2);
	inv_moi[5] = inv_moi_mat(2, 0);
}

void RigidObjectMotion3D::set_angle(
	// in radius
	double _x_ang,
	double _y_ang,
	double _z_ang)
{
	x_ang = _x_ang;
	trim_to_pi(x_ang);
	y_ang = _y_ang;
	trim_to_pi(y_ang);
	z_ang = _z_ang;
	trim_to_pi(z_ang);
	ix.x = 1.0, ix.y = 0.0, ix.z = 0.0;
	iy.x = 0.0, iy.y = 1.0, iy.z = 0.0;
	iz.x = 0.0, iz.y = 0.0, iz.z = 1.0;
	rotate_axses_by_angle(pos_ang, ix, iy, iz);
}

void RigidObjectMotion3D::set_translation_velocity_bc(
	double _vx,
	double _vy,
	double _vz)
{
	set_vx_bc(_vx); set_vy_bc(_vy);	set_vz_bc(_vz);
	set_vx_ang_bc(0.0);	set_vy_ang_bc(0.0);	set_vz_ang_bc(0.0);
}

void RigidObjectMotion3D::update_motion(double dt) noexcept
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
	cur_time += dt;
	union
	{
		struct { double adj_vx_bc, adj_vy_bc, adj_vz_bc; };
		struct { size_t ivx_bc, ivy_bc, ivz_bc; };
	};
	adj_vx_bc = vx_bc * (*pvx_bc_ratio)(cur_time);
	adj_vy_bc = vy_bc * (*pvy_bc_ratio)(cur_time);
	adj_vz_bc = vz_bc * (*pvz_bc_ratio)(cur_time);
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
	vec_tmp2[0] = vec_tmp1[2] * vy_ang - vec_tmp1[1] * vz_ang;
	vec_tmp2[1] = vec_tmp1[0] * vz_ang - vec_tmp1[2] * vx_ang;
	vec_tmp2[2] = vec_tmp1[1] * vx_ang - vec_tmp1[0] * vy_ang;
	// m - w cross I*w
	vec_tmp2[0] = (mx_ext + mx_cont) * 1.0 - vec_tmp2[0];
	vec_tmp2[1] = (my_ext + my_cont) * 1.0 - vec_tmp2[1];
	vec_tmp2[2] = (mz_ext + mz_cont) * 1.0 - vec_tmp2[2];
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
	//x_ang += vx_ang * dt;
	trim_to_pi(x_ang);
	//y_ang += vy_ang * dt;
	trim_to_pi(y_ang);
	//z_ang += vz_ang * dt;
	trim_to_pi(z_ang);
	// update ix, iy, iz
	ix.x = 1.0, ix.y = 0.0, ix.z = 0.0;
	iy.x = 0.0, iy.y = 1.0, iy.z = 0.0;
	iz.x = 0.0, iz.y = 0.0, iz.z = 1.0;
	rotate_axses_by_angle(pos_ang, ix, iy, iz);
	ix.normalize();
	iy.normalize();
	iz.normalize();
}
