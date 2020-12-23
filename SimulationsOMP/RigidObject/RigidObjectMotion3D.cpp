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
	double _moi_data[6]
	)
{
	x_ori = _x;	y_ori = _y;	z_ori = _z;
	x = _x; y = _y;	z = _z;
	m = _m;
	moi[0] = _moi_data[0];	moi[1] = _moi_data[1];
	moi[2] = _moi_data[2];	moi[3] = _moi_data[3];
	moi[4] = _moi_data[4];	moi[5] = _moi_data[5];

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
	double _z_ang
	)
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
	double _vz
	)
{
	set_vx_bc(_vx);
	set_vy_bc(_vy);
	set_vz_bc(_vz);
	set_vx_ang_bc(0.0);
	set_vy_ang_bc(0.0);
	set_vz_ang_bc(0.0);
}
