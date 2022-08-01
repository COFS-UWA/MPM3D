#ifndef __Rigid_Object_By_T3D_Mesh_h__
#define __Rigid_Object_By_T3D_Mesh_h__

#include "RigidObjectMotion3D.h"
#include "RigidMeshT3D.h"

class RigidObjectByT3DMesh :
	public RigidObjectMotion3D,
	public RigidMeshT3D
{
protected:
	double density;

public:
	explicit RigidObjectByT3DMesh();
	~RigidObjectByT3DMesh();

	int init(double _density, const char* file_name,
		double dx, double dy, double dz,
		double dx_ang, double dy_ang, double dz_ang, // in degree
		double ghx,	double ghy, double ghz);

	inline double get_density() const noexcept { return density; }
	inline void get_bbox(Cube& bbox) const noexcept
	{
		const double hx = (grid.xu - grid.xl) * 0.5;
		const double hy = (grid.yu - grid.yl) * 0.5;
		const double hz = (grid.zu - grid.zl) * 0.5;
		const double rx = abs(ix.x * hx) + abs(iy.x * hy) + abs(iz.x * hz);
		const double ry = abs(ix.y * hx) + abs(iy.y * hy) + abs(iz.y * hz);
		const double rz = abs(ix.z * hx) + abs(iy.z * hy) + abs(iz.z * hz);
		bbox.xl = x - rx;
		bbox.xu = x + rx;
		bbox.yl = y - ry;
		bbox.yu = y + ry;
		bbox.zl = z - rz;
		bbox.zu = z + rz;
	}

	inline void set_position(double _x, double _y, double _z) noexcept
	{
		x_ori = _x; y_ori = _y; z_ori = _z;
		ux = 0.0; uy = 0.0; uz = 0.0;
		x = _x; y = _y; z = _z;
	}
	inline void set_pos_ang(double _x_ang, double _y_ang, double _z_ang) noexcept
	{
		x_ang = _x_ang; y_ang = _y_ang; z_ang = _z_ang;
	}
	inline void set_velocity(double _vx, double _vy, double _vz) noexcept
	{ vx = _vx; vy = _vy; vz = _vz; }
	inline void set_v_ang(double _vx_ang, double _vy_ang, double _vz_ang) noexcept
	{ vx_ang = _vx_ang; vy_ang = _vy_ang; vz_ang = _vz_ang; }

	inline void set_cont_force(const Force3D& cf) noexcept
	{
		fx_cont = cf.fx; fy_cont = cf.fy; fz_cont = cf.fz;
		mx_cont = cf.mx; my_cont = cf.my; mz_cont = cf.mz;
	}

	inline void reset_cont_force() noexcept
	{
		fx_cont = 0.0; fy_cont = 0.0; fz_cont = 0.0;
		mx_cont = 0.0; my_cont = 0.0; mz_cont = 0.0;
	}

	using RigidObjectMotion3D::get_local_point;
	using RigidObjectMotion3D::get_global_point;
	using RigidObjectMotion3D::get_local_vector;
	using RigidObjectMotion3D::get_global_vector;
	using RigidMeshT3D::init_max_dist;

	bool detect_collision_with_point(
		double p_x, double p_y, double p_z, double p_r,
		double& dist, Vector3D& lnorm, Point3D& lcontpos
		) const noexcept;

	void init_from_hdf5_res_file(
		const double _density, const double _m,
		const double _moi[6], const double _inv_moi[6],
		const double _T_mat[9],
		const Vector3D &acc, const Vector3D &acc_ang,
		const Vector3D &vec, const Vector3D &vec_ang,
		const Point3D &pos, const Vector3D &pos_ang,
		const Force3D &force_cont, const Force3D &force_ext,
		const double grid_xl, const double grid_yl, const double grid_zl,
		const double grid_hx, const double grid_hy, const double grid_hz,
		const unsigned long long grid_x_num, const unsigned long long grid_y_num,
		const unsigned long long grid_z_num,
		const unsigned long long _face_in_grid_list_len,
		const unsigned long long _face_num,
		unsigned char *&grid_pos_type, unsigned long long *&face_in_grid_range,
		unsigned long long *&face_in_grid_list, PointToTriangleDistance *&pt_tri_dist);
};

#endif