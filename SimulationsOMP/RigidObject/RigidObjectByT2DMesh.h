#ifndef __Rigid_Object_By_T2D_Mesh_h__
#define __Rigid_Object_By_T2D_Mesh_h__

#include "RigidObjectMotion2D.h"
#include "RigidMeshT2D.h"

class RigidObjectByT2DMesh :
	public RigidObjectMotion2D,
	public RigidMeshT2D
{
protected:
	double density;

public:
	explicit RigidObjectByT2DMesh();
	~RigidObjectByT2DMesh();

	int init(double _density, const char* file_name,
		double dx, double dy, double dx_ang, // in degree
		double ghx,	double ghy);

	inline double get_density() const noexcept { return density; }
	inline void get_bbox(Rect &bbox) const noexcept
	{
		const double hx = (grid.xu - grid.xl) * 0.5;
		const double hy = (grid.yu - grid.yl) * 0.5;
		const double rx = abs(ix.x * hx) + abs(iy.x * hy);
		const double ry = abs(ix.y * hx) + abs(iy.y * hy);
		bbox.xl = x - rx;
		bbox.xu = x + rx;
		bbox.yl = y - ry;
		bbox.yu = y + ry;
	}

	inline void set_position(double _x, double _y) noexcept
	{
		x_ori = _x; y_ori = _y;
		ux = 0.0; uy = 0.0;
		x = _x; y = _y;
	}
	inline void set_pos_ang(double _ang) noexcept { ang = _ang; }
	inline void set_velocity(double _vx, double _vy, double _vang) noexcept
	{ vx = _vx; vy = _vy; v_ang = _vang; }

	inline void set_ext_force(const Force2D& cf) noexcept
	{ fx_ext = cf.fx; fy_ext = cf.fy;	m_cont = cf.m; }
	inline void set_ext_force(double fx, double fy, double m) noexcept
	{ fx_ext = fx; fy_ext = fy; m_cont = m; }
	inline void set_cont_force(const Force2D& cf) noexcept
	{ fx_cont = cf.fx; fy_cont = cf.fy; m_cont = cf.m; }
	inline void reset_cont_force() noexcept
	{ fx_cont = 0.0; fy_cont = 0.0; m_cont = 0.0; }

	using RigidObjectMotion2D::get_local_point;
	using RigidObjectMotion2D::get_global_point;
	using RigidObjectMotion2D::get_local_vector;
	using RigidObjectMotion2D::get_global_vector;
	using RigidMeshT2D::init_max_dist;

	bool detect_collision_with_point(
		double p_x, double p_y, double p_r,
		double& dist, Vector2D& lnorm, Point2D& lcontpos
		) const noexcept;

	void init_from_hdf5_res_file(
		const double _density, const double _m,
		const double _moi, const double _inv_moi,
		const double _T_mat[4],
		const Vector2D &acc, const double acc_ang,
		const Vector2D &vec, const double vec_ang,
		const Point2D &pos, const double ang,
		const Force2D &force_cont, const Force2D &force_ext,
		const double grid_xl, const double grid_yl,
		const double grid_hx, const double grid_hy,
		const unsigned long long grid_x_num, const unsigned long long grid_y_num,
		const unsigned long long _line_in_grid_list_len,
		const unsigned long long _line_num,
		unsigned char *&grid_pos_type, unsigned long long *&line_in_grid_range,
		unsigned long long *&line_in_grid_list, PointToLineDistance *&pt_ln_dist);
};

#endif