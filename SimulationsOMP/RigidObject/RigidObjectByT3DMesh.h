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