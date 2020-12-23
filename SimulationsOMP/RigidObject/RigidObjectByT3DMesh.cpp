#include "SimulationsOMP_pcp.h"

#include "GeometryUtils.h"
#include "TetrahedronMesh.h"
#include "RigidObjectByT3DMesh.h"

RigidObjectByT3DMesh::RigidObjectByT3DMesh() :
	RigidObjectMotion3D(), RigidMeshT3D() {}

RigidObjectByT3DMesh::~RigidObjectByT3DMesh() {}

int RigidObjectByT3DMesh::init(
	double _density,
	const char* filename,
	double dx,
	double dy,
	double dz,
	// in degree
	double dx_ang,
	double dy_ang,
	double dz_ang,
	double ghx,
	double ghy,
	double ghz
	)
{
	TetrahedronMesh tmesh;
	int res = tmesh.load_mesh_from_hdf5(filename);
	if (res)
		return res;

	const Point3D& cen = tmesh.get_centre();
	const double mh_x = cen.x + dx;
	const double mh_y = cen.y + dy;
	const double mh_z = cen.z + dz;

	tmesh.translate_mesh(-cen.x, -cen.y, -cen.z);
	
	density = _density;
	// moment of intertia
	size_t elem_num = tmesh.get_elem_num();
	auto* nodes = tmesh.get_nodes();
	auto* elems = tmesh.get_elems();
	double elem_moi_data[6];
	double moi_data[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		auto &e = elems[e_id];
		auto &n1 = nodes[e.n1];
		auto &n2 = nodes[e.n2];
		auto &n3 = nodes[e.n3];
		auto &n4 = nodes[e.n4];
		cal_tetrahedron_moi(0.0, 0.0, 0.0,
			n1, n2, n3, n4, e.vol, elem_moi_data);
		moi_data[0] += density * elem_moi_data[0];
		moi_data[1] += density * elem_moi_data[1];
		moi_data[2] += density * elem_moi_data[2];
		moi_data[3] += density * elem_moi_data[3];
		moi_data[5] += density * elem_moi_data[5];
		moi_data[4] += density * elem_moi_data[4];
	}
	RigidObjectMotion3D::init(mh_x, mh_y, mh_z,
		tmesh.get_vol() * density, moi_data);

	dx_ang = deg_to_rad(dx_ang);
	dy_ang = deg_to_rad(dy_ang);
	dz_ang = deg_to_rad(dz_ang);
	RigidObjectMotion3D::set_angle(dx_ang, dy_ang, dz_ang);
	
	RigidMeshT3D::init_from_mesh(tmesh, ghx, ghy, ghz);
	return 0;
}

bool RigidObjectByT3DMesh::detect_collision_with_point(
	double p_x,
	double p_y,
	double p_z,
	double p_r,
	double& dist,
	Vector3D& lnorm,
	Point3D& lcontpos
	) const noexcept
{
	Point3D gp(p_x, p_y, p_z);
	get_local_point(gp, lcontpos);
	if (RigidMeshT3D::detect_collision_with_point(lcontpos, p_r, dist, lnorm))
		return true;
	return false;
}

void RigidObjectByT3DMesh::init_from_hdf5_res_file(
	const double _density,
	const double _m,
	const double _moi[6],
	const double _inv_moi[6],
	const double _T_mat[9],
	const Vector3D& acc,
	const Vector3D& acc_ang,
	const Vector3D& vec,
	const Vector3D& vec_ang,
	const Point3D& pos,
	const Vector3D& pos_ang,
	const Force3D& force_cont,
	const Force3D& force_ext,
	const double grid_xl,
	const double grid_yl,
	const double grid_zl,
	const double grid_hx,
	const double grid_hy,
	const double grid_hz,
	const unsigned long long grid_x_num,
	const unsigned long long grid_y_num,
	const unsigned long long grid_z_num,
	const unsigned long long _face_in_grid_list_len,
	const unsigned long long _face_num,
	unsigned char*& _grid_pos_type,
	unsigned long long*& _face_in_grid_range,
	unsigned long long*& _face_in_grid_list,
	PointToTriangleDistance*& _pt_tri_dist)
{
	clear();
	// rigid object motion
	density = _density; m = _m;
	moi[0] = _moi[0]; moi[1] = _moi[1]; moi[2] = _moi[2];
	moi[3] = _moi[3]; moi[4] = _moi[4]; moi[5] = _moi[5];
	inv_moi[0] = _inv_moi[0]; inv_moi[1] = _inv_moi[1];
	inv_moi[2] = _inv_moi[2]; inv_moi[3] = _inv_moi[3];
	inv_moi[4] = _inv_moi[4]; inv_moi[5] = _inv_moi[5];
	T_mat_data[0] = _T_mat[0]; T_mat_data[1] = _T_mat[1];
	T_mat_data[1] = _T_mat[1]; T_mat_data[2] = _T_mat[2];
	T_mat_data[3] = _T_mat[3]; T_mat_data[4] = _T_mat[4];
	T_mat_data[5] = _T_mat[5]; T_mat_data[6] = _T_mat[6];
	T_mat_data[7] = _T_mat[7]; T_mat_data[8] = _T_mat[8];
	ax = acc.x; ay = acc.y; az = acc.z;
	ax_ang = acc_ang.x; ay_ang = acc_ang.y; az_ang = acc_ang.z;
	vx = vec.x; vy = vec.y; vz = vec.z;
	vx_ang = vec_ang.x; vy_ang = vec_ang.y; vz_ang = vec_ang.z;
	x_ori = pos.x; y_ori = pos.y; z_ori = pos.z;
	ux = 0.0; uy = 0.0; uz = 0.0;
	x = x_ori; y = y_ori; z = z_ori;
	x_ang = pos_ang.x; y_ang = pos_ang.y; z_ang = pos_ang.z;
	fx_ext = force_ext.fx; fy_ext = force_ext.fy;
	fz_ext = force_ext.fz; mx_ext = force_ext.mx;
	my_ext = force_ext.my; mz_ext = force_ext.mz;
	fx_cont = force_cont.fx; fy_cont = force_cont.fy;
	fz_cont = force_cont.fz; mx_cont = force_cont.mx;
	my_cont = force_cont.my; mz_cont = force_cont.mz;
	rotate_axses_by_angle(pos_ang, ix, iy, iz);
	// bg grid
	grid.alloc_grid(grid_xl, grid_yl, grid_zl,
		grid_hx, grid_hy, grid_hz,
		grid_x_num, grid_y_num, grid_z_num);
	face_in_grid_list_len = _face_in_grid_list_len;
	face_num = _face_num;
	grid_pos_type = new GridPosType[grid.num];
	face_in_grid_range = new size_t[grid.num + 1];
	face_in_grid_list = new size_t[face_in_grid_list_len];
	pt_tri_dist = new PointToTriangleDistance[face_num];
	_grid_pos_type = reinterpret_cast<unsigned char *>(grid_pos_type);
	_face_in_grid_range = reinterpret_cast<size_t *>(face_in_grid_range);
	_face_in_grid_list = reinterpret_cast<size_t *>(face_in_grid_list);
	_pt_tri_dist = pt_tri_dist;
}
