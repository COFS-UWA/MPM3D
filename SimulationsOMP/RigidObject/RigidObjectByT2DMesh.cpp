#include "SimulationsOMP_pcp.h"

#include "GeometryUtils.h"
#include "TriangleUtils.h"
#include "TriangleMesh.h"
#include "RigidObjectByT2DMesh.h"

RigidObjectByT2DMesh::RigidObjectByT2DMesh() :
	RigidObjectMotion2D(), RigidMeshT2D() {}

RigidObjectByT2DMesh::~RigidObjectByT2DMesh() {}

int RigidObjectByT2DMesh::init(
	double _density,
	const char* filename,
	double dx,
	double dy,
	double dang, // in degree
	double ghx,
	double ghy)
{
	TriangleMesh tmesh;
	int res = tmesh.load_mesh_from_hdf5(filename);
	if (res)
		return res;

	dang = deg_to_rad(dang);
	tmesh.rotate_mesh(dang);

	const Point2D& cen = tmesh.get_centre();
	const double mh_x = cen.x + dx;
	const double mh_y = cen.y + dy;
	tmesh.translate_mesh(-cen.x, -cen.y);
	
	density = _density;
	// moment of intertia
	const size_t elem_num = tmesh.get_elem_num();
	auto* nodes = tmesh.get_nodes();
	auto* elems = tmesh.get_elems();
	double elem_moi, moi = 0.0;
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		auto &e = elems[e_id];
		auto &n1 = nodes[e.n1];
		auto &n2 = nodes[e.n2];
		auto &n3 = nodes[e.n3];
		cal_triangle_moi(0.0, 0.0, n1, n2, n3, e.area, elem_moi);
		moi += density * elem_moi;
	}
	RigidObjectMotion2D::init(mh_x, mh_y, tmesh.get_area() * density, moi);
	
	RigidMeshT2D::init_from_mesh(tmesh, ghx, ghy);
	return 0;
}

bool RigidObjectByT2DMesh::detect_collision_with_point(
	double p_x,
	double p_y,
	double p_r,
	double& dist,
	Vector2D& lnorm,
	Point2D& lcontpos
	) const noexcept
{
	Point2D gp(p_x, p_y);
	get_local_point(gp, lcontpos);
	if (RigidMeshT2D::detect_collision_with_point(lcontpos, p_r, dist, lnorm))
		return true;
	return false;
}

void RigidObjectByT2DMesh::init_from_hdf5_res_file(
	const double _density,
	const double _m,
	const double _moi,
	const double _inv_moi,
	const double _T_mat[4],
	const Vector2D &acc,
	const double acc_ang,
	const Vector2D &vec,
	const double vec_ang,
	const Point2D &pos,
	const double pos_ang,
	const Force2D& force_cont,
	const Force2D& force_ext,
	const double grid_xl,
	const double grid_yl,
	const double grid_hx,
	const double grid_hy,
	const unsigned long long grid_x_num,
	const unsigned long long grid_y_num,
	const unsigned long long _edge_in_grid_list_len,
	const unsigned long long _edge_num,
	unsigned char*& _grid_pos_type,
	unsigned long long*& _edge_in_grid_range,
	unsigned long long*& _edge_in_grid_list,
	PointToLineDistance*& _pt_ln_dist)
{
	clear();
	// rigid object motion
	density = _density;
	m = _m;
	moi = _moi;
	inv_moi = _inv_moi;
	T_mat_data[0] = _T_mat[0];
	T_mat_data[1] = _T_mat[1];
	T_mat_data[2] = _T_mat[2];
	T_mat_data[3] = _T_mat[3];
	ax = acc.x; ay = acc.y;	a_ang = acc_ang;
	vx = vec.x; vy = vec.y; v_ang = vec_ang;
	x_ori = pos.x; y_ori = pos.y; 
	ux = 0.0; uy = 0.0;
	x = x_ori; y = y_ori; ang = pos_ang;
	fx_ext = force_ext.fx; fy_ext = force_ext.fy; m_ext = force_ext.m;
	fx_cont = force_cont.fx; fy_cont = force_cont.fy; m_cont = force_cont.m;
	rotate_axses_by_angle(pos_ang, ix, iy);
	// bg grid
	grid.alloc_grid(grid_xl, grid_yl, grid_hx, grid_hy, grid_x_num, grid_y_num);
	edge_in_grid_list_len = _edge_in_grid_list_len;
	edge_num = _edge_num;
	grid_pos_type = new GridPosType[grid.num];
	edge_in_grid_range = new size_t[grid.num + 1];
	edge_in_grid_list = new size_t[edge_in_grid_list_len];
	pt_ln_dist = new PointToLineDistance[edge_num];
	_grid_pos_type = reinterpret_cast<unsigned char *>(grid_pos_type);
	_edge_in_grid_range = reinterpret_cast<size_t *>(edge_in_grid_range);
	_edge_in_grid_list = reinterpret_cast<size_t *>(edge_in_grid_list);
	_pt_ln_dist = pt_ln_dist;
}
