#include "SimulationsOMP_pcp.h"

#include <float.h>
#include <iostream>

#include "SearchGrid3DTetrahedron.hpp"
#include "ExtractSurfaceFromT3DMesh.h"
#include "SearchGrid3DTriangle.hpp"
#include "SimulationsOMPUtils.h"
#include "RigidMeshT3D.h"

RigidMeshT3D::RigidMeshT3D() :
	grid_pos_type(nullptr),
	face_in_grid_range(nullptr),
	face_in_grid_list(nullptr),
	pt_tri_dist(nullptr),
	max_dist(DBL_MAX) {}

RigidMeshT3D::~RigidMeshT3D() { clear(); }

int RigidMeshT3D::init_from_mesh(
	TetrahedronMesh& tmesh,
	double ghx,
	double ghy,
	double ghz
	)
{
	clear();

	const Cube& mh_bbox = tmesh.get_bounding_box();
	grid.alloc_grid(
		mh_bbox.xl, mh_bbox.yl, mh_bbox.zl,
		mh_bbox.xu, mh_bbox.yu, mh_bbox.zu,
		ghx, ghy, ghz, 0.01);
	grid_xu_min_hx = grid.xu - grid.hx;
	grid_yu_min_hy = grid.yu - grid.hy;
	grid_zu_min_hz = grid.zu - grid.hz;
	grid_max_x_id = grid.x_num - 1;
	grid_max_y_id = grid.y_num - 1;
	grid_max_z_id = grid.z_num - 1;

	// init grid_pos_type
	grid_pos_type = new GridPosType[grid.num];
	for (size_t g_id = 0; g_id < grid.num; ++g_id)
		grid_pos_type[g_id] = GridPosType::Outside;

	// apply tetrahedron to mesh
	SearchGrid3DTetrahedron<TetrahedronMesh::Element> teh_grid;
	teh_grid.alloc_grid<Grid3D<>>(grid);
	teh_grid.apply_tetrahedrons<TetrahedronMesh::Node>(
		tmesh.get_nodes(), tmesh.get_elems(), tmesh.get_elem_num());
	const SearchGrid3DTetrahedron<TetrahedronMesh::Element>::Grid *teh_grids = teh_grid.get_grids();
	for (size_t g_id = 0; g_id < grid.num; ++g_id)
	{
		auto& g = teh_grids[g_id];
		if (g.ptehs)
			grid_pos_type[g_id] = GridPosType::Inside;
		//std::cout << g_id << ", ";
		//for (auto pts = g.ptehs; pts; pts = pts->next)
		//	std::cout << pts->pteh - tmesh.get_elems() << ", ";
		//std::cout << "\n";
	}
	teh_grid.clear(); // release memory

	// boundary surface
	ExtractSurfaceFromT3DMesh surf_extractor;
	surf_extractor.init_from_tetrahedron(
		tmesh.get_elems(), tmesh.get_elem_num());
	face_num = surf_extractor.get_face_num();
	const ExtractSurfaceFromT3DMesh::Face* faces = surf_extractor.get_faces();
	const TetrahedronMesh::Node* nodes = tmesh.get_nodes();
	pt_tri_dist = new PointToTriangleDistance[face_num];
	for (size_t f_id = 0; f_id < face_num; ++f_id)
	{
		const ExtractSurfaceFromT3DMesh::Face& f = faces[f_id];
		pt_tri_dist[f_id].init_triangle(nodes[f.n1], nodes[f.n2], nodes[f.n3]);
		//std::cout << f_id << ": " << f.n1 << ", " << f.n2 << ", " << f.n3 << "\n";
	}

	// apply boundary surface
	SearchGrid3DTriangle<ExtractSurfaceFromT3DMesh::Face> surf_grid;
	surf_grid.alloc_grid<Grid3D<>>(grid);
	surf_grid.apply_triangles<TetrahedronMesh::Node>(
		tmesh.get_nodes(),
		surf_extractor.get_faces(),
		surf_extractor.get_face_num());
	const SearchGrid3DTriangle<ExtractSurfaceFromT3DMesh::Face>::Grid
		*surf_grids = surf_grid.get_grids();
	face_in_grid_range = new size_t[grid.num + 1];
	face_in_grid_range[0] = 0;
	size_t surf_num = 0;
	for (size_t g_id = 0; g_id < grid.num; ++g_id)
	{
		const auto &g = surf_grids[g_id];
		if (g.ptris)
			grid_pos_type[g_id] = GridPosType::AtBoundary;
		for (auto* ptri = g.ptris; ptri; ptri = ptri->next)
			++surf_num;
		face_in_grid_range[g_id + 1] = surf_num;
	}
	face_in_grid_list_len = surf_num;
	face_in_grid_list = new size_t[surf_num];
	surf_num = 0;
	for (size_t g_id = 0; g_id < grid.num; ++g_id)
	{
		const size_t gf_start_id = surf_num;
		const auto& g = surf_grids[g_id];
		for (auto* ptri = g.ptris; ptri; ptri = ptri->next)
		{
			face_in_grid_list[surf_num] = surf_extractor.get_face_id(*(ptri->ptri));
			++surf_num;
		}
		SimulationsOMP::swap_sort_acc<size_t>(
			face_in_grid_list + gf_start_id,
			surf_num - gf_start_id);
		//std::cout << g_id << ": ";
		//for (size_t f_id = gf_start_id; f_id < surf_num; f_id++)
		//	std::cout << face_in_grid_list[f_id] << ", ";
		//std::cout << "\n";
	}
	surf_grid.clear(); // release memory
	surf_extractor.clear();

	return 0;
}

void RigidMeshT3D::init_max_dist(double _dist) noexcept
{
	GridPosType* cur_g_pos;
	const size_t gn_x_num = grid.x_num + 1;
	const size_t gn_y_num = grid.y_num + 1;
	const size_t gn_z_num = grid.z_num + 1;
	
	if (max_dist != DBL_MAX)
	{
		// reset previous max_dist result
		cur_g_pos = grid_pos_type;
		for (size_t z_id = 0; z_id < grid.z_num; ++z_id)
			for (size_t y_id = 0; y_id < grid.y_num; ++y_id)
				for (size_t x_id = 0; x_id < grid.x_num; ++x_id)
				{
					*reinterpret_cast<unsigned char *>(cur_g_pos) &= 0100;
					++cur_g_pos;
				}
	}

	max_dist = _dist;

	const size_t gn_xy_num = gn_x_num * gn_y_num;
	double *gns = new double[gn_xy_num * gn_z_num];
	double *cur_gn = gns;
	double gn_r = grid.hx < grid.hy ? grid.hx : grid.hy;
	gn_r = gn_r < grid.hz ? gn_r : grid.hz;
	gn_r *= 0.5;
	Vector3D norm;
	Point3D gn_pt;
	gn_pt.z = grid.zl;
	for (size_t z_id = 0; z_id < gn_z_num; ++z_id)
	{
		gn_pt.y = grid.yl;
		for (size_t y_id = 0; y_id < gn_y_num; ++y_id)
		{
			gn_pt.x = grid.xl;
			for (size_t x_id = 0; x_id < gn_x_num; ++x_id)
			{
				detect_collision_with_point(gn_pt, gn_r, *cur_gn, norm);
				++cur_gn;
				gn_pt.x += grid.hx;
			}
			gn_pt.y += grid.hy;
		}
		gn_pt.z += grid.hz;
	}
	
	cur_g_pos = grid_pos_type;
	for (size_t z_id = 0; z_id < grid.z_num; ++z_id)
		for (size_t y_id = 0; y_id < grid.y_num; ++y_id)
			for (size_t x_id = 0; x_id < grid.x_num; ++x_id)
			{
				cur_gn = gns + (z_id * gn_xy_num + y_id * gn_x_num + x_id);
				const double n1_dist = *cur_gn;
				const double n2_dist = *(cur_gn + 1);
				const double n3_dist = *(cur_gn + gn_x_num);
				const double n4_dist = *(cur_gn + gn_x_num + 1);
				const double n5_dist = *(cur_gn + gn_xy_num);
				const double n6_dist = *(cur_gn + gn_xy_num + 1);
				const double n7_dist = *(cur_gn + gn_xy_num + gn_x_num);
				const double n8_dist = *(cur_gn + gn_xy_num + gn_x_num + 1);
				if (n1_dist <= -max_dist && n2_dist <= -max_dist &&
					n3_dist <= -max_dist && n4_dist <= -max_dist &&
					n5_dist <= -max_dist && n6_dist <= -max_dist &&
					n7_dist <= -max_dist && n8_dist <= -max_dist)
					*cur_g_pos = GridPosType::FarOutside;
				if (n1_dist >= max_dist && n2_dist >= max_dist &&
					n3_dist >= max_dist && n4_dist >= max_dist &&
					n5_dist >= max_dist && n6_dist >= max_dist &&
					n7_dist >= max_dist && n8_dist >= max_dist)
					*cur_g_pos = GridPosType::FarInside;
				++cur_g_pos;
			}

	delete[] gns;
}

bool RigidMeshT3D::detect_collision_with_point(
	const Point3D &pt,
	const double p_r,
	double &dist,
	Vector3D &norm
	) const noexcept
{
	const double p_xl = pt.x - p_r, p_xu = pt.x + p_r;
	const double p_yl = pt.y - p_r, p_yu = pt.y + p_r;
	const double p_zl = pt.z - p_r, p_zu = pt.z + p_r;
	if (p_xu < grid.xl || p_xl > grid.xu ||
		p_yu < grid.yl || p_yl > grid.yu ||
		p_zu < grid.zl || p_zl > grid.zu)
		return false;

	const size_t p_x_id = pt.x < 0.0 ? 0 : (pt.x > grid_xu_min_hx ? grid_max_x_id : size_t((pt.x - grid.xl) / grid.hx));
	const size_t p_y_id = pt.y < 0.0 ? 0 : (pt.y > grid_yu_min_hy ? grid_max_y_id : size_t((pt.y - grid.yl) / grid.hy));
	const size_t p_z_id = pt.z < 0.0 ? 0 : (pt.z > grid_zu_min_hz ? grid_max_z_id : size_t((pt.z - grid.zl) / grid.hz));
	const size_t p_off = offset_from_xyz_id(p_x_id, p_y_id, p_z_id);
	if (unsigned char(grid_pos_type[p_off]) > 2) // far inside or outside mesh
		return false;
	
	IdRange id_range;
	id_range.xl_id = p_xl < 0.0 ? 0 : size_t((p_xl - grid.xl) / grid.hx);
	id_range.yl_id = p_yl < 0.0 ? 0 : size_t((p_yl - grid.yl) / grid.hy);
	id_range.zl_id = p_zl < 0.0 ? 0 : size_t((p_zl - grid.zl) / grid.hz);
	id_range.xu_id = p_xu > grid_xu_min_hx ? grid.x_num : size_t(ceil((p_xu - grid.xl) / grid.hx));
	id_range.yu_id = p_yu > grid_yu_min_hy ? grid.y_num : size_t(ceil((p_yu - grid.yl) / grid.hy));
	id_range.zu_id = p_zu > grid_zu_min_hz ? grid.z_num : size_t(ceil((p_zu - grid.zl) / grid.hz));

	// cal stride
	size_t id_len = (id_range.xu_id - id_range.xl_id) > (id_range.yu_id - id_range.yl_id) ?
					(id_range.xu_id - id_range.xl_id) : (id_range.yu_id - id_range.yl_id);
	id_len = id_len > (id_range.zu_id - id_range.zl_id) ?
			 id_len : (id_range.zu_id - id_range.zl_id);
	size_t stride = 1;
	while (stride < id_len)
		stride <<= 1;

	// cal dist
	ClosestFaceRes closest_face;
	closest_face.face_id = SIZE_MAX;
	closest_face.distance = max_dist;
	search_closest_face(closest_face, pt, stride, id_range);

	// cal normal
	if (closest_face.face_id != SIZE_MAX)
	{
		pt_tri_dist[closest_face.face_id].cal_normal_to_point(
			pt, closest_face.normal_type, norm);
		if (closest_face.distance >= 0.0) // inside object
			norm.reverse();
		dist = closest_face.distance;
		return true;
	}
	return false;
}

void RigidMeshT3D::search_closest_face(
	ClosestFaceRes& closest_face,
	const Point3D &pt,
	size_t stride,
	const IdRange &id_range
	) const noexcept
{
	if (stride > 1)
	{
		stride >>= 1;
		const size_t xm_id = id_range.xl_id + stride;
		const size_t ym_id = id_range.yl_id + stride;
		const size_t zm_id = id_range.zl_id + stride;
		const double xl = grid.xl + double(id_range.xl_id) * grid.hx;
		const double yl = grid.yl + double(id_range.yl_id) * grid.hy;
		const double zl = grid.zl + double(id_range.zl_id) * grid.hz;
		const double xm = grid.xl + double(xm_id) * grid.hx;
		const double ym = grid.yl + double(ym_id) * grid.hy;
		const double zm = grid.zl + double(zm_id) * grid.hz;
		const double xu = grid.xl + double(id_range.xu_id) * grid.hx;
		const double yu = grid.yl + double(id_range.yu_id) * grid.hy;
		const double zu = grid.zl + double(id_range.zu_id) * grid.hz;
		IdDist child_id_dists[8];
		IdRange child_id_ranges[8];
		Cube child_box;
#define Set_Id_Range(cid, xld, xud, yld, yud, zld, zud) \
		IdRange &cir##cid = child_id_ranges[cid]; \
		cir##cid.xl_id = xld; \
		cir##cid.xu_id = xud; \
		cir##cid.yl_id = yld; \
		cir##cid.yu_id = yud; \
		cir##cid.zl_id = zld; \
		cir##cid.zu_id = zud
#define Set_Id_Dist(_cid, _xl, _xu, _yl, _yu, _zl, _zu) \
		child_box.xl = _xl; \
		child_box.xu = _xu; \
		child_box.yl = _yl; \
		child_box.yu = _yu; \
		child_box.zl = _zl; \
		child_box.zu = _zu; \
		IdDist &cid##_cid = child_id_dists[_cid]; \
		cid##_cid.id = _cid; \
		cid##_cid.dist = cal_cube_point_distance(child_box, pt)
#define Search_Closest_Face(child_id) \
		if (cid##child_id.dist >= abs(closest_face.distance)) \
			return; \
		search_closest_face(closest_face, pt, \
			stride, child_id_ranges[cid##child_id.id])
		if (xm_id < id_range.xu_id)
		{
			if (ym_id < id_range.yu_id)
			{
				if (zm_id < id_range.zu_id)
				{
					// 0
					Set_Id_Range(0,
						id_range.xl_id, xm_id,
						id_range.yl_id, ym_id,
						id_range.zl_id, zm_id);
					Set_Id_Dist(0, xl, xm, yl, ym, zl, zm);
					// 1
					Set_Id_Range(1,
						xm_id, id_range.xu_id,
						id_range.yl_id, ym_id,
						id_range.zl_id, zm_id);
					Set_Id_Dist(1, xm, xu, yl, ym, zl, zm);
					// 2
					Set_Id_Range(2,
						id_range.xl_id, xm_id,
						ym_id, id_range.yu_id,
						id_range.zl_id, zm_id);
					Set_Id_Dist(2, xl, xm, ym, yu, zl, zm);
					// 3
					Set_Id_Range(3,
						xm_id, id_range.xu_id,
						ym_id, id_range.yu_id,
						id_range.zl_id, zm_id);
					Set_Id_Dist(3, xm, xu, ym, yu, zl, zm);
					// 4
					Set_Id_Range(4,
						id_range.xl_id, xm_id,
						id_range.yl_id, ym_id,
						zm_id, id_range.zu_id);
					Set_Id_Dist(4, xl, xm, yl, ym, zm, zu);
					// 5
					Set_Id_Range(5,
						xm_id, id_range.xu_id,
						id_range.yl_id, ym_id,
						zm_id, id_range.zu_id);
					Set_Id_Dist(5, xm, xu, yl, ym, zm, zu);
					// 6
					Set_Id_Range(6,
						id_range.xl_id, xm_id,
						ym_id, id_range.yu_id,
						zm_id, id_range.zu_id);
					Set_Id_Dist(6, xl, xm, ym, yu, zm, zu);
					// 7
					Set_Id_Range(7,
						xm_id, id_range.xu_id,
						ym_id, id_range.yu_id,
						zm_id, id_range.zu_id);
					Set_Id_Dist(7, xm, xu, ym, yu, zm, zu);
					//
					sort_acc_id_dist_pairs_8(child_id_dists);
					Search_Closest_Face(0);
					Search_Closest_Face(1);
					Search_Closest_Face(2);
					Search_Closest_Face(3);
					Search_Closest_Face(4);
					Search_Closest_Face(5);
					Search_Closest_Face(6);
					Search_Closest_Face(7);
				}
				else // zm_id > zu_id
				{
					// 0
					Set_Id_Range(0,
						id_range.xl_id, xm_id,
						id_range.yl_id, ym_id,
						id_range.zl_id, id_range.zu_id);
					Set_Id_Dist(0, xl, xm, yl, ym, zl, zu);
					// 1
					Set_Id_Range(1,
						xm_id, id_range.xu_id,
						id_range.yl_id, ym_id,
						id_range.zl_id, id_range.zu_id);
					Set_Id_Dist(1, xm, xu, yl, ym, zl, zu);
					// 2
					Set_Id_Range(2,
						id_range.xl_id, xm_id,
						ym_id, id_range.yu_id,
						id_range.zl_id, id_range.zu_id);
					Set_Id_Dist(2, xl, xm, ym, yu, zl, zu);
					// 3
					Set_Id_Range(3,
						xm_id, id_range.xu_id,
						ym_id, id_range.yu_id,
						id_range.zl_id, id_range.zu_id);
					Set_Id_Dist(3, xm, xu, ym, yu, zl, zu);
					//
					sort_acc_id_dist_pairs_4(child_id_dists);
					Search_Closest_Face(0);
					Search_Closest_Face(1);
					Search_Closest_Face(2);
					Search_Closest_Face(3);
				}
			}
			else // ym_id >= yu_id
			{
				if (zm_id < id_range.zu_id)
				{
					// 0
					Set_Id_Range(0,
						id_range.xl_id, xm_id,
						id_range.yl_id, id_range.yu_id,
						id_range.zl_id, zm_id);
					Set_Id_Dist(0, xl, xm, yl, yu, zl, zm);
					// 1
					Set_Id_Range(1,
						xm_id, id_range.xu_id,
						id_range.yl_id, id_range.yu_id,
						id_range.zl_id, zm_id);
					Set_Id_Dist(1, xm, xu, yl, yu, zl, zm);
					// 2
					Set_Id_Range(2,
						id_range.xl_id, xm_id,
						id_range.yl_id, id_range.yu_id,
						zm_id, id_range.zu_id);
					Set_Id_Dist(2, xl, xm, yl, yu, zm, zu);
					// 3
					Set_Id_Range(3,
						xm_id, id_range.xu_id,
						id_range.yl_id, id_range.yu_id,
						zm_id, id_range.zu_id);
					Set_Id_Dist(3, xm, xu, yl, yu, zm, zu);
					//
					sort_acc_id_dist_pairs_4(child_id_dists);
					Search_Closest_Face(0);
					Search_Closest_Face(1);
					Search_Closest_Face(2);
					Search_Closest_Face(3);
				}
				else // ym_id >= yu_id && zm_id >= zu_id
				{
					// 0
					Set_Id_Range(0,
						id_range.xl_id, xm_id,
						id_range.yl_id, id_range.yu_id,
						id_range.zl_id, id_range.zu_id);
					Set_Id_Dist(0, xl, xm, yl, yu, zl, zu);
					// 1
					Set_Id_Range(1,
						xm_id, id_range.xu_id,
						id_range.yl_id, id_range.yu_id,
						id_range.zl_id, id_range.zu_id);
					Set_Id_Dist(1, xm, xu, yl, yu, zl, zu);
					//
					sort_acc_id_dist_pairs_2(child_id_dists);
					Search_Closest_Face(0);
					Search_Closest_Face(1);
				}
			}
		}
		else // xm_id >= xu_id
		{
			if (ym_id < id_range.yu_id)
			{
				if (zm_id < id_range.zu_id)
				{
					// 0
					Set_Id_Range(0,
						id_range.xl_id, id_range.xu_id,
						id_range.yl_id, ym_id,
						id_range.zl_id, zm_id);
					Set_Id_Dist(0, xl, xu, yl, ym, zl, zm);
					// 1
					Set_Id_Range(1,
						id_range.xl_id, id_range.xu_id,
						ym_id, id_range.yu_id,
						id_range.zl_id, zm_id);
					Set_Id_Dist(1, xl, xu, ym, yu, zl, zm);
					// 2
					Set_Id_Range(2,
						id_range.xl_id, id_range.xu_id,
						id_range.yl_id, ym_id,
						zm_id, id_range.zu_id);
					Set_Id_Dist(2, xl, xu, yl, ym, zm, zu);
					// 3
					Set_Id_Range(3,
						id_range.xl_id, id_range.xu_id,
						ym_id, id_range.yu_id,
						zm_id, id_range.zu_id);
					Set_Id_Dist(3, xl, xu, ym, yu, zm, zu);
					//
					sort_acc_id_dist_pairs_4(child_id_dists);
					Search_Closest_Face(0);
					Search_Closest_Face(1);
					Search_Closest_Face(2);
					Search_Closest_Face(3);
				}
				else // xm_id >= xu_id && zm >= zu_id
				{
					// 0
					Set_Id_Range(0,
						id_range.xl_id, id_range.xu_id,
						id_range.yl_id, ym_id,
						id_range.zl_id, id_range.zu_id);
					Set_Id_Dist(0, xl, xu, yl, ym, zl, zu);
					// 1
					Set_Id_Range(1,
						id_range.xl_id, id_range.xu_id,
						ym_id, id_range.yu_id,
						id_range.zl_id, id_range.zu_id);
					Set_Id_Dist(1, xl, xu, ym, yu, zl, zu);
					//
					sort_acc_id_dist_pairs_2(child_id_dists);
					Search_Closest_Face(0);
					Search_Closest_Face(1);
				}
			}
			else // ym_id >= yu_id && xm_id >= xu_id
			{
				if (zm_id < id_range.zu_id)
				{
					// 0
					Set_Id_Range(0,
						id_range.xl_id, id_range.xu_id,
						id_range.yl_id, id_range.yu_id,
						id_range.zl_id, zm_id);
					Set_Id_Dist(0, xl, xu, yl, yu, zl, zm);
					// 1
					Set_Id_Range(1,
						id_range.xl_id, id_range.xu_id,
						id_range.yl_id, id_range.yu_id,
						zm_id, id_range.zu_id);
					Set_Id_Dist(1, xl, xu, yl, yu, zm, zu);
					//
					sort_acc_id_dist_pairs_2(child_id_dists);
					Search_Closest_Face(0);
					Search_Closest_Face(1);
				}
				else // xm_id >= xu_id &&
				// ym_id >= yu_id && zm_id >= zu_id
				{
					Set_Id_Range(0,
						id_range.xl_id, id_range.xu_id,
						id_range.yl_id, id_range.yu_id,
						id_range.zl_id, id_range.zu_id);
					Set_Id_Dist(0, xl, xu, yl, yu, zl, zu);
					Search_Closest_Face(0);
				}
			}
		}
		return;
	}

	double dist_tmp;
	unsigned char norm_type;
	const size_t g_off = offset_from_xyz_id(id_range.xl_id, id_range.yl_id, id_range.zl_id);
	const size_t end_f_id = face_in_grid_range[g_off + 1];
	for (size_t f_id = face_in_grid_range[g_off];
		 f_id < end_f_id; ++f_id)
	{
		const PointToTriangleDistance& ptd = pt_tri_dist[face_in_grid_list[f_id]];
		norm_type = ptd.cal_distance_to_point(pt, dist_tmp);
		if (abs(closest_face.distance) > abs(dist_tmp))
		{
			closest_face.face_id = face_in_grid_list[f_id];
			closest_face.distance = dist_tmp;
			closest_face.normal_type = norm_type;
		}
	}
}
