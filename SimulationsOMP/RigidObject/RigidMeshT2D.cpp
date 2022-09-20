#include "SimulationsOMP_pcp.h"

#include <float.h>
#include <iostream>

#include "SearchGrid2DLine.hpp"
#include "SearchGrid2DTriangle.hpp"
#include "ExtractEdgeFromT2DMesh.h"
#include "SimulationsOMPUtils.h"
#include "RigidMeshT2D.h"

RigidMeshT2D::RigidMeshT2D() :
	grid_pos_type(nullptr),
	edge_in_grid_range(nullptr),
	edge_in_grid_list(nullptr),
	pt_ln_dist(nullptr) {}

RigidMeshT2D::~RigidMeshT2D() { clear(); }

int RigidMeshT2D::init_from_mesh(TriangleMesh& tmesh, double ghx, double ghy)
{
	clear();

	const Rect& mh_bbox = tmesh.get_bounding_box();
	grid.alloc_grid(
		mh_bbox.xl, mh_bbox.yl,
		mh_bbox.xu, mh_bbox.yu,
		ghx, ghy, 0.01);
	
	// init max_dist and max_stride
	max_dist = DBL_MAX;
	max_stride = get_max_stride();

	// init grid_pos_type
	grid_pos_type = new GridPosType[grid.num];
	for (size_t g_id = 0; g_id < grid.num; ++g_id)
		grid_pos_type[g_id] = GridPosType::Outside;

	// apply tetrahedron to mesh
	SearchGrid2DTriangle<TriangleMesh::Element> tri_grid;
	tri_grid.alloc_grid<Grid2D<>>(grid);
	tri_grid.apply_triangles<TriangleMesh::Node>(
		tmesh.get_nodes(), tmesh.get_elems(), tmesh.get_elem_num());
	const SearchGrid2DTriangle<TriangleMesh::Element>::Grid *tri_grids = tri_grid.get_grids();
	for (size_t g_id = 0; g_id < grid.num; ++g_id)
	{
		auto& g = tri_grids[g_id];
		if (g.ptris)
			grid_pos_type[g_id] = GridPosType::Inside;
	}
	tri_grid.clear(); // release memory

	// boundary surface
	ExtractEdgeFromT2DMesh edge_extractor;
	edge_extractor.init_from_triangle(
		tmesh.get_elems(), tmesh.get_elem_num());
	edge_num = edge_extractor.get_boundary_edge_num();
	const ExtractEdgeFromT2DMesh::Edge* edges = edge_extractor.get_boundary_edges();
	const TriangleMesh::Node* nodes = tmesh.get_nodes();
	pt_ln_dist = new PointToLineDistance[edge_num];
	for (size_t e_id = 0; e_id < edge_num; ++e_id)
	{
		const ExtractEdgeFromT2DMesh::Edge& e = edges[e_id];
		pt_ln_dist[e_id].init_line(nodes[e.n1], nodes[e.n2]);
	}

	// apply boundary surface
	SearchGrid2DLine<ExtractEdgeFromT2DMesh::Edge> edge_grid;
	edge_grid.alloc_grid<Grid2D<>>(grid);
	edge_grid.apply_lines<TriangleMesh::Node>(
		tmesh.get_nodes(),
		edge_extractor.get_boundary_edges(),
		edge_extractor.get_boundary_edge_num());
	const SearchGrid2DLine<ExtractEdgeFromT2DMesh::Edge>::Grid
		*edge_grids = edge_grid.get_grids();
	edge_in_grid_range = new size_t[grid.num + 1];
	edge_in_grid_range[0] = 0;
	size_t edge_num = 0;
	for (size_t g_id = 0; g_id < grid.num; ++g_id)
	{
		const auto &g = edge_grids[g_id];
		if (g.plns)
			grid_pos_type[g_id] = GridPosType::AtBoundary;
		for (auto* pln = g.plns; pln; pln = pln->next)
			++edge_num;
		edge_in_grid_range[g_id + 1] = edge_num;
	}
	edge_in_grid_list_len = edge_num;
	edge_in_grid_list = new size_t[edge_num];
	edge_num = 0;
	for (size_t g_id = 0; g_id < grid.num; ++g_id)
	{
		const size_t ge_start_id = edge_num;
		const auto& g = edge_grids[g_id];
		for (auto* pln = g.plns; pln; pln = pln->next)
		{
			edge_in_grid_list[edge_num] = edge_extractor.get_edge_id(*(pln->pln));
			++edge_num;
		}
		SimulationsOMP::swap_sort_acc<size_t>(
			edge_in_grid_list + ge_start_id, edge_num - ge_start_id);
	}
	edge_grid.clear(); // release memory
	edge_extractor.clear();

	for (size_t y_id = 0; y_id < grid.y_num; y_id++)
	{
		for (size_t x_id = 0; x_id < grid.x_num; x_id++)
		{
			if (grid_pos_type[x_id + y_id * grid.x_num] == GridPosType::AtBoundary)
				std::cout << "#";
			else
				std::cout << "@";
		}
		std::cout << "\n";
	}

	return 0;
}

void RigidMeshT2D::init_max_dist(double _dist) noexcept
{
	GridPosType* cur_g_pos;
	const size_t gn_x_num = grid.x_num + 1;
	const size_t gn_y_num = grid.y_num + 1;
	
	// reset previous max_dist result
	if (max_dist != DBL_MAX)
	{
		cur_g_pos = grid_pos_type;
		for (size_t y_id = 0; y_id < grid.y_num; ++y_id)
			for (size_t x_id = 0; x_id < grid.x_num; ++x_id)
			{
				// init grid type
				*reinterpret_cast<unsigned char *>(cur_g_pos) &= 0100;
				++cur_g_pos;
			}
	}

	// gns is the array of distance from point to rigid object
	double *gns = new double[gn_x_num * gn_y_num];
	double *cur_gn = gns;
	double g_h = (grid.hx < grid.hy ? grid.hx : grid.hy) * 0.5;
	Vector2D norm;
	Point2D gn_pt;
	gn_pt.y = grid.yl;
	for (size_t y_id = 0; y_id < gn_y_num; ++y_id)
	{
		gn_pt.x = grid.xl;
		for (size_t x_id = 0; x_id < gn_x_num; ++x_id)
		{
			detect_collision_with_point(gn_pt, g_h, *cur_gn, norm);
			*cur_gn -= g_h;
			++cur_gn;
			gn_pt.x += grid.hx;
		}
		gn_pt.y += grid.hy;
	}
	
	cur_g_pos = grid_pos_type;
	for (size_t y_id = 0; y_id < grid.y_num; ++y_id)
		for (size_t x_id = 0; x_id < grid.x_num; ++x_id)
		{
			cur_gn = gns + (y_id * gn_x_num + x_id);
			const double n1_dist = *cur_gn;
			const double n2_dist = *(cur_gn + 1);
			const double n3_dist = *(cur_gn + gn_x_num);
			const double n4_dist = *(cur_gn + gn_x_num + 1);
			if (n1_dist <= -_dist && n2_dist <= -_dist &&
				n3_dist <= -_dist && n4_dist <= -_dist)
				*cur_g_pos = GridPosType::FarOutside;
			else if (n1_dist >= _dist && n2_dist >= _dist &&
					 n3_dist >= _dist && n4_dist >= _dist)
				*cur_g_pos = GridPosType::FarInside;
			++cur_g_pos;
		}
	delete[] gns;

	// update max_dist and max_stride
	max_dist = _dist;
	max_stride = size_t(ceil(_dist / grid.hx));
	size_t max_stride_tmp = size_t(ceil(_dist / grid.hy));
	if (max_stride < max_stride_tmp)
		max_stride = max_stride_tmp;
}

bool RigidMeshT2D::detect_collision_with_point(
	const Point2D &pt,
	const double p_r,
	double &dist,
	Vector2D &norm
	) const noexcept
{
	if ((pt.x + p_r) < grid.xl || (pt.x - p_r) > grid.xu ||
		(pt.y + p_r) < grid.yl || (pt.y - p_r) > grid.yu)
		return false;

	const size_t p_x_id = grid.get_x_id(pt.x);
	const size_t p_y_id = grid.get_y_id(pt.y);
	// far inside or outside mesh
	if (unsigned char(grid_pos_type[offset_from_xy_id(p_x_id, p_y_id)]) > 2)
		return false;
	
	// cal dist
	IdRange id_range;
	id_range.xl_id = p_x_id > max_stride ? (p_x_id - max_stride) : 0;
	id_range.yl_id = p_y_id > max_stride ? (p_y_id - max_stride) : 0;
	size_t id_tmp;
	id_tmp = p_x_id + max_stride + 1;
	id_range.xu_id = id_tmp < grid.x_num ? id_tmp : grid.x_num;
	id_tmp = p_y_id + max_stride + 1;
	id_range.yu_id = id_tmp < grid.y_num ? id_tmp : grid.y_num;

	ClosestEdgeRes closest_edge;
	closest_edge.edge_id = SIZE_MAX;
	closest_edge.distance = max_dist;
	search_closest_edge(closest_edge, pt, id_range);

	// cal normal
	if (closest_edge.edge_id != SIZE_MAX)
	{
		pt_ln_dist[closest_edge.edge_id].cal_normal_to_point(
			pt, closest_edge.normal_type, norm);
		if (closest_edge.distance >= 0.0) // inside object
			norm.reverse();
		dist = closest_edge.distance + p_r;
		return dist >= 0.0;
	}
	return false;
}

void RigidMeshT2D::search_closest_edge(
	ClosestEdgeRes& closest_edge, const Point2D &pt,
	const IdRange &id_range) const noexcept
{
	size_t stride = id_range.xu_id - id_range.xl_id;
	size_t stride_tmp = id_range.yu_id - id_range.yl_id;
	if (stride < stride_tmp)
		stride = stride_tmp;
	if (stride > 1)
	{
		stride >>= 1;
		const size_t xm_id = id_range.xl_id + stride;
		const size_t ym_id = id_range.yl_id + stride;
		const double xl = grid.xl + double(id_range.xl_id) * grid.hx;
		const double yl = grid.yl + double(id_range.yl_id) * grid.hy;
		const double xm = grid.xl + double(xm_id) * grid.hx;
		const double ym = grid.yl + double(ym_id) * grid.hy;
		const double xu = grid.xl + double(id_range.xu_id) * grid.hx;
		const double yu = grid.yl + double(id_range.yu_id) * grid.hy;
		IdDist child_id_dists[4];
		IdRange child_id_ranges[4];
		Rect child_box;
		//
#define Set_Id_Range(cid, xld, xud, yld, yud) \
		IdRange &cir##cid = child_id_ranges[cid]; \
		cir##cid.xl_id = xld; \
		cir##cid.xu_id = xud; \
		cir##cid.yl_id = yld; \
		cir##cid.yu_id = yud;
		//
#define Set_Id_Dist(_cid, _xl, _xu, _yl, _yu) \
		child_box.xl = _xl; \
		child_box.xu = _xu; \
		child_box.yl = _yl; \
		child_box.yu = _yu; \
		IdDist &cid##_cid = child_id_dists[_cid]; \
		cid##_cid.id = _cid; \
		cid##_cid.dist = cal_rect_point_distance(child_box, pt)
		//
#define Search_Closest_Face(child_id) \
		if (cid##child_id.dist >= abs(closest_edge.distance)) \
			return; \
		search_closest_edge(closest_edge, pt, child_id_ranges[cid##child_id.id])
		//
		if (xm_id < id_range.xu_id)
		{
			if (ym_id < id_range.yu_id)
			{
				// 0
				Set_Id_Range(0, id_range.xl_id, xm_id, id_range.yl_id, ym_id);
				Set_Id_Dist(0, xl, xm, yl, ym);
				// 1
				Set_Id_Range(1, xm_id, id_range.xu_id, id_range.yl_id, ym_id);
				Set_Id_Dist(1, xm, xu, yl, ym);
				// 2
				Set_Id_Range(2, id_range.xl_id, xm_id, ym_id, id_range.yu_id);
				Set_Id_Dist(2, xl, xm, ym, yu);
				// 3
				Set_Id_Range(3, xm_id, id_range.xu_id, ym_id, id_range.yu_id);
				Set_Id_Dist(3, xm, xu, ym, yu);
				//
				sort_acc_id_dist_pairs_4(child_id_dists);
				Search_Closest_Face(0);
				Search_Closest_Face(1);
				Search_Closest_Face(2);
				Search_Closest_Face(3);
			}
			else // ym_id >= yu_id
			{
				// 0
				Set_Id_Range(0,id_range.xl_id, xm_id, id_range.yl_id, id_range.yu_id);
				Set_Id_Dist(0, xl, xm, yl, yu);
				// 1
				Set_Id_Range(1, xm_id, id_range.xu_id, id_range.yl_id, id_range.yu_id);
				Set_Id_Dist(1, xm, xu, yl, yu);
				//
				sort_acc_id_dist_pairs_2(child_id_dists);
				Search_Closest_Face(0);
				Search_Closest_Face(1);
			}
		}
		else // xm_id >= xu_id
		{
			if (ym_id < id_range.yu_id)
			{
				// 0
				Set_Id_Range(0, id_range.xl_id, id_range.xu_id, id_range.yl_id, ym_id);
				Set_Id_Dist(0, xl, xu, yl, ym);
				// 1
				Set_Id_Range(1, id_range.xl_id, id_range.xu_id, ym_id, id_range.yu_id);
				Set_Id_Dist(1, xl, xu, ym, yu);
				//
				sort_acc_id_dist_pairs_2(child_id_dists);
				Search_Closest_Face(0);
				Search_Closest_Face(1);
			}
			else // ym_id >= yu_id && xm_id >= xu_id
			{
				Set_Id_Range(0, id_range.xl_id, id_range.xu_id, id_range.yl_id, id_range.yu_id);
				Set_Id_Dist(0, xl, xu, yl, yu);
				Search_Closest_Face(0);
			}
		}
		return;
	}

	double dist_tmp;
	unsigned char norm_type;
	const size_t g_off = offset_from_xy_id(id_range.xl_id, id_range.yl_id);
	const size_t end_e_id = edge_in_grid_range[g_off + 1];
	for (size_t e_id = edge_in_grid_range[g_off]; e_id < end_e_id; ++e_id)
	{
		const PointToLineDistance& pld = pt_ln_dist[edge_in_grid_list[e_id]];
		norm_type = pld.cal_distance_to_point(pt, dist_tmp);
		if (abs(dist_tmp) < abs(closest_edge.distance))
		{
			closest_edge.edge_id = edge_in_grid_list[e_id];
			closest_edge.distance = dist_tmp;
			closest_edge.normal_type = norm_type;
		}
	}
}
