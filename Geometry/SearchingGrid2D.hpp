#ifndef __Searching_Grid_2D_hpp__
#define __Searching_Grid_2D_hpp__

#include "ItemBuffer.hpp"
#include "Geometry.h"

// Accelerate spatial searching of Triangle mesh
// Assumptions:
//	1. Triangle mesh has is_in_triangle(MeshElement &e, double x, double y)
//	2. Triangle mesh has get_elems() and get_elem_num()
//  3. Triangle mesh has get_bounding_box()
template <typename TriangleMesh>
class SearchingGrid2D
{
public:
	typedef typename TriangleMesh::Element MeshElement;
	typedef typename TriangleMesh::Node MeshNode;
	
	struct ElemPointer
	{
		MeshElement *e;
		ElemPointer *next;
	};

	struct Grid
	{
		size_t x_id, y_id;
		ElemPointer *pelems;
	};

protected:
	double hx, hy;
	double x_min, x_max;
	double y_min, y_max;
	size_t x_num, y_num, num;
	Grid *grids;

	TriangleMesh *mesh;
	size_t node_num;
	MeshNode *nodes;
	size_t elem_num;
	MeshElement *elems;

	MemoryUtils::ItemBuffer<ElemPointer> pe_buffer;

public:
	SearchingGrid2D() :
		x_num(0), y_num(0), num(0),
		grids(nullptr), mesh(nullptr) {}
	~SearchingGrid2D() { clear(); }

	inline double get_hx() { return hx; }
	inline double get_hy() { return hy; }
	inline double get_x_min() { return x_min; }
	inline double get_x_max() { return x_max; }
	inline double get_y_min() { return y_min; }
	inline double get_y_max() { return y_max; }
	inline size_t get_grid_num() { return num; }
	inline Grid *get_grids() { return grids; }

	int init(TriangleMesh &_mesh, double hx, double hy)
	{
		Rect &mesh_bbox = _mesh.get_bounding_box();
		if (alloc_grids(mesh_bbox, hx, hy) < 0)
			return -1;

		set_mesh_info(_mesh);
		
		for (size_t e_id = 0; e_id < elem_num; ++e_id)
			add_elem_to_grids(elems[e_id]);

		return 0;
	}

	void clear()
	{
		if (grids)
		{
			delete[] grids;
			grids = nullptr;
		}
		x_num = 0;
		y_num = 0;
		num = 0;
		mesh = nullptr;
		node_num = 0;
		nodes = nullptr;
		elem_num = 0;
		elems = nullptr;
		pe_buffer.clear();
	}

	template <typename Point2D>
	inline MeshElement *find_in_which_element(Point2D &point)
	{
		if (point.x < x_min || point.x > x_max ||
			point.y < y_min || point.y > y_max)
			return nullptr;

		size_t x_id = size_t((point.x - x_min) / hx);
		size_t y_id = size_t((point.y - y_min) / hy);
		Grid &g = get_grid(x_id, y_id);
		MeshElement *elem;
		for (ElemPointer *pelem = g.pelems; pelem; pelem = pelem->next)
		{
			elem = pelem->e;
			if (mesh->is_in_triangle(*elem, point))
				return elem;
		}
		return nullptr;
	}

	int alloc_grids(Rect box, double _hx, double _hy)
	{
		clear();
		double x_len = box.xu - box.xl;
		double y_len = box.yu - box.yl;
		if (x_len <= 0.0 || y_len <= 0.0)
			return -1;

		hx = _hx;
		hy = _hy;
		x_num = size_t(ceil(x_len / hx));
		y_num = size_t(ceil(y_len / hy));
		double x_pad = (double(x_num) * hx - x_len) * 0.5;
		double y_pad = (double(y_num) * hy - y_len) * 0.5;
		x_min = box.xl - x_pad;
		x_max = box.xu + x_pad;
		y_min = box.yl - y_pad;
		y_max = box.yu + y_pad;
		num = x_num * y_num;
		grids = new Grid[num];
		Grid *pg = grids;
		for (size_t y_id = 0; y_id < y_num; ++y_id)
			for (size_t x_id = 0; x_id < x_num; ++x_id)
			{
				pg->x_id = x_id;
				pg->y_id = y_id;
				pg->pelems = nullptr;
				++pg;
			}
		return 0;
	}

	inline void set_mesh_info(TriangleMesh &_mesh)
	{
		mesh = &_mesh;
		node_num = mesh->get_node_num();
		nodes = mesh->get_nodes();
		elem_num = mesh->get_elem_num();
		elems = mesh->get_elems();
	}

	void add_elem_to_grids(MeshElement &elem)
	{
		Rect elem_bbox;
		get_elem_bbox(elem, elem_bbox);

		size_t xl_id = size_t(floor((elem_bbox.xl - x_min) / hx));
		size_t xu_id = size_t(ceil( (elem_bbox.xu - x_min) / hx));
		size_t yl_id = size_t(floor((elem_bbox.yl - y_min) / hy));
		size_t yu_id = size_t(ceil( (elem_bbox.yu - y_min) / hy));
		Rect grid_box;
		for (size_t y_id = yl_id; y_id < yu_id; ++y_id)
		{
			grid_box.yl = y_min + double(y_id) * hy;
			grid_box.yu = grid_box.yl + hy;
			for (size_t x_id = xl_id; x_id < xu_id; ++x_id)
			{
				grid_box.xl = x_min + double(x_id) * hx;
				grid_box.xu = grid_box.xl + hx;
				if (detect_AABB_triangle_collision(grid_box, elem))
				{
					Grid &cur_grid = get_grid(x_id, y_id);
					add_elem_to_grid(cur_grid, elem);
				}
			}
		}
	}

	inline Grid &get_grid(size_t x_id, size_t y_id)
	{
		return grids[x_num * y_id + x_id];
	}

protected:
	inline void add_elem_to_grid(Grid &g, MeshElement &e)
	{
		ElemPointer *pe = pe_buffer.alloc();
		pe->e = &e;
		pe->next = g.pelems;
		g.pelems = pe;
	}

	inline void clear_elems_in_grid(Grid &g)
	{
		ElemPointer *pe = g.pelems;
		ElemPointer *pe_tmp;
		while (pe)
		{
			pe_tmp = pe;
			pe = pe->next;
			pe_buffer.del(pe_tmp);
		}
	}

	inline void get_elem_bbox(MeshElement& elem, Rect &elem_bbox)
	{
		MeshNode &n1 = nodes[elem.n1];
		MeshNode &n2 = nodes[elem.n2];
		MeshNode &n3 = nodes[elem.n3];

		double xl, xu, yl, yu;
		xl = n1.x;
		xu = xl;
		if (xl > n2.x)
			xl = n2.x;
		if (xu < n2.x)
			xu = n2.x;
		if (xl > n3.x)
			xl = n3.x;
		if (xu < n3.x)
			xu = n3.x;

		yl = n1.y;
		yu = yl;
		if (yl > n2.y)
			yl = n2.y;
		if (yu < n2.y)
			yu = n2.y;
		if (yl > n3.y)
			yl = n3.y;
		if (yu < n3.y)
			yu = n3.y;

		elem_bbox.xl = xl;
		elem_bbox.xu = xu;
		elem_bbox.yl = yl;
		elem_bbox.yu = yu;
	}

	inline void swap(double &a, double &b) { double c = a; a = b; b = c; }
public:
	// test if aligned-axis bounding box intersects triangle
	bool detect_AABB_triangle_collision(Rect &aabb, MeshElement &e)
	{
		MeshNode &n1 = nodes[e.n1];
		MeshNode &n2 = nodes[e.n2];
		MeshNode &n3 = nodes[e.n3];

		// whether the triangle nodes locate in box
		// efficient when triangle is much smaller than grid
		if (aabb.is_in_box(n1) ||
			aabb.is_in_box(n2) || 
			aabb.is_in_box(n3))
			return true;

		// whether box corners locate in triangle
		// efficient when grid is much smaller than triangle
		if (mesh->is_in_triangle(e, aabb.xl, aabb.yl) ||
			mesh->is_in_triangle(e, aabb.xl, aabb.yu) ||
			mesh->is_in_triangle(e, aabb.xu, aabb.yl) ||
			mesh->is_in_triangle(e, aabb.xu, aabb.yu))
			return true;

		// take grid centre as origin
		double aabb_xc = (aabb.xl + aabb.xu) * 0.5;
		double aabb_yc = (aabb.yl + aabb.yu) * 0.5;
		double n_x0 = n1.x - aabb_xc;
		double n_y0 = n1.y - aabb_yc;
		double n_x1 = n2.x - aabb_xc;
		double n_y1 = n2.y - aabb_yc;
		double n_x2 = n3.x - aabb_xc;
		double n_y2 = n3.y - aabb_yc;

		double r, v_min, v_max;
		// a31
		r = (hx * abs(n_y1 - n_y0) + hy * abs(n_x1 - n_x0)) * 0.5;
		v_min = n_x1 * n_y0 - n_x0 * n_y1;
		v_max = (n_y0 - n_y1) * n_x2 + (n_x1 - n_x0) * n_y2;
		if (v_min > v_max)
			swap(v_min, v_max);
		if (v_min > r || v_max < -r)
			return false;
		// a32
		r = (hx * abs(n_y2 - n_y1) + hy * abs(n_x2 - n_x1)) * 0.5;
		v_min = n_x2 * n_y1 - n_x1 * n_y2;
		v_max = (n_y1 - n_y2) * n_x0 + (n_x2 - n_x1) * n_y0;
		if (v_min > v_max)
			swap(v_min, v_max);
		if (v_min > r || v_max < -r)
			return false;
		// a33
		r = (hx * abs(n_y0 - n_y2) + hy * abs(n_x0 - n_x2)) * 0.5;
		v_min = n_x0 * n_y2 - n_x2 * n_y0;
		v_max = (n_y2 - n_y0) * n_x1 + (n_x0 - n_x2) * n_y1;
		if (v_min > v_max)
			swap(v_min, v_max);
		if (v_min > r || v_max < -r)
			return false;
		// no seperating axis
		return true;
	}
};

#endif