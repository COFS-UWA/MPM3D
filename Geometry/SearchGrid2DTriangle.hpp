#ifndef __Search_Grid_2D_Triangle_hpp__
#define __Search_Grid_2D_Triangle_hpp__

#include "ItemBuffer.hpp"
#include "Geometry2D.h"
#include "Grid2D.hpp"
#include "TriangleUtils.h"
#include "DetectCollisionSAT.hpp"

// Tri has n1, n2, n3, n4
template <typename Tri>
class SearchGrid2DTriangle
{
public:
	struct TriPointer
	{
		Tri *ptri;
		TriPointer* next;
	};

	struct Grid { TriPointer *ptris; };
	
	inline SearchGrid2DTriangle() {}
	~SearchGrid2DTriangle() { clear(); }

	inline double get_xl() const noexcept { return grid.xl; }
	inline double get_yl() const noexcept { return grid.yl; }
	inline double get_xu() const noexcept { return grid.xu; }
	inline double get_yu() const noexcept { return grid.yu; }
	inline double get_hx() const noexcept { return grid.hx; }
	inline double get_hy() const noexcept { return grid.hy; }
	inline size_t get_x_num() const noexcept { return grid.x_num; }
	inline size_t get_y_num() const noexcept { return grid.y_num; }
	inline size_t get_num() const noexcept { return grid.num; }
	inline const Grid* get_grids() const noexcept { return grid.grids; }
	inline Grid* get_grids() noexcept { return grid.grids; }
	
	inline void clear()
	{
		grid.clear();
		tri_pt_buffer.clear();
	}
	
	inline int alloc_grid(double xl, double yl,
		double hx, double hy, size_t x_num, size_t y_num)
	{
		if (grid.alloc_grid(xl, yl, hx, hy, x_num, y_num) < 0)
			return -1;
		for (size_t g_id = 0; g_id < grid.num; ++g_id)
			grid.grids[g_id].ptris = nullptr;
		return 0;
	}

	inline int alloc_grid(
		double xl, double yl,
		double xu, double yu,
		double hx, double hy)
	{
		if (grid.alloc_grid(xl, yl, xu, yu, hx, hy) < 0)
			return -1;
		for (size_t g_id = 0; g_id < grid.num; ++g_id)
			grid.grids[g_id].ptris = nullptr;
		return 0;
	}

	template <typename Grid2>
	inline int alloc_grid(const Grid2 &other)
	{
		if (grid.alloc_grid(other) < 0)
			return -1;
		for (size_t g_id = 0; g_id < grid.num; ++g_id)
			grid.grids[g_id].ptris = nullptr;
		return 0;
	}

	template <typename Node>
	void apply_triangles(const Node *nodes, const Tri *tris, size_t tri_num)
	{
		size_t buf_size = tri_num * 8;
		if (buf_size < 16) buf_size = 16;
		tri_pt_buffer.set_page_size(buf_size);
		for (size_t t_id = 0; t_id < tri_num; ++t_id)
			add_tri_to_grids<Node>(tris[t_id], nodes);
	}
	
protected:
	Grid2D<Grid> grid;
	MemoryUtils::ItemBuffer<TriPointer, 2> tri_pt_buffer;

	inline void add_tri_to_grid(Grid &g, const Tri &tri)
	{
		TriPointer *tri_pt = tri_pt_buffer.alloc();
		tri_pt->ptri = const_cast<Tri *>(&tri);
		tri_pt->next = g.ptris;
		g.ptris = tri_pt;
	}
	
	template <typename Node>
	void add_tri_to_grids(const Tri& tri, const Node* nodes) noexcept
	{
		const Node& n1 = nodes[tri.n1];
		const Node& n2 = nodes[tri.n2];
		const Node& n3 = nodes[tri.n3];
		
		DetectTriangleAABBCollisionSAT tri_aabb_collision;
		tri_aabb_collision.init_triangle(n1, n2, n3);

		Rect tri_bbox;
		tri_aabb_collision.get_tri_bbox(tri_bbox);
		const size_t xl_id = grid.get_x_id(tri_bbox.xl);
		const size_t xu_id = grid.get_x_id(tri_bbox.xu) + 1;
		const size_t yl_id = grid.get_y_id(tri_bbox.yl);
		const size_t yu_id = grid.get_y_id(tri_bbox.yu) + 1; 
		const double grid_box_x_min = grid.xl + double(xl_id) * grid.hx;
		Rect grid_box;
		grid_box.yl = grid.yl + double(yl_id) * grid.hy;
		for (size_t y_id = yl_id; y_id < yu_id; ++y_id)
		{
			grid_box.yu = grid_box.yl + grid.hy;
			grid_box.xl = grid_box_x_min;
			for (size_t x_id = xl_id; x_id < xu_id; ++x_id)
			{
				grid_box.xu = grid_box.xl + grid.hx;
				if (tri_aabb_collision.detect(grid_box))
				{
					Grid &g = grid.grid_by_xy_id(x_id, y_id);
					add_tri_to_grid(g, tri);
				}
				grid_box.xl = grid_box.xu;
			}
			grid_box.yl = grid_box.yu;
		}
	}
};

#endif