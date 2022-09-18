#ifndef __Search_Grid_2D_Line_hpp__
#define __Search_Grid_2D_Line_hpp__

#include "ItemBuffer.hpp"
#include "Geometry2D.h"
#include "Grid2D.hpp"
#include "DetectCollisionSAT.hpp"

template <typename Line>
class SearchGrid2DLine
{
public:
	struct LinePointer
	{
		Line *pln;
		LinePointer* next;
	};

	struct Grid { LinePointer* plns; };

	inline SearchGrid2DLine() {}
	~SearchGrid2DLine() { clear(); }

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
		ln_pt_buffer.clear();
	}

	inline int alloc_grid(double xl, double yl,
		double hx, double hy, size_t x_num, size_t y_num)
	{
		if (grid.alloc_grid(xl, yl, hx, hy, x_num, y_num) < 0)
			return -1;
		for (size_t g_id = 0; g_id < grid.num; ++g_id)
			grid.grids[g_id].plns = nullptr;
		return 0;
	}

	inline int alloc_grid(double xl, double yl,
		double xu, double yu, double hx, double hy)
	{
		if (grid.alloc_grid(xl, yl, xu, yu, hx, hy) < 0)
			return -1;
		for (size_t g_id = 0; g_id < grid.num; ++g_id)
			grid.grids[g_id].plns = nullptr;
		return 0;
	}

	template <typename Grid2>
	inline int alloc_grid(const Grid2 &other)
	{
		if (grid.alloc_grid(other) < 0)
			return -1;
		for (size_t g_id = 0; g_id < grid.num; ++g_id)
			grid.grids[g_id].plns = nullptr;
		return 0;
	}
	
	template <typename Node>
	void apply_lines(const Node* nodes, const Line *lns, size_t ln_num) noexcept
	{
		ln_pt_buffer.set_page_size(ln_num * 4);
		for (size_t l_id = 0; l_id < ln_num; ++l_id)
			add_line_to_grids<Node>(lns[l_id], nodes);
	}

protected:
	Grid2D<Grid> grid;
	MemoryUtils::ItemBuffer<LinePointer, 2> ln_pt_buffer;

	inline void add_line_to_grid(Grid& g, const Line &ln)
	{
		LinePointer* ln_pt = ln_pt_buffer.alloc();
		ln_pt->pln = const_cast<Line *>(&ln);
		ln_pt->next = g.plns;
		g.plns = ln_pt;
	}

	template <typename Node>
	void add_line_to_grids(const Line &ln, const Node *nodes) noexcept
	{
		const Node &n1 = nodes[ln.n1];
		const Node &n2 = nodes[ln.n2];
				
		DetectLineAABBCollisionSAT ln_aabb_collision;
		ln_aabb_collision.init_line(n1, n2);

		Rect ln_bbox;
		ln_aabb_collision.get_ln_bbox(ln_bbox);
		const size_t xl_id = grid.get_x_id(ln_bbox.xl);
		const size_t xu_id = grid.get_x_id(ln_bbox.xu) + 1;
		const size_t yl_id = grid.get_y_id(ln_bbox.yl);
		const size_t yu_id = grid.get_y_id(ln_bbox.yu) + 1;
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
				if (ln_aabb_collision.detect(grid_box))
				{
					Grid& g = grid.grid_by_xy_id(x_id, y_id);
					add_line_to_grid(g, ln);
				}
				grid_box.xl = grid_box.xu;
			}
			grid_box.yl = grid_box.yu;
		}
	}
};

#endif