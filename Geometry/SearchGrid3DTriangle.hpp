#ifndef __Search_Grid_3D_Triangle_hpp__
#define __Search_Grid_3D_Triangle_hpp__

#include "ItemBuffer.hpp"
#include "Geometry3D.h"
#include "DetectCollisionSAT.hpp"
#include "Grid3D.hpp"

template <typename Tri>
class SearchGrid3DTriangle
{
public:
	struct TriPointer
	{
		Tri *ptri;
		TriPointer* next;
	};

	struct Grid { TriPointer* ptris; };

	inline SearchGrid3DTriangle() {}
	~SearchGrid3DTriangle() { clear(); }

	inline double get_xl() const noexcept { return grid.xl; }
	inline double get_yl() const noexcept { return grid.yl; }
	inline double get_zl() const noexcept { return grid.zl; }
	inline double get_xu() const noexcept { return grid.xu; }
	inline double get_yu() const noexcept { return grid.yu; }
	inline double get_zu() const noexcept { return grid.zu; }
	inline double get_hx() const noexcept { return grid.hx; }
	inline double get_hy() const noexcept { return grid.hy; }
	inline double get_hz() const noexcept { return grid.hz; }
	inline size_t get_x_num() const noexcept { return grid.x_num; }
	inline size_t get_y_num() const noexcept { return grid.y_num; }
	inline size_t get_z_num() const noexcept { return grid.z_num; }
	inline size_t get_xy_num() const noexcept { return grid.xy_num; }
	inline size_t get_num() const noexcept { return grid.num; }
	inline const Grid* get_grids() const noexcept { return grid.grids; }
	inline Grid* get_grids() noexcept { return grid.grids; }

	inline void clear()
	{
		grid.clear();
		tri_pt_buffer.clear();
	}

	inline int alloc_grid(
		double xl, double yl, double zl,
		double hx, double hy, double hz,
		size_t x_num, size_t y_num, size_t z_num)
	{ return grid.alloc_grid(xl, yl, zl, hx, hy, hz, x_num, y_num, z_num); }

	inline int alloc_grid(
		double xl, double yl, double zl,
		double xu, double yu, double zu,
		double hx, double hy, double hz)
	{ return grid.alloc_grid(xl, yl, zl, xu, yu, zu, hx, hy, hz); }

	template <typename Grid2>
	inline int alloc_grid(const Grid2 &other)
	{ return grid.alloc_grid(other); }
	
	template <typename Node>
	void apply_triangles(
		const Node* nodes,
		const Tri *tris,
		size_t tri_num
		) noexcept
	{
		tri_pt_buffer.set_page_size(tri_num);
		for (size_t t_id = 0; t_id < tri_num; ++t_id)
			add_triangle_to_grids<Node>(tris[t_id], nodes);
	}

protected:
	Grid3D<Grid> grid;
	MemoryUtils::ItemBuffer<TriPointer, 2> tri_pt_buffer;

	// seperating axises:
	// 4 tri normal + 3 * 6 edge cross products
	Vector3D seperating_axes[22];

	inline void add_triangle_to_grid(Grid& g, const Tri &tri)
	{
		TriPointer* tri_pt = tri_pt_buffer.alloc();
		tri_pt->ptri = const_cast<Tri *>(&tri);
		tri_pt->next = g.ptris;
		g.ptris = tri_pt;
	}

	template <typename Node>
	void add_triangle_to_grids(
		const Tri &tri,
		const Node *nodes
		) noexcept
	{
		const Node &n1 = nodes[tri.n1];
		const Node &n2 = nodes[tri.n2];
		const Node &n3 = nodes[tri.n3];
		
		Detect3DTriangleAABBCollisionSAT tri_aabb_collision;
		tri_aabb_collision.init_triangle(n1, n2, n3);

		Cube tri_bbox;
		tri_aabb_collision.get_tri_bbox(tri_bbox);
		const size_t xl_id = size_t(floor((tri_bbox.xl - grid.xl) / grid.hx));
		const size_t xu_id = size_t(ceil((tri_bbox.xu - grid.xl) / grid.hx));
		const size_t yl_id = size_t(floor((tri_bbox.yl - grid.yl) / grid.hy));
		const size_t yu_id = size_t(ceil((tri_bbox.yu - grid.yl) / grid.hy));
		const size_t zl_id = size_t(floor((tri_bbox.zl - grid.zl) / grid.hz));
		const size_t zu_id = size_t(ceil((tri_bbox.zu - grid.zl) / grid.hz));
		const double grid_box_x_min = grid.xl + double(xl_id) * grid.hx;
		const double grid_box_y_min = grid.yl + double(yl_id) * grid.hy;
		Cube grid_box;
		grid_box.zl = grid.zl + double(zl_id) * grid.hz;
		for (size_t z_id = zl_id; z_id < zu_id; ++z_id)
		{
			grid_box.zu = grid_box.zl + grid.hz;
			grid_box.yl = grid_box_y_min;
			for (size_t y_id = yl_id; y_id < yu_id; ++y_id)
			{
				grid_box.yu = grid_box.yl + grid.hy;
				grid_box.xl = grid_box_x_min;
				for (size_t x_id = xl_id; x_id < xu_id; ++x_id)
				{
					grid_box.xu = grid_box.xl + grid.hx;
					if (tri_aabb_collision.detect(grid_box))
					{
						Grid& g = grid.grid_by_xyz_id(x_id, y_id, z_id);
						add_triangle_to_grid(g, tri);
					}
					grid_box.xl = grid_box.xu;
				}
				grid_box.yl = grid_box.yu;
			}
			grid_box.zl = grid_box.zu;
		}
	}
};

#endif