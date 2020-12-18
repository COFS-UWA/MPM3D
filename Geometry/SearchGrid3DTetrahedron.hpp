#ifndef __Search_Grid_3D_Tetrahedron_hpp__
#define __Search_Grid_3D_Tetrahedron_hpp__

#include "ItemBuffer.hpp"
#include "Geometry3D.h"
#include "TetrahedronUtils.h"
#include "DetectCollisionSAT.hpp"
#include "Grid3D.hpp"

// Teh has n1, n2, n3, n4
template <typename Teh>
class SearchGrid3DTetrahedron
{
public:
	struct TehPointer
	{
		Teh *pteh;
		TehPointer* next;
	};

	struct Grid { TehPointer *ptehs; };
	
	inline SearchGrid3DTetrahedron() {}
	~SearchGrid3DTetrahedron() { clear(); }

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
		teh_pt_buffer.clear();
	}
	
	inline int alloc_grid(double xl, double yl, double zl,
		double hx, double hy, double hz,
		size_t x_num, size_t y_num, size_t z_num)
	{ return grid.alloc_grid(xl, yl, zl, hx, hy, hz, x_num, y_num, z_num); }

	inline int alloc_grid(double xl, double yl, double zl,
		double xu, double yu, double zu,
		double hx, double hy, double hz)
	{ return grid.alloc_grid(xl, yl, zl, xu, yu, zu, hx, hy, hz); }

	template <typename Grid2>
	inline int alloc_grid(const Grid2 &other)
	{ return grid.alloc_grid(other); }

	template <typename Node>
	void apply_tetrahedrons(const Node *nodes, const Teh *tehs, size_t teh_num)
	{
		size_t buf_size = teh_num / 16;
		if (buf_size < 16) buf_size = 16;
		teh_pt_buffer.set_page_size(buf_size);
		for (size_t t_id = 0; t_id < teh_num; ++t_id)
			add_teh_to_grids(tehs[t_id], nodes);
	}
	
protected:
	Grid3D<Grid> grid;
	MemoryUtils::ItemBuffer<TehPointer, 2> teh_pt_buffer;

	// 22 seperating axises:
	// 4 face normal + 3 * 6 edge cross products
	Vector3D seperating_axes[22];

	inline void add_teh_to_grid(Grid &g, const Teh &teh)
	{
		TehPointer *teh_pt = teh_pt_buffer.alloc();
		teh_pt->pteh = const_cast<Teh *>(&teh);
		teh_pt->next = g.ptehs;
		g.ptehs = teh_pt;
	}
	
	template <typename Node>
	void add_teh_to_grids(const Teh& teh, const Node* nodes) noexcept
	{
		const Node& n1 = nodes[teh.n1];
		const Node& n2 = nodes[teh.n2];
		const Node& n3 = nodes[teh.n3];
		const Node& n4 = nodes[teh.n4];
		
		DetectTetrahedronAABBCollisionSAT teh_aabb_collision;
		teh_aabb_collision.init_tetrahedron(n1, n2, n3, n4);

		Cube teh_bbox;
		teh_aabb_collision.get_teh_bbox(teh_bbox);
		const size_t xl_id = size_t(floor((teh_bbox.xl - grid.xl) / grid.hx));
		const size_t xu_id = size_t(ceil((teh_bbox.xu - grid.xl) / grid.hx));
		const size_t yl_id = size_t(floor((teh_bbox.yl - grid.yl) / grid.hy));
		const size_t yu_id = size_t(ceil((teh_bbox.yu - grid.yl) / grid.hy));
		const size_t zl_id = size_t(floor((teh_bbox.zl - grid.zl) / grid.hz));
		const size_t zu_id = size_t(ceil((teh_bbox.zu - grid.zl) / grid.hz));
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
					if (teh_aabb_collision.detect(grid_box))
					{
						Grid &g = grid.grid_by_xyz_id(x_id, y_id, z_id);
						add_teh_to_grid(g, teh);
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