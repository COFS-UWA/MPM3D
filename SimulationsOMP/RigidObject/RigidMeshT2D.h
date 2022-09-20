#ifndef __Rigid_Mesh_T2D_h__
#define __Rigid_Mesh_T2D_h__

#include "Geometry2D.h"
#include "Grid2D.hpp"
#include "TriangleUtils.h"
#include "TriangleMesh.h"

class RigidMeshT2D
{
public:
	enum class GridPosType : unsigned char
	{
		AtBoundary = 0, // 000
		Inside = 1, // 001
		Outside = 2, // 010
		FarInside = 5, // 101
		FarOutside = 6 // 110
	};

protected:
	Grid2D<> grid;
	GridPosType* grid_pos_type;
	size_t *edge_in_grid_range;
	size_t edge_in_grid_list_len;
	size_t *edge_in_grid_list;
	size_t edge_num;
	PointToLineDistance*pt_ln_dist;

	inline size_t offset_from_xy_id(size_t x_id, size_t y_id) const noexcept
	{ return grid.offset_from_xy_id(x_id, y_id); }

	double max_dist;
	size_t max_stride;
	inline size_t get_max_stride() const noexcept
	{
		return grid.x_num > grid.y_num ? grid.x_num : grid.y_num;
	}

	struct IdRange
	{
		size_t xl_id, xu_id;
		size_t yl_id, yu_id;
	};
	struct IdDist
	{
		size_t id;
		double dist;
	};
	inline static void sort_acc_id_dist_pairs_4(IdDist id_dists[4]) noexcept
	{
		for (size_t i = 0; i < 3; ++i)
		{
			size_t min_id = i;
			for (size_t j = i + 1; j < 4; ++j)
				if (id_dists[min_id].dist > id_dists[j].dist)
					min_id = j;
			if (min_id != i)
			{
				IdDist& idi = id_dists[i];
				const size_t id_tmp = idi.id;
				const double dist_tmp = idi.dist;
				IdDist& idmin = id_dists[min_id];
				idi.id = idmin.id;
				idi.dist = idmin.dist;
				idmin.id = id_tmp;
				idmin.dist = dist_tmp;
			}
		}
	}
	inline static void sort_acc_id_dist_pairs_2(IdDist id_dists[2]) noexcept
	{
		if (id_dists[0].dist > id_dists[1].dist)
		{
			IdDist& id0 = id_dists[0];
			const size_t id_tmp = id0.id;
			const double dist_tmp = id0.dist;
			IdDist& id1 = id_dists[1];
			id0.id = id1.id;
			id0.dist = id1.dist;
			id1.id = id_tmp;
			id1.dist = dist_tmp;
		}
	}

	struct ClosestEdgeRes
	{
		size_t edge_id;
		double distance;
		unsigned char normal_type;
	};

	void search_closest_edge(ClosestEdgeRes& cloeset_face,
		const Point2D& pt, const IdRange& id_range) const noexcept;
	
public:
	RigidMeshT2D();
	~RigidMeshT2D();

	inline void clear() noexcept
	{
		grid.clear();
		edge_in_grid_list_len = 0;
		edge_num = 0;
		if (grid_pos_type)
		{
			delete[] grid_pos_type;
			grid_pos_type = nullptr;
		}
		if (edge_in_grid_range)
		{
			delete[] edge_in_grid_range;
			edge_in_grid_range = nullptr;
		}
		if (edge_in_grid_list)
		{
			delete[] edge_in_grid_list;
			edge_in_grid_list = nullptr;
		}
		if (pt_ln_dist)
		{
			delete[] pt_ln_dist;
			pt_ln_dist = nullptr;
		}
	}

	inline const double get_grid_xl() const noexcept { return grid.xl; }
	inline const double get_grid_yl() const noexcept { return grid.yl; }
	inline const double get_grid_hx() const noexcept { return grid.hx; }
	inline const double get_grid_hy() const noexcept { return grid.hy; }
	inline const double get_grid_xu() const noexcept { return grid.xu; }
	inline const double get_grid_yu() const noexcept { return grid.yu; }
	inline const size_t get_grid_x_num() const noexcept { return grid.x_num; }
	inline const size_t get_grid_y_num() const noexcept { return grid.y_num; }
	inline const size_t get_grid_num() const noexcept { return grid.num; }
	inline void get_bbox(Rect& bbox) const noexcept
	{
		bbox.xl = grid.xl; bbox.yl = grid.yl;
		bbox.xu = grid.xu; bbox.yu = grid.yu;
	}

	// grid_num
	inline const GridPosType* get_grid_pos_type() const noexcept { return grid_pos_type; }
	// grid_num + 1
	inline const size_t* get_edge_in_grid_range() const noexcept { return edge_in_grid_range; }
	inline size_t get_edge_in_grid_list_len() const noexcept { return edge_in_grid_list_len; }
	inline const size_t* get_edge_in_grid_list() const noexcept { return edge_in_grid_list; }
	inline size_t get_edge_num() const noexcept { return edge_num; }
	inline const PointToLineDistance *get_pt_ln_dist() const noexcept { return pt_ln_dist; }

	int init_from_mesh(TriangleMesh &tmesh, double ghx, double ghy);
	// Speed up searching
	// must be called after init_from_mesh()
	void init_max_dist(double _dist) noexcept;
	
	bool detect_collision_with_point(const Point2D& pt,
		const double p_r, double& dist, Vector2D& norm) const noexcept;
};

#endif