#ifndef __Rigid_Mesh_T3D_h__
#define __Rigid_Mesh_T3D_h__

#include "Geometry3D.h"
#include "TetrahedronUtils.h"
#include "Grid3D.hpp"
#include "TetrahedronMesh.h"

class RigidMeshT3D
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
	Grid3D<> grid;
	GridPosType* grid_pos_type;
	size_t *face_in_grid_range;
	size_t face_in_grid_list_len;
	size_t *face_in_grid_list;
	size_t face_num;
	PointToTriangleDistance *pt_tri_dist;

	inline size_t offset_from_xyz_id(size_t x_id, size_t y_id, size_t z_id) const noexcept
	{ return grid.offset_from_xyz_id(x_id, y_id, z_id); }

	double max_dist;
	size_t max_stride;
	inline size_t get_max_stride() const noexcept
	{
		const size_t tmp = grid.x_num > grid.y_num ? grid.x_num : grid.y_num;
		return tmp > grid.z_num ? tmp : grid.z_num;
	}

	struct IdRange
	{
		size_t xl_id, xu_id;
		size_t yl_id, yu_id;
		size_t zl_id, zu_id;
	};
	struct IdDist
	{
		size_t id;
		double dist;
	};
	inline static void sort_acc_id_dist_pairs_8(IdDist id_dists[8]) noexcept
	{
		for (size_t i = 0; i < 7; ++i)
		{
			size_t min_id = i;
			for (size_t j = i + 1; j < 8; ++j)
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

	struct ClosestFaceRes
	{
		size_t face_id;
		double distance;
		unsigned char normal_type;
	};

	void search_closest_face(ClosestFaceRes& cloeset_face,
		const Point3D& pt, const IdRange& id_range) const noexcept;
	
public:
	RigidMeshT3D();
	~RigidMeshT3D();

	inline void clear() noexcept
	{
		grid.clear();
		face_in_grid_list_len = 0;
		face_num = 0;
		if (grid_pos_type)
		{
			delete[] grid_pos_type;
			grid_pos_type = nullptr;
		}
		if (face_in_grid_range)
		{
			delete[] face_in_grid_range;
			face_in_grid_range = nullptr;
		}
		if (face_in_grid_list)
		{
			delete[] face_in_grid_list;
			face_in_grid_list = nullptr;
		}
		if (pt_tri_dist)
		{
			delete[] pt_tri_dist;
			pt_tri_dist = nullptr;
		}
	}

	inline const size_t get_grid_xl() const noexcept { return grid.xl; }
	inline const size_t get_grid_yl() const noexcept { return grid.yl; }
	inline const size_t get_grid_zl() const noexcept { return grid.zl; }
	inline const size_t get_grid_hx() const noexcept { return grid.hx; }
	inline const size_t get_grid_hy() const noexcept { return grid.hy; }
	inline const size_t get_grid_hz() const noexcept { return grid.hz; }
	inline const size_t get_grid_xu() const noexcept { return grid.xu; }
	inline const size_t get_grid_yu() const noexcept { return grid.yu; }
	inline const size_t get_grid_zu() const noexcept { return grid.zu; }
	inline const size_t get_grid_x_num() const noexcept { return grid.x_num; }
	inline const size_t get_grid_y_num() const noexcept { return grid.y_num; }
	inline const size_t get_grid_z_num() const noexcept { return grid.z_num; }
	inline const size_t get_grid_num() const noexcept { return grid.num; }
	inline void get_bbox(Cube& bbox) const noexcept
	{
		bbox.xl = grid.xl; bbox.yl = grid.yl; bbox.zl = grid.zl;
		bbox.xu = grid.xu; bbox.yu = grid.yu; bbox.zu = grid.zu;
	}

	// grid_num
	inline const GridPosType* get_grid_pos_type() const noexcept { return grid_pos_type; }
	// grid_num + 1
	inline const size_t* get_face_in_grid_range() const noexcept { return face_in_grid_range; }
	inline size_t get_face_in_grid_list_len() const noexcept { return face_in_grid_list_len; }
	inline const size_t* get_face_in_grid_list() const noexcept { return face_in_grid_list; }
	inline size_t get_face_num() const noexcept { return face_num; }
	inline const PointToTriangleDistance* get_pt_tri_dist() const noexcept { return pt_tri_dist; }

	int init_from_mesh(TetrahedronMesh &tmesh,
		double ghx, double ghy, double ghz);
	// Speed up searching
	// must be called after init_from_mesh()
	void init_max_dist(double _dist) noexcept;
	
	bool detect_collision_with_point(const Point3D& pt,
		const double p_r, double& dist, Vector3D& norm) const noexcept;
};

#endif