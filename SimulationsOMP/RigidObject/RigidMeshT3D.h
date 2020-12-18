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
		AtBoundary = 0,
		Inside = 1,
		Outside = 2,
		FarInside = 5,
		FarOutside = 6
	};

protected:
	Grid3D<> grid;
	GridPosType* grid_pos_type;
	size_t *face_in_grid_range;
	size_t *face_in_grid_list;
	PointToTriangleDistance* pt_tri_dist;
	
	// grid_xu_min_hx = grid_xu - grid_hx
	double grid_xu_min_hx, grid_yu_min_hy, grid_zu_min_hz;
	// grid_max_x_id = grid_x_num - 1
	size_t grid_max_x_id, grid_max_y_id, grid_max_z_id;

	inline size_t offset_from_xyz_id(size_t x_id, size_t y_id, size_t z_id) const noexcept
	{ return grid.offset_from_xyz_id(x_id, y_id, z_id); }

	double max_dist;
	
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
	void search_closest_face(ClosestFaceRes &cloeset_face,
		const Point3D &pt, size_t stride,
		const IdRange &id_range) const noexcept;

public:
	RigidMeshT3D();
	~RigidMeshT3D();

	inline void clear() noexcept
	{
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

	int init_from_mesh(TetrahedronMesh &tmesh,
		double ghx, double ghy, double ghz);
	//  speed up searching, called after init_from_hdf5
	void init_max_dist(double _dist);
	
	bool detect_collision_with_point(const Point3D& pt,
		const double p_r, double& dist, Vector3D& norm) const noexcept;
};

#endif