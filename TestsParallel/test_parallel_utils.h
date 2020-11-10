#ifndef __Test_Parallel_Utils_h__
#define __Test_Parallel_Utils_h__

#include "ItemArray.hpp"
#include "Geometry2D.h"
#include "Geometry3D.h"

typedef MemoryUtils::ItemArray<size_t> IndexArray;
typedef MemoryUtils::ItemArray<Point2D> Point2DArray;
typedef MemoryUtils::ItemArray<Point3D> Point3DArray;

template <typename Model>
void find_2d_nodes_on_x_line(Model& md, IndexArray& id_array,
	float x, bool need_reset_array = true, float tol = 1.0e-3)
{
	if (need_reset_array)
		id_array.reset();

	tol = abs(x) < 1.0 ? tol : abs(x) * tol;
	const float xl = x - tol;
	const float xu = x + tol;
	size_t node_num = md.get_node_num();
	typedef Model::NodePos NodePos;
	const NodePos *node_pos = md.get_node_pos();
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		const NodePos& np = node_pos[n_id];
		if (np.x > xl && np.x < xu)
			id_array.add(n_id);
	}
}

template <typename Model>
void find_2d_nodes_on_y_line(Model& md, IndexArray& id_array,
	float y, bool need_reset_array = true, float tol = 1.0e-3)
{
	if (need_reset_array)
		id_array.reset();

	tol = abs(y) < 1.0 ? tol : abs(y) * tol;
	const float yl = y - tol;
	const float yu = y + tol;
	size_t node_num = md.get_node_num();
	typedef Model::NodePos NodePos;
	const NodePos* node_pos = md.get_node_pos();
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		const NodePos& np = node_pos[n_id];
		if (np.y > yl && np.y < yu)
			id_array.add(n_id);
	}
}

template <typename Model>
void find_2d_pcls(Model& md, IndexArray& pt_array,
	Rect& range, bool need_reset_array = true)
{
	if (need_reset_array)
		pt_array.reset();

	size_t pcl_num = md.get_pcl_num();
	typedef Model::PclPos PclPos;
	const PclPos *pcl_pos = md.get_pcl_pos();
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		const PclPos& pp = pcl_pos[p_id];
		if (pp.x >= range.xl && pp.x <= range.xu &&
			pp.y >= range.yl && pp.y <= range.yu)
			pt_array.add(p_id);
	}
}

template <typename Model>
void find_3d_nodes_on_x_plane(Model& md, IndexArray& id_array,
	double x, bool need_reset_array = true, double tol = 1.0e-3)
{
	typedef typename Model::Node Node;

	if (need_reset_array)
		id_array.reset();

	size_t node_num = md.get_node_num();
	typedef Model::NodePos NodePos;
	const NodePos* node_pos = md.get_node_pos();
	tol = abs(x) < 1.0 ? tol : abs(x) * tol;
	double xl = x - tol;
	double xu = x + tol;
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node& n = node_pos[n_id];
		if (n.x > xl && n.x < xu)
			id_array.add(n.id);
	}
}

template <typename Model>
void find_3d_nodes_on_y_plane(Model& md, IndexArray& id_array,
	double y, bool need_reset_array = true, double tol = 1.0e-3)
{
	typedef typename Model::Node Node;

	if (need_reset_array)
		id_array.reset();

	size_t node_num = md.get_node_num();
	typedef Model::NodePos NodePos;
	const NodePos* node_pos = md.get_node_pos();
	tol = abs(y) < 1.0 ? tol : abs(y) * tol;
	double yl = y - tol;
	double yu = y + tol;
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node& n = node_pos[n_id];
		if (n.y > yl && n.y < yu)
			id_array.add(n.id);
	}
}

template <typename Model>
void find_3d_nodes_on_z_plane(Model& md, IndexArray& id_array,
	double z, bool need_reset_array = true, double tol = 1.0e-3)
{
	typedef typename Model::Node Node;

	if (need_reset_array)
		id_array.reset();

	size_t node_num = md.get_node_num();
	typedef Model::NodePos NodePos;
	const NodePos* node_pos = md.get_node_pos();
	tol = abs(z) < 1.0 ? tol : abs(z) * tol;
	double zl = z - tol;
	double zu = z + tol;
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node& n = node_pos[n_id];
		if (n.z > zl && n.z < zu)
			id_array.add(n.id);
	}
}

template <typename Model>
void find_3d_pcls(Model& md, IndexArray& pt_array,
	Cube& range, bool need_reset_array = true)
{
	if (need_reset_array)
		pt_array.reset();

	size_t pcl_num = md.get_pcl_num();
	typedef Model::Position Position;
	const Position* pcl_pos = md.get_pcl_pos();
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		const Position& pp = pcl_pos[p_id];
		if (pp.x >= range.xl && pp.x <= range.xu &&
			pp.y >= range.yl && pp.y <= range.yu &&
			pp.z >= range.zl && pp.z <= range.zu)
			pt_array.add(p_id);
	}
}

#endif