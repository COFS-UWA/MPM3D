#ifndef __Test_Parallel_Utils_h__
#define __Test_Parallel_Utils_h__

#include "ItemArray.hpp"
#include "Geometry2D.h"

typedef MemoryUtils::ItemArray<size_t> IndexArray;
typedef MemoryUtils::ItemArray<Point2D> Point2DArray;

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
	const float rxl = float(range.xl);
	const float rxu = float(range.xu);
	const float ryl = float(range.yl);
	const float ryu = float(range.yu);
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		const PclPos& pp = pcl_pos[p_id];
		if (pp.x >= rxl && pp.x <= rxu &&
			pp.y >= ryl && pp.y <= ryu)
			pt_array.add(p_id);
	}
}

#endif