#ifndef __Test_Parallel_Utils2_h__
#define __Test_Parallel_Utils2_h__

#include "ItemArray.hpp"
#include "Geometry2D.h"
#include "Geometry3D.h"

typedef MemoryUtils::ItemArray<size_t> IndexArray;

template <typename Model>
void find_2d_nodes_on_x_line(Model& md, IndexArray& id_array,
	double x, bool need_reset_array = true, double tol = 1.0e-3)
{
	typedef typename Model::Node Node;

	if (need_reset_array)
		id_array.reset();

	size_t node_num = md.get_node_num();
	Node* nodes = md.get_nodes();
	tol = abs(x) < 1.0 ? tol : abs(x) * tol;
	double xl = x - tol;
	double xu = x + tol;
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node& n = nodes[n_id];
		if (n.x > xl && n.x < xu)
			id_array.add(n.id);
	}
}

template <typename Model>
void find_2d_nodes_on_y_line(Model& md, IndexArray& id_array,
	double y, bool need_reset_array = true, double tol = 1.0e-3)
{
	typedef typename Model::Node Node;

	if (need_reset_array)
		id_array.reset();

	size_t node_num = md.get_node_num();
	Node* nodes = md.get_nodes();
	tol = abs(y) < 1.0 ? tol : abs(y) * tol;
	double yl = y - tol;
	double yu = y + tol;
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node& n = nodes[n_id];
		if (n.y > yl && n.y < yu)
			id_array.add(n.id);
	}
}

template <typename Model>
void find_2d_pcls(Model& md, IndexArray& pt_array,
	Rect& range, bool need_reset_array = true)
{
	typedef typename Model::Particle Particle;

	if (need_reset_array)
		pt_array.reset();

	size_t pcl_num = md.get_pcl_num();
	Particle* pcls = md.get_pcls();
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		Particle& pcl = pcls[p_id];
		if (pcl.x >= range.xl && pcl.x <= range.xu &&
			pcl.y >= range.yl && pcl.y <= range.yu)
			pt_array.add(pcl.id);
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
	Node* nodes = md.get_nodes();
	tol = abs(x) < 1.0 ? tol : abs(x) * tol;
	double xl = x - tol;
	double xu = x + tol;
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node& n = nodes[n_id];
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
	Node* nodes = md.get_nodes();
	tol = abs(y) < 1.0 ? tol : abs(y) * tol;
	double yl = y - tol;
	double yu = y + tol;
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node& n = nodes[n_id];
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
	Node* nodes = md.get_nodes();
	tol = abs(z) < 1.0 ? tol : abs(z) * tol;
	double zl = z - tol;
	double zu = z + tol;
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node& n = nodes[n_id];
		if (n.z > zl && n.z < zu)
			id_array.add(n.id);
	}
}

template <typename Model>
void find_3d_pcls(Model& md, IndexArray& pt_array,
	Cube& range, bool need_reset_array = true)
{
	typedef typename Model::Particle Particle;

	if (need_reset_array)
		pt_array.reset();

	size_t pcl_num = md.get_pcl_num();
	Particle* pcls = md.get_pcls();
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		Particle& pcl = pcls[p_id];
		if (pcl.x >= range.xl && pcl.x <= range.xu &&
			pcl.y >= range.yl && pcl.y <= range.yu &&
			pcl.z >= range.zl && pcl.z <= range.zu)
			pt_array.add(pcl.id);
	}
}

#endif