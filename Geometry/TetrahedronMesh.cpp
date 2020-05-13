#include "Geometry_pcp.h"

#include <iostream>

#include "hdf5.h"
#include "ItemArrayFast.hpp"
#include "Geometry.h"
#include "NumPairHashTable.hpp"

#include "TetrahedronMesh.h"

inline bool TetrahedronMesh::is_in_tetrahedron(Element &elem, Point3D &p)
{
	Node &n1 = nodes[elem.n1];
	Node &n2 = nodes[elem.n2];
	Node &n3 = nodes[elem.n3];
	Node &n4 = nodes[elem.n4];

	double vol1, vol2, vol3, vol4;
	vol1 = cal_tetrahedron_vol<Node, Point3D>(n1, n2, n3, p);
	vol2 = cal_tetrahedron_vol<Node, Point3D>(n1, n4, n2, p);
	vol3 = cal_tetrahedron_vol<Node, Point3D>(n2, n4, n3, p);
	vol4 = cal_tetrahedron_vol<Node, Point3D>(n1, n3, n4, p);

	if (vol1 >= 0.0 && vol1 <= elem.vol && vol2 >= 0.0 && vol2 <= elem.vol &&
		vol3 >= 0.0 && vol3 <= elem.vol && vol4 >= 0.0 && vol4 <= elem.vol)
		return true;
	return false;
}

void TetrahedronMesh::compress_node_and_elem_indices()
{
	NumPairHashTable<size_t> table(node_num);
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node &n = nodes[n_id];
		table.add_key(n.id, n_id);
		n.id = n_id;
	}

	size_t new_id;
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		Element &e = elems[e_id];
		e.id = e_id;
		// n1
		table.find_key(e.n1, new_id);
		e.n1 = new_id;
		// n2
		table.find_key(e.n2, new_id);
		e.n2 = new_id;
		// n3
		table.find_key(e.n3, new_id);
		e.n3 = new_id;
		// n4
		table.find_key(e.n4, new_id);
		e.n4 = new_id;
	}
}

void TetrahedronMesh::cal_area_and_reorder_node()
{
	size_t elem_nd;
	volume = 0.0;
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		Element &elem = elems[e_id];
		Node &n1 = nodes[elem.n1];
		Node &n2 = nodes[elem.n2];
		Node &n3 = nodes[elem.n3];
		Node &n4 = nodes[elem.n4];
		elem.vol = cal_tetrahedron_vol<Node>(n1, n2, n3, n4);
		if (elem.vol < 0.0)
		{
			elem_nd = elem.n2;
			elem.n2 = elem.n3;
			elem.n3 = elem_nd;
			elem.vol = -elem.vol;
		}
		volume += elem.vol;
	}
}

int TetrahedronMesh::init_mesh(
	double *node_coords,
	size_t node_num,
	size_t *elem_indices,
	size_t elem_num
	)
{
	if (node_num == 0 || elem_num == 0)
		return -1;

	clear_nodes();
	alloc_nodes(node_num);
	double *pnode_coord = node_coords;
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		nodes[n_id].id = n_id;
		nodes[n_id].x = *pnode_coord;
		++pnode_coord;
		nodes[n_id].y = *pnode_coord;
		++pnode_coord;
		nodes[n_id].z = *pnode_coord;
		++pnode_coord;
	}

	clear_elements();
	alloc_elements(elem_num);
	size_t *pelem_index = elem_indices;
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		elems[e_id].id = e_id;
		elems[e_id].n1 = *pelem_index;
		++pelem_index;
		elems[e_id].n2 = *pelem_index;
		++pelem_index;
		elems[e_id].n3 = *pelem_index;
		++pelem_index;
		elems[e_id].n4 = *pelem_index;
		++pelem_index;
	}

	return 0;
}

namespace
{
	struct NodeData
	{
		long long id;
		double x;
		double y;
		double z;
	};
	inline hid_t get_node_dt_id(void)
	{
		hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(NodeData));
		H5Tinsert(res, "index", HOFFSET(NodeData, id), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "x", HOFFSET(NodeData, x), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "y", HOFFSET(NodeData, y), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "z", HOFFSET(NodeData, z), H5T_NATIVE_DOUBLE);
		return res;
	}

	struct ElemData
	{
		long long id;
		long long n1;
		long long n2;
		long long n3;
		long long n4;
	};
	inline hid_t get_elem_dt_id(void)
	{
		hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(ElemData));
		H5Tinsert(res, "index", HOFFSET(ElemData, id), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "n1", HOFFSET(ElemData, n1), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "n2", HOFFSET(ElemData, n2), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "n3", HOFFSET(ElemData, n3), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "n4", HOFFSET(ElemData, n4), H5T_NATIVE_ULLONG);
		return res;
	}
};

int TetrahedronMesh::load_mesh_from_hdf5(const char *file_name)
{
	clear_nodes();
	clear_elements();
	hid_t mesh_file = H5Fopen(file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (mesh_file < 0)
		return -1;

	herr_t res;
	hsize_t dim;
	hid_t space_id, dtype_id;
	MemoryUtils::ItemArrayFast<char> mem_buf;
	hid_t mesh_grp = H5Gopen(mesh_file, "Mesh", H5P_DEFAULT);
	// node
	hid_t node_dset = H5Dopen(mesh_grp, "Nodes", H5P_DEFAULT);
	space_id = H5Dget_space(node_dset);
	H5Sget_simple_extent_dims(space_id, &dim, nullptr);
	node_num = dim;
	dtype_id = get_node_dt_id();
	mem_buf.resize(sizeof(NodeData) * node_num);
	NodeData *node_data_buf = (NodeData *)mem_buf.get_mem();
	res = H5Dread(
		node_dset,
		dtype_id,
		H5S_ALL,
		H5S_ALL,
		H5P_DEFAULT,
		node_data_buf
		);
	H5Sclose(space_id);
	H5Tclose(dtype_id);
	H5Dclose(node_dset);
	alloc_nodes(node_num);
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		NodeData &nd = node_data_buf[n_id];
		Node &n = nodes[n_id];
		n.id = nd.id;
		n.x = nd.x;
		n.y = nd.y;
		n.z = nd.z;
	}
	// element
	hid_t elem_dset = H5Dopen(mesh_grp, "Elements", H5P_DEFAULT);
	space_id = H5Dget_space(elem_dset);
	H5Sget_simple_extent_dims(space_id, &dim, nullptr);
	elem_num = dim;
	dtype_id = get_elem_dt_id();
	mem_buf.resize(sizeof(ElemData) * elem_num);
	ElemData *elem_data_buf = (ElemData *)mem_buf.get_mem();
	res = H5Dread(
		elem_dset,
		dtype_id,
		H5S_ALL,
		H5S_ALL,
		H5P_DEFAULT,
		elem_data_buf
		);
	H5Sclose(space_id);
	H5Tclose(dtype_id);
	H5Dclose(elem_dset);
	alloc_elements(elem_num);
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		ElemData &ed = elem_data_buf[e_id];
		Element &e = elems[e_id];
		e.id = ed.id;
		e.n1 = ed.n1;
		e.n2 = ed.n2;
		e.n3 = ed.n3;
		e.n4 = ed.n4;
	}

	H5Gclose(mesh_grp);
	H5Fclose(mesh_file);

	compress_node_and_elem_indices();
	cal_area_and_reorder_node();

	return 0;
}


int TetrahedronMesh::init_bounding_box()
{
	if (node_num == 0)
		return -1;
	bounding_box.xl = nodes[0].x;
	bounding_box.xu = bounding_box.xl;
	bounding_box.yl = nodes[0].y;
	bounding_box.yu = bounding_box.yl;
	bounding_box.zl = nodes[0].z;
	bounding_box.zu = bounding_box.zl;
	for (size_t n_id = 1; n_id < node_num; ++n_id)
	{
		Node &n = nodes[n_id];
		if (bounding_box.xl > n.x)
			bounding_box.xl = n.x;
		if (bounding_box.xu < n.x)
			bounding_box.xu = n.x;
		if (bounding_box.yl > n.y)
			bounding_box.yl = n.y;
		if (bounding_box.yu < n.y)
			bounding_box.yu = n.y;
		if (bounding_box.zl > n.z)
			bounding_box.zl = n.z;
		if (bounding_box.zu < n.z)
			bounding_box.zu = n.z;
	}
	return 0;
}

namespace
{
	void sort_acc(size_t ids[4])
	{
		size_t tmp, min_id;
		min_id = 0;
		if (ids[0] > ids[1])
			min_id = 1;
		if (ids[0] > ids[2])
			min_id = 2;
		if (ids[0] > ids[3])
			min_id = 3;
		if (min_id != 0)
		{
			tmp = ids[0];
			ids[0] = ids[min_id];
			ids[min_id] = tmp;
		}
		min_id = 1;
		if (ids[1] > ids[2])
			min_id = 2;
		if (ids[1] > ids[3])
			min_id = 3;
		if (min_id != 1)
		{
			tmp = ids[1];
			ids[1] = ids[min_id];
			ids[min_id] = tmp;
		}
		if (ids[2] > ids[3])
		{
			tmp = ids[2];
			ids[2] = ids[3];
			ids[3] = tmp;
		}
	}
};

int TetrahedronMesh::init_edges()
{
	NumPairHashTable<size_t> table(node_num*2);
	union
	{
		struct { size_t n1_id, n2_id, n3_id, n4_id; };
		size_t n_ids[4];
	};
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		Element &e = elems[e_id];
		n1_id = e.n1;
		n2_id = e.n2;
		n3_id = e.n3;
		n4_id = e.n4;
		sort_acc(n_ids); // n1_id < n2_id < n3_id < n4_id
		table.add_pair(n1_id, n2_id);
		table.add_pair(n1_id, n3_id);
		table.add_pair(n1_id, n4_id);
		table.add_pair(n2_id, n3_id);
		table.add_pair(n2_id, n4_id);
		table.add_pair(n3_id, n4_id);
	}
	alloc_edges(table.get_pair_num());
	table.output_pairs((size_t *)edges);
	return 0;
}
