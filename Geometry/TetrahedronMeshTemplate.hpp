#ifndef __Tetrahedron_Mesh_Template_hpp__
#define __Tetrahedron_Mesh_Template_hpp__

#include "ItemArray.hpp"
#include "ItemArrayFast.hpp"
#include "ItemBuffer.hpp"

#include "hdf5.h"

#include "Geometry.h"
#include "NumPairHashTable.hpp"

/*============================= 
Assumptions on type:
struct Node
{
	size_t id;
	double x, y, z;
};
struct Element
{
	size_t id;
	size_t n1, n2, n3, n4;
	double vol;
};
struct Edge { size_t n1, n2; };
 ==============================*/

template <typename _Node, typename _Element, typename _Edge>
struct TetrahedronMeshTemplate
{
public:
	typedef _Node Node;
	typedef _Element Element;
	typedef _Edge Edge;

protected:
	size_t node_num;
	Node *nodes;
	size_t elem_num;
	Element *elems;
	
	size_t edge_num;
	Edge *edges;

	// geometric properties
	double volume;
	Point3D centre;
	Cube bounding_box;

public:
	TetrahedronMeshTemplate() :
		nodes(nullptr), node_num(0),
		elems(nullptr), elem_num(0),
		edges(nullptr), edge_num(0) {}
	~TetrahedronMeshTemplate() { clear(); }

	inline void clear()
	{
		clear_nodes();
		clear_elements();
		clear_edges();
	}

	inline void clear_nodes()
	{
		if (nodes)
		{
			delete[] nodes;
			nodes = nullptr;
		}
		node_num = 0;
	}

	inline void clear_elements()
	{
		if (elems)
		{
			delete[] elems;
			elems = nullptr;
		}
		elem_num = 0;
	}

	inline void clear_edges()
	{
		if (edges)
		{
			delete[] edges;
			edges = nullptr;
		}
		edge_num = 0;
	}

	inline void alloc_nodes(size_t num)
	{
		clear_nodes();
		if (num == 0)
			return;
		nodes = new Node[num];
		node_num = num;
	}

	inline void alloc_elements(size_t num)
	{
		clear_elements();
		if (num == 0)
			return;
		elems = new Element[num];
		elem_num = num;
	}

	inline void alloc_edges(size_t num)
	{
		clear_edges();
		if (num == 0)
			return;
		edges = new Edge[num];
		edge_num = num;
	}

	inline size_t get_node_num() const noexcept { return node_num; }
	inline Node *get_nodes() const noexcept { return nodes; }
	inline size_t get_elem_num() const noexcept { return elem_num; }
	inline Element *get_elems() const noexcept { return elems; }
	inline size_t get_edge_num() const noexcept { return edge_num; }
	inline Edge *get_edges() const noexcept { return edges; }

	inline double get_vol() { return volume; }
	inline Point3D get_centre() { return centre; }
	inline Cube get_bounding_box() { return bounding_box; }

	inline bool is_in_tetrahedron(Element &elem, double x, double y, double z)
	{
		Point3D p = { x, y, z };
		return is_in_tetrahedron<Point3D>(elem, p);
	}

	template <typename Point3D>
	inline bool is_in_tetrahedron(Element &elem, Point3D &p)
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

protected: // helper for load_mesh_from_hdf5()
	struct NodeData
	{
		long long id;
		double x;
		double y;
		double z;
	};
	inline static hid_t get_node_dt_id()
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
	inline static hid_t get_elem_dt_id()
	{
		hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(ElemData));
		H5Tinsert(res, "index", HOFFSET(ElemData, id), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "n1", HOFFSET(ElemData, n1), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "n2", HOFFSET(ElemData, n2), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "n3", HOFFSET(ElemData, n3), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "n4", HOFFSET(ElemData, n4), H5T_NATIVE_ULLONG);
		return res;
	}

	// compress node and elem indices
	void compress_node_and_elem_indices()
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

	// reorder node to counter-clockwise order
	void cal_area_and_reorder_node()
	{
		size_t elem_nd;
		volume = 0.0;
		centre.x = 0.0;
		centre.y = 0.0;
		centre.z = 0.0;
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
			centre.x += (n1.x + n2.x + n3.x + n4.x) * 0.25 * elem.vol;
			centre.y += (n1.y + n2.y + n3.y + n4.y) * 0.25 * elem.vol;
			centre.z += (n1.z + n2.z + n3.z + n4.z) * 0.25 * elem.vol;
		}
		centre.x /= volume;
		centre.y /= volume;
		centre.z /= volume;
	}

	// get bounding box
	void cal_bounding_box()
	{
		if (node_num == 0) return;
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
	}

public:
	int load_mesh_from_hdf5(const char *file_name)
	{
		clear();
		hid_t mesh_file = H5Fopen(file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
		if (mesh_file < 0)
			return -1;
		hid_t mesh_grp = H5Gopen(mesh_file, "Mesh", H5P_DEFAULT);
		if (mesh_grp < 0)
			return -1;
		herr_t res;
		hsize_t dim;
		hid_t space_id, dtype_id;
		MemoryUtils::ItemArrayFast<char> mem_buf;
		// read nodes
		hid_t node_dset = H5Dopen(mesh_grp, "Nodes", H5P_DEFAULT);
		space_id = H5Dget_space(node_dset);
		H5Sget_simple_extent_dims(space_id, &dim, nullptr);
		node_num = dim;
		dtype_id = get_node_dt_id();
		mem_buf.resize(sizeof(NodeData) * node_num);
		NodeData *node_data_buf = (NodeData *)mem_buf.get_mem();
		res = H5Dread(node_dset, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, node_data_buf);
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
		// read elements
		hid_t elem_dset = H5Dopen(mesh_grp, "Elements", H5P_DEFAULT);
		space_id = H5Dget_space(elem_dset);
		H5Sget_simple_extent_dims(space_id, &dim, nullptr);
		elem_num = dim;
		dtype_id = get_elem_dt_id();
		mem_buf.resize(sizeof(ElemData) * elem_num);
		ElemData *elem_data_buf = (ElemData *)mem_buf.get_mem();
		res = H5Dread(elem_dset, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, elem_data_buf);
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
		// finish reading
		H5Gclose(mesh_grp);
		H5Fclose(mesh_file);
		// pre process data
		compress_node_and_elem_indices();
		cal_area_and_reorder_node();
		// init geometric properties
		cal_bounding_box();
		// init edges
		init_edges();
		return res;
	}

	int init_mesh(double* node_coords, size_t node_num, size_t* elem_indices, size_t elem_num)
	{
		if (node_num == 0 || elem_num == 0)
			return -1;
		clear();
		// init nodes
		alloc_nodes(node_num);
		double* pnode_coord = node_coords;
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
		// init elements
		alloc_elements(elem_num);
		size_t* pelem_index = elem_indices;
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
		// pre process data
		compress_node_and_elem_indices();
		cal_area_and_reorder_node();
		// init geometric properties
		cal_bounding_box();
		// init edges
		init_edges();
		return 0;
	}

	// for restarting calculation from file
	// initialize properties after loading
	// nodes and elements from file  
	void init_mesh_properties_after_loading()
	{
		cal_area_and_reorder_node();
		cal_bounding_box();
		init_edges();
	}

protected: // helpers for init edges
	// sort ids[4] in accending order
	inline static void sort_acc(size_t ids[4])
	{
		size_t tmp, min_id;
		for (size_t i = 0; i < 3; ++i)
		{
			min_id = i;
			for (size_t j = i + 1; j < 4; ++j)
			{
				if (ids[j] < ids[min_id])
					min_id = j;
			}
			if (min_id != i)
			{
				tmp = ids[min_id];
				ids[min_id] = ids[i];
				ids[i] = tmp;
			}
		}
	}

public:
	int init_edges()
	{
		clear_edges();

		NumPairHashTable<size_t> table(node_num * 2);
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
};

#endif