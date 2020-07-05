#ifndef __Triangle_Mesh_Template_hpp__
#define __Triangle_Mesh_Template_hpp__

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
	double x, y;
};
struct Element
{
	size_t id;
	size_t n1, n2, n3;
	double area;
};
struct Edge { size_t n1, n2; };
 ==============================*/

template <typename _Node, typename _Element, typename _Edge>
struct TriangleMeshTemplate
{
public:
	typedef _Node Node;
	typedef _Element Element;
	typedef _Edge Edge;

protected:
	size_t node_num;
	Node* nodes;
	size_t elem_num;
	Element* elems;

	size_t edge_num;
	Edge* edges;

	// geometric properties
	double area;
	Point2D centre;
	Rect bounding_box;

public:
	TriangleMeshTemplate() :
		nodes(nullptr), node_num(0),
		elems(nullptr), elem_num(0),
		edges(nullptr), edge_num(0) {}
	~TriangleMeshTemplate() { clear(); }

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
		if (num == 0) return;
		clear_nodes();
		nodes = new Node[num];
		node_num = num;
	}

	inline void alloc_elements(size_t num)
	{
		if (num == 0) return;
		clear_elements();
		elems = new Element[num];
		elem_num = num;
	}

	inline void alloc_edges(size_t num)
	{
		if (num == 0) return;
		clear_edges();
		edges = new Edge[num];
		edge_num = num;
	}

	inline size_t get_node_num() const noexcept { return node_num; }
	inline Node* get_nodes() const noexcept { return nodes; }
	inline size_t get_elem_num() const noexcept { return elem_num; }
	inline Element* get_elems() const noexcept { return elems; }
	inline size_t get_edge_num() const noexcept { return edge_num; }
	inline Edge* get_edges() const noexcept { return edges; }

	inline double get_area() { return area; }
	inline Point2D get_centre() { return centre; }
	inline Rect get_bounding_box() { return bounding_box; }

	inline bool is_in_triangle(Element& elem, double x, double y)
	{
		Point2D p = { x, y };
		return is_in_triangle<Point2D>(elem, p);
	}

	template <typename Point3D>
	inline bool is_in_triangle(Element& elem, Point3D& p)
	{
		Node& n1 = nodes[elem.n1];
		Node& n2 = nodes[elem.n2];
		Node& n3 = nodes[elem.n3];

		double area1, area2, area3;
		area1 = cal_triangle_area<Node, Point2D>(n1, n2, p);
		area2 = cal_triangle_area<Node, Point2D>(n2, n3, p);
		area3 = cal_triangle_area<Node, Point2D>(n3, n1, p);

		if (area1 >= 0.0 && area1 <= elem.area &&
			area2 >= 0.0 && area2 <= elem.area &&
			area3 >= 0.0 && area3 <= elem.area)
			return true;
		return false;
	}

public:
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

protected: // helper for load_mesh_from_hdf5()
	struct NodeData
	{
		long long id;
		double x;
		double y;
	};
	inline static hid_t get_node_dt_id()
	{
		hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(NodeData));
		H5Tinsert(res, "index", HOFFSET(NodeData, id), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "x", HOFFSET(NodeData, x), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "y", HOFFSET(NodeData, y), H5T_NATIVE_DOUBLE);
		return res;
	}

	struct ElemData
	{
		long long id;
		long long n1;
		long long n2;
		long long n3;
	};
	inline static hid_t get_elem_dt_id()
	{
		hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(ElemData));
		H5Tinsert(res, "index", HOFFSET(ElemData, id), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "n1", HOFFSET(ElemData, n1), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "n2", HOFFSET(ElemData, n2), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "n3", HOFFSET(ElemData, n3), H5T_NATIVE_ULLONG);
		return res;
	}

	// compress node and elem indices
	void compress_node_and_elem_indices()
	{
		NumPairHashTable<size_t> table(node_num);
		for (size_t n_id = 0; n_id < node_num; ++n_id)
		{
			Node& n = nodes[n_id];
			table.add_key(n.id, n_id);
			n.id = n_id;
		}

		size_t new_id;
		for (size_t e_id = 0; e_id < elem_num; ++e_id)
		{
			Element& e = elems[e_id];
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
		}
	}

	// reorder node to counter-clockwise order
	void cal_area_and_reorder_node()
	{
		size_t elem_nd;
		area = 0.0;
		centre.x = 0.0;
		centre.y = 0.0;
		for (size_t e_id = 0; e_id < elem_num; ++e_id)
		{
			Element& elem = elems[e_id];
			Node& n1 = nodes[elem.n1];
			Node& n2 = nodes[elem.n2];
			Node& n3 = nodes[elem.n3];
			elem.area = cal_triangle_area<Node>(n1, n2, n3);
			if (elem.area < 0.0)
			{
				elem_nd = elem.n2;
				elem.n2 = elem.n3;
				elem.n3 = elem_nd;
				elem.area = -elem.area;
			}
			area += elem.area;
			centre.x += (n1.x + n2.x + n3.x) / 3.0 * elem.area;
			centre.y += (n1.y + n2.y + n3.y) / 3.0 * elem.area;
		}
		centre.x /= area;
		centre.y /= area;
	}

	// get bounding box
	void cal_bounding_box()
	{
		if (node_num == 0) return;
		bounding_box.xl = nodes[0].x;
		bounding_box.xu = bounding_box.xl;
		bounding_box.yl = nodes[0].y;
		bounding_box.yu = bounding_box.yl;
		for (size_t n_id = 1; n_id < node_num; ++n_id)
		{
			Node& n = nodes[n_id];
			if (bounding_box.xl > n.x)
				bounding_box.xl = n.x;
			if (bounding_box.xu < n.x)
				bounding_box.xu = n.x;
			if (bounding_box.yl > n.y)
				bounding_box.yl = n.y;
			if (bounding_box.yu < n.y)
				bounding_box.yu = n.y;
		}
	}

public:
	int load_mesh_from_hdf5(const char* file_name)
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
		NodeData* node_data_buf = (NodeData*)mem_buf.get_mem();
		res = H5Dread(node_dset, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, node_data_buf);
		H5Sclose(space_id);
		H5Tclose(dtype_id);
		H5Dclose(node_dset);
		alloc_nodes(node_num);
		for (size_t n_id = 0; n_id < node_num; ++n_id)
		{
			NodeData& nd = node_data_buf[n_id];
			Node& n = nodes[n_id];
			n.id = nd.id;
			n.x = nd.x;
			n.y = nd.y;
		}
		// read elements
		hid_t elem_dset = H5Dopen(mesh_grp, "Elements", H5P_DEFAULT);
		space_id = H5Dget_space(elem_dset);
		H5Sget_simple_extent_dims(space_id, &dim, nullptr);
		elem_num = dim;
		dtype_id = get_elem_dt_id();
		mem_buf.resize(sizeof(ElemData) * elem_num);
		ElemData* elem_data_buf = (ElemData*)mem_buf.get_mem();
		res = H5Dread(elem_dset, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, elem_data_buf);
		H5Sclose(space_id);
		H5Tclose(dtype_id);
		H5Dclose(elem_dset);
		alloc_elements(elem_num);
		for (size_t e_id = 0; e_id < elem_num; ++e_id)
		{
			ElemData& ed = elem_data_buf[e_id];
			Element& e = elems[e_id];
			e.id = ed.id;
			e.n1 = ed.n1;
			e.n2 = ed.n2;
			e.n3 = ed.n3;
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

protected: // helpers for init edges
	// sort ids[3] in accending order
	inline static void sort_acc(size_t ids[3])
	{
		size_t tmp, min_id;
		min_id = 0;
		if (ids[0] > ids[1])
			min_id = 1;
		if (ids[0] > ids[2])
			min_id = 2;
		if (min_id != 0)
		{
			tmp = ids[0];
			ids[0] = ids[min_id];
			ids[min_id] = tmp;
		}
		if (ids[1] > ids[2])
		{
			tmp = ids[1];
			ids[1] = ids[min_id];
			ids[min_id] = tmp;
		}
	}

public:
	int init_edges()
	{
		clear_edges();

		NumPairHashTable<size_t> table(node_num * 2);
		union
		{
			struct { size_t n1_id, n2_id, n3_id; };
			size_t n_ids[3];
		};

		for (size_t e_id = 0; e_id < elem_num; ++e_id)
		{
			Element& e = elems[e_id];
			n1_id = e.n1;
			n2_id = e.n2;
			n3_id = e.n3;
			sort_acc(n_ids); // n1_id < n2_id < n3_id
			table.add_pair(n1_id, n2_id);
			table.add_pair(n1_id, n3_id);
			table.add_pair(n2_id, n3_id);
		}
		alloc_edges(table.get_pair_num());
		table.output_pairs((size_t*)edges);
		return 0;
	}
};

#endif