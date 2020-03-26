#ifndef __Tetrahedron_Mesh_H__
#define __Tetrahedron_Mesh_H__

#include "ItemArray.hpp"
#include "ItemBuffer.hpp"

#include "Geometry.h"

struct TetrahedronMesh
{
public:
	//template <typename Item>
	//struct ItemPointer
	//{
	//	Item *item;
	//	ItemPointer *next;
	//};

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

	struct Edge
	{
		size_t n1, n2;
	};

protected:
	size_t node_num;
	Node *nodes;
	size_t elem_num;
	Element *elems;

	// geometric properties
	Cube bounding_box;

	Point3D centre;
	
public:
	TetrahedronMesh() :
		nodes(nullptr), node_num(0),
		elems(nullptr), elem_num(0) {}
	~TetrahedronMesh() { clear(); }

	inline void clear(void)
	{
		clear_nodes();
		clear_elements();
	}

	inline void clear_nodes(void)
	{
		if (nodes) delete[] nodes;
		nodes = nullptr;
		node_num = 0;
	}

	inline void clear_elements(void)
	{
		if (elems) delete[] elems;
		elems = nullptr;
		elem_num = 0;
	}

	inline void alloc_nodes(size_t num)
	{
		if (num == 0) return;
		clear_nodes();
		nodes = new Node[num];
		node_num = num;
	}

	void alloc_elements(size_t num)
	{
		if (num == 0) return;
		clear_elements();
		elems = new Element[num];
		elem_num = num;
	}

	inline size_t get_node_num() const noexcept { return node_num; }
	inline Node *get_nodes(void) const noexcept { return nodes; }
	inline size_t get_elem_num(void) const noexcept { return elem_num; }
	inline Element *get_elems(void) const noexcept { return elems; }

	bool is_in_tetrahedron(Element &elem, Point3D &p);

protected: // helper functions for initializing mesh 
	void compress_node_and_elem_indices();
	// reorder node to counter-clockwise order
	void cal_area_and_reorder_node();

public:
	// read mesh
	int init_mesh(double *node_coords, size_t node_num, size_t *elem_indices, size_t elem_num);
	int load_mesh_hdf5(const char *file_name);
	
	int init_bounding_box();
	int init_edges();
};

inline double distance(TetrahedronMesh::Node &n1, TetrahedronMesh::Node &n2) noexcept
{
	double x_diff = n1.x - n2.x;
	double y_diff = n1.y - n2.y;
	return sqrt(x_diff * x_diff + y_diff * y_diff);
}

#endif