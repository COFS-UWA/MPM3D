#include "QtPostProcessor_pcp.h"

#include <fstream>
#include "hdf5.h"

#include "TetrahedronMesh.h"

namespace // helper funcs
{
template<typename _Type>
inline void swap(_Type &a, _Type &b) noexcept { _Type tmp = a; a = b, b = tmp; }
}

bool TetrahedronMesh::is_in_tetrahedron(Element &elem, Point3D &p)
{
	Node &n1 = nodes[elem.n1];
	Node &n2 = nodes[elem.n2];
	Node &n3 = nodes[elem.n3];
	Node &n4 = nodes[elem.n4];
	//double a = (n2.x - n1.x) * (p.y - n1.y) - (p.x - n1.x) * (n2.y - n1.y);
	//double b = (p.x - n1.x) * (n3.y - n1.y) - (n3.x - n1.x) * (p.y - n1.y);
	//double c = elem.area - a - b;
	//return 0.0 <= a && a <= elem.area
	//	&& 0.0 <= b && b <= elem.area 
	//	&& 0.0 <= c && c <= elem.area;
	return false;
}

void TetrahedronMesh::compress_node_and_elem_indices()
{
	//// "compress" nodes and elements index
	//HashTable node_id_map(node_num);
	//for (size_t i = 0; i < node_num; ++i)
	//{
	//	node_id_map.add_pair(nodes[i].id, i);
	//	nodes[i].id = i;
	//}
	//for (size_t i = 0; i < elem_num; ++i)
	//{
	//	size_t n_id;
	//	Element &elem = elems[i];
	//	node_id_map.get_pair(elem.n1, n_id);
	//	elem.n1 = n_id;
	//	node_id_map.get_pair(elem.n2, n_id);
	//	elem.n2 = n_id;
	//	node_id_map.get_pair(elem.n3, n_id);
	//	elem.n3 = n_id;
	//	elem.id = i;
	//}
}

void TetrahedronMesh::cal_area_and_reorder_node(void)
{
	//// Calculate area of triangle.
	//area = 0.0;
	//x_mc = 0.0;
	//y_mc = 0.0;
	//moi_area = 0.0;
	//double elem_area;
	//for (size_t i = 0; i < elem_num; ++i)
	//{
	//	Element &elem = elems[i];
	//	Node &n1 = nodes[elem.n1];
	//	Node &n2 = nodes[elem.n2];
	//	Node &n3 = nodes[elem.n3];
	//	double x1, y1, x2, y2, x3, y3;
	//	x1 = n1.x;
	//	y1 = n1.y;
	//	x2 = n2.x;
	//	y2 = n2.y;
	//	x3 = n3.x;
	//	y3 = n3.y;
	//	elem.area = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
	//	// Ensure that nodes index of element is in counter-clockwise sequence.
	//	if (elem.area < 0.0)
	//	{
	//		elem.area = -elem.area;
	//		size_t n_tmp = elem.n2;
	//		elem.n2 = elem.n3;
	//		elem.n3 = n_tmp;
	//	}

	//	elem_area = elem.area / 2.0;
	//	area += elem_area;
	//	x_mc += elem_area * (x1 + x2 + x3) / 3.0;
	//	y_mc += elem_area * (y1 + y2 + y3) / 3.0;
	//	moi_area += (x1 * x1 + x2 * x2 + x3 * x3
	//		+ x1 * x2 + x2 * x3 + x3 * x1
	//		+ y1 * y1 + y2 * y2 + y3 * y3
	//		+ y1 * y2 + y2 * y3 + y3 * y1) * elem_area / 6.0;
	//}
	//x_mc /= area;
	//y_mc /= area;
	//moi_area -= area * (x_mc * x_mc + y_mc * y_mc);
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

int TetrahedronMesh::load_mesh_hdf5(const char *file_name)
{
	hid_t mesh_file = H5Fopen(file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
	H5Dclose(mesh_file);
	return 0;
}


int TetrahedronMesh::init_bounding_box()
{
	if (node_num == 0)
		return -1;

	// Init bounding box
	bounding_box.xl = nodes[0].x;
	bounding_box.xu = bounding_box.xl;
	bounding_box.yl = nodes[0].y;
	bounding_box.yu = bounding_box.yl;
	bounding_box.zl = nodes[0].z;
	bounding_box.zu = bounding_box.zl;
	for (size_t n_id = 1; n_id < node_num; ++n_id)
	{
		if (bounding_box.xl > nodes[n_id].x)
			bounding_box.xl = nodes[n_id].x;
		if (bounding_box.xu < nodes[n_id].x)
			bounding_box.xu = nodes[n_id].x;
		if (bounding_box.yl > nodes[n_id].y)
			bounding_box.yl = nodes[n_id].y;
		if (bounding_box.yu < nodes[n_id].y)
			bounding_box.yu = nodes[n_id].y;
		if (bounding_box.zl > nodes[n_id].z)
			bounding_box.zl = nodes[n_id].z;
		if (bounding_box.zu < nodes[n_id].z)
			bounding_box.zu = nodes[n_id].z;
	}

	return 0;
}
