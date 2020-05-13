#include "Tests_pcp.h"

#include "TetrahedronMesh.h"

#include "test_geometry.h"

void test_tetrahedron_mesh()
{
	TetrahedronMesh mesh;
	
	//double node_cds[] = {
	//	-1.0, -1.0, -1.0,
	//	 1.0, -1.0, -1.0,
	//	 1.0,  1.0, -1.0,
	//	-1.0,  1.0, -1.0,
	//	-1.0, -1.0,  1.0,
	//	 1.0, -1.0,  1.0,
	//	 1.0,  1.0,  1.0,
	//	-1.0,  1.0,  1.0
	//};

	//size_t elem_ids[] = {
	//	0, 1, 3, 4,
	//	6, 3, 2, 1,
	//	1, 5, 4, 3,
	//	6, 3, 5, 1,
	//	4, 5, 7, 3,
	//	5, 6, 7, 3
	//};

	//mesh.init_mesh(
	//	node_cds,
	//	sizeof(node_cds) / (sizeof(node_cds[0])*3),
	//	elem_ids,
	//	sizeof(elem_ids) / (sizeof(elem_ids[0])*4)
	//	);

	//mesh.init_edges();
	//size_t edge_num = mesh.get_edge_num();
	//TetrahedronMesh::Edge *edges = mesh.get_edges();
	//std::cout << edge_num << "\n";
	//// edges[19]:
	//// 0, 1; 0, 3; 0, 4; 
	//// 1, 2; 1, 3; 1, 4; 1, 5; 1, 6;
	//// 2, 3; 2, 6;
	//// 3, 4; 3, 5; 3, 6; 3, 7;
	//// 4, 5; 4, 7; 4, 3;
	//// 5, 6; 5, 7; 6, 7;
	//for (size_t i = 0; i < edge_num; i++)
	//{
	//	TetrahedronMesh::Edge &edge = edges[i];
	//	std::cout << edge.n1 << " " << edge.n2 << "\n";
	//}

	// test hdf5
	mesh.load_mesh_from_hdf5("..\\..\\Asset\\brick_mesh.h5");
	std::cout << mesh.get_vol() << "\n";
	// print nodes
	size_t node_num = mesh.get_node_num();
	TetrahedronMesh::Node *nodes = mesh.get_nodes();
	std::cout << node_num << " nodes:\n";
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		TetrahedronMesh::Node &n = nodes[n_id];
		std::cout << n.id << " " << n.x << " " << n.y << " " << n.z << "\n";
	}
	// print elements
	size_t elem_num = mesh.get_elem_num();
	TetrahedronMesh::Element *elems = mesh.get_elems();
	std::cout << elem_num << " elements:\n";
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		TetrahedronMesh::Element &e = elems[e_id];
		std::cout << e.id << " " << e.n1 << " " << e.n2 << " "
				  << e.n3 << " " << e.n4 << "\n";
	}
	
	Cube bbox = mesh.get_bounding_box();
	std::cout << "bbox: " 
			  << bbox.xl << " " << bbox.xu << " "
			  << bbox.yl << " " << bbox.yu << "\n";
}
