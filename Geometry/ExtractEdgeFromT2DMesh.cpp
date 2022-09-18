#include "Geometry_pcp.h"

#include "ExtractEdgeFromT2DMesh.h"

ExtractEdgeFromT2DMesh::ExtractEdgeFromT2DMesh() :
	boundary_edges(nullptr), boundary_edge_num(0) {}

ExtractEdgeFromT2DMesh::~ExtractEdgeFromT2DMesh() { clear(); }

inline static void sort_acc_node_index(size_t ns[2], size_t ids[2]) noexcept
{
	size_t tmp;
	if (ns[1] < ns[0])
	{
		tmp = ns[0];
		ns[0] = ns[1];
		ns[1] = tmp;
		ids[0] = 1;
		ids[1] = 0;
		return;
	}
	//
	ids[0] = 0;
	ids[1] = 1;
}

size_t ExtractEdgeFromT2DMesh::try_adding_tri_to_map(EdgeMap& edge_map, Edge& edge)
{
	size_t ids[2];
	sort_acc_node_index(reinterpret_cast<size_t *>(&edge), ids);
	// record the order of node index
	unsigned char content = ids[1] << 4 | ids[0] << 2;
	auto res = edge_map.emplace(edge, content);
	if (res.second == false)
	{
		// face already in map
		res.first->second |= 1;
		return 1;
	}
	return 0;
}
