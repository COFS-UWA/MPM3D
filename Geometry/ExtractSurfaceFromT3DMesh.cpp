#include "Geometry_pcp.h"

#include "ExtractSurfaceFromT3DMesh.h"

ExtractSurfaceFromT3DMesh::ExtractSurfaceFromT3DMesh() :
	faces(nullptr), face_num(0) {}

ExtractSurfaceFromT3DMesh::~ExtractSurfaceFromT3DMesh() { clear(); }

inline static void sort_acc_node_index(
	size_t* const ns, size_t* const ids) noexcept
{
	size_t tmp;
	size_t min_id = 0;
	if (ns[1] < ns[0])
		min_id = 1;
	if (ns[2] < ns[min_id])
		min_id = 2;
	if (min_id)
	{
		tmp = ns[0];
		ns[0] = ns[min_id];
		ns[min_id] = tmp;
		tmp = ids[0];
		ids[0] = ids[min_id];
		ids[min_id] = tmp;
	}
	if (ns[2] < ns[1])
	{
		tmp = ns[1];
		ns[1] = ns[2];
		ns[2] = tmp;
		tmp = ids[1];
		ids[1] = ids[2];
		ids[2] = tmp;
	}
}

size_t ExtractSurfaceFromT3DMesh::try_adding_teh_to_map(
	FaceMap& face_map,
	Face& face)
{
	size_t ids[3] = { 0, 1, 2 };
	sort_acc_node_index(reinterpret_cast<size_t *>(&face), ids);
	unsigned char content = ids[2] << 6 | ids[1] << 4 | ids[0] << 2;
	auto res = face_map.emplace(face, content);
	if (res.second == false)
	{
		// face already in map
		res.first->second |= 1;
		return 1;
	}
	return 0;
}
