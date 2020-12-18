#include "Geometry_pcp.h"

#include "TetrahedronMesh.h"

int TetrahedronMesh::init_search_grid(
	double _hx,
	double _hy,
	double _hz
	)
{
	return search_bg_grid.init(*this, _hx, _hy, _hz);
}

void TetrahedronMesh::translate_mesh(
	double dx,
	double dy,
	double dz
	) noexcept
{
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node& n = nodes[n_id];
		n.x += dx;
		n.y += dy;
		n.z += dz;
	}

	centre.x += dx;
	centre.y += dy;
	centre.z += dz;

	bounding_box.xl += dx;
	bounding_box.xu += dx;
	bounding_box.yl += dy;
	bounding_box.yu += dy;
	bounding_box.zl += dz;
	bounding_box.zu += dz;
}
