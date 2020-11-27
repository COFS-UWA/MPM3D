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
