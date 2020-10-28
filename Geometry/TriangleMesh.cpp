#include "Geometry_pcp.h"

#include "TriangleMesh.h"

int TriangleMesh::init_search_grid(double _hx, double _hy)
{
	return search_bg_grid.init(*this, _hx, _hy);
}
