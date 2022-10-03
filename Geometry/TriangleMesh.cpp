#include "Geometry_pcp.h"

#include "GeometryUtils.h"
#include "TriangleMesh.h"

int TriangleMesh::init_search_grid(double _hx, double _hy)
{
	return search_bg_grid.init(*this, _hx, _hy);
}

void TriangleMesh::rotate_mesh(double dang) noexcept
{
	trim_to_pi(dang);
	Vector2D ix, iy;
	rotate_axses_by_angle(dang, ix, iy);

	double n_x, n_y;
	bounding_box.xl = DBL_MAX;
	bounding_box.xu = -DBL_MAX;
	bounding_box.yl = DBL_MAX;
	bounding_box.yu = -DBL_MAX;
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node& n = nodes[n_id];
		n_x = n.x;
		n_y = n.y;
		n.x = ix.x * n_x + ix.y * n_y;
		n.y = iy.x * n_x + iy.y * n_y;
		if (bounding_box.xl > n.x)
			bounding_box.xl = n.x;
		if (bounding_box.xu < n.x)
			bounding_box.xu = n.x;
		if (bounding_box.yl > n.y)
			bounding_box.yl = n.y;
		if (bounding_box.yu < n.y)
			bounding_box.yu = n.y;
	}

	n_x = centre.x;
	n_y = centre.y;
	centre.x = ix.x * n_x + ix.y * n_y;
	centre.y = iy.x * n_x + iy.y * n_y;
}

void TriangleMesh::translate_mesh(double dx, double dy) noexcept
{
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node& n = nodes[n_id];
		n.x += dx;
		n.y += dy;
	}
	centre.x += dx;
	centre.y += dy;
	bounding_box.xl += dx;
	bounding_box.xu += dx;
	bounding_box.yl += dy;
	bounding_box.yu += dy;
}

