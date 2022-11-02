#include "Geometry_pcp.h"

#include <float.h>
#include "GeometryUtils.h"
#include "TetrahedronMesh.h"

int TetrahedronMesh::init_search_grid(double _hx, double _hy, double _hz)
{
	return search_bg_grid.init(*this, _hx, _hy, _hz);
}

void TetrahedronMesh::rotate_mesh(
	double dx_ang,
	double dy_ang,
	double dz_ang
	) noexcept
{
	trim_to_pi(dx_ang);
	trim_to_pi(dy_ang);
	trim_to_pi(dz_ang);
	const Vector3D pos_ang(dx_ang, dy_ang, dz_ang);
	Vector3D ix(1.0, 0.0, 0.0);
	Vector3D iy(0.0, 1.0, 0.0);
	Vector3D iz(0.0, 0.0, 1.0);
	rotate_axses_by_angle(pos_ang, ix, iy, iz);

	double n_x, n_y, n_z;
	bounding_box.xl = DBL_MAX;
	bounding_box.xu = -DBL_MAX;
	bounding_box.yl = DBL_MAX;
	bounding_box.yu = -DBL_MAX;
	bounding_box.zl = DBL_MAX;
	bounding_box.zu = -DBL_MAX;
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node& n = nodes[n_id];
		n_x = n.x;
		n_y = n.y;
		n_z = n.z;
		n.x = ix.x * n_x + ix.y * n_y + ix.z * n_z;
		n.y = iy.x * n_x + iy.y * n_y + iy.z * n_z;
		n.z = iz.x * n_x + iz.y * n_y + iz.z * n_z;
		if (bounding_box.xl > n.x)
			bounding_box.xl = n.x;
		if (bounding_box.xu < n.x)
			bounding_box.xu = n.x;
		if (bounding_box.yl > n.y)
			bounding_box.yl = n.y;
		if (bounding_box.yu < n.y)
			bounding_box.yu = n.y;
		if (bounding_box.zl > n.z)
			bounding_box.zl = n.z;
		if (bounding_box.zu < n.z)
			bounding_box.zu = n.z;
	}

	n_x = centre.x; n_y = centre.y; n_z = centre.z;
	centre.x = ix.x * n_x + ix.y * n_y + ix.z * n_z;
	centre.y = iy.x * n_x + iy.y * n_y + iy.z * n_z;
	centre.z = iz.x * n_x + iz.y * n_y + iz.z * n_z;
}

void TetrahedronMesh::translate_mesh(double dx, double dy, double dz) noexcept
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
