#include "SimulationsOMP_pcp.h"

#include "TetrahedronMesh.h"
#include "RigidObjectByT3DMesh.h"

RigidObjectByT3DMesh::RigidObjectByT3DMesh() :
	RigidObjectMotion3D(), RigidMeshT3D() {}

RigidObjectByT3DMesh::~RigidObjectByT3DMesh() {}

int RigidObjectByT3DMesh::init(
	double _density,
	const char* filename,
	double dx,
	double dy,
	double dz,
	double dx_ang,
	double dy_ang,
	double dz_ang
	)
{
	TetrahedronMesh tmesh;
	int res = tmesh.load_mesh_from_hdf5(filename);
	if (res)
		return res;

	const Point3D& cen = tmesh.get_centre();
	double mh_x = cen.x + dx;
	double mh_y = cen.y + dy;
	double mh_z = cen.z + dz;

	tmesh.translate_mesh(-cen.x, -cen.y, -cen.z);
	
	density = _density;
	// moment of intertia
	size_t elem_num = tmesh.get_elem_num();
	auto* nodes = tmesh.get_nodes();
	auto* elems = tmesh.get_elems();
	double elem_moi_data[6];
	double moi_data[6] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		auto &e = elems[e_id];
		auto &n1 = nodes[e.n1];
		auto &n2 = nodes[e.n2];
		auto &n3 = nodes[e.n3];
		auto &n4 = nodes[e.n4];
		cal_tetrahedron_moi(0.0, 0.0, 0.0,
			n1, n2, n3, n4, e.vol, elem_moi_data);
		moi_data[0] += density * elem_moi_data[0];
		moi_data[1] += density * elem_moi_data[1];
		moi_data[2] += density * elem_moi_data[2];
		moi_data[3] += density * elem_moi_data[3];
		moi_data[5] += density * elem_moi_data[5];
		moi_data[4] += density * elem_moi_data[4];
	}
	RigidObjectMotion3D::init(
		mh_x, mh_y, mh_z,
		tmesh.get_vol() * density,
		moi_data
		);
	RigidObjectMotion3D::set_angle(dx_ang, dy_ang, dz_ang);
	
	//RigidMeshT3D::init_from_mesh(tmesh);
	return 0;
}

bool RigidObjectByT3DMesh::detect_collision_with_point(
	double p_x,
	double p_y,
	double p_z,
	double p_r,
	double& dist,
	Vector3D& lnorm,
	Point3D& lcontpos
	) const noexcept
{
	Point3D gp(p_x, p_y, p_z);
	get_local_point(gp, lcontpos);
	if (RigidMeshT3D::detect_collision_with_point(lcontpos, p_r, dist, lnorm))
		return true;
	return false;
}
