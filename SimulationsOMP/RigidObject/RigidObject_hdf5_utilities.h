#ifndef __Rigid_Object_hdf5_utilities_h__
#define __Rigid_Object_hdf5_utilities_h__

#include "RigidCircle.h"
#include "RigidCube.h"
#include "RigidCylinder.h"
#include "RigidCone.h"
#include "RigidObjectByT3DMesh.h"
#include "ResultFile_hdf5.h"

namespace RigidObject_hdf5_utilities
{
	inline hid_t get_pt_to_tri_dist_dt_id()
	{
		hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(PointToTriangleDistance));
		// n1, n2, n3
		H5Tinsert(res, "n1_x", HOFFSET(PointToTriangleDistance, n1) + HOFFSET(Point3D, x), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "n1_y", HOFFSET(PointToTriangleDistance, n1) + HOFFSET(Point3D, y), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "n1_z", HOFFSET(PointToTriangleDistance, n1) + HOFFSET(Point3D, z), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "n2_x", HOFFSET(PointToTriangleDistance, n2) + HOFFSET(Point3D, x), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "n2_y", HOFFSET(PointToTriangleDistance, n2) + HOFFSET(Point3D, y), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "n2_z", HOFFSET(PointToTriangleDistance, n2) + HOFFSET(Point3D, z), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "n3_x", HOFFSET(PointToTriangleDistance, n3) + HOFFSET(Point3D, x), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "n3_y", HOFFSET(PointToTriangleDistance, n3) + HOFFSET(Point3D, y), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "n3_z", HOFFSET(PointToTriangleDistance, n3) + HOFFSET(Point3D, z), H5T_NATIVE_DOUBLE);
		// ix1, iy1, iz1
		H5Tinsert(res, "ix1_x", HOFFSET(PointToTriangleDistance, ix1) + HOFFSET(Vector3D, x), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "ix1_y", HOFFSET(PointToTriangleDistance, ix1) + HOFFSET(Vector3D, y), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "ix1_z", HOFFSET(PointToTriangleDistance, ix1) + HOFFSET(Vector3D, z), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "iy1_x", HOFFSET(PointToTriangleDistance, iy1) + HOFFSET(Vector3D, x), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "iy1_y", HOFFSET(PointToTriangleDistance, iy1) + HOFFSET(Vector3D, y), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "iy1_z", HOFFSET(PointToTriangleDistance, iy1) + HOFFSET(Vector3D, z), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "iz1_x", HOFFSET(PointToTriangleDistance, iz1) + HOFFSET(Vector3D, x), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "iz1_y", HOFFSET(PointToTriangleDistance, iz1) + HOFFSET(Vector3D, y), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "iz1_z", HOFFSET(PointToTriangleDistance, iz1) + HOFFSET(Vector3D, z), H5T_NATIVE_DOUBLE);
		// ix2, iy2, iz2
		H5Tinsert(res, "ix2_x", HOFFSET(PointToTriangleDistance, ix2) + HOFFSET(Vector3D, x), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "ix2_y", HOFFSET(PointToTriangleDistance, ix2) + HOFFSET(Vector3D, y), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "ix2_z", HOFFSET(PointToTriangleDistance, ix2) + HOFFSET(Vector3D, z), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "iy2_x", HOFFSET(PointToTriangleDistance, iy2) + HOFFSET(Vector3D, x), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "iy2_y", HOFFSET(PointToTriangleDistance, iy2) + HOFFSET(Vector3D, y), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "iy2_z", HOFFSET(PointToTriangleDistance, iy2) + HOFFSET(Vector3D, z), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "iz2_x", HOFFSET(PointToTriangleDistance, iz2) + HOFFSET(Vector3D, x), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "iz2_y", HOFFSET(PointToTriangleDistance, iz2) + HOFFSET(Vector3D, y), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "iz2_z", HOFFSET(PointToTriangleDistance, iz2) + HOFFSET(Vector3D, z), H5T_NATIVE_DOUBLE);
		// ix3, iy3, iz3
		H5Tinsert(res, "ix3_x", HOFFSET(PointToTriangleDistance, ix3) + HOFFSET(Vector3D, x), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "ix3_y", HOFFSET(PointToTriangleDistance, ix3) + HOFFSET(Vector3D, y), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "ix3_z", HOFFSET(PointToTriangleDistance, ix3) + HOFFSET(Vector3D, z), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "iy3_x", HOFFSET(PointToTriangleDistance, iy3) + HOFFSET(Vector3D, x), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "iy3_y", HOFFSET(PointToTriangleDistance, iy3) + HOFFSET(Vector3D, y), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "iy3_z", HOFFSET(PointToTriangleDistance, iy3) + HOFFSET(Vector3D, z), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "iz3_x", HOFFSET(PointToTriangleDistance, iz3) + HOFFSET(Vector3D, x), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "iz3_y", HOFFSET(PointToTriangleDistance, iz3) + HOFFSET(Vector3D, y), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "iz3_z", HOFFSET(PointToTriangleDistance, iz3) + HOFFSET(Vector3D, z), H5T_NATIVE_DOUBLE);
		// a1, a2, a3
		H5Tinsert(res, "a1", HOFFSET(PointToTriangleDistance, a1), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "a2", HOFFSET(PointToTriangleDistance, a2), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "a3", HOFFSET(PointToTriangleDistance, a3), H5T_NATIVE_DOUBLE);
		return res;
	}

	int output_rigid_circle_to_hdf5_file(RigidObject::RigidCircle& rc, ResultFile_hdf5& rf, hid_t rc_grp_id);
	int load_rigid_circle_from_hdf5_file(RigidObject::RigidCircle& rc, ResultFile_hdf5& rf, hid_t rc_grp_id);

	int output_rigid_cube_to_hdf5_file(RigidCube& rc, ResultFile_hdf5& rf, hid_t rc_grp_id);
	int load_rigid_cube_from_hdf5_file(RigidCube& rc, ResultFile_hdf5& rf, hid_t rc_grp_id);

	int output_rigid_cylinder_to_hdf5_file(RigidCylinder &rc, ResultFile_hdf5 &rf, hid_t rc_grp_id);
	int load_rigid_cylinder_from_hdf5_file(RigidCylinder &rc, ResultFile_hdf5 &rf, hid_t rc_grp_id);
	
	int output_rigid_cone_to_hdf5_file(RigidCone &rc, ResultFile_hdf5& rf, hid_t rc_grp_id);
	int load_rigid_cone_from_hdf5_file(RigidCone &rc, ResultFile_hdf5& rf, hid_t rc_grp_id);

	int output_rigid_object_by_3dmesh_to_hdf5_file(RigidObjectByT3DMesh &rb, ResultFile_hdf5 &rf, hid_t rb_grp_id);
	int output_rigid_object_by_3dmesh_state_to_hdf5_file(RigidObjectByT3DMesh& rb, ResultFile_hdf5& rf, hid_t rb_grp_id);
	int load_rigid_object_by_3dmesh_from_hdf5_file(RigidObjectByT3DMesh& rb, ResultFile_hdf5& rf, hid_t rb_grp_id);
	int load_rigid_object_by_3dmesh_state_from_hdf5_file(RigidObjectByT3DMesh& rb, ResultFile_hdf5& rf, hid_t rb_grp_id);
}

#endif