#ifndef __Rigid_Body_hdf5_utilities_h__
#define __Rigid_Body_hdf5_utilities_h__

#include "RigidBody/RigidCircle.h"
#include "RigidBody/RigidRect.h"
#include "RigidBody/RigidTetrahedronMesh.h"
#include "ResultFile_hdf5.h"

namespace RigidBody_hdf5_utilities
{
	// rigid circle
	int output_rigid_circle_to_hdf5_file(
		RigidCircle& rc, ResultFile_hdf5& rf, hid_t rc_grp_id);
	int load_rigid_circle_from_hdf5_file(
		RigidCircle& rc, ResultFile_hdf5& rf, hid_t rc_grp_id);

	// rigid rect
	int output_rigid_rect_to_hdf5_file(
		RigidRect& rr, ResultFile_hdf5& rf, hid_t rr_grp_id);
	int load_rigid_rect_from_hdf5_file(
		RigidRect& rr, ResultFile_hdf5& rf, hid_t rr_grp_id);

	struct RigidTehMeshNodeData
	{
		unsigned long long id;
		double x, y, z;
		inline void from_node(const RigidTetrahedronMesh::Node& n)
		{
			id = n.id; x = n.x; y = n.y; z = n.z;
		}
		inline void to_node(RigidTetrahedronMesh::Node& n)
		{
			n.id = id; n.x = x; n.y = y; n.z = z;
		}
	};

	inline hid_t get_rigid_teh_mesh_node_dt_id()
	{
		hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(RigidTehMeshNodeData));
		H5Tinsert(res, "id", HOFFSET(RigidTehMeshNodeData, id), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "x", HOFFSET(RigidTehMeshNodeData, x), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "y", HOFFSET(RigidTehMeshNodeData, y), H5T_NATIVE_DOUBLE);
		H5Tinsert(res, "z", HOFFSET(RigidTehMeshNodeData, z), H5T_NATIVE_DOUBLE);
		return res;
	}

	struct RigidTehMeshElemData
	{
		unsigned long long id;
		unsigned long long n1, n2, n3, n4;
		inline void from_elem(const RigidTetrahedronMesh::Element& e)
		{
			id = e.id; n1 = e.n1; n2 = e.n2; n3 = e.n3; n4 = e.n4;
		}
		inline void to_elem(RigidTetrahedronMesh::Element& e)
		{
			e.id = id; e.n1 = n1; e.n2 = n2; e.n3 = n3; e.n4 = n4;
		}
	};

	inline hid_t get_rigid_teh_mesh_elem_dt_id()
	{
		hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(RigidTehMeshElemData));
		H5Tinsert(res, "id", HOFFSET(RigidTehMeshElemData, id), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "n1", HOFFSET(RigidTehMeshElemData, n1), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "n2", HOFFSET(RigidTehMeshElemData, n2), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "n3", HOFFSET(RigidTehMeshElemData, n3), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "n4", HOFFSET(RigidTehMeshElemData, n4), H5T_NATIVE_ULLONG);
		return res;
	}

	struct RigidTehMeshFaceData
	{
		unsigned long long id;
		unsigned long long n1, n2, n3;
		inline void from_face(const RigidTetrahedronMesh::Face& f)
		{
			id = f.id; n1 = f.n1; n2 = f.n2; n3 = f.n3;
		}
		inline void to_face(RigidTetrahedronMesh::Face& f)
		{
			f.id = id; f.n1 = n1; f.n2 = n2; f.n3 = n3;
		}
	};

	inline hid_t get_rigid_teh_mesh_face_dt_id()
	{
		hid_t res = H5Tcreate(H5T_COMPOUND, sizeof(RigidTehMeshFaceData));
		H5Tinsert(res, "id", HOFFSET(RigidTehMeshFaceData, id), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "n1", HOFFSET(RigidTehMeshFaceData, n1), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "n2", HOFFSET(RigidTehMeshFaceData, n2), H5T_NATIVE_ULLONG);
		H5Tinsert(res, "n3", HOFFSET(RigidTehMeshFaceData, n3), H5T_NATIVE_ULLONG);
		return res;
	}

	// rigid tetrahedron mesh
	// only output state info, no mesh info
	int output_rigid_tetrahedron_mesh_state_to_hdf5_file(
		RigidTetrahedronMesh& rtm, ResultFile_hdf5& rf, hid_t rtm_grp_id);
	int load_rigid_tetrahedron_mesh_state_from_hdf5_file(
		RigidTetrahedronMesh& rtm, ResultFile_hdf5& rf, hid_t rtm_grp_id);
	int output_rigid_tetrahedron_mesh_to_hdf5_file(
		RigidTetrahedronMesh& rtm, ResultFile_hdf5& rf, hid_t rtm_grp_id);
	int load_rigid_tetrahedron_mesh_from_hdf5_file(
		RigidTetrahedronMesh& rtm, ResultFile_hdf5& rf, hid_t rtm_grp_id);
}

#endif