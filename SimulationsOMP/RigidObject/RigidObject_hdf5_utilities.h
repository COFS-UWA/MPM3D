#ifndef __Rigid_Object_hdf5_utilities_h__
#define __Rigid_Object_hdf5_utilities_h__

#include "ResultFile_hdf5.h"
#include "RigidCylinder.h"
#include "RigidCone.h"

namespace RigidObject_hdf5_utilities
{
	int output_rigid_cylinder_to_hdf5_file(RigidCylinder &rc, ResultFile_hdf5 &rf, hid_t rc_grp_id);
	int load_rigid_cylinder_from_hdf5_file(RigidCylinder &rc, ResultFile_hdf5 &rf, hid_t rc_grp_id);
	int output_rigid_cone_to_hdf5_file(RigidCone &rc, ResultFile_hdf5& rf, hid_t rc_grp_id);
	int load_rigid_cone_from_hdf5_file(RigidCone &rc, ResultFile_hdf5& rf, hid_t rc_grp_id);
}

#endif