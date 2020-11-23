#include "SimulationsOMP_pcp.h"

#include "RigidObject_hdf5_utilities.h"

namespace RigidObject_hdf5_utilities
{
	int output_rigid_cylinder_to_hdf5_file(
		RigidCylinder& rc,
		ResultFile_hdf5& rf,
		hid_t rc_grp_id
		)
	{
		if (rc_grp_id < 0)
			return -1;

		rf.write_attribute(rc_grp_id, "h", rc.get_h());
		rf.write_attribute(rc_grp_id, "r", rc.get_r());
		const Point3D &cen = rc.get_centre();
		rf.write_attribute(rc_grp_id, "x", cen.x);
		rf.write_attribute(rc_grp_id, "y", cen.y);
		rf.write_attribute(rc_grp_id, "z", cen.z);
		const Vector3D& v = rc.get_velocity();
		rf.write_attribute(rc_grp_id, "vx", v.x);
		rf.write_attribute(rc_grp_id, "vy", v.y);
		rf.write_attribute(rc_grp_id, "vz", v.z);
		const ContactForce3D &cf = rc.get_cont_force();
		rf.write_attribute(rc_grp_id, "fx", cf.fx);
		rf.write_attribute(rc_grp_id, "fy", cf.fy);
		rf.write_attribute(rc_grp_id, "fz", cf.fz);
		rf.write_attribute(rc_grp_id, "mx", cf.mx);
		rf.write_attribute(rc_grp_id, "my", cf.my);
		rf.write_attribute(rc_grp_id, "mz", cf.mz);
		return 0;
	}
	
	int load_rigid_cylinder_from_hdf5_file(
		RigidCylinder& rc,
		ResultFile_hdf5& rf,
		hid_t rc_grp_id
		)
	{
		if (rc_grp_id < 0)
			return -1;
		
		double h, r;
		double x, y, z, vx, vy, vz;
		double fx, fy, fz, mx, my, mz;
		rf.read_attribute(rc_grp_id, "h", h);
		rf.read_attribute(rc_grp_id, "r", r);
		rf.read_attribute(rc_grp_id, "x", x);
		rf.read_attribute(rc_grp_id, "y", y);
		rf.read_attribute(rc_grp_id, "z", z);
		rf.read_attribute(rc_grp_id, "vx", vx);
		rf.read_attribute(rc_grp_id, "vy", vy);
		rf.read_attribute(rc_grp_id, "vz", vz);
		rf.read_attribute(rc_grp_id, "fx", fx);
		rf.read_attribute(rc_grp_id, "fy", fy);
		rf.read_attribute(rc_grp_id, "fz", fz);
		rf.read_attribute(rc_grp_id, "mx", mx);
		rf.read_attribute(rc_grp_id, "my", my);
		rf.read_attribute(rc_grp_id, "mz", mz);
		rc.init(x, y, z, h, r);
		rc.set_vbc(vx, vy, vz);
		rc.set_cont_force(fx, fy, fz, mx, my, mz);
		return 0;
	}

	int output_rigid_cone_to_hdf5_file(
		RigidCone& rc,
		ResultFile_hdf5& rf,
		hid_t rc_grp_id
		)
	{
		if (rc_grp_id < 0)
			return -1;

		rf.write_attribute(rc_grp_id, "r", rc.get_r());
		rf.write_attribute(rc_grp_id, "h_tip", rc.get_h_tip());
		rf.write_attribute(rc_grp_id, "h_shaft", rc.get_h_shaft());
		const Point3D& cen = rc.get_centre();
		rf.write_attribute(rc_grp_id, "x", cen.x);
		rf.write_attribute(rc_grp_id, "y", cen.y);
		rf.write_attribute(rc_grp_id, "z", cen.z);
		const Vector3D& v = rc.get_velocity();
		rf.write_attribute(rc_grp_id, "vx", v.x);
		rf.write_attribute(rc_grp_id, "vy", v.y);
		rf.write_attribute(rc_grp_id, "vz", v.z);
		const ContactForce3D& cf = rc.get_cont_force();
		rf.write_attribute(rc_grp_id, "fx", cf.fx);
		rf.write_attribute(rc_grp_id, "fy", cf.fy);
		rf.write_attribute(rc_grp_id, "fz", cf.fz);
		rf.write_attribute(rc_grp_id, "mx", cf.mx);
		rf.write_attribute(rc_grp_id, "my", cf.my);
		rf.write_attribute(rc_grp_id, "mz", cf.mz);
		return 0;
	}

	int load_rigid_cone_from_hdf5_file(
		RigidCone& rc,
		ResultFile_hdf5& rf,
		hid_t rc_grp_id
		)
	{
		if (rc_grp_id < 0)
			return -1;

		double r, h_tip, h_shaft;
		double x, y, z, vx, vy, vz;
		double fx, fy, fz, mx, my, mz;
		rf.read_attribute(rc_grp_id, "r", r);
		rf.read_attribute(rc_grp_id, "h_tip", h_tip);
		rf.read_attribute(rc_grp_id, "h_shaft", h_shaft);
		rf.read_attribute(rc_grp_id, "x", x);
		rf.read_attribute(rc_grp_id, "y", y);
		rf.read_attribute(rc_grp_id, "z", z);
		rf.read_attribute(rc_grp_id, "vx", vx);
		rf.read_attribute(rc_grp_id, "vy", vy);
		rf.read_attribute(rc_grp_id, "vz", vz);
		rf.read_attribute(rc_grp_id, "fx", fx);
		rf.read_attribute(rc_grp_id, "fy", fy);
		rf.read_attribute(rc_grp_id, "fz", fz);
		rf.read_attribute(rc_grp_id, "mx", mx);
		rf.read_attribute(rc_grp_id, "my", my);
		rf.read_attribute(rc_grp_id, "mz", mz);
		rc.init(x, y, z, r, h_tip, h_shaft);
		rc.set_vbc(vx, vy, vz);
		rc.set_cont_force(fx, fy, fz, mx, my, mz);
		return 0;
	}
}
