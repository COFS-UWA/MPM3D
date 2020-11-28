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
		const Force3D &cf = rc.get_cont_force();
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
		const Force3D& cf = rc.get_cont_force();
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

	int output_rigid_cube_to_hdf5_file(
		RigidCube& rc,
		ResultFile_hdf5& rf,
		hid_t rc_grp_id
		)
	{
		if (rc_grp_id < 0)
			return -1;

		rf.write_attribute(rc_grp_id, "hx", rc.get_hx());
		rf.write_attribute(rc_grp_id, "hy", rc.get_hy());
		rf.write_attribute(rc_grp_id, "hz", rc.get_hz());
		rf.write_attribute(rc_grp_id, "density", rc.get_density());
		const Point3D& cen = rc.get_centre();
		rf.write_attribute(rc_grp_id, "x", cen.x);
		rf.write_attribute(rc_grp_id, "y", cen.y);
		rf.write_attribute(rc_grp_id, "z", cen.z);
		const Vector3D& a = rc.get_acceleration();
		rf.write_attribute(rc_grp_id, "ax", a.x);
		rf.write_attribute(rc_grp_id, "ay", a.y);
		rf.write_attribute(rc_grp_id, "az", a.z);
		const Vector3D& v = rc.get_velocity();
		rf.write_attribute(rc_grp_id, "vx", v.x);
		rf.write_attribute(rc_grp_id, "vy", v.y);
		rf.write_attribute(rc_grp_id, "vz", v.z);
		const Force3D& cf = rc.get_cont_force();
		rf.write_attribute(rc_grp_id, "fx_cont", cf.fx);
		rf.write_attribute(rc_grp_id, "fy_cont", cf.fy);
		rf.write_attribute(rc_grp_id, "fz_cont", cf.fz);
		rf.write_attribute(rc_grp_id, "mx_cont", cf.mx);
		rf.write_attribute(rc_grp_id, "my_cont", cf.my);
		rf.write_attribute(rc_grp_id, "mz_cont", cf.mz);
		const Force3D& ef = rc.get_ext_force();
		rf.write_attribute(rc_grp_id, "fx_ext", ef.fx);
		rf.write_attribute(rc_grp_id, "fy_ext", ef.fy);
		rf.write_attribute(rc_grp_id, "fz_ext", ef.fz);
		rf.write_attribute(rc_grp_id, "mx_ext", ef.mx);
		rf.write_attribute(rc_grp_id, "my_ext", ef.my);
		rf.write_attribute(rc_grp_id, "mz_ext", ef.mz);
		return 0;
	}

	int load_rigid_cube_from_hdf5_file(
		RigidCube& rc,
		ResultFile_hdf5& rf,
		hid_t rc_grp_id
		)
	{
		if (rc_grp_id < 0)
			return -1;

		double hx, hy, hz, density;
		double x, y, z, ax, ay, az, vx, vy, vz;
		double fx_cont, fy_cont, fz_cont, mx_cont, my_cont, mz_cont;
		double fx_ext, fy_ext, fz_ext, mx_ext, my_ext, mz_ext;
		rf.read_attribute(rc_grp_id, "hx", hx);
		rf.read_attribute(rc_grp_id, "hy", hy);
		rf.read_attribute(rc_grp_id, "hz", hz);
		rf.read_attribute(rc_grp_id, "density", density);
		rf.read_attribute(rc_grp_id, "x", x);
		rf.read_attribute(rc_grp_id, "y", y);
		rf.read_attribute(rc_grp_id, "z", z);
		rf.read_attribute(rc_grp_id, "x", ax);
		rf.read_attribute(rc_grp_id, "y", ay);
		rf.read_attribute(rc_grp_id, "z", az);
		rf.read_attribute(rc_grp_id, "vx", vx);
		rf.read_attribute(rc_grp_id, "vy", vy);
		rf.read_attribute(rc_grp_id, "vz", vz);
		rf.read_attribute(rc_grp_id, "fx_cont", fx_cont);
		rf.read_attribute(rc_grp_id, "fy_cont", fy_cont);
		rf.read_attribute(rc_grp_id, "fz_cont", fz_cont);
		rf.read_attribute(rc_grp_id, "mx_cont", mx_cont);
		rf.read_attribute(rc_grp_id, "my_cont", my_cont);
		rf.read_attribute(rc_grp_id, "mz_cont", mz_cont);
		rf.read_attribute(rc_grp_id, "fx_ext", fx_ext);
		rf.read_attribute(rc_grp_id, "fy_ext", fy_ext);
		rf.read_attribute(rc_grp_id, "fz_ext", fz_ext);
		rf.read_attribute(rc_grp_id, "mx_ext", mx_ext);
		rf.read_attribute(rc_grp_id, "my_ext", my_ext);
		rf.read_attribute(rc_grp_id, "mz_exxt", mz_ext);
		rc.init(x, y, z, hx, hy, hz, density);
		rc.set_acceleration(ax, ay, az);
		rc.set_velocity(vx, vy, vz);
		rc.set_cont_force(fx_cont, fy_cont, fz_cont, mx_cont, my_cont, mz_cont);
		rc.set_ext_force(fx_ext, fy_ext, fz_ext, mx_ext, my_ext, mz_ext);
		return 0;
	}
}
