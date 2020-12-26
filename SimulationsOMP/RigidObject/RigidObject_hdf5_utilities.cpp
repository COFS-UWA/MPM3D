#include "SimulationsOMP_pcp.h"

#include "RigidObject_hdf5_utilities.h"

namespace RigidObject_hdf5_utilities
{
	// rigid circle
	int output_rigid_circle_to_hdf5_file(
		RigidObject::RigidCircle& rc,
		ResultFile_hdf5& rf,
		hid_t rc_grp_id
		)
	{
		if (rc_grp_id < 0)
			return -1;

		rf.write_attribute(rc_grp_id, "radius", rc.get_radius());
		rf.write_attribute(rc_grp_id, "density", rc.get_density());
		rf.write_attribute(rc_grp_id, "ax", rc.get_ax());
		rf.write_attribute(rc_grp_id, "ay", rc.get_ay());
		rf.write_attribute(rc_grp_id, "a_angle", rc.get_a_ang());
		rf.write_attribute(rc_grp_id, "vx", rc.get_vx());
		rf.write_attribute(rc_grp_id, "vy", rc.get_vy());
		rf.write_attribute(rc_grp_id, "v_angle", rc.get_v_ang());
		rf.write_attribute(rc_grp_id, "x", rc.get_x());
		rf.write_attribute(rc_grp_id, "y", rc.get_y());
		rf.write_attribute(rc_grp_id, "angle", rc.get_ang());
		const Force2D& f_cont = rc.get_cont_force();
		rf.write_attribute(rc_grp_id, "fx_cont", f_cont.fx);
		rf.write_attribute(rc_grp_id, "fy_cont", f_cont.fy);
		rf.write_attribute(rc_grp_id, "m_cont", f_cont.m);
		const Force2D& f_ext = rc.get_ext_force();
		rf.write_attribute(rc_grp_id, "fx_ext", f_ext.fx);
		rf.write_attribute(rc_grp_id, "fy_ext", f_ext.fy);
		rf.write_attribute(rc_grp_id, "m_ext", f_ext.m);

		if (rc.has_ax_bc())
			rf.write_attribute(rc_grp_id, "ax_bc", rc.get_ax_bc());
		if (rc.has_ay_bc())
			rf.write_attribute(rc_grp_id, "ay_bc", rc.get_ay_bc());
		if (rc.has_a_ang_bc())
			rf.write_attribute(rc_grp_id, "a_ang_bc", rc.get_a_ang_bc());
		if (rc.has_vx_bc())
			rf.write_attribute(rc_grp_id, "vx_bc", rc.get_vx_bc());
		if (rc.has_vy_bc())
			rf.write_attribute(rc_grp_id, "vy_bc", rc.get_vy_bc());
		if (rc.has_v_ang_bc())
			rf.write_attribute(rc_grp_id, "v_ang_bc", rc.get_v_ang_bc());

		return 0;
	}

	int load_rigid_circle_from_hdf5_file(
		RigidObject::RigidCircle& rc,
		ResultFile_hdf5& rf,
		hid_t rc_grp_id
		)
	{
		if (rc_grp_id < 0)
			return -1;

		double rc_radius, rc_density;
		double rc_ax, rc_ay, rc_a_ang;
		double rc_vx, rc_vy, rc_v_ang;
		double rc_x, rc_y, rc_ang;
		rf.read_attribute(rc_grp_id, "radius", rc_radius);
		rf.read_attribute(rc_grp_id, "density", rc_density);
		rf.read_attribute(rc_grp_id, "ax", rc_ax);
		rf.read_attribute(rc_grp_id, "ay", rc_ay);
		rf.read_attribute(rc_grp_id, "a_angle", rc_a_ang);
		rf.read_attribute(rc_grp_id, "vx", rc_vx);
		rf.read_attribute(rc_grp_id, "vy", rc_vy);
		rf.read_attribute(rc_grp_id, "v_angle", rc_v_ang);
		rf.read_attribute(rc_grp_id, "x", rc_x);
		rf.read_attribute(rc_grp_id, "y", rc_y);
		rf.read_attribute(rc_grp_id, "angle", rc_ang);

		Force2D f_cont, f_ext;
		rf.read_attribute(rc_grp_id, "fx_cont", f_cont.fx);
		rf.read_attribute(rc_grp_id, "fy_cont", f_cont.fy);
		rf.read_attribute(rc_grp_id, "m_cont", f_cont.m);
		rf.read_attribute(rc_grp_id, "fx_ext", f_ext.fx);
		rf.read_attribute(rc_grp_id, "fy_ext", f_ext.fy);
		rf.read_attribute(rc_grp_id, "m_ext", f_ext.m);

		rc.set_init_state(
			rc_radius, rc_density,
			rc_ax, rc_ay, rc_a_ang,
			rc_vx, rc_vy, rc_v_ang,
			rc_x, rc_y, rc_ang,
			f_cont, f_ext
			);

		// boundary conditions
		double bc_value;
		if (rf.has_attribute(rc_grp_id, "ax_bc"))
		{
			rf.read_attribute(rc_grp_id, "ax_bc", bc_value);
			rc.set_ax_bc(bc_value);
		}
		if (rf.has_attribute(rc_grp_id, "ay_bc"))
		{
			rf.read_attribute(rc_grp_id, "ay_bc", bc_value);
			rc.set_ay_bc(bc_value);
		}
		if (rf.has_attribute(rc_grp_id, "a_ang_bc"))
		{
			rf.read_attribute(rc_grp_id, "a_ang_bc", bc_value);
			rc.set_a_ang_bc(bc_value);
		}
		if (rf.has_attribute(rc_grp_id, "vx_bc"))
		{
			rf.read_attribute(rc_grp_id, "vx_bc", bc_value);
			rc.set_vx_bc(bc_value);
		}
		if (rf.has_attribute(rc_grp_id, "vy_bc"))
		{
			rf.read_attribute(rc_grp_id, "vy_bc", bc_value);
			rc.set_vy_bc(bc_value);
		}
		if (rf.has_attribute(rc_grp_id, "v_ang_bc"))
		{
			rf.read_attribute(rc_grp_id, "v_ang_bc", bc_value);
			rc.set_v_ang_bc(bc_value);
		}

		return 0;
	}
	
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

	int output_rigid_object_by_3dmesh_to_hdf5_file(
		RigidObjectByT3DMesh& rb,
		ResultFile_hdf5& rf,
		hid_t rb_grp_id
		)
	{
		if (rb_grp_id < 0)
			return -1;
		
		rf.write_attribute(rb_grp_id, "density", rb.get_density());

		// rigid object motion
		rf.write_attribute(rb_grp_id, "m", rb.get_m());
		rf.write_attribute(rb_grp_id, "moi", 6, rb.get_moi());
		rf.write_attribute(rb_grp_id, "inv_moi", 6, rb.get_inv_moi());
		rf.write_attribute(rb_grp_id, "T_mat", 9, rb.get_T_mat());
		const Vector3D& a = rb.get_a();
		rf.write_attribute(rb_grp_id, "ax", a.x);
		rf.write_attribute(rb_grp_id, "ay", a.y);
		rf.write_attribute(rb_grp_id, "az", a.z);
		const Vector3D& a_ang = rb.get_a_ang();
		rf.write_attribute(rb_grp_id, "ax_ang", a_ang.x);
		rf.write_attribute(rb_grp_id, "ay_ang", a_ang.y);
		rf.write_attribute(rb_grp_id, "az_ang", a_ang.z);
		const Vector3D& v = rb.get_v();
		rf.write_attribute(rb_grp_id, "vx", v.x);
		rf.write_attribute(rb_grp_id, "vy", v.y);
		rf.write_attribute(rb_grp_id, "vz", v.z);
		const Vector3D& v_ang = rb.get_v_ang();
		rf.write_attribute(rb_grp_id, "vx_ang", v_ang.x);
		rf.write_attribute(rb_grp_id, "vy_ang", v_ang.y);
		rf.write_attribute(rb_grp_id, "vz_ang", v_ang.z);
		const Point3D &pos = rb.get_pos();
		rf.write_attribute(rb_grp_id, "x", pos.x);
		rf.write_attribute(rb_grp_id, "y", pos.y);
		rf.write_attribute(rb_grp_id, "z", pos.z);
		const Vector3D& pos_ang = rb.get_pos_ang();
		rf.write_attribute(rb_grp_id, "x_ang", pos_ang.x);
		rf.write_attribute(rb_grp_id, "y_ang", pos_ang.y);
		rf.write_attribute(rb_grp_id, "z_ang", pos_ang.z);
		const Force3D& cf = rb.get_force_contact();
		rf.write_attribute(rb_grp_id, "fx_cont", cf.fx);
		rf.write_attribute(rb_grp_id, "fy_cont", cf.fy);
		rf.write_attribute(rb_grp_id, "fz_cont", cf.fz);
		rf.write_attribute(rb_grp_id, "mx_cont", cf.mx);
		rf.write_attribute(rb_grp_id, "my_cont", cf.my);
		rf.write_attribute(rb_grp_id, "mz_cont", cf.mz);
		const Force3D& ef = rb.get_force_ext();
		rf.write_attribute(rb_grp_id, "fx_ext", ef.fx);
		rf.write_attribute(rb_grp_id, "fy_ext", ef.fy);
		rf.write_attribute(rb_grp_id, "fz_ext", ef.fz);
		rf.write_attribute(rb_grp_id, "mx_ext", ef.mx);
		rf.write_attribute(rb_grp_id, "my_ext", ef.my);
		rf.write_attribute(rb_grp_id, "mz_ext", ef.mz);
		if (rb.has_ax_bc())
			rf.write_attribute(rb_grp_id, "ax_bc", rb.get_ax_bc());
		if (rb.has_ay_bc())
			rf.write_attribute(rb_grp_id, "ay_bc", rb.get_ay_bc());
		if (rb.has_az_bc())
			rf.write_attribute(rb_grp_id, "az_bc", rb.get_az_bc());
		if (rb.has_ax_ang_bc())
			rf.write_attribute(rb_grp_id, "ax_ang_bc", rb.get_ax_ang_bc());
		if (rb.has_ay_ang_bc())
			rf.write_attribute(rb_grp_id, "ay_ang_bc", rb.get_ay_ang_bc());
		if (rb.has_az_ang_bc())
			rf.write_attribute(rb_grp_id, "az_ang_bc", rb.get_az_ang_bc());
		if (rb.has_vx_bc())
			rf.write_attribute(rb_grp_id, "vx_bc", rb.get_vx_bc());
		if (rb.has_vy_bc())
			rf.write_attribute(rb_grp_id, "vy_bc", rb.get_vy_bc());
		if (rb.has_vz_bc())
			rf.write_attribute(rb_grp_id, "vz_bc", rb.get_vz_bc());
		if (rb.has_vx_ang_bc())
			rf.write_attribute(rb_grp_id, "vx_ang_bc", rb.get_vx_ang_bc());
		if (rb.has_vy_ang_bc())
			rf.write_attribute(rb_grp_id, "vy_ang_bc", rb.get_vy_ang_bc());
		if (rb.has_vz_ang_bc())
			rf.write_attribute(rb_grp_id, "vz_ang_bc", rb.get_vz_ang_bc());

		// rigid mesh
		rf.write_attribute(rb_grp_id, "grid_xl", rb.get_grid_xl());
		rf.write_attribute(rb_grp_id, "grid_yl", rb.get_grid_yl());
		rf.write_attribute(rb_grp_id, "grid_zl", rb.get_grid_zl());
		rf.write_attribute(rb_grp_id, "grid_hx", rb.get_grid_hx());
		rf.write_attribute(rb_grp_id, "grid_hy", rb.get_grid_hy());
		rf.write_attribute(rb_grp_id, "grid_hz", rb.get_grid_hz());
		rf.write_attribute(rb_grp_id, "grid_x_num", rb.get_grid_x_num());
		rf.write_attribute(rb_grp_id, "grid_y_num", rb.get_grid_y_num());
		rf.write_attribute(rb_grp_id, "grid_z_num", rb.get_grid_z_num());
		rf.write_attribute(rb_grp_id, "face_grid_list_len", rb.get_face_in_grid_list_len());
		rf.write_attribute(rb_grp_id, "face_num", rb.get_face_num());

		rf.write_dataset(rb_grp_id,
			"grid_pos_type",
			rb.get_grid_num(),
			(const unsigned char *)rb.get_grid_pos_type());
		rf.write_dataset(rb_grp_id,
			"face_in_grid_range",
			rb.get_grid_num() + 1,
			(unsigned long long*)rb.get_face_in_grid_range());
		rf.write_dataset(rb_grp_id,
			"face_in_grid_list",
			rb.get_face_in_grid_list_len(),
			(unsigned long long*)rb.get_face_in_grid_list());
		
		hid_t pt2tri_dt = get_pt_to_tri_dist_dt_id();
		rf.write_dataset(rb_grp_id,
			"pt_tri_dist",
			rb.get_face_num(),
			rb.get_pt_tri_dist(),
			pt2tri_dt
			);
		H5Tclose(pt2tri_dt);

		return 0;
	}

	int output_rigid_object_by_3dmesh_state_to_hdf5_file(
		RigidObjectByT3DMesh& rb,
		ResultFile_hdf5& rf,
		hid_t rb_grp_id
		)
	{
		if (rb_grp_id < 0)
			return -1;

		const Vector3D& a = rb.get_a();
		rf.write_attribute(rb_grp_id, "ax", a.x);
		rf.write_attribute(rb_grp_id, "ay", a.y);
		rf.write_attribute(rb_grp_id, "az", a.z);
		const Vector3D& a_ang = rb.get_a_ang();
		rf.write_attribute(rb_grp_id, "ax_ang", a_ang.x);
		rf.write_attribute(rb_grp_id, "ay_ang", a_ang.y);
		rf.write_attribute(rb_grp_id, "az_ang", a_ang.z);
		const Vector3D& v = rb.get_v();
		rf.write_attribute(rb_grp_id, "vx", v.x);
		rf.write_attribute(rb_grp_id, "vy", v.y);
		rf.write_attribute(rb_grp_id, "vz", v.z);
		const Vector3D& v_ang = rb.get_v_ang();
		rf.write_attribute(rb_grp_id, "vx_ang", v_ang.x);
		rf.write_attribute(rb_grp_id, "vy_ang", v_ang.y);
		rf.write_attribute(rb_grp_id, "vz_ang", v_ang.z);
		const Point3D& pos = rb.get_pos();
		rf.write_attribute(rb_grp_id, "x", pos.x);
		rf.write_attribute(rb_grp_id, "y", pos.y);
		rf.write_attribute(rb_grp_id, "z", pos.z);
		const Vector3D& pos_ang = rb.get_pos_ang();
		rf.write_attribute(rb_grp_id, "x_ang", pos_ang.x);
		rf.write_attribute(rb_grp_id, "y_ang", pos_ang.y);
		rf.write_attribute(rb_grp_id, "z_ang", pos_ang.z);
		const Force3D& cf = rb.get_force_contact();
		rf.write_attribute(rb_grp_id, "fx_cont", cf.fx);
		rf.write_attribute(rb_grp_id, "fy_cont", cf.fy);
		rf.write_attribute(rb_grp_id, "fz_cont", cf.fz);
		rf.write_attribute(rb_grp_id, "mx_cont", cf.mx);
		rf.write_attribute(rb_grp_id, "my_cont", cf.my);
		rf.write_attribute(rb_grp_id, "mz_cont", cf.mz);
		const Force3D& ef = rb.get_force_ext();
		rf.write_attribute(rb_grp_id, "fx_ext", ef.fx);
		rf.write_attribute(rb_grp_id, "fy_ext", ef.fy);
		rf.write_attribute(rb_grp_id, "fz_ext", ef.fz);
		rf.write_attribute(rb_grp_id, "mx_ext", ef.mx);
		rf.write_attribute(rb_grp_id, "my_ext", ef.my);
		rf.write_attribute(rb_grp_id, "mz_ext", ef.mz);
		if (rb.has_ax_bc())
			rf.write_attribute(rb_grp_id, "ax_bc", rb.get_ax_bc());
		if (rb.has_ay_bc())
			rf.write_attribute(rb_grp_id, "ay_bc", rb.get_ay_bc());
		if (rb.has_az_bc())
			rf.write_attribute(rb_grp_id, "az_bc", rb.get_az_bc());
		if (rb.has_ax_ang_bc())
			rf.write_attribute(rb_grp_id, "ax_ang_bc", rb.get_ax_ang_bc());
		if (rb.has_ay_ang_bc())
			rf.write_attribute(rb_grp_id, "ay_ang_bc", rb.get_ay_ang_bc());
		if (rb.has_az_ang_bc())
			rf.write_attribute(rb_grp_id, "az_ang_bc", rb.get_az_ang_bc());
		if (rb.has_vx_bc())
			rf.write_attribute(rb_grp_id, "vx_bc", rb.get_vx_bc());
		if (rb.has_vy_bc())
			rf.write_attribute(rb_grp_id, "vy_bc", rb.get_vy_bc());
		if (rb.has_vz_bc())
			rf.write_attribute(rb_grp_id, "vz_bc", rb.get_vz_bc());
		if (rb.has_vx_ang_bc())
			rf.write_attribute(rb_grp_id, "vx_ang_bc", rb.get_vx_ang_bc());
		if (rb.has_vy_ang_bc())
			rf.write_attribute(rb_grp_id, "vy_ang_bc", rb.get_vy_ang_bc());
		if (rb.has_vz_ang_bc())
			rf.write_attribute(rb_grp_id, "vz_ang_bc", rb.get_vz_ang_bc());

		return 0;
	}
	
	int load_rigid_object_by_3dmesh_from_hdf5_file(RigidObjectByT3DMesh& rb, ResultFile_hdf5& rf, hid_t rb_grp_id)
	{
		// rigid object motion
		double density, m, moi[6], inv_moi[6], T_mat[9];
		rf.read_attribute(rb_grp_id, "density", density);
		rf.read_attribute(rb_grp_id, "m", m);
		rf.read_attribute(rb_grp_id, "moi", 6, moi);
		rf.read_attribute(rb_grp_id, "inv_moi", 6, inv_moi);
		rf.read_attribute(rb_grp_id, "T_mat", 9, T_mat);

		Vector3D acc, acc_ang;
		rf.read_attribute(rb_grp_id, "ax", acc.x);
		rf.read_attribute(rb_grp_id, "ay", acc.y);
		rf.read_attribute(rb_grp_id, "az", acc.z);
		rf.read_attribute(rb_grp_id, "ax_ang", acc_ang.x);
		rf.read_attribute(rb_grp_id, "ay_ang", acc_ang.y);
		rf.read_attribute(rb_grp_id, "az_ang", acc_ang.z);
		Vector3D vec, vec_ang;
		rf.read_attribute(rb_grp_id, "vx", vec.x);
		rf.read_attribute(rb_grp_id, "vy", vec.y);
		rf.read_attribute(rb_grp_id, "vz", vec.z);
		rf.read_attribute(rb_grp_id, "vx_ang", vec_ang.x);
		rf.read_attribute(rb_grp_id, "vy_ang", vec_ang.y);
		rf.read_attribute(rb_grp_id, "vz_ang", vec_ang.z);
		Point3D pos;
		Vector3D pos_ang;
		rf.read_attribute(rb_grp_id, "x", pos.x);
		rf.read_attribute(rb_grp_id, "y", pos.y);
		rf.read_attribute(rb_grp_id, "z", pos.z);
		rf.read_attribute(rb_grp_id, "x_ang", pos_ang.x);
		rf.read_attribute(rb_grp_id, "y_ang", pos_ang.y);
		rf.read_attribute(rb_grp_id, "z_ang", pos_ang.z);
		Force3D force_cont, force_ext;
		rf.read_attribute(rb_grp_id, "fx_cont", force_cont.fx);
		rf.read_attribute(rb_grp_id, "fy_cont", force_cont.fy);
		rf.read_attribute(rb_grp_id, "fz_cont", force_cont.fz);
		rf.read_attribute(rb_grp_id, "mx_cont", force_cont.mx);
		rf.read_attribute(rb_grp_id, "my_cont", force_cont.my);
		rf.read_attribute(rb_grp_id, "mz_cont", force_cont.mz);
		rf.read_attribute(rb_grp_id, "fx_ext", force_ext.fx);
		rf.read_attribute(rb_grp_id, "fy_ext", force_ext.fy);
		rf.read_attribute(rb_grp_id, "fz_ext", force_ext.fz);
		rf.read_attribute(rb_grp_id, "mx_ext", force_ext.mx);
		rf.read_attribute(rb_grp_id, "my_ext", force_ext.my);
		rf.read_attribute(rb_grp_id, "mz_ext", force_ext.mz);
		double bc_value;
		if (rf.has_attribute(rb_grp_id, "ax_bc"))
		{
			rf.read_attribute(rb_grp_id, "ax_bc", bc_value);
			rb.set_ax_bc(bc_value);
		}
		if (rf.has_attribute(rb_grp_id, "ay_bc"))
		{
			rf.read_attribute(rb_grp_id, "ay_bc", bc_value);
			rb.set_ay_bc(bc_value);
		}
		if (rf.has_attribute(rb_grp_id, "az_bc"))
		{
			rf.read_attribute(rb_grp_id, "az_bc", bc_value);
			rb.set_az_bc(bc_value);
		}
		if (rf.has_attribute(rb_grp_id, "ax_ang_bc"))
		{
			rf.read_attribute(rb_grp_id, "ax_ang_bc", bc_value);
			rb.set_ax_ang_bc(bc_value);
		}
		if (rf.has_attribute(rb_grp_id, "ay_ang_bc"))
		{
			rf.read_attribute(rb_grp_id, "ay_ang_bc", bc_value);
			rb.set_ay_ang_bc(bc_value);
		}
		if (rf.has_attribute(rb_grp_id, "az_ang_bc"))
		{
			rf.read_attribute(rb_grp_id, "az_ang_bc", bc_value);
			rb.set_az_ang_bc(bc_value);
		}
		if (rf.has_attribute(rb_grp_id, "vx_bc"))
		{
			rf.read_attribute(rb_grp_id, "vx_bc", bc_value);
			rb.set_vx_bc(bc_value);
		}
		if (rf.has_attribute(rb_grp_id, "vy_bc"))
		{
			rf.read_attribute(rb_grp_id, "vy_bc", bc_value);
			rb.set_vy_bc(bc_value);
		}
		if (rf.has_attribute(rb_grp_id, "vz_bc"))
		{
			rf.read_attribute(rb_grp_id, "vz_bc", bc_value);
			rb.set_vz_bc(bc_value);
		}
		if (rf.has_attribute(rb_grp_id, "vx_ang_bc"))
		{
			rf.read_attribute(rb_grp_id, "vx_ang_bc", bc_value);
			rb.set_vx_ang_bc(bc_value);
		}
		if (rf.has_attribute(rb_grp_id, "vy_ang_bc"))
		{
			rf.read_attribute(rb_grp_id, "vy_ang_bc", bc_value);
			rb.set_vy_ang_bc(bc_value);
		}
		if (rf.has_attribute(rb_grp_id, "vz_ang_bc"))
		{
			rf.read_attribute(rb_grp_id, "vz_ang_bc", bc_value);
			rb.set_vz_ang_bc(bc_value);
		}

		// rigid mesh
		double grid_xl, grid_yl, grid_zl;
		double grid_xu, grid_yu, grid_zu;
		double grid_hx, grid_hy, grid_hz;
		unsigned long long grid_x_num, grid_y_num, grid_z_num;
		unsigned long long face_in_grid_list_len, face_num;
		rf.read_attribute(rb_grp_id, "grid_xl", grid_xl);
		rf.read_attribute(rb_grp_id, "grid_yl", grid_yl);
		rf.read_attribute(rb_grp_id, "grid_zl", grid_zl);
		rf.read_attribute(rb_grp_id, "grid_hx", grid_hx);
		rf.read_attribute(rb_grp_id, "grid_hy", grid_hy);
		rf.read_attribute(rb_grp_id, "grid_hz", grid_hz);
		rf.read_attribute(rb_grp_id, "grid_x_num", grid_x_num);
		rf.read_attribute(rb_grp_id, "grid_y_num", grid_y_num);
		rf.read_attribute(rb_grp_id, "grid_z_num", grid_z_num);
		rf.read_attribute(rb_grp_id, "face_grid_list_len", face_in_grid_list_len);
		rf.read_attribute(rb_grp_id, "face_num", face_num);

		unsigned char *grid_pos_type;
		unsigned long long* face_in_grid_range;
		unsigned long long* face_in_grid_list;
		PointToTriangleDistance* pt_tri_dist;
		rb.init_from_hdf5_res_file(
			density, m, moi, inv_moi, T_mat,
			acc, acc_ang, vec, vec_ang, pos, pos_ang,
			force_cont, force_ext,
			grid_xl, grid_yl, grid_zl,
			grid_hx, grid_hy, grid_hz,
			grid_x_num, grid_y_num, grid_z_num,
			face_in_grid_list_len, face_num,
			grid_pos_type, face_in_grid_range,
			face_in_grid_list, pt_tri_dist);
		size_t grid_num = grid_x_num * grid_y_num * grid_z_num;
		rf.read_dataset(
			rb_grp_id,
			"grid_pos_type",
			grid_num,
			grid_pos_type);
		rf.read_dataset(
			rb_grp_id,
			"face_in_grid_range",
			grid_num + 1,
			face_in_grid_range);
		rf.read_dataset(
			rb_grp_id,
			"face_in_grid_list",
			face_in_grid_list_len,
			face_in_grid_list);

		hid_t pt2tri_dt = get_pt_to_tri_dist_dt_id();
		rf.read_dataset(
			rb_grp_id,
			"pt_tri_dist",
			face_num,
			pt_tri_dist,
			pt2tri_dt
			);
		H5Tclose(pt2tri_dt);

		return 0;
	}
}
