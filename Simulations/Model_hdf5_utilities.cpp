#include "Simulations_pcp.h"

#include "Model_hdf5_utilities.h"

namespace Model_hdf5_utilities
{
	// rigid circle
	int output_rigid_circle_to_hdf5_file(
		RigidCircle& rc,
		ResultFile_hdf5& rf,
		hid_t rc_grp_id
		)
	{
		if (rc_grp_id < 0)
			return -1;

		rf.write_attribute(rc_grp_id, "radius", rc.get_radius());
		rf.write_attribute(rc_grp_id, "density", rc.get_density());
		rf.write_attribute(rc_grp_id, "rfx", rc.get_rfx());
		rf.write_attribute(rc_grp_id, "rfy", rc.get_rfy());
		rf.write_attribute(rc_grp_id, "rm", rc.get_rm());
		rf.write_attribute(rc_grp_id, "ax", rc.get_ax());
		rf.write_attribute(rc_grp_id, "ay", rc.get_ay());
		rf.write_attribute(rc_grp_id, "a_angle", rc.get_a_ang());
		rf.write_attribute(rc_grp_id, "vx", rc.get_vx());
		rf.write_attribute(rc_grp_id, "vy", rc.get_vy());
		rf.write_attribute(rc_grp_id, "v_angle", rc.get_v_ang());
		rf.write_attribute(rc_grp_id, "x", rc.get_x());
		rf.write_attribute(rc_grp_id, "y", rc.get_y());
		rf.write_attribute(rc_grp_id, "angle", rc.get_ang());

		// boundary conditions
		rf.write_attribute(rc_grp_id, "rfx_bc", rc.get_rfx_bc());
		rf.write_attribute(rc_grp_id, "rfy_bc", rc.get_rfy_bc());
		rf.write_attribute(rc_grp_id, "rm_bc", rc.get_rm_bc());
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
			rf.write_attribute(rc_grp_id, "v_ang_bc", rc.get_v_ang());

		return 0;
	}

	int load_rigid_circle_from_hdf5_file(
		RigidCircle& rc,
		ResultFile_hdf5& rf,
		hid_t rc_grp_id
	)
	{
		double rc_radius, rc_density;
		double rc_rfx, rc_rfy, rc_rm;
		double rc_ax, rc_ay, rc_a_ang;
		double rc_vx, rc_vy, rc_v_ang;
		double rc_x, rc_y, rc_ang;

		rf.read_attribute(rc_grp_id, "radius", rc_radius);
		rf.read_attribute(rc_grp_id, "density", rc_density);
		rf.read_attribute(rc_grp_id, "rfx", rc_rfx);
		rf.read_attribute(rc_grp_id, "rfy", rc_rfy);
		rf.read_attribute(rc_grp_id, "rm", rc_rm);
		rf.read_attribute(rc_grp_id, "ax", rc_ax);
		rf.read_attribute(rc_grp_id, "ay", rc_ay);
		rf.read_attribute(rc_grp_id, "a_angle", rc_a_ang);
		rf.read_attribute(rc_grp_id, "vx", rc_vx);
		rf.read_attribute(rc_grp_id, "vy", rc_vy);
		rf.read_attribute(rc_grp_id, "v_angle", rc_v_ang);
		rf.read_attribute(rc_grp_id, "x", rc_x);
		rf.read_attribute(rc_grp_id, "y", rc_y);
		rf.read_attribute(rc_grp_id, "angle", rc_ang);

		rc.set_init_state(
			rc_radius, rc_density,
			rc_rfx, rc_rfy, rc_rm,
			rc_ax, rc_ay, rc_a_ang,
			rc_vx, rc_vy, rc_v_ang,
			rc_x, rc_y, rc_ang
		);

		// boundary conditions
		double bc_value;
		rf.read_attribute(rc_grp_id, "rfx_bc", bc_value);
		rc.add_rfx_bc(bc_value);
		rf.read_attribute(rc_grp_id, "rfy_bc", bc_value);
		rc.add_rfy_bc(bc_value);
		rf.read_attribute(rc_grp_id, "rm_bc", bc_value);
		rc.add_rm_bc(bc_value);
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
			rc.set_a_ang_bc(bc_value);
		}
		if (rf.has_attribute(rc_grp_id, "vy_bc"))
		{
			rf.read_attribute(rc_grp_id, "vy_bc", bc_value);
			rc.set_a_ang_bc(bc_value);
		}
		if (rf.has_attribute(rc_grp_id, "v_ang_bc"))
		{
			rf.read_attribute(rc_grp_id, "v_ang_bc", bc_value);
			rc.set_a_ang_bc(bc_value);
		}

		return 0;
	}

}