#include "Simulations_pcp.h"

#include "Model_hdf5_utilities.h"

namespace Model_hdf5_utilities
{
	// material model container
	int output_material_model_container_to_hdf5_file(
		MatModel::MatModelContainer& mc,
		ResultFile_hdf5& rf,
		hid_t mc_grp_id
		)
	{
		if (mc_grp_id < 0)
			return -1;

		size_t mm_id, mm_num;

		// linear elasticity
		mm_num = mc.get_num_LinearElasticity();
		if (mm_num)
		{
			rf.write_attribute(mc_grp_id, "LinearElasticity_num", mm_num);

			LinearElasticityStateData* mm_data = new LinearElasticityStateData[mm_num];
			mm_id = 0;
			for (MatModel::LinearElasticity* iter = mc.first_LinearElasticity();
				 mc.is_not_end_LinearElasticity(iter);
				 iter = mc.next_LinearElasticity(iter))
			{
				mm_data[mm_id].from_mm(*iter);
				++mm_id;
			}
			hid_t le_dt_id = get_le_hdf5_dt_id();
			rf.write_dataset(mc_grp_id, "LinearElasticity", mm_num, mm_data, le_dt_id);
			H5Tclose(le_dt_id);
			delete[] mm_data;
		}

		// modified cam clay
		mm_num = mc.get_num_ModifiedCamClay();
		if (mm_num)
		{
			rf.write_attribute(mc_grp_id, "ModifiedCamClay_num", mm_num);

			ModifiedCamClayStateData* mm_data = new ModifiedCamClayStateData[mm_num];
			mm_id = 0;
			for (MatModel::ModifiedCamClay* iter = mc.first_ModifiedCamClay();
				 mc.is_not_end_ModifiedCamClay(iter);
				 iter = mc.next_ModifiedCamClay(iter))
			{
				mm_data[mm_id].from_mm(*iter);
				++mm_id;
			}
			hid_t mcc_dt_id = get_mcc_hdf5_dt_id();
			rf.write_dataset(mc_grp_id, "ModifiedCamClay", mm_num, mm_data, mcc_dt_id);
			H5Tclose(mcc_dt_id);
			delete[] mm_data;
		}

		// undrained modified cam clay
		mm_num = mc.get_num_UndrainedModifiedCamClay();
		if (mm_num)
		{
			rf.write_attribute(mc_grp_id, "UndrainedModifiedCamClay_num", mm_num);

			UndrainedModifiedCamClayStateData* mm_data
				= new UndrainedModifiedCamClayStateData[mm_num];
			mm_id = 0;
			for (MatModel::UndrainedModifiedCamClay* iter = mc.first_UndrainedModifiedCamClay();
				 mc.is_not_end_UndrainedModifiedCamClay(iter);
				 iter = mc.next_UndrainedModifiedCamClay(iter))
			{
				mm_data[mm_id].from_mm(*iter);
				++mm_id;
			}
			hid_t mcc_dt_id = get_undrained_mcc_hdf5_dt_id();
			rf.write_dataset(mc_grp_id, "UndrainedModifiedCamClay", mm_num, mm_data, mcc_dt_id);
			H5Tclose(mcc_dt_id);
			delete[] mm_data;
		}

		return 0;
	}

	int load_material_model_container_from_hdf5_file(
		MatModel::MatModelContainer& mc,
		ResultFile_hdf5& rf,
		hid_t mc_grp_id
		)
	{
		if (mc_grp_id < 0)
			return -1;

		hid_t mm_dset_id;
		size_t mm_num;

		// linear elasticity
		if (rf.has_dataset(mc_grp_id, "LinearElasticity"))
		{
			rf.read_attribute(mc_grp_id, "LinearElasticity_num", mm_num);

			LinearElasticityStateData* mm_data = new LinearElasticityStateData[mm_num];
			hid_t le_dt_id = get_le_hdf5_dt_id();
			rf.read_dataset(
				mc_grp_id,
				"LinearElasticity",
				mm_num,
				mm_data,
				le_dt_id
			);
			H5Tclose(le_dt_id);
			MatModel::LinearElasticity* mms = mc.add_LinearElasticity(mm_num);
			for (size_t mm_id = 0; mm_id < mm_num; ++mm_id)
			{
				LinearElasticityStateData& mmd = mm_data[mm_id];
				MatModel::LinearElasticity& mm = mms[mm_id];
				mmd.to_mm(mm);
			}
			delete[] mm_data;
		}

		// modified cam clay
		if (rf.has_dataset(mc_grp_id, "ModifiedCamClay"))
		{
			rf.read_attribute(mc_grp_id, "ModifiedCamClay_num", mm_num);

			// get data
			ModifiedCamClayStateData* mm_data = new ModifiedCamClayStateData[mm_num];
			hid_t mcc_dt_id = get_mcc_hdf5_dt_id();
			rf.read_dataset(
				mc_grp_id,
				"ModifiedCamClay",
				mm_num,
				mm_data,
				mcc_dt_id
			);
			H5Tclose(mcc_dt_id);
			MatModel::ModifiedCamClay* mms = mc.add_ModifiedCamClay(mm_num);
			for (size_t mm_id = 0; mm_id < mm_num; ++mm_id)
			{
				ModifiedCamClayStateData& mmd = mm_data[mm_id];
				MatModel::ModifiedCamClay& mm = mms[mm_id];
				mmd.to_mm(mm);
			}
			delete[] mm_data;
		}

		// undrained modified cam clay
		if (rf.has_dataset(mc_grp_id, "UndrainedModifiedCamClay"))
		{
			rf.read_attribute(mc_grp_id, "UndrainedModifiedCamClay_num", mm_num);

			// get data
			UndrainedModifiedCamClayStateData* mm_data
				= new UndrainedModifiedCamClayStateData[mm_num];
			hid_t umcc_dt_id = get_undrained_mcc_hdf5_dt_id();
			rf.read_dataset(
				mc_grp_id,
				"UndrainedModifiedCamClay",
				mm_num,
				mm_data,
				umcc_dt_id
			);
			H5Tclose(umcc_dt_id);
			MatModel::UndrainedModifiedCamClay* mms
				= mc.add_UndrainedModifiedCamClay(mm_num);
			for (size_t mm_id = 0; mm_id < mm_num; ++mm_id)
			{
				UndrainedModifiedCamClayStateData& mmd = mm_data[mm_id];
				MatModel::UndrainedModifiedCamClay& mm = mms[mm_id];
				mmd.to_mm(mm);
			}
			delete[] mm_data;
		}

		return 0;
	}

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

	// only output state info, no mesh info
	int output_rigid_tetrahedron_mesh_state_to_hdf5_file(
		RigidTetrahedronMesh& rtm,
		ResultFile_hdf5& rf,
		hid_t rtm_grp_id
		)
	{
		rf.write_attribute(rtm_grp_id, "density", rtm.get_density());

		rf.write_attribute(rtm_grp_id, "ax", rtm.get_ax());
		rf.write_attribute(rtm_grp_id, "ay", rtm.get_ay());
		rf.write_attribute(rtm_grp_id, "az", rtm.get_az());
		rf.write_attribute(rtm_grp_id, "vx", rtm.get_vx());
		rf.write_attribute(rtm_grp_id, "vy", rtm.get_vy());
		rf.write_attribute(rtm_grp_id, "vz", rtm.get_vz());
		rf.write_attribute(rtm_grp_id, "x", rtm.get_x());
		rf.write_attribute(rtm_grp_id, "y", rtm.get_y());
		rf.write_attribute(rtm_grp_id, "z", rtm.get_z());

		rf.write_attribute(rtm_grp_id, "fx_ext", rtm.get_fx_ext());
		rf.write_attribute(rtm_grp_id, "fy_ext", rtm.get_fy_ext());
		rf.write_attribute(rtm_grp_id, "fz_ext", rtm.get_fz_ext());
		rf.write_attribute(rtm_grp_id, "fx_contact", rtm.get_fx_contact());
		rf.write_attribute(rtm_grp_id, "fy_contact", rtm.get_fy_contact());
		rf.write_attribute(rtm_grp_id, "fz_contact", rtm.get_fz_contact());

		rf.write_attribute(rtm_grp_id, "ax_ang", rtm.get_ax_ang());
		rf.write_attribute(rtm_grp_id, "ay_ang", rtm.get_ay_ang());
		rf.write_attribute(rtm_grp_id, "az_ang", rtm.get_az_ang());
		rf.write_attribute(rtm_grp_id, "vx_ang", rtm.get_vx_ang());
		rf.write_attribute(rtm_grp_id, "vy_ang", rtm.get_vy_ang());
		rf.write_attribute(rtm_grp_id, "vz_ang", rtm.get_vz_ang());
		rf.write_attribute(rtm_grp_id, "x_ang", rtm.get_x_ang());
		rf.write_attribute(rtm_grp_id, "y_ang", rtm.get_y_ang());
		rf.write_attribute(rtm_grp_id, "z_ang", rtm.get_z_ang());

		rf.write_attribute(rtm_grp_id, "mx_ext", rtm.get_mx_ext());
		rf.write_attribute(rtm_grp_id, "my_ext", rtm.get_my_ext());
		rf.write_attribute(rtm_grp_id, "mz_ext", rtm.get_mz_ext());
		rf.write_attribute(rtm_grp_id, "mx_contact", rtm.get_mx_contact());
		rf.write_attribute(rtm_grp_id, "my_contact", rtm.get_my_contact());
		rf.write_attribute(rtm_grp_id, "mz_contact", rtm.get_mz_contact());

		if (rtm.has_ax_bc())
			rf.write_attribute(rtm_grp_id, "ax_bc", rtm.get_ax_bc());
		if (rtm.has_ay_bc())
			rf.write_attribute(rtm_grp_id, "ay_bc", rtm.get_ay_bc());
		if (rtm.has_az_bc())
			rf.write_attribute(rtm_grp_id, "az_bc", rtm.get_az_bc());

		if (rtm.has_vx_bc())
			rf.write_attribute(rtm_grp_id, "vx_bc", rtm.get_vx_bc());
		if (rtm.has_vy_bc())
			rf.write_attribute(rtm_grp_id, "vy_bc", rtm.get_vy_bc());
		if (rtm.has_vz_bc())
			rf.write_attribute(rtm_grp_id, "vz_bc", rtm.get_vz_bc());

		if (rtm.has_ax_ang_bc())
			rf.write_attribute(rtm_grp_id, "ax_ang_bc", rtm.get_ax_ang_bc());
		if (rtm.has_ay_ang_bc())
			rf.write_attribute(rtm_grp_id, "ay_ang_bc", rtm.get_ay_ang_bc());
		if (rtm.has_az_ang_bc())
			rf.write_attribute(rtm_grp_id, "az_ang_bc", rtm.get_az_ang_bc());
		
		if (rtm.has_vx_ang_bc())
			rf.write_attribute(rtm_grp_id, "vx_ang_bc", rtm.get_vx_ang_bc());
		if (rtm.has_vy_ang_bc())
			rf.write_attribute(rtm_grp_id, "vy_ang_bc", rtm.get_vy_ang_bc());
		if (rtm.has_vz_ang_bc())
			rf.write_attribute(rtm_grp_id, "vz_ang_bc", rtm.get_vz_ang_bc());
		
		return 0;
	}

	int load_rigid_tetrahedron_mesh_state_from_hdf5_file(
		RigidTetrahedronMesh& rtm,
		ResultFile_hdf5& rf,
		hid_t rtm_grp_id
		)
	{
		if (rtm_grp_id < 0)
			return -1;

		double rtm_density;
		double rtm_fx_contact, rtm_fy_contact, rtm_fz_contact;
		double rtm_ax, rtm_ay, rtm_az;
		double rtm_vx, rtm_vy, rtm_vz;
		double rtm_x, rtm_y, rtm_z;
		double rtm_mx_contact, rtm_my_contact, rtm_mz_contact;
		double rtm_ax_ang, rtm_ay_ang, rtm_az_ang;
		double rtm_vx_ang, rtm_vy_ang, rtm_vz_ang;
		double rtm_x_ang, rtm_y_ang, rtm_z_ang;

		rf.read_attribute(rtm_grp_id, "density", rtm_density);

		rf.read_attribute(rtm_grp_id, "ax", rtm_ax);
		rf.read_attribute(rtm_grp_id, "ay", rtm_ay);
		rf.read_attribute(rtm_grp_id, "az", rtm_az);
		rf.read_attribute(rtm_grp_id, "vx", rtm_vx);
		rf.read_attribute(rtm_grp_id, "vy", rtm_vy);
		rf.read_attribute(rtm_grp_id, "vz", rtm_vz);
		rf.read_attribute(rtm_grp_id, "x", rtm_x);
		rf.read_attribute(rtm_grp_id, "y", rtm_y);
		rf.read_attribute(rtm_grp_id, "z", rtm_z);
		rf.read_attribute(rtm_grp_id, "fx_contact", rtm_fx_contact);
		rf.read_attribute(rtm_grp_id, "fy_contact", rtm_fy_contact);
		rf.read_attribute(rtm_grp_id, "fz_contact", rtm_fz_contact);

		rf.read_attribute(rtm_grp_id, "ax_ang", rtm_ax_ang);
		rf.read_attribute(rtm_grp_id, "ay_ang", rtm_ay_ang);
		rf.read_attribute(rtm_grp_id, "az_ang", rtm_az_ang);
		rf.read_attribute(rtm_grp_id, "vx_ang", rtm_vx_ang);
		rf.read_attribute(rtm_grp_id, "vy_ang", rtm_vy_ang);
		rf.read_attribute(rtm_grp_id, "vz_ang", rtm_vz_ang);
		rf.read_attribute(rtm_grp_id, "x_ang", rtm_x_ang);
		rf.read_attribute(rtm_grp_id, "y_ang", rtm_y_ang);
		rf.read_attribute(rtm_grp_id, "z_ang", rtm_z_ang);
		rf.read_attribute(rtm_grp_id, "mx_contact", rtm_mx_contact);
		rf.read_attribute(rtm_grp_id, "my_contact", rtm_my_contact);
		rf.read_attribute(rtm_grp_id, "mz_contact", rtm_mz_contact);

		rtm.set_init_state(
			rtm_density,
			rtm_fx_contact, rtm_fy_contact, rtm_fz_contact,
			rtm_ax, rtm_ay, rtm_az,
			rtm_vx, rtm_vy, rtm_vz,
			rtm_x, rtm_y, rtm_z,
			rtm_mx_contact, rtm_my_contact, rtm_mz_contact,
			rtm_ax_ang, rtm_ay_ang, rtm_az_ang,
			rtm_vx_ang, rtm_vy_ang, rtm_vz_ang,
			rtm_x_ang, rtm_y_ang, rtm_z_ang
			);

		double bc_value;
		rf.read_attribute(rtm_grp_id, "fx_ext", bc_value);
		rtm.add_fx_ext(bc_value);
		rf.read_attribute(rtm_grp_id, "fy_ext", bc_value);
		rtm.add_fy_ext(bc_value);
		rf.read_attribute(rtm_grp_id, "fz_ext", bc_value);
		rtm.add_fz_ext(bc_value);
		if (rf.has_attribute(rtm_grp_id, "ax_bc"))
		{
			rf.read_attribute(rtm_grp_id, "ax_bc", bc_value);
			rtm.set_ax_bc(bc_value);
		}
		if (rf.has_attribute(rtm_grp_id, "ay_bc"))
		{
			rf.read_attribute(rtm_grp_id, "ay_bc", bc_value);
			rtm.set_ay_bc(bc_value);
		}
		if (rf.has_attribute(rtm_grp_id, "az_bc"))
		{
			rf.read_attribute(rtm_grp_id, "az_bc", bc_value);
			rtm.set_az_bc(bc_value);
		}
		if (rf.has_attribute(rtm_grp_id, "vx_bc"))
		{
			rf.read_attribute(rtm_grp_id, "vx_bc", bc_value);
			rtm.set_vx_bc(bc_value);
		}
		if (rf.has_attribute(rtm_grp_id, "vy_bc"))
		{
			rf.read_attribute(rtm_grp_id, "vy_bc", bc_value);
			rtm.set_vy_bc(bc_value);
		}
		if (rf.has_attribute(rtm_grp_id, "vz_bc"))
		{
			rf.read_attribute(rtm_grp_id, "vz_bc", bc_value);
			rtm.set_vz_bc(bc_value);
		}

		rf.read_attribute(rtm_grp_id, "mx_ext", bc_value);
		rtm.add_mx_ext(bc_value);
		rf.read_attribute(rtm_grp_id, "my_ext", bc_value);
		rtm.add_my_ext(bc_value);
		rf.read_attribute(rtm_grp_id, "mz_ext", bc_value);
		rtm.add_mz_ext(bc_value);
		if (rf.has_attribute(rtm_grp_id, "ax_ang_bc"))
		{
			rf.read_attribute(rtm_grp_id, "ax_ang_bc", bc_value);
			rtm.set_ax_bc(bc_value);
		}
		if (rf.has_attribute(rtm_grp_id, "ay_ang_bc"))
		{
			rf.read_attribute(rtm_grp_id, "ay_ang_bc", bc_value);
			rtm.set_ay_bc(bc_value);
		}
		if (rf.has_attribute(rtm_grp_id, "az_ang_bc"))
		{
			rf.read_attribute(rtm_grp_id, "az_ang_bc", bc_value);
			rtm.set_az_bc(bc_value);
		}
		if (rf.has_attribute(rtm_grp_id, "vx_ang_bc"))
		{
			rf.read_attribute(rtm_grp_id, "vx_ang_bc", bc_value);
			rtm.set_vx_bc(bc_value);
		}
		if (rf.has_attribute(rtm_grp_id, "vy_ang_bc"))
		{
			rf.read_attribute(rtm_grp_id, "vy_ang_bc", bc_value);
			rtm.set_vy_bc(bc_value);
		}
		if (rf.has_attribute(rtm_grp_id, "vz_ang_bc"))
		{
			rf.read_attribute(rtm_grp_id, "vz_ang_bc", bc_value);
			rtm.set_vz_bc(bc_value);
		}
		
		return 0;
	}

	int output_rigid_tetrahedron_mesh_to_hdf5_file(
		RigidTetrahedronMesh& rtm,
		ResultFile_hdf5& rf,
		hid_t rtm_grp_id
		)
	{
		if (rtm_grp_id < 0)
			return -1;
		
		int res;
		size_t node_num = rtm.get_node_num();
		rf.write_attribute(rtm_grp_id, "node_num", node_num);
		RigidTetrahedronMesh::Node* nodes = rtm.get_nodes();
		RigidTehMeshNodeData* node_data = new RigidTehMeshNodeData[node_num];
		for (size_t n_id = 0; n_id < node_num; ++n_id)
		{
			RigidTehMeshNodeData &nd = node_data[n_id];
			RigidTetrahedronMesh::Node &n = nodes[n_id];
			nd.from_node(n);
		}
		hid_t node_dt_id = get_rigid_teh_mesh_node_dt_id();
		res = rf.write_dataset(
			rtm_grp_id,
			"NodeData",
			node_num,
			node_data,
			node_dt_id
		);
		H5Tclose(node_dt_id);
		delete[] node_data;

		size_t elem_num = rtm.get_elem_num();
		rf.write_attribute(rtm_grp_id, "elem_num", elem_num);
		RigidTetrahedronMesh::Element* elems = rtm.get_elems();
		RigidTehMeshElemData* elem_data = new RigidTehMeshElemData[elem_num];
		for (size_t e_id = 0; e_id < elem_num; ++e_id)
		{
			RigidTehMeshElemData& ed = elem_data[e_id];
			RigidTetrahedronMesh::Element &e = elems[e_id];
			ed.from_elem(e);
		}
		hid_t elem_dt_id = get_rigid_teh_mesh_elem_dt_id();
		res = rf.write_dataset(
			rtm_grp_id,
			"ElementData",
			elem_num,
			elem_data,
			elem_dt_id
			);
		H5Tclose(elem_dt_id);
		delete[] elem_data;

		//size_t bface_num = rtm.get_bface_num();
		//rf.write_attribute(rtm_grp_id, "bface_num", bface_num);
		//RigidTetrahedronMesh::Face* bfaces = rtm.get_bfaces();
		//RigidTehMeshFaceData* bface_data = new RigidTehMeshFaceData[bface_num];
		//for (size_t f_id = 0; f_id < bface_num; ++f_id)
		//{
		//	RigidTehMeshFaceData &fd = bface_data[f_id];
		//	RigidTetrahedronMesh::Face &f = bfaces[f_id];
		//	fd.from_face(f);
		//}
		//hid_t face_dt_id = get_rigid_teh_mesh_face_dt_id();
		//res = rf.write_dataset(
		//	rtm_grp_id,
		//	"BoundaryFaceData",
		//	bface_num,
		//	bface_data,
		//	face_dt_id
		//);
		//H5Tclose(face_dt_id);
		//delete[] bface_data;
		
		return 0;
	}
	
	int load_rigid_tetrahedron_mesh_from_hdf5_file(
		RigidTetrahedronMesh& rtm,
		ResultFile_hdf5& rf,
		hid_t rtm_grp_id
		)
	{
		if (rtm_grp_id < 0)
			return -1;

		int res;
		size_t node_num;
		rf.read_attribute(rtm_grp_id, "node_num", node_num);
		RigidTehMeshNodeData* node_data = new RigidTehMeshNodeData[node_num];
		hid_t node_dt_id = get_rigid_teh_mesh_node_dt_id();
		res = rf.read_dataset(
			rtm_grp_id,
			"NodeData",
			node_num,
			node_data,
			node_dt_id
			);
		H5Tclose(node_dt_id);
		RigidTetrahedronMesh::Node* nodes = rtm.alloc_nodes(node_num);
		for (size_t n_id = 0; n_id < node_num; ++n_id)
		{
			RigidTehMeshNodeData& nd = node_data[n_id];
			RigidTetrahedronMesh::Node& n = nodes[n_id];
			nd.to_node(n);
		}
		delete[] node_data;

		size_t elem_num;
		rf.read_attribute(rtm_grp_id, "elem_num", elem_num);
		RigidTehMeshElemData* elem_data = new RigidTehMeshElemData[elem_num];
		hid_t elem_dt_id = get_rigid_teh_mesh_elem_dt_id();
		res = rf.read_dataset(
			rtm_grp_id,
			"ElementData",
			elem_num,
			elem_data,
			elem_dt_id
			);
		H5Tclose(elem_dt_id);
		RigidTetrahedronMesh::Element* elems = rtm.alloc_elements(elem_num);
		for (size_t e_id = 0; e_id < elem_num; ++e_id)
		{
			RigidTehMeshElemData& ed = elem_data[e_id];
			RigidTetrahedronMesh::Element& e = elems[e_id];
			ed.to_elem(e);
		}
		delete[] elem_data;

		rtm.init_mesh_properties_after_loading();

		//size_t bface_num;
		//rf.read_attribute(rtm_grp_id, "bface_num", bface_num);
		//RigidTehMeshFaceData* bface_data = new RigidTehMeshFaceData[bface_num];
		//hid_t face_dt_id = get_rigid_teh_mesh_face_dt_id();
		//res = rf.read_dataset(
		//	rtm_grp_id,
		//	"BoundaryFaceData",
		//	bface_num,
		//	bface_data,
		//	face_dt_id
		//);
		//H5Tclose(face_dt_id);
		//RigidTetrahedronMesh::Face* bfaces = rtm.alloc_bfaces(bface_num);
		//for (size_t f_id = 0; f_id < bface_num; ++f_id)
		//{
		//	RigidTehMeshFaceData& fd = bface_data[f_id];
		//	RigidTetrahedronMesh::Face& f = bfaces[f_id];
		//	fd.to_face(f);
		//}
		//delete[] bface_data;

		return 0;
	}
}
