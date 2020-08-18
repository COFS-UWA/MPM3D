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

}
