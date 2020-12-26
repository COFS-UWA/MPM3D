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

		// von mises
		mm_num = mc.get_num_VonMises();
		if (mm_num)
		{
			rf.write_attribute(mc_grp_id, "VonMises_num", mm_num);

			VonMisesStateData *mm_data = new VonMisesStateData[mm_num];
			mm_id = 0;
			for (MatModel::VonMises* iter = mc.first_VonMises();
				mc.is_not_end_VonMises(iter); iter = mc.next_VonMises(iter))
			{
				mm_data[mm_id].from_mm(*iter);
				++mm_id;
			}
			hid_t vm_dt_id = get_von_mises_hdf5_dt_id();
			rf.write_dataset(mc_grp_id, "VonMises", mm_num, mm_data, vm_dt_id);
			H5Tclose(vm_dt_id);
			delete[] mm_data;
		}

		// tresca
		mm_num = mc.get_num_Tresca();
		if (mm_num)
		{
			rf.write_attribute(mc_grp_id, "Tresca_num", mm_num);

			TrescaStateData* mm_data = new TrescaStateData[mm_num];
			mm_id = 0;
			for (MatModel::Tresca* iter = mc.first_Tresca();
				 mc.is_not_end_Tresca(iter); iter = mc.next_Tresca(iter))
			{
				mm_data[mm_id].from_mm(*iter);
				++mm_id;
			}
			hid_t tc_dt_id = get_tresca_hdf5_dt_id();
			rf.write_dataset(mc_grp_id, "Tresca", mm_num, mm_data, tc_dt_id);
			H5Tclose(tc_dt_id);
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

		// von mises
		if (rf.has_dataset(mc_grp_id, "VonMises"))
		{
			rf.read_attribute(mc_grp_id, "VonMises_num", mm_num);

			// get data
			VonMisesStateData* mm_data = new VonMisesStateData[mm_num];
			hid_t vm_dt_id = get_von_mises_hdf5_dt_id();
			rf.read_dataset(
				mc_grp_id,
				"VonMises",
				mm_num,
				mm_data,
				vm_dt_id
				);
			H5Tclose(vm_dt_id);
			MatModel::VonMises* mms = mc.add_VonMises(mm_num);
			for (size_t mm_id = 0; mm_id < mm_num; ++mm_id)
			{
				VonMisesStateData& mmd = mm_data[mm_id];
				MatModel::VonMises &mm = mms[mm_id];
				mmd.to_mm(mm);
			}
			delete[] mm_data;
		}

		// tresca
		if (rf.has_dataset(mc_grp_id, "Tresca"))
		{
			rf.read_attribute(mc_grp_id, "Tresca_num", mm_num);

			// get data
			TrescaStateData* mm_data = new TrescaStateData[mm_num];
			hid_t tc_dt_id = get_tresca_hdf5_dt_id();
			rf.read_dataset(
				mc_grp_id,
				"Tresca",
				mm_num,
				mm_data,
				tc_dt_id
			);
			H5Tclose(tc_dt_id);
			MatModel::Tresca* mms = mc.add_Tresca(mm_num);
			for (size_t mm_id = 0; mm_id < mm_num; ++mm_id)
			{
				TrescaStateData& mmd = mm_data[mm_id];
				MatModel::Tresca& mm = mms[mm_id];
				mmd.to_mm(mm);
			}
			delete[] mm_data;
		}
		
		return 0;
	}
}
