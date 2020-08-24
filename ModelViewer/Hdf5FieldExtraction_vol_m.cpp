#include "ModelViewer_pcp.h"

#include "Hdf5FieldExtraction_vol_m.h"

bool Hdf5FieldExtraction_vol_m::validate_data_type()
{
	hid_t pcl_dt_id = data_loader->get_pcl_data_type();
	int mem_num = H5Tget_nmembers(pcl_dt_id);
	int tag = 0;
	for (size_t mem_id = 0; mem_id < mem_num; ++mem_id)
	{
		const char* mem_name = H5Tget_member_name(pcl_dt_id, mem_id);
		if (strcmp(mem_name, "m_s") == 0)
		{
			m_s_offset = H5Tget_member_offset(pcl_dt_id, mem_id);
			++tag;
		}
		else if (strcmp(mem_name, "density_s") == 0)
		{
			density_s_offset = H5Tget_member_offset(pcl_dt_id, mem_id);
			++tag;
		}
		else if (strcmp(mem_name, "n") == 0)
		{
			n_offset = H5Tget_member_offset(pcl_dt_id, mem_id);
			++tag;
		}

		if (tag == 3)
		{
			all_fld_is_found = true;
			break;
		}
	}
	return all_fld_is_found;
}

int Hdf5FieldExtraction_vol_m::extract_pcl_fld_data(double* pcl_fld_data)
{
	if (!all_fld_is_found)
		return -1;

	Hdf5DataLoader& loader = *data_loader;
	size_t pcl_size = loader.get_pcl_size();
	size_t pcl_num = loader.get_pcl_num();
	char *cur_pcl = loader.get_pcl_field_data();
	double pcl_m_s, pcl_density_s, pcl_n;
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		pcl_m_s = *(double *)(cur_pcl + m_s_offset);
		pcl_density_s = *(double *)(cur_pcl + density_s_offset);
		pcl_n = *(double *)(cur_pcl + n_offset);
		pcl_fld_data[pcl_id] = pcl_m_s / (pcl_density_s * (1.0 - pcl_n));
		cur_pcl += pcl_size;
	}
	return 0;
}
