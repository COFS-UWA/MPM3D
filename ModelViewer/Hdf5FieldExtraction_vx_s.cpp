#include "ModelViewer_pcp.h"

#include "Hdf5FieldExtraction_vx_s.h"

bool Hdf5FieldExtraction_vx_s::validate_data_type()
{
	hid_t pcl_dt_id = data_loader->get_pcl_data_type();
	int mem_num = H5Tget_nmembers(pcl_dt_id);
	for (size_t mem_id = 0; mem_id < mem_num; ++mem_id)
	{
		const char* mem_name = H5Tget_member_name(pcl_dt_id, mem_id);
		if (strcmp(mem_name, "vx_s") == 0)
		{
			vx_s_offset = H5Tget_member_offset(pcl_dt_id, mem_id);
			vx_s_fld_is_found = true;
			break;
		}
	}
	return vx_s_fld_is_found;
}

int Hdf5FieldExtraction_vx_s::extract_pcl_fld_data(double* pcl_fld_data)
{
	if (!vx_s_fld_is_found)
		return -1;

	Hdf5DataLoader& loader = *data_loader;
	size_t pcl_size = loader.get_pcl_size();
	size_t pcl_num = loader.get_pcl_num();
	char *cur_pcl = loader.get_pcl_field_data();
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		pcl_fld_data[pcl_id] = *(double *)(cur_pcl + vx_s_offset);
		cur_pcl += pcl_size;
	}
	return 0;
}
