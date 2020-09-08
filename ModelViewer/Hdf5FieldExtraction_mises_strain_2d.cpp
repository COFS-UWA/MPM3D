#include "ModelViewer_pcp.h"

#include "Hdf5FieldExtraction_mises_strain_2d.h"

bool Hdf5FieldExtraction_mises_strain_2d::validate_data_type()
{
	hid_t pcl_dt_id = data_loader->get_pcl_data_type();
	int mem_num = H5Tget_nmembers(pcl_dt_id);
	int tag = 0;
	for (size_t mem_id = 0; mem_id < mem_num; ++mem_id)
	{
		const char* mem_name = H5Tget_member_name(pcl_dt_id, mem_id);
		if (strcmp(mem_name, "e11") == 0)
		{
			e11_offset = H5Tget_member_offset(pcl_dt_id, mem_id);
			++tag;
		}
		else if (strcmp(mem_name, "e22") == 0)
		{
			e22_offset = H5Tget_member_offset(pcl_dt_id, mem_id);
			++tag;
		}
		else if (strcmp(mem_name, "e12") == 0)
		{
			e12_offset = H5Tget_member_offset(pcl_dt_id, mem_id);
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

int Hdf5FieldExtraction_mises_strain_2d::extract_pcl_fld_data(double* pcl_fld_data)
{
	if (!all_fld_is_found)
		return -1;

	Hdf5DataLoader& loader = *data_loader;
	size_t pcl_size = loader.get_pcl_size();
	size_t pcl_num = loader.get_pcl_num();
	char *cur_pcl = loader.get_pcl_field_data();
	double pcl_e11, pcl_e22, pcl_e12;
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		pcl_e11 = *(double *)(cur_pcl + e11_offset);
		pcl_e22 = *(double *)(cur_pcl + e22_offset);
		pcl_e12 = *(double *)(cur_pcl + e12_offset);
		pcl_fld_data[pcl_id] = 2.0 / 3.0 * sqrt(0.5 *
			((pcl_e11 - pcl_e22) * (pcl_e11 - pcl_e22)
			+ pcl_e11 * pcl_e11 + pcl_e22 * pcl_e22	+ 6.0 * pcl_e12 * pcl_e12));
		cur_pcl += pcl_size;
	}
	return 0;
}
