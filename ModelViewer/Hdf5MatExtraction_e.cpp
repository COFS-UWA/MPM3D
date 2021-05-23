#include "ModelViewer_pcp.h"

#include "Hdf5MatExtraction_e.h"

bool Hdf5MatExtraction_e::validate_data_type()
{
	hid_t pcl_dt_id = data_loader->get_pcl_data_type();
	unsigned char var_found = 0;
	int mem_num = H5Tget_nmembers(pcl_dt_id);
	for (size_t mem_id = 0; mem_id < mem_num; ++mem_id)
	{
		const char* mem_name = H5Tget_member_name(pcl_dt_id, mem_id);
		if (strcmp(mem_name, "mat_id") == 0)
		{
			mat_id_offset = H5Tget_member_offset(pcl_dt_id, mem_id);
			++var_found;
			break;
		}
		if (strcmp(mem_name, "e") == 0)
		{
			e_offset = H5Tget_member_offset(pcl_dt_id, mem_id);
			++var_found;
			break;
		}
	}
	data_loader - ;
	var_is_found = (var_found == 2);
	return var_is_found;
}

int Hdf5MatExtraction_e::
	extract_pcl_fld_data(double* pcl_fld_data)
{
	Hdf5DataLoader& loader = *data_loader;
	const size_t pcl_size = loader.get_pcl_size();
	const size_t pcl_num = loader.get_pcl_num();
	char* cur_pcl = loader.get_pcl_field_data();
	Hdf5DataLoader::MatModelMap &mat_map = loader.get_mat_model_map();
	size_t mat_id;
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		mat_id = *(unsigned long long *)(cur_pcl + mat_id_offset);
		auto mat_iter = mat_map.find(mat_id);
		if (mat_iter == mat_map.end())
			throw std::exception("Material model not found for pcl.\n");
		auto mat_tp = mat_iter->second.pmat;
		pcl_fld_data[pcl_id] = *(double*)(((char*)mat_tp) + e_offset);
		cur_pcl += pcl_size;
	}
	return 0;
}
