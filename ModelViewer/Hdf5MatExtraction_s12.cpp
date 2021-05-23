#include "ModelViewer_pcp.h"

#include "Hdf5MatExtraction_s12.h"

bool Hdf5MatExtraction_s12::validate_data_type()
{
	hid_t pcl_dt_id = data_loader->get_pcl_data_type();
	int mem_num = H5Tget_nmembers(pcl_dt_id);
	for (size_t mem_id = 0; mem_id < mem_num; ++mem_id)
	{
		const char* mem_name = H5Tget_member_name(pcl_dt_id, mem_id);
		if (strcmp(mem_name, "mat_id") == 0)
		{
			mat_id_offset = H5Tget_member_offset(pcl_dt_id, mem_id);
			mat_id_is_found = true;
			break;
		}
	}
	return mat_id_is_found;
}

int Hdf5MatExtraction_s12::extract_pcl_fld_data(double* pcl_fld_data)
{
	Hdf5DataLoader& loader = *data_loader;
	size_t pcl_size = loader.get_pcl_size();
	size_t pcl_num = loader.get_pcl_num();
	char* cur_pcl = loader.get_pcl_field_data();
	Hdf5DataLoader::MatModelMap &mat_map = loader.get_mat_model_map();
	size_t mat_id;
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		mat_id = *(unsigned long long *)(cur_pcl + mat_id_offset);
		auto mat_iter = mat_map.find(mat_id);
		if (mat_iter == mat_map.end())
			throw std::exception("Material model not found for pcl.\n");
		auto mat_type = mat_iter->second.type;
		auto &mat_info = loader.get_mat_model_info(mat_type);
		auto mat_tp = mat_iter->second.pmat;
		pcl_fld_data[pcl_id] = *(double*)(((char*)mat_tp) + mat_info.s12_off);
		cur_pcl += pcl_size;
	}
	return 0;
}
