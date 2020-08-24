#include "ModelViewer_pcp.h"

#include <exception>

#include "Hdf5DataLoader.h"

Hdf5DataLoader::Hdf5DataLoader() :
	res_file(nullptr), th_id(-1), pcl_dt_id(-1), pcl_size(0),
	frame_id(std::numeric_limits<size_t>::max()), pcl_num(0)
{

}

Hdf5DataLoader::~Hdf5DataLoader()
{
	close_res_file();
}

void Hdf5DataLoader::close_res_file()
{
	if (!res_file)
		return;
	ResultFile_hdf5& rf = *res_file;
	if (th_id >= 0)
	{
		rf.close_group(th_id);
		th_id = -1;
	}
	if (pcl_dt_id >= 0)
	{
		H5Tclose(pcl_dt_id);
		pcl_dt_id = -1;
	}
	res_file = nullptr;
}

int Hdf5DataLoader::set_time_history(
	ResultFile_hdf5& rf,
	const char* th_name
	)
{
	char exception_msg[250];
	close_res_file();

	res_file = &rf;

	hid_t th_grp_id = rf.get_time_history_grp_id();
	// check if the time history exists
	if (!rf.has_group(th_grp_id, th_name))
	{
		snprintf(exception_msg, sizeof(exception_msg),
			"hdf5 file has no time history %s.", th_name);
		throw std::exception(exception_msg);
		close_res_file();
		return -1;
	}
	th_id = rf.open_group(th_grp_id, th_name);

	if (!rf.has_group(th_id, "frame_0"))
	{
		snprintf(exception_msg, sizeof(exception_msg),
			"hdf5 file doesn't have frame_0 in time history %s.", th_name);
		throw std::exception(exception_msg);
		close_res_file();
		return -2;
	}
	
	hid_t frame_grp_id = rf.open_group(th_id, "frame_0");
	hid_t pcl_grp_id = rf.open_group(frame_grp_id, "ParticleData");
	hid_t pcl_dset_id = rf.open_dataset(pcl_grp_id, "field");
	// get field data type
	pcl_dt_id = H5Dget_type(pcl_dset_id);
	pcl_size = H5Tget_size(pcl_dt_id);
	rf.close_dataset(pcl_dset_id);
	rf.close_group(pcl_grp_id);
	rf.close_group(frame_grp_id);

	return 0;
}

int Hdf5DataLoader::load_frame_data(size_t fm_id)
{
	if (frame_id == fm_id)
		return 0;

	int res;
	ResultFile_hdf5& rf = *res_file;
	char frame_name[50];
	snprintf(frame_name, sizeof(frame_name), "frame_%zu", fm_id);
	hid_t frame_grp_id = rf.open_group(th_id, frame_name);
	hid_t pcl_data_id = rf.open_group(frame_grp_id, "ParticleData");

	// pcl num
	res = rf.read_attribute(pcl_data_id, "pcl_num", pcl_num);
	if (res) return res;

	// pcl data
	pcl_fld_data_mem.reserve(pcl_size * pcl_num);
	char* pcls_data = pcl_fld_data_mem.get_mem();
	res = rf.read_dataset(pcl_data_id, "field", pcl_num, (void*)pcls_data, pcl_dt_id);
	if (res) return res;

	rf.close_group(pcl_data_id);
	rf.close_group(frame_grp_id);
	frame_id = fm_id;
	return 0;
}
