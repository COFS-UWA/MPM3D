#ifndef __Hdf5_Data_Loader_h__
#define __Hdf5_Data_Loader_h__

#include "hdf5.h"
#include "ItemArray.hpp"
#include "ResultFile_hdf5.h"

class Hdf5DataLoader
{
protected:
	ResultFile_hdf5 *res_file;
	hid_t th_id;

	hid_t pcl_dt_id;
	size_t pcl_size;

	size_t frame_id;
	size_t pcl_num;
	MemoryUtils::ItemArray<char> pcl_fld_data_mem;

public:
	Hdf5DataLoader();
	~Hdf5DataLoader();
	void close_res_file();
	
	int set_time_history(ResultFile_hdf5& rf, const char* th_name);

	inline hid_t get_time_history_id() const noexcept { return th_id; }
	inline hid_t get_pcl_data_type() const noexcept { return pcl_dt_id; }
	inline size_t get_pcl_size() const noexcept { return pcl_size; }

	int load_frame_data(size_t fm_id);

	inline size_t get_pcl_num() const noexcept { return pcl_num; }
	inline char* get_pcl_field_data() const noexcept
	{ return pcl_fld_data_mem.get_mem(); }
};

#endif