#ifndef __hdf5_utils_H__
#define __hdf5_utils_H__

#include "hdf5.h"

inline void close_hdf5_dataset(hid_t dset_id)
{
	if (dset_id >= 0)
	{
		H5Dclose(dset_id);
		dset_id = -1;
	}
}

inline void close_hdf5_group(hid_t grp_id)
{
	if (grp_id >= 0)
	{
		H5Gclose(grp_id);
		grp_id = -1;
	}
}

inline void close_hdf5_file(hid_t file_id)
{
	if (file_id >= 0)
	{
		H5Fclose(file_id);
		file_id = -1;
	}
}

#endif