#include "SimulationCore_pcp.h"

#include "ResultFile_hdf5.h"

ResultFile_hdf5::ResultFile_hdf5() : 
	ResultFile(ResultFileType::Hdf5),
	file_id(-1), md_grp_id(-1), th_grp_id(-1) {}

ResultFile_hdf5::~ResultFile_hdf5() { close(); }

int ResultFile_hdf5::create(const char *file_name, bool over_write)
{
	if (over_write)
		file_id = H5Fcreate(file_name, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	else
		file_id = H5Fcreate(file_name, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
	return file_id < 0 ? -1 : 0;
}

int ResultFile_hdf5::open(const char *file_name, bool read_only)
{
	if (read_only)
	{
		file_id = H5Fopen(file_name, H5F_ACC_RDONLY, H5P_DEFAULT);
	}
	else
	{
		file_id = H5Fopen(file_name, H5F_ACC_RDWR, H5P_DEFAULT);
	}
	return file_id < 0 ? -1 : 0;
}

void ResultFile_hdf5::close(void)
{
	if (md_grp_id >= 0)
	{
		H5Gclose(md_grp_id);
		md_grp_id = -1;
	}
	if (th_grp_id >= 0)
	{
		H5Gclose(th_grp_id);
		th_grp_id = -1;
	}
	if (file_id >= 0)
	{
		herr_t status = H5Fclose(file_id);
		file_id = -1;
	}
}

// =========================== group operation ==========================
hid_t ResultFile_hdf5::create_group(hid_t parent_id, const char *name)
{
	hid_t cpl_id = H5Pcreate(H5P_FILE_CREATE);
	// can be iterated according to creation order
	H5Pset_link_creation_order(cpl_id, H5P_CRT_ORDER_TRACKED | H5P_CRT_ORDER_INDEXED);
	hid_t grp_id = H5Gcreate(parent_id, name, H5P_DEFAULT, cpl_id, H5P_DEFAULT);
	H5Pclose(cpl_id);
	return grp_id;
}

hid_t ResultFile_hdf5::open_group(hid_t parent_id, const char *name)
{
	return H5Gopen(parent_id, name, H5P_DEFAULT);
}

void ResultFile_hdf5::close_group(hid_t id) { if (id >= 0) H5Gclose(id); }

bool ResultFile_hdf5::has_group(hid_t parent_id, const char *name)
{
	herr_t res = H5Lexists(parent_id, name, H5P_DEFAULT);
	return res > 0 ? true : false;
}

// ============================== dataset operation ========================
hid_t ResultFile_hdf5::open_dataset(hid_t parent_id, const char *name)
{
	return H5Dopen(parent_id, name, H5P_DEFAULT);
}

void ResultFile_hdf5::close_dataset(hid_t id) { if (id >= 0) H5Dclose(id); }

bool ResultFile_hdf5::has_dataset(hid_t parent_id, const char *name)
{
	herr_t res = H5Lexists(parent_id, name, H5P_DEFAULT);
	return res > 0 ? true : false;
}

// =========================================================================
hid_t ResultFile_hdf5::get_model_data_grp_id(void)
{
	if (md_grp_id < 0)
	{
		if (has_group(file_id, "ModelData"))
			md_grp_id = open_group(file_id, "ModelData");
		else
			md_grp_id = create_group(file_id, "ModelData");
	}
	return md_grp_id;
}

hid_t ResultFile_hdf5::get_time_history_grp_id(void)
{
	if (th_grp_id < 0)
	{
		if (has_group(file_id, "TimeHistory"))
			th_grp_id = open_group(file_id, "TimeHistory");
		else
			th_grp_id = create_group(file_id, "TimeHistory");
	}
	return th_grp_id;
}

// ======================== routines for dataset =======================
int ResultFile_hdf5::write_dataset(
	hid_t grp_id,
	const char *dset_name, 
	size_t num,
	double *data
	)
{
	hid_t dataspace_id, dset_id;
	dataspace_id = H5Screate_simple(1, &num, nullptr);
	dset_id = H5Dcreate(grp_id, dset_name, H5T_NATIVE_DOUBLE, dataspace_id,
						H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (dset_id < 0) return -1;
	herr_t res = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, dataspace_id,
						  dataspace_id, H5P_DEFAULT, data);
	H5Dclose(dset_id);
	H5Sclose(dataspace_id);
	return res < 0 ? -2 : 0;
}

int ResultFile_hdf5::write_dataset(
	hid_t grp_id,
	const char *dset_name,
	size_t row_num,
	size_t col_num,
	double *data
	)
{
	hid_t dataspace_id, dset_id;
	hsize_t dims[2] = { row_num, col_num };
	dataspace_id = H5Screate_simple(2, dims, nullptr);
	dset_id = H5Dcreate(grp_id, dset_name, H5T_NATIVE_DOUBLE, dataspace_id,
						H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (dset_id < 0) return -1;
	herr_t res = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, dataspace_id, dataspace_id, H5P_DEFAULT, data);
	H5Dclose(dset_id);
	H5Sclose(dataspace_id);
	return res < 0 ? -2 : 0;
}

int ResultFile_hdf5::write_dataset(
	hid_t grp_id,
	const char *dset_name,
	size_t row_num,
	size_t col_num,
	unsigned long long *data
	)
{
	hid_t dataspace_id, dset_id;
	hsize_t dims[2] = { row_num, col_num };
	dataspace_id = H5Screate_simple(2, dims, nullptr);
	dset_id = H5Dcreate(grp_id, dset_name, H5T_NATIVE_ULLONG, dataspace_id,
						H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (dset_id < 0) return -1;
	herr_t res = H5Dwrite(dset_id, H5T_NATIVE_ULLONG, dataspace_id,
						  dataspace_id, H5P_DEFAULT, data);
	H5Dclose(dset_id);
	H5Sclose(dataspace_id);
	return res < 0 ? -2 : 0;
}


int ResultFile_hdf5::read_dataset(
	hid_t grp_id,
	const char *dset_name,
	size_t num,
	double *data
	)
{
	hid_t dataspace_id, dset_id;
	dataspace_id = H5Screate_simple(1, &num, nullptr);
	dset_id = H5Dopen(grp_id, dset_name, H5P_DEFAULT);
	if (dset_id < 0) return -1;
	herr_t res = H5Dread(dset_id, H5T_NATIVE_DOUBLE, dataspace_id,
						 dataspace_id, H5P_DEFAULT, data);
	H5Dclose(dset_id);
	H5Sclose(dataspace_id);
	return res < 0 ? -2 : 0;
}

int ResultFile_hdf5::read_dataset(
	hid_t grp_id,
	const char *dset_name,
	size_t row_num,
	size_t col_num,
	double *data
	)
{
	hid_t dataspace_id, dset_id;
	hsize_t dims[2] = { row_num, col_num };
	dataspace_id = H5Screate_simple(2, dims, nullptr);
	dset_id = H5Dopen(grp_id, dset_name, H5P_DEFAULT);
	if (dset_id < 0) return -1;
	herr_t res = H5Dread(dset_id, H5T_NATIVE_DOUBLE, dataspace_id,
						 dataspace_id, H5P_DEFAULT, data);
	H5Dclose(dset_id);
	H5Sclose(dataspace_id);
	return res < 0 ? -2 : 0;
}

int ResultFile_hdf5::read_dataset(
	hid_t grp_id, 
	const char *dset_name,
	size_t row_num,
	size_t col_num,
	unsigned long long *data
	)
{
	hid_t dataspace_id, dset_id;
	hsize_t dims[2] = { row_num, col_num };
	dataspace_id = H5Screate_simple(2, dims, nullptr);
	dset_id = H5Dopen(grp_id, dset_name, H5P_DEFAULT);
	if (dset_id < 0) return -1;
	herr_t res = H5Dread(dset_id, H5T_NATIVE_ULLONG, dataspace_id,
						 dataspace_id, H5P_DEFAULT, data);
	H5Dclose(dset_id);
	H5Sclose(dataspace_id);
	return res < 0 ? -2 : 0;
}


// self-defined datatype
int ResultFile_hdf5::write_dataset(
	hid_t grp_id,
	const char *dset_name,
	size_t num,
	void *data,
	hid_t datatype_id
	)
{
	hid_t dataspace_id, dset_id;
	dataspace_id = H5Screate_simple(1, &num, nullptr);
	dset_id = H5Dcreate(grp_id, dset_name, datatype_id, dataspace_id,
		H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (dset_id < 0) return -1;
	herr_t res = H5Dwrite(dset_id, datatype_id, dataspace_id,
		dataspace_id, H5P_DEFAULT, data);
	H5Dclose(dset_id);
	H5Sclose(dataspace_id);
	return res < 0 ? -2 : 0;
}

int ResultFile_hdf5::read_dataset(hid_t grp_id, const char *dset_name,
	size_t num, void *data, hid_t datatype_id)
{
	hid_t dataspace_id, dset_id;
	dataspace_id = H5Screate_simple(1, &num, nullptr);
	dset_id = H5Dopen(grp_id, dset_name, H5P_DEFAULT);
	if (dset_id < 0) return -1;
	herr_t res = H5Dread(dset_id, datatype_id, dataspace_id,
						 dataspace_id, H5P_DEFAULT, data);
	H5Dclose(dset_id);
	H5Sclose(dataspace_id);
	return res < 0 ? -2 : 0;
}

// ============================ routines for attributes ========================
int ResultFile_hdf5::write_attribute(
	hid_t grp_id,
	const char *name,
	double value
	)
{
	hid_t dataspace_id = H5Screate(H5S_SCALAR);
	hid_t attr_id = H5Acreate(grp_id, name, H5T_NATIVE_DOUBLE,
							  dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attr_id, H5T_NATIVE_DOUBLE, &value);
	H5Aclose(attr_id);
	H5Sclose(dataspace_id);
	return 0;
}

int ResultFile_hdf5::write_attribute(
	hid_t grp_id,
	const char *name,
	size_t value
	)
{
	unsigned long long _value = (unsigned long long)value;
	hid_t dataspace_id = H5Screate(H5S_SCALAR);
	hid_t attr_id = H5Acreate(grp_id, name, H5T_NATIVE_ULLONG,
							  dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attr_id, H5T_NATIVE_ULLONG, &_value);
	H5Aclose(attr_id);
	H5Sclose(dataspace_id);
	return 0;
}

int ResultFile_hdf5::write_attribute(
	hid_t grp_id, 
	const char *name,
	size_t num,
	const char *str
	)
{
	hid_t dataspace_id = H5Screate_simple(1, &num, nullptr);
	hid_t attr_id = H5Acreate(grp_id, name, H5T_NATIVE_CHAR,
							  dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
	H5Awrite(attr_id, H5T_NATIVE_CHAR, str);
	H5Aclose(attr_id);
	H5Sclose(dataspace_id);
	return 0;
}

int ResultFile_hdf5::read_attribute(
	hid_t grp_id,
	const char *name,
	double &value
	)
{
	hid_t dataspace_id = H5Screate(H5S_SCALAR);
	hid_t attr_id = H5Aopen(grp_id, name, H5P_DEFAULT);
	if (attr_id < 0) return -1;
	H5Aread(attr_id, H5T_NATIVE_DOUBLE, &value);
	H5Aclose(attr_id);
	H5Sclose(dataspace_id);
	return 0;
}

int ResultFile_hdf5::read_attribute(
	hid_t grp_id,
	const char *name,
	size_t &value
	)
{
	unsigned long long _value;
	hid_t dataspace_id = H5Screate(H5S_SCALAR);
	hid_t attr_id = H5Aopen(grp_id, name, H5P_DEFAULT);
	if (attr_id < 0) return -1;
	H5Aread(attr_id, H5T_NATIVE_ULLONG, &_value);
	H5Aclose(attr_id);
	H5Sclose(dataspace_id);
	value = (size_t)_value;
	return 0;
}

int ResultFile_hdf5::read_attribute(
	hid_t grp_id,
	const char *name,
	size_t num,
	const char *str
	)
{
	hid_t dataspace_id = H5Screate_simple(1, &num, nullptr);
	hid_t attr_id = H5Aopen(grp_id, name, H5P_DEFAULT);
	if (attr_id < 0) return -1;
	H5Awrite(attr_id, H5T_NATIVE_CHAR, str);
	H5Aclose(attr_id);
	H5Sclose(dataspace_id);
	return 0;
}
