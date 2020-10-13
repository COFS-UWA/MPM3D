#ifndef __Hdf5_Field_Extraction_h__
#define __Hdf5_Field_Extraction_h__

#include "hdf5.h"
#include "Hdf5DataLoader.h"

class Hdf5FieldExtraction
{
protected:
	Hdf5DataLoader *data_loader;

public:
	Hdf5FieldExtraction() : data_loader(nullptr) {}
	Hdf5FieldExtraction(Hdf5DataLoader& loader) : data_loader(&loader) {}
	virtual ~Hdf5FieldExtraction() {}
	
	inline void set_data_loader(Hdf5DataLoader& loader) noexcept
	{ data_loader = &loader; }

	virtual bool validate_data_type() = 0;
	virtual int extract_pcl_fld_data(double *pcl_fld_data) = 0;
	virtual int extract_pcl_fld_data_f(float* pcl_fld_data) { return 0; }
};

#endif