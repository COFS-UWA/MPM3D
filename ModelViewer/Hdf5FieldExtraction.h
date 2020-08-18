#ifndef __Hdf5_Field_Extraction_h__
#define __Hdf5_Field_Extraction_h__

#include "hdf5.h"
#include "Hdf5DataLoader.h"

class Hdf5FieldExtraction
{
public:
	Hdf5FieldExtraction() {}
	virtual ~Hdf5FieldExtraction() {}
	virtual bool validate_data_type(hid_t dt) = 0;
	virtual int extract_data(void *hdf5_data, double *fld_data) = 0;
};

#endif