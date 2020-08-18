#ifndef __Hdf5_Field_Extraction_x_h__
#define __Hdf5_Field_Extraction_x_h__

#include "Hdf5FieldExtraction.h";

class Hdf5FieldExtraction_x : public Hdf5FieldExtraction
{
public:
	Hdf5FieldExtraction_x();
	~Hdf5FieldExtraction_x();
	bool validate_data_type(hid_t dt);
	int extract_data(double *fld_data);
};

#endif