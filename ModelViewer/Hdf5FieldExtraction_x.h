#ifndef __Hdf5_Field_Extraction_x_h__
#define __Hdf5_Field_Extraction_x_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_x : public Hdf5FieldExtraction
{
protected:
	bool x_fld_is_found;
	size_t x_offset;

public:
	Hdf5FieldExtraction_x() : x_fld_is_found(false) {}
	Hdf5FieldExtraction_x(Hdf5DataLoader& loader) :
		x_fld_is_found(false), Hdf5FieldExtraction(loader) {}
	~Hdf5FieldExtraction_x() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double* pcl_fld_data) override;
};

#endif