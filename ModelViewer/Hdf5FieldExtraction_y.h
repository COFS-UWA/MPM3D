#ifndef __Hdf5_Field_Extraction_y_h__
#define __Hdf5_Field_Extraction_y_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_y : public Hdf5FieldExtraction
{
protected:
	bool y_fld_is_found;
	size_t y_offset;

public:
	Hdf5FieldExtraction_y() : y_fld_is_found(false) {}
	Hdf5FieldExtraction_y(Hdf5DataLoader& loader) :
		Hdf5FieldExtraction(loader), y_fld_is_found(false) {}
	~Hdf5FieldExtraction_y() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double* pcl_fld_data) override;
};

#endif