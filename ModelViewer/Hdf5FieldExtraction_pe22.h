#ifndef __Hdf5_Field_Extraction_pe22_h__
#define __Hdf5_Field_Extraction_pe22_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_pe22 : public Hdf5FieldExtraction
{
protected:
	bool pe22_fld_is_found;
	size_t pe22_offset;

public:
	Hdf5FieldExtraction_pe22() : pe22_fld_is_found(false) {}
	Hdf5FieldExtraction_pe22(Hdf5DataLoader& loader) :
		pe22_fld_is_found(false), Hdf5FieldExtraction(loader) {}
	~Hdf5FieldExtraction_pe22() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double* pcl_fld_data) override;
};

#endif