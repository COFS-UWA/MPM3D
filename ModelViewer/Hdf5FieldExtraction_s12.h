#ifndef __Hdf5_Field_Extraction_s12_h__
#define __Hdf5_Field_Extraction_s12_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_s12 : public Hdf5FieldExtraction
{
protected:
	bool s12_fld_is_found;
	size_t s12_offset;

public:
	Hdf5FieldExtraction_s12() : s12_fld_is_found(false) {}
	Hdf5FieldExtraction_s12(Hdf5DataLoader& loader) :
		s12_fld_is_found(false), Hdf5FieldExtraction(loader) {}
	~Hdf5FieldExtraction_s12() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double* pcl_fld_data) override;
	int extract_pcl_fld_data_f(float* pcl_fld_data) override;
};

#endif