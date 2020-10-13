#ifndef __Hdf5_Field_Extraction_s33_h__
#define __Hdf5_Field_Extraction_s33_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_s33 : public Hdf5FieldExtraction
{
protected:
	bool s33_fld_is_found;
	size_t s33_offset;

public:
	Hdf5FieldExtraction_s33() : s33_fld_is_found(false) {}
	Hdf5FieldExtraction_s33(Hdf5DataLoader& loader) :
		s33_fld_is_found(false), Hdf5FieldExtraction(loader) {}
	~Hdf5FieldExtraction_s33() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double* pcl_fld_data) override;
	int extract_pcl_fld_data_f(float* pcl_fld_data) override;
};

#endif