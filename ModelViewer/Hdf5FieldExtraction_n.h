#ifndef __Hdf5_Field_Extraction_n_h__
#define __Hdf5_Field_Extraction_n_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_n : public Hdf5FieldExtraction
{
protected:
	bool n_fld_is_found;
	size_t n_offset;

public:
	Hdf5FieldExtraction_n() : n_fld_is_found(false) {}
	Hdf5FieldExtraction_n(Hdf5DataLoader& loader) :
		n_fld_is_found(false), Hdf5FieldExtraction(loader) {}
	~Hdf5FieldExtraction_n() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double* pcl_fld_data) override;
	//int extract_pcl_fld_data_f(float* pcl_fld_data) override;
};

#endif