#ifndef __Hdf5_Field_Extraction_density_h__
#define __Hdf5_Field_Extraction_density_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_density : public Hdf5FieldExtraction
{
protected:
	bool density_fld_is_found;
	size_t density_offset;

public:
	Hdf5FieldExtraction_density() : density_fld_is_found(false) {}
	Hdf5FieldExtraction_density(Hdf5DataLoader& loader) :
		density_fld_is_found(false), Hdf5FieldExtraction(loader) {}
	~Hdf5FieldExtraction_density() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double* pcl_fld_data) override;
	int extract_pcl_fld_data_f(float* pcl_fld_data) override;
};

#endif