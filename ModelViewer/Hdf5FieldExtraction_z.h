#ifndef __Hdf5_Field_Extraction_z_h__
#define __Hdf5_Field_Extraction_z_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_z : public Hdf5FieldExtraction
{
protected:
	bool z_fld_is_found;
	size_t z_offset;

public:
	Hdf5FieldExtraction_z() : z_fld_is_found(false) {}
	Hdf5FieldExtraction_z(Hdf5DataLoader& loader) :
		Hdf5FieldExtraction(loader), z_fld_is_found(false) {}
	~Hdf5FieldExtraction_z() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double* pcl_fld_data) override;
	int extract_pcl_fld_data_f(float *pcl_fld_data) override;
};

#endif