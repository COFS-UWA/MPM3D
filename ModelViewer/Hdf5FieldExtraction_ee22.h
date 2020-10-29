#ifndef __Hdf5_Field_Extraction_ee22_h__
#define __Hdf5_Field_Extraction_ee22_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_ee22 : public Hdf5FieldExtraction
{
protected:
	bool ee22_fld_is_found;
	size_t ee22_offset;

public:
	Hdf5FieldExtraction_ee22() : ee22_fld_is_found(false) {}
	Hdf5FieldExtraction_ee22(Hdf5DataLoader& loader) :
		ee22_fld_is_found(false), Hdf5FieldExtraction(loader) {}
	~Hdf5FieldExtraction_ee22() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double* pcl_fld_data) override;
};

#endif