#ifndef __Hdf5_Field_Extraction_s22_h__
#define __Hdf5_Field_Extraction_s22_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_s22 : public Hdf5FieldExtraction
{
protected:
	bool s22_fld_is_found;
	size_t s22_offset;

public:
	Hdf5FieldExtraction_s22() : s22_fld_is_found(false) {}
	Hdf5FieldExtraction_s22(Hdf5DataLoader& loader) :
		s22_fld_is_found(false), Hdf5FieldExtraction(loader) {}
	~Hdf5FieldExtraction_s22() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double* pcl_fld_data) override;
};

#endif