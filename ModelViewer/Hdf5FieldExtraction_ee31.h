#ifndef __Hdf5_Field_Extraction_ee31_h__
#define __Hdf5_Field_Extraction_ee31_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_ee31 : public Hdf5FieldExtraction
{
protected:
	bool ee31_fld_is_found;
	size_t ee31_offset;

public:
	Hdf5FieldExtraction_ee31() : ee31_fld_is_found(false) {}
	Hdf5FieldExtraction_ee31(Hdf5DataLoader& loader) :
		ee31_fld_is_found(false), Hdf5FieldExtraction(loader) {}
	~Hdf5FieldExtraction_ee31() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double* pcl_fld_data) override;
};

#endif