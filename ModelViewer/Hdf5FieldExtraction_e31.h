#ifndef __Hdf5_Field_Extraction_e31_h__
#define __Hdf5_Field_Extraction_e31_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_e31 : public Hdf5FieldExtraction
{
protected:
	bool e31_fld_is_found;
	size_t e31_offset;

public:
	Hdf5FieldExtraction_e31() : e31_fld_is_found(false) {}
	Hdf5FieldExtraction_e31(Hdf5DataLoader& loader) :
		e31_fld_is_found(false), Hdf5FieldExtraction(loader) {}
	~Hdf5FieldExtraction_e31() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double* pcl_fld_data) override;
};

#endif