#ifndef __Hdf5_Field_Extraction_pe31_h__
#define __Hdf5_Field_Extraction_pe31_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_pe31 : public Hdf5FieldExtraction
{
protected:
	bool pe31_fld_is_found;
	size_t pe31_offset;

public:
	Hdf5FieldExtraction_pe31() : pe31_fld_is_found(false) {}
	Hdf5FieldExtraction_pe31(Hdf5DataLoader& loader) :
		pe31_fld_is_found(false), Hdf5FieldExtraction(loader) {}
	~Hdf5FieldExtraction_pe31() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double* pcl_fld_data) override;
};

#endif