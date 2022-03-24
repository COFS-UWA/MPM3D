#ifndef __Hdf5_Field_Extraction_is_cavitated_h__
#define __Hdf5_Field_Extraction_is_cavitated_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_is_cavitated : public Hdf5FieldExtraction
{
protected:
	bool is_cavitated_fld_is_found;
	size_t is_cavitated_offset;

public:
	Hdf5FieldExtraction_is_cavitated() : is_cavitated_fld_is_found(false) {}
	Hdf5FieldExtraction_is_cavitated(Hdf5DataLoader& loader) :
		Hdf5FieldExtraction(loader), is_cavitated_fld_is_found(false) {}
	~Hdf5FieldExtraction_is_cavitated() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double *pcl_fld_data) override;
	int extract_pcl_fld_data_f(float *pcl_fld_data) override;
};

#endif