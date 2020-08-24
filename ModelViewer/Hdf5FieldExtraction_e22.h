#ifndef __Hdf5_Field_Extraction_e22_h__
#define __Hdf5_Field_Extraction_e22_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_e22 : public Hdf5FieldExtraction
{
protected:
	bool e22_fld_is_found;
	size_t e22_offset;

public:
	Hdf5FieldExtraction_e22() : e22_fld_is_found(false) {}
	Hdf5FieldExtraction_e22(Hdf5DataLoader& loader) :
		e22_fld_is_found(false), Hdf5FieldExtraction(loader) {}
	~Hdf5FieldExtraction_e22() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double* pcl_fld_data) override;
};

#endif