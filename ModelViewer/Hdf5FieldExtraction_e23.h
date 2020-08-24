#ifndef __Hdf5_Field_Extraction_e23_h__
#define __Hdf5_Field_Extraction_e23_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_e23 : public Hdf5FieldExtraction
{
protected:
	bool e23_fld_is_found;
	size_t e23_offset;

public:
	Hdf5FieldExtraction_e23() : e23_fld_is_found(false) {}
	Hdf5FieldExtraction_e23(Hdf5DataLoader& loader) :
		e23_fld_is_found(false), Hdf5FieldExtraction(loader) {}
	~Hdf5FieldExtraction_e23() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double* pcl_fld_data) override;
};

#endif