#ifndef __Hdf5_Field_Extraction_pe11_h__
#define __Hdf5_Field_Extraction_pe11_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_pe11 : public Hdf5FieldExtraction
{
protected:
	bool pe11_fld_is_found;
	size_t pe11_offset;

public:
	Hdf5FieldExtraction_pe11() : pe11_fld_is_found(false) {}
	Hdf5FieldExtraction_pe11(Hdf5DataLoader& loader) :
		pe11_fld_is_found(false), Hdf5FieldExtraction(loader) {}
	~Hdf5FieldExtraction_pe11() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double* pcl_fld_data) override;
};

#endif