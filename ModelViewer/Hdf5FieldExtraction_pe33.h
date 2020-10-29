#ifndef __Hdf5_Field_Extraction_pe33_h__
#define __Hdf5_Field_Extraction_pe33_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_pe33 : public Hdf5FieldExtraction
{
protected:
	bool pe33_fld_is_found;
	size_t pe33_offset;

public:
	Hdf5FieldExtraction_pe33() : pe33_fld_is_found(false) {}
	Hdf5FieldExtraction_pe33(Hdf5DataLoader& loader) :
		pe33_fld_is_found(false), Hdf5FieldExtraction(loader) {}
	~Hdf5FieldExtraction_pe33() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double* pcl_fld_data) override;
};

#endif