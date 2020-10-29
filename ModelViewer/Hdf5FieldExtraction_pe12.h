#ifndef __Hdf5_Field_Extraction_pe12_h__
#define __Hdf5_Field_Extraction_pe12_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_pe12 : public Hdf5FieldExtraction
{
protected:
	bool pe12_fld_is_found;
	size_t pe12_offset;

public:
	Hdf5FieldExtraction_pe12() : pe12_fld_is_found(false) {}
	Hdf5FieldExtraction_pe12(Hdf5DataLoader& loader) :
		pe12_fld_is_found(false), Hdf5FieldExtraction(loader) {}
	~Hdf5FieldExtraction_pe12() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double* pcl_fld_data) override;
};

#endif