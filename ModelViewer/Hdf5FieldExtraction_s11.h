#ifndef __Hdf5_Field_Extraction_s11_h__
#define __Hdf5_Field_Extraction_s11_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_s11 : public Hdf5FieldExtraction
{
protected:
	bool s11_fld_is_found;
	size_t s11_offset;

public:
	Hdf5FieldExtraction_s11() : s11_fld_is_found(false) {}
	Hdf5FieldExtraction_s11(Hdf5DataLoader& loader) :
		s11_fld_is_found(false), Hdf5FieldExtraction(loader) {}
	~Hdf5FieldExtraction_s11() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double* pcl_fld_data) override;
};

#endif