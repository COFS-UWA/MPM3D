#ifndef __Hdf5_Field_Extraction_s31_h__
#define __Hdf5_Field_Extraction_s31_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_s31 : public Hdf5FieldExtraction
{
protected:
	bool s31_fld_is_found;
	size_t s31_offset;

public:
	Hdf5FieldExtraction_s31() : s31_fld_is_found(false) {}
	Hdf5FieldExtraction_s31(Hdf5DataLoader& loader) :
		s31_fld_is_found(false), Hdf5FieldExtraction(loader) {}
	~Hdf5FieldExtraction_s31() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double* pcl_fld_data) override;
};

#endif