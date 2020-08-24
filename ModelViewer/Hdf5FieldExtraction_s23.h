#ifndef __Hdf5_Field_Extraction_s23_h__
#define __Hdf5_Field_Extraction_s23_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_s23 : public Hdf5FieldExtraction
{
protected:
	bool s23_fld_is_found;
	size_t s23_offset;

public:
	Hdf5FieldExtraction_s23() : s23_fld_is_found(false) {}
	Hdf5FieldExtraction_s23(Hdf5DataLoader& loader) :
		s23_fld_is_found(false), Hdf5FieldExtraction(loader) {}
	~Hdf5FieldExtraction_s23() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double* pcl_fld_data) override;
};

#endif