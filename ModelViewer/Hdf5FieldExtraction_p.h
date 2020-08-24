#ifndef __Hdf5_Field_Extraction_p_h__
#define __Hdf5_Field_Extraction_p_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_p : public Hdf5FieldExtraction
{
protected:
	bool p_fld_is_found;
	size_t p_offset;

public:
	Hdf5FieldExtraction_p() : p_fld_is_found(false) {}
	Hdf5FieldExtraction_p(Hdf5DataLoader& loader) :
		p_fld_is_found(false), Hdf5FieldExtraction(loader) {}
	~Hdf5FieldExtraction_p() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double* pcl_fld_data) override;
};

#endif