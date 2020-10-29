#ifndef __Hdf5_Field_Extraction_ee33_h__
#define __Hdf5_Field_Extraction_ee33_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_ee33 : public Hdf5FieldExtraction
{
protected:
	bool ee33_fld_is_found;
	size_t ee33_offset;

public:
	Hdf5FieldExtraction_ee33() : ee33_fld_is_found(false) {}
	Hdf5FieldExtraction_ee33(Hdf5DataLoader& loader) :
		ee33_fld_is_found(false), Hdf5FieldExtraction(loader) {}
	~Hdf5FieldExtraction_ee33() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double* pcl_fld_data) override;
};

#endif