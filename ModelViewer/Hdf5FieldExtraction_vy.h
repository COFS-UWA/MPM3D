#ifndef __Hdf5_Field_Extraction_vy_h__
#define __Hdf5_Field_Extraction_vy_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_vy : public Hdf5FieldExtraction
{
protected:
	bool vy_fld_is_found;
	size_t vy_offset;

public:
	Hdf5FieldExtraction_vy() : vy_fld_is_found(false) {}
	Hdf5FieldExtraction_vy(Hdf5DataLoader& loader) :
		vy_fld_is_found(false), Hdf5FieldExtraction(loader) {}
	~Hdf5FieldExtraction_vy() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double* pcl_fld_data) override;
};

#endif