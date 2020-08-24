#ifndef __Hdf5_Field_Extraction_vy_f_h__
#define __Hdf5_Field_Extraction_vy_f_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_vy_f : public Hdf5FieldExtraction
{
protected:
	bool vy_f_fld_is_found;
	size_t vy_f_offset;

public:
	Hdf5FieldExtraction_vy_f() : vy_f_fld_is_found(false) {}
	Hdf5FieldExtraction_vy_f(Hdf5DataLoader& loader) :
		vy_f_fld_is_found(false), Hdf5FieldExtraction(loader) {}
	~Hdf5FieldExtraction_vy_f() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double* pcl_fld_data) override;
};

#endif