#ifndef __Hdf5_Field_Extraction_vx_f_h__
#define __Hdf5_Field_Extraction_vx_f_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_vx_f : public Hdf5FieldExtraction
{
protected:
	bool vx_f_fld_is_found;
	size_t vx_f_offset;

public:
	Hdf5FieldExtraction_vx_f() : vx_f_fld_is_found(false) {}
	Hdf5FieldExtraction_vx_f(Hdf5DataLoader& loader) :
		vx_f_fld_is_found(false), Hdf5FieldExtraction(loader) {}
	~Hdf5FieldExtraction_vx_f() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double* pcl_fld_data) override;
};

#endif