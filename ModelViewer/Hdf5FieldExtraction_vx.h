#ifndef __Hdf5_Field_Extraction_vx_h__
#define __Hdf5_Field_Extraction_vx_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_vx : public Hdf5FieldExtraction
{
protected:
	bool vx_fld_is_found;
	size_t vx_offset;

public:
	Hdf5FieldExtraction_vx() : vx_fld_is_found(false) {}
	Hdf5FieldExtraction_vx(Hdf5DataLoader& loader) :
		vx_fld_is_found(false), Hdf5FieldExtraction(loader) {}
	~Hdf5FieldExtraction_vx() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double* pcl_fld_data) override;
};

#endif