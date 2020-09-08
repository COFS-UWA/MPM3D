#ifndef __Hdf5_Field_Extraction_mises_strain_2d_h__
#define __Hdf5_Field_Extraction_mises_strain_2d_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_mises_strain_2d : public Hdf5FieldExtraction
{
protected:
	bool all_fld_is_found;
	size_t e11_offset, e22_offset, e12_offset;

public:
	Hdf5FieldExtraction_mises_strain_2d() : all_fld_is_found(false) {}
	Hdf5FieldExtraction_mises_strain_2d(Hdf5DataLoader& loader) :
		all_fld_is_found(false), Hdf5FieldExtraction(loader) {}
	~Hdf5FieldExtraction_mises_strain_2d() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double* pcl_fld_data) override;
};

#endif