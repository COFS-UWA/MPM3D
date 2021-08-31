#ifndef __Hdf5_Field_Extraction_mises_strain_3d_h__
#define __Hdf5_Field_Extraction_mises_strain_3d_h__

#include "Hdf5FieldExtraction.h"

// mises / shear strain
class Hdf5FieldExtraction_mises_strain_3d : public Hdf5FieldExtraction
{
protected:
	bool all_fld_is_found;
	size_t e11_offset, e22_offset, e33_offset;
	size_t e12_offset, e23_offset, e31_offset;

public:
	Hdf5FieldExtraction_mises_strain_3d() : all_fld_is_found(false) {}
	Hdf5FieldExtraction_mises_strain_3d(Hdf5DataLoader& loader) :
		all_fld_is_found(false), Hdf5FieldExtraction(loader) {}
	~Hdf5FieldExtraction_mises_strain_3d() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double* pcl_fld_data) override;
};

#endif