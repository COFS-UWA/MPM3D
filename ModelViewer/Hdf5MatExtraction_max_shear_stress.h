#ifndef __Hdf5_Mat_Extraction_max_shear_stress_h__
#define __Hdf5_Mat_Extraction_max_shear_stress_h__

#include "Hdf5FieldExtraction.h"

class Hdf5MatExtraction_max_shear_stress :
	public Hdf5FieldExtraction
{
protected:
	bool mat_id_is_found;
	size_t mat_id_offset;

public:
	Hdf5MatExtraction_max_shear_stress() :
		mat_id_is_found(false) {}
	Hdf5MatExtraction_max_shear_stress(Hdf5DataLoader& loader) :
		Hdf5FieldExtraction(loader), mat_id_is_found(false) {}
	~Hdf5MatExtraction_max_shear_stress() {}

	bool validate_data_type() override;
	bool need_mat_model_data() override { return true; }
	int extract_pcl_fld_data(double* pcl_fld_data) override;
};

#endif