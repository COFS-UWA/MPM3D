#ifndef __Hdf5_Mat_Extraction_e_h__
#define __Hdf5_Mat_Extraction_e_h__

#include "Hdf5FieldExtraction.h"

// extract e (void ratio) from constitutive model
class Hdf5MatExtraction_e : public Hdf5FieldExtraction
{
protected:
	bool var_is_found;
	size_t mat_id_offset;
	size_t e_offset;

public:
	Hdf5MatExtraction_e() : var_is_found(false) {}
	Hdf5MatExtraction_e(Hdf5DataLoader& loader) :
		Hdf5FieldExtraction(loader), var_is_found(false) {}
	~Hdf5MatExtraction_e() {}

	bool validate_data_type() override;
	bool need_mat_model_data() override { return true; }
	int extract_pcl_fld_data(double* pcl_fld_data) override;
};

#endif