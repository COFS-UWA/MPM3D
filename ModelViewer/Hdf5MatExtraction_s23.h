#ifndef __Hdf5_Mat_Extraction_s23_h__
#define __Hdf5_Mat_Extraction_s23_h__

#include "Hdf5FieldExtraction.h"

class Hdf5MatExtraction_s23 : public Hdf5FieldExtraction
{
protected:
	bool mat_id_is_found;
	size_t mat_id_offset;

public:
	Hdf5MatExtraction_s23() : mat_id_is_found(false) {}
	Hdf5MatExtraction_s23(Hdf5DataLoader& loader) :
		Hdf5FieldExtraction(loader), mat_id_is_found(false) {}
	~Hdf5MatExtraction_s23() {}

	bool validate_data_type() override;
	bool need_mat_model_data() override { return true; }
	int extract_pcl_fld_data(double* pcl_fld_data) override;
};

#endif