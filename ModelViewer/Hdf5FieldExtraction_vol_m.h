#ifndef __Hdf5_Field_Extraction_vol_m_h__
#define __Hdf5_Field_Extraction_vol_m_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_vol_m : public Hdf5FieldExtraction
{
protected:
	bool all_fld_is_found;
	size_t m_s_offset, density_s_offset, n_offset;

public:
	Hdf5FieldExtraction_vol_m() : all_fld_is_found(false) {}
	Hdf5FieldExtraction_vol_m(Hdf5DataLoader& loader) :
		all_fld_is_found(false), Hdf5FieldExtraction(loader) {}
	~Hdf5FieldExtraction_vol_m() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double* pcl_fld_data) override;
};

#endif