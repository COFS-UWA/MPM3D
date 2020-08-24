#ifndef __Hdf5_Field_Extraction_m_f_h__
#define __Hdf5_Field_Extraction_m_f_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_m_f : public Hdf5FieldExtraction
{
protected:
	bool all_fld_is_found;
	size_t m_s_offset, density_s_offset, n_offset;
	size_t density_f_offset;

public:
	Hdf5FieldExtraction_m_f() : all_fld_is_found(false) {}
	Hdf5FieldExtraction_m_f(Hdf5DataLoader& loader) :
		all_fld_is_found(false), Hdf5FieldExtraction(loader) {}
	~Hdf5FieldExtraction_m_f() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double* pcl_fld_data) override;
};

#endif