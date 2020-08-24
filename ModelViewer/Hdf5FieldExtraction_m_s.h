#ifndef __Hdf5_Field_Extraction_m_s_h__
#define __Hdf5_Field_Extraction_m_s_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_m_s : public Hdf5FieldExtraction
{
protected:
	bool m_s_fld_is_found;
	size_t m_s_offset;

public:
	Hdf5FieldExtraction_m_s() : m_s_fld_is_found(false) {}
	Hdf5FieldExtraction_m_s(Hdf5DataLoader& loader) :
		m_s_fld_is_found(false), Hdf5FieldExtraction(loader) {}
	~Hdf5FieldExtraction_m_s() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double* pcl_fld_data) override;
};

#endif