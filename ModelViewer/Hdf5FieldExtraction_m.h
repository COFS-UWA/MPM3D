#ifndef __Hdf5_Field_Extraction_m_h__
#define __Hdf5_Field_Extraction_m_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_m : public Hdf5FieldExtraction
{
protected:
	bool m_fld_is_found;
	size_t m_offset;

public:
	Hdf5FieldExtraction_m() : m_fld_is_found(false) {}
	Hdf5FieldExtraction_m(Hdf5DataLoader& loader) :
		m_fld_is_found(false), Hdf5FieldExtraction(loader) {}
	~Hdf5FieldExtraction_m() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double* pcl_fld_data) override;
};

#endif