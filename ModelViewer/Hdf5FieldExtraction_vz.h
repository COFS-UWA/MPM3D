#ifndef __Hdf5_Field_Extraction_vz_h__
#define __Hdf5_Field_Extraction_vz_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_vz : public Hdf5FieldExtraction
{
protected:
	bool vz_fld_is_found;
	size_t vz_offset;

public:
	Hdf5FieldExtraction_vz() : vz_fld_is_found(false) {}
	Hdf5FieldExtraction_vz(Hdf5DataLoader& loader) :
		vz_fld_is_found(false), Hdf5FieldExtraction(loader) {}
	~Hdf5FieldExtraction_vz() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double* pcl_fld_data) override;
};

#endif