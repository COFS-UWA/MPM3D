#ifndef __Hdf5_Field_Extraction_u_cav_h__
#define __Hdf5_Field_Extraction_u_cav_h__

#include "Hdf5FieldExtraction.h"

class Hdf5FieldExtraction_u_cav : public Hdf5FieldExtraction
{
protected:
	bool u_cav_fld_is_found;
	size_t u_cav_offset;

public:
	Hdf5FieldExtraction_u_cav() : u_cav_fld_is_found(false) {}
	Hdf5FieldExtraction_u_cav(Hdf5DataLoader& loader) :
		Hdf5FieldExtraction(loader), u_cav_fld_is_found(false) {}
	~Hdf5FieldExtraction_u_cav() {}

	bool validate_data_type() override;
	int extract_pcl_fld_data(double *pcl_fld_data) override;
	int extract_pcl_fld_data_f(float *pcl_fld_data) override;
};

#endif