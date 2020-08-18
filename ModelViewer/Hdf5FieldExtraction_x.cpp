#include "ModelViewer_pcp.h"

#include "Hdf5FieldExtraction_x.h"

Hdf5FieldExtraction_x::Hdf5FieldExtraction_x()
{

}

Hdf5FieldExtraction_x::~Hdf5FieldExtraction_x()
{

}

bool Hdf5FieldExtraction_x::validate_data_type(hid_t dt)
{

	return true;
}

int Hdf5FieldExtraction_x::extract_data(double* fld_data)
{

	return 0;
}
