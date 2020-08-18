#ifndef __Hdf5_Data_Loader_h__
#define __Hdf5_Data_Loader_h__

#include "hdf5.h"

class Hdf5DataLoader
{
protected:


public:
	Hdf5DataLoader();
	~Hdf5DataLoader();

	void* get_field_data();
};

#endif