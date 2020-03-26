#ifndef __Result_File_Type_H__
#define __Result_File_Type_H__

enum class ResultFileType : unsigned short
{
	Invalid = 0,
	PlainBin = 1,
	XML = 2,
	Hdf5 = 3
};

#endif