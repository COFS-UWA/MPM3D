#ifndef __Result_File_h__
#define __Result_File_h__

#include "ResultFileType.h"

class ResultFile
{
protected:
	const char *type;

public:
	ResultFile(const char *type_name = "ResultFile") : type(type_name) {}
	~ResultFile() {}
	inline const char *get_type(void) const { return type; }
};

#endif