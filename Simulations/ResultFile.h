#ifndef __Result_File_H__
#define __Result_File_H__

#include "ResultFileType.h"

class ResultFile
{
protected:
	ResultFileType type;

public:
	ResultFile(ResultFileType _tp) : type(_tp) {}
	~ResultFile() {}
	inline ResultFileType get_type(void) { return type; }
};

#endif