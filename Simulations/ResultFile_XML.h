#ifndef __Result_File_XML_h__
#define __Result_File_XML_h__

#include <fstream>
#include "ResultFile.h"

class ResultFile_XML : public ResultFile
{
protected:
	std::fstream file;

public:
	ResultFile_XML() : ResultFile("ResultFile_XML") {}
	~ResultFile_XML() { finalize(); }
	int init(const char *file_name)
	{
		file.open(file_name, std::ios::binary | std::ios::out);
		if (!file.is_open())
			return -1;

		// xml version info
		const char *file_header = "<?xml version=\"1.0\" encoding=\"ascii\"?>\n"
			"<ResultFile>\n";
		file.write(file_header, strlen(file_header));

		return 0;
	}
	void finalize()
	{
		const char *file_ending = "</ResultFile>\n";
		file.write(file_ending, strlen(file_ending));
		file.close();
	}

	std::fstream &get_file() noexcept { return file; }
};

#endif