#ifndef __RESULT_FILE_XML_H__
#define __RESULT_FILE_XML_H__

#include "SimulationCore_pcp.h"

#include <fstream>
#include "ResultFile.h"

class ResultFile_XML : public ResultFile
{
protected:
	std::fstream file;

public:
	ResultFile_XML() : ResultFile(ResultFileType::XML) {}
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
	void finalize(void)
	{
		const char *file_ending = "</ResultFile>\n";
		file.write(file_ending, strlen(file_ending));
		file.close();
	}

	std::fstream &get_file(void) noexcept { return file; }
};

#endif