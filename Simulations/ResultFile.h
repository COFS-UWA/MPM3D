#ifndef __Result_File_h__
#define __Result_File_h__

class ResultFile
{
protected:
	const char *type;

public:
	ResultFile(const char *type_name = "ResultFile") : type(type_name) {}
	~ResultFile() {}
	inline const char *get_type() const { return type; }

private: // non-copyable
	ResultFile(ResultFile& other);
};

#endif