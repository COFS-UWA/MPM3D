#ifndef __Model_Data_T2D_ME_p_h__
#define __Model_Data_T2D_ME_p_h__

#include <string>

#include "ResultFile_XML.h"
#include "ResultFile_hdf5.h"

#include "ModelData.h"

int model_data_output_func_t2d_me_p_to_xml_res_file(ModelData &_self);
int model_data_output_func_t2d_me_p_to_hdf5_res_file(ModelData &_self);

/*=============================================================
Class ModelData_T2D_ME_p
==============================================================*/
class ModelData_T2D_ME_p : public ModelData
{
public:
	ModelData_T2D_ME_p(const char *_name = "ModelData") :
		ModelData(_name, "ModelData_T2D_ME_p") {}
	~ModelData_T2D_ME_p() {}

	friend int model_data_output_func_t2d_me_p_to_xml_res_file(ModelData &_self);
	inline void set_res_file(ResultFile_XML &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &model_data_output_func_t2d_me_p_to_xml_res_file;
	}

	friend int model_data_output_func_t2d_me_p_to_hdf5_res_file(ModelData &_self);
	inline void set_res_file(ResultFile_hdf5 &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &model_data_output_func_t2d_me_p_to_hdf5_res_file;
	}

public:
	template <typename ResultFileType>
	inline int output_model(Model& md, ResultFileType& rf)
	{
		set_model(md);
		set_res_file(rf);
		return output();
	}
};

#endif