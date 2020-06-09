#ifndef __Model_Data_T3D_CHM_s_H__
#define __Model_Data_T3D_CHM_s_H__

#include <string>

#include "ModelData.h"

#include "ResultFile_XML.h"
#include "ResultFile_hdf5.h"

int model_data_output_func_t3d_chm_s_to_xml_res_file(ModelData &_self);
int model_data_output_func_t3d_chm_s_to_hdf5_res_file(ModelData &_self);

/*=============================================================
Class ModelData_T3D_CHM_s
==============================================================*/
class ModelData_T3D_CHM_s : public ModelData
{
public:
	ModelData_T3D_CHM_s(const char *_name = "ModelData") :
		ModelData(_name, "ModelData_T3D_CHM_s") {}
	~ModelData_T3D_CHM_s() {}

	friend int model_data_output_func_t3d_chm_s_to_xml_res_file(ModelData &_self);
	inline void set_res_file(ResultFile_XML &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &model_data_output_func_t3d_chm_s_to_xml_res_file;
	}

	friend int model_data_output_func_t3d_chm_s_to_hdf5_res_file(ModelData &_self);
	inline void set_res_file(ResultFile_hdf5 &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &model_data_output_func_t3d_chm_s_to_hdf5_res_file;
	}

public:
	template <typename ResultFileType>
	inline int output_model(Model &md, ResultFileType &rf)
	{
		set_model(md);
		set_res_file(rf);
		return output();
	}
};

#endif