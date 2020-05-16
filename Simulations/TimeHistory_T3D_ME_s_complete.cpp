#include "Simulations_pcp.h"

#include "Model_T3D_ME_s.h"
#include "Step_T3D_ME_s.h"

#include "Model_T3D_ME_s_hdf5_utilities.h"

#include "TimeHistory_T3D_ME_s_complete.h"

int TimeHistory_T3D_ME_s_complete::init()
{
	if (is_init) return 0;
	
	const char *res_file_type = res_file->get_type();
	if (!strcmp(res_file_type, "ResultFile_XML"))
	{

	}
	else if (!strcmp(res_file_type, "ResultFile_hdf5"))
	{
		ResultFile_hdf5 &rf = *static_cast<ResultFile_hdf5 *>(res_file);
		hid_t th_grp_id = rf.get_time_history_grp_id();
		th_id = rf.create_group(th_grp_id, name.c_str());
	}
	else // unknown file type
	{
		return -1;
	}

	is_init = true;
	return 0;
}

void TimeHistory_T3D_ME_s_complete::close()
{
	const char *res_file_type = res_file->get_type();
	if (!strcmp(res_file_type, "ResultFile_XML"))
	{

	}
	else if (!strcmp(res_file_type, "ResultFile_hdf5"))
	{
		ResultFile_hdf5 &rf = *static_cast<ResultFile_hdf5 *>(res_file);
		rf.write_attribute(th_id, "output_num", output_id);
		if (th_id > 0)
		{
			rf.close_group(th_id);
			th_id = -1;
		}
	}

	is_init = false;
}

int time_history_output_func_t3d_me_s_to_xml_res_file(TimeHistory &_self)
{
	TimeHistory_T3D_ME_s_complete &th = static_cast<TimeHistory_T3D_ME_s_complete &>(_self);
	ResultFile_XML &rf = static_cast<ResultFile_XML &>(th.get_res_file());
	std::fstream &file = rf.get_file();

	char str_buffer[512];
#define str_buffer_len (sizeof(str_buffer) / sizeof(str_buffer[0]))

	// time history
	Step_T3D_ME_s &step = static_cast<Step_T3D_ME_s &>(th.get_step());
	const char *time_history_info = ""
		"<TimeHistory>\n"
		"    <substep_index> %zu </substep_index>\n"
		"    <total_substep_index> %zu </total_substep_index>\n"
		"    <current_time> %16.10e </current_time>\n"
		"    <total_time> %16.10e </total_time>\n";
	snprintf(str_buffer, str_buffer_len, time_history_info,
		step.get_substep_index(),  step.get_total_substep_index(),
		step.get_current_time(), step.get_total_time());
	file.write(str_buffer, strlen(str_buffer));

	// output material points data
	Model_T3D_ME_s &model = static_cast<Model_T3D_ME_s &>(th.get_model());
	const char *material_point_info = ""
		"    <MaterialPoints num=%zu>\n";
	snprintf(str_buffer, str_buffer_len, material_point_info, model.pcl_num);
	file.write(str_buffer, strlen(str_buffer));
	// field data: x, y, z, vol, s11, s22, s33, s12, s23, s31
	file << "        <field_data>\n"
			"        <!-- x, y, z, vol, s11, s22, s33, s12, s23, s31 -->\n";
	const char *field_data_format = "        %16.10e, %16.10e, %16.10e, %16.10e, %16.10e, %16.10e, %16.10e, %16.10e, %16.10e, %16.10e\n";
	for (size_t pcl_id = 0; pcl_id < model.pcl_num; ++pcl_id)
	{
		Model_T3D_ME_s::Particle &pcl = model.pcls[pcl_id];
		snprintf(str_buffer, str_buffer_len, field_data_format,
				 pcl.x, pcl.y, pcl.z, pcl.m / pcl.density,
				 pcl.s11, pcl.s22, pcl.s33, pcl.s12, pcl.s23, pcl.s31);
	}
	file << "        </field_data>\n"
			"    </MaterialPointObject>\n";

	// ending
	file << "</TimeHistory>\n";

	return 0;
}

int time_history_output_func_t3d_me_s_to_hdf5_res_file(TimeHistory &_self)
{
	TimeHistory_T3D_ME_s_complete &th = static_cast<TimeHistory_T3D_ME_s_complete &>(_self);
	Step_T3D_ME_s &step = static_cast<Step_T3D_ME_s &>(th.get_step());
	Model_T3D_ME_s &md = static_cast<Model_T3D_ME_s &>(step.get_model());
	ResultFile_hdf5 &rf = static_cast<ResultFile_hdf5 &>(*th.res_file);

	char frame_name[30];
	snprintf(frame_name, 30, "frame_%zu", th.output_id);
	hid_t frame_grp_id = rf.create_group(th.th_id, frame_name);
	
	rf.write_attribute(frame_grp_id, "substep_index", step.get_substep_index());
	rf.write_attribute(frame_grp_id, "total_substep_index", step.get_total_substep_index());
	rf.write_attribute(frame_grp_id, "current_time", step.get_current_time());
	rf.write_attribute(frame_grp_id, "total_time", step.get_total_time());
	
	// output particle data
	using Model_T3D_ME_s_hdf5_utilities::output_pcl_data_to_hdf5_file;
	output_pcl_data_to_hdf5_file(md, rf, frame_grp_id);
	// output consititutive model
	using Model_T3D_ME_s_hdf5_utilities::output_material_model_to_hdf5_file;
	output_material_model_to_hdf5_file(md, rf, frame_grp_id);
	
	rf.close_group(frame_grp_id);

	++th.output_id;
	return 0;
}
