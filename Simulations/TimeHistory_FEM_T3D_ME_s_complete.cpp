#include "Simulations_pcp.h"

#include "Model_FEM_T3D_ME_s.h"
#include "Step_FEM_T3D_ME_s.h"

#include "Model_FEM_T3D_ME_s_hdf5_utilities.h"

#include "TimeHistory_FEM_T3D_ME_s_complete.h"

int TimeHistory_FEM_T3D_ME_s_complete::init()
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

void TimeHistory_FEM_T3D_ME_s_complete::close()
{
	const char *res_file_type = res_file->get_type();
	if (!strcmp(res_file_type, "ResultFile_XML"))
	{

	}
	else if (!strcmp(res_file_type, "ResultFile_hdf5"))
	{
		ResultFile_hdf5 &rf = *static_cast<ResultFile_hdf5 *>(res_file);
		if (is_init)
		{
			rf.write_attribute(th_id, "output_num", output_id);
			rf.close_group(th_id);
		}
	}

	is_init = false;
}

int time_history_output_func_fem_t3d_me_s_complete_to_xml_res_file(TimeHistory &_self)
{
	TimeHistory_FEM_T3D_ME_s_complete &th = static_cast<TimeHistory_FEM_T3D_ME_s_complete &>(_self);
	ResultFile_XML &rf = static_cast<ResultFile_XML &>(th.get_res_file());
	std::fstream &file = rf.get_file();

	char str_buffer[512];
#define str_buffer_len (sizeof(str_buffer) / sizeof(str_buffer[0]))

	// time history
	Step_FEM_T3D_ME_s &step = static_cast<Step_FEM_T3D_ME_s &>(th.get_step());
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

	// ending
	file << "</TimeHistory>\n";

	return 0;
}

int time_history_output_func_fem_t3d_me_s_complete_to_hdf5_res_file(TimeHistory &_self)
{
	TimeHistory_FEM_T3D_ME_s_complete &th = static_cast<TimeHistory_FEM_T3D_ME_s_complete &>(_self);
	Step_FEM_T3D_ME_s &step = static_cast<Step_FEM_T3D_ME_s &>(th.get_step());
	Model_FEM_T3D_ME_s &md = static_cast<Model_FEM_T3D_ME_s &>(step.get_model());
	ResultFile_hdf5 &rf = static_cast<ResultFile_hdf5 &>(th.get_res_file());

	char frame_name[30];
	snprintf(frame_name, 30, "frame_%zu", th.output_id);
	hid_t frame_grp_id = rf.create_group(th.th_id, frame_name);
	
	rf.write_attribute(frame_grp_id, "substep_index", step.get_substep_index());
	rf.write_attribute(frame_grp_id, "total_substep_index", step.get_total_substep_index());
	rf.write_attribute(frame_grp_id, "current_time", step.get_current_time());
	rf.write_attribute(frame_grp_id, "total_time", step.get_total_time());
	
	using Model_FEM_T3D_ME_s_hdf5_utilities::time_history_complete_output_to_hdf5_file;
	time_history_complete_output_to_hdf5_file(md, rf, frame_grp_id);
	
	rf.close_group(frame_grp_id);

	++th.output_id;
	return 0;
}
