#ifndef __Time_History_T2D_ME_s_complete_H__
#define __Time_History_T2D_ME_s_complete_H__

#include "ResultFile_XML.h"
#include "ResultFile_hdf5.h"

#include "TimeHistory.h"

int time_history_output_func_t2d_me_s_to_xml_res_file(TimeHistory &_self);
int time_history_output_func_t2d_me_s_to_hdf5_res_file(TimeHistory &_self);

/* ===========================================================
Class TimeHistory_T2D_ME_s_complete
  =========================================================== */
class TimeHistory_T2D_ME_s_complete : public TimeHistory
{
protected:
	size_t output_id;
	bool is_init;
	int init();
	void close();

public:
	TimeHistory_T2D_ME_s_complete(const char *_name) :
		TimeHistory(_name, "TimeHistory_T2D_ME_s_complete"),
		output_id(0), is_init(false),
		th_id(-1) {}
	~TimeHistory_T2D_ME_s_complete() { close(); }

	int init_per_step() { return init(); }

	friend int time_history_output_func_t2d_me_s_to_xml_res_file(TimeHistory &_self);
	inline void set_res_file(ResultFile_XML &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &time_history_output_func_t2d_me_s_to_xml_res_file;
	}

protected:
	hid_t th_id;
public:
	friend int time_history_output_func_t2d_me_s_to_hdf5_res_file(TimeHistory &_self);
	inline void set_res_file(ResultFile_hdf5 &_res_file) noexcept
	{
		res_file = &_res_file;
		output_func = &time_history_output_func_t2d_me_s_to_hdf5_res_file;
	}
};

#endif