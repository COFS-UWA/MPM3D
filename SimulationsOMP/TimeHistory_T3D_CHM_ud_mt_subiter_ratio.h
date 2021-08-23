#ifndef __Time_History_T3D_CHM_ud_mt_subiter_ratio_h__
#define __Time_History_T3D_CHM_ud_mt_subiter_ratio_h__

#include <fstream>
#include "TimeHistory.h"

int time_history_output_func_T3D_CHM_ud_mt_subiter_ratio_to_text_res_file(TimeHistory &_self);

/* ===========================================================
Class TimeHistory_T3D_CHM_ud_mt_subiter_ratio
  =========================================================== */
class TimeHistory_T3D_CHM_ud_mt_subiter_ratio : public TimeHistory
{
protected:
	std::fstream res_file;
	size_t output_id;
	bool is_init;
	void close();

public:
	TimeHistory_T3D_CHM_ud_mt_subiter_ratio(const char *_name) :
		TimeHistory(_name, "TimeHistory_T3D_CHM_ud_mt_subiter_ratio"),
		output_id(0), is_init(false) {}
	~TimeHistory_T3D_CHM_ud_mt_subiter_ratio() { close(); }

	friend int time_history_output_func_T3D_CHM_ud_mt_subiter_ratio_to_text_res_file(TimeHistory &_self);
	inline void set_res_file(const char *filename) noexcept
	{
		if (is_init)
			close();
		is_init = true;
		res_file.open(filename, std::ios::binary | std::ios::out);
		res_file << "substep index, time, e_kin, max_e_kin, subiter index\n";
		output_func = &time_history_output_func_T3D_CHM_ud_mt_subiter_ratio_to_text_res_file;
	}
};

#endif