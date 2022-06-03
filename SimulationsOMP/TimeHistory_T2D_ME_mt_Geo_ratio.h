#ifndef __Time_History_T2D_ME_mt_Geo_ratio_h__
#define __Time_History_T2D_ME_mt_Geo_ratio_h__

#include <fstream>
#include "TimeHistory.h"

int time_history_output_func_T2D_ME_mt_geo_ratio_to_text_res_file(TimeHistory &_self);

/* ===========================================================
Class TimeHistory_T2D_ME_mt_Geo_ratio
  =========================================================== */
class TimeHistory_T2D_ME_mt_Geo_ratio : public TimeHistory
{
protected:
	std::fstream res_file;
	size_t output_id;
	bool is_init;
	void close();

public:
	TimeHistory_T2D_ME_mt_Geo_ratio(const char *_name) :
		TimeHistory(_name, "TimeHistory_T2D_ME_mt_Geo_ratio"),
		output_id(0), is_init(false) {}
	~TimeHistory_T2D_ME_mt_Geo_ratio() { close(); }

	friend int time_history_output_func_T2D_ME_mt_geo_ratio_to_text_res_file(TimeHistory &_self);
	inline void set_res_file(const char *filename) noexcept
	{
		if (is_init)
			close();
		is_init = true;
		res_file.open(filename, std::ios::binary | std::ios::out);
		res_file << "substep index, time, f_ub, e_kin, f_ub_ratio, e_kin_ratio\n";
		output_func = &time_history_output_func_T2D_ME_mt_geo_ratio_to_text_res_file;
	}
};

#endif