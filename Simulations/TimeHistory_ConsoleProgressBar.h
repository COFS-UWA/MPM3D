#ifndef __Time_History_Output_Func_Console_Progress_Bar_h__
#define __Time_History_Output_Func_Console_Progress_Bar_h__

#include <chrono>

#include "TimeHistory.h"

int time_history_output_func_console_progress_bar(TimeHistory &_self);

/* ===========================================================
Class TimeHistory_ConsoleProgressBar
=========================================================== */
class TimeHistory_ConsoleProgressBar : public TimeHistory
{
protected:
	int width; // width must <= 200
	size_t prev_pos, cur_pos;
	float width_div_100;
	
	std::chrono::system_clock::time_point start_time;
	std::chrono::system_clock::time_point end_time;
	
public:
	TimeHistory_ConsoleProgressBar();
	~TimeHistory_ConsoleProgressBar();

	// Initialize each steps
	int init_per_step(void) override;
	// output function
	friend int time_history_output_func_console_progress_bar(TimeHistory &_self);
	// Finalize each steps
	void finalize_per_step(void) override;

	void set_output_init_state(bool _need = true) noexcept override {}

	inline void set_width(int wd) noexcept
	{
		if (wd > 0)
		{
			width = wd < 200 ? wd : 200;
			width_div_100 = (float)width / 100.0f;
		}
	}

protected:
	void print_progress(void);
};


#endif