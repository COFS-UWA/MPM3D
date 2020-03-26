#ifndef __TIME_HISTORY_OUTPUT_H__
#define __TIME_HISTORY_OUTPUT_H__

#include <string>

#include "ResultFile.h"

class Model;
class Step;
class TimeHistoryOutput;
typedef int (*TimeHistoryOutputFunc)(TimeHistoryOutput &_self);
int time_history_output_func_null(TimeHistoryOutput &_self);

/*=============================================================
Class TimeHistory
 ==============================================================*/
class TimeHistoryOutput
{
	friend Step;
protected:
	const char *type;
	std::string name;
	size_t interval_num; // number of output this step:
	double interval_time;
	double next_time;
	bool need_output_init_state; // true if output the initial state
	
	Model *model;
	Step *step;
	ResultFile *res_file;
	TimeHistoryOutputFunc output_func;

public:
	TimeHistoryOutput(const char *_name,
		const char *_type = "TimeHistory",
		TimeHistoryOutputFunc _output_func = &time_history_output_func_null) :
		name(_name), type(_type),
		interval_num(1), need_output_init_state(false),
		model(nullptr), step(nullptr), res_file(nullptr),
		output_func(_output_func), next(nullptr) {}
	~TimeHistoryOutput() {}
	
	inline const char *get_type(void) const { return type; }
	inline void set_name(const char *_name) { name = _name; }
	inline const char *get_name(void) const { return name.c_str(); }
	inline void set_interval_num(size_t num) { interval_num = num; }
	inline size_t get_interval_num(void) const { return interval_num; }
	inline void set_output_init_state(bool _need = true) noexcept { need_output_init_state = _need; }
	
	inline Model &get_model(void) const noexcept { return *model; }
	inline Step &get_step(void) const noexcept { return *step; }
	inline ResultFile &get_res_file(void) const noexcept { return *res_file; }

public:
	// Initialize before each steps
	virtual int init_per_step(void) { return 0; }
	// Finalize after each steps
	virtual void finalize_per_step(void) {}
	// Output substep
	inline int output(void) { return (*output_func)(*this); }

protected:
	TimeHistoryOutput *next; // used by step
};

#endif