#ifndef __Time_History_h__
#define __Time_History_h__

#include <string>
#include "LinkList.hpp"

#include "ResultFile.h"

class Model;
class Step;
class TimeHistory;

typedef int (*TimeHistoryFunc)(TimeHistory &_self);
int time_history_output_func_null(TimeHistory &_self);

/*=============================================================
Class TimeHistory
 ==============================================================*/
class TimeHistory
{
	friend Step;
protected:
	std::string name;
	const char *type;

	size_t interval_num; // number of output this step:
	bool need_output_init_state; // true if output the initial state
	bool need_output_final_state; // true if output the final state

	Model *model;
	Step *step;
	ResultFile *res_file;

public:
	TimeHistory(const char *_name,
		const char *_type = "TimeHistory",
		TimeHistoryFunc _output_func = &time_history_output_func_null) :
		name(_name), type(_type), interval_num(1),
		need_output_init_state(false),
		need_output_final_state(false),
		model(nullptr), step(nullptr), res_file(nullptr),
		output_func(_output_func) {}
	~TimeHistory() {}
	
	inline void set_name(const char *_name) { name = _name; }
	inline void set_interval_num(size_t num) { interval_num = num; }
	// overrided by class without initial output
	virtual void set_output_init_state(bool _need = true) noexcept { need_output_init_state = _need; }
	virtual void set_output_final_state(bool _need = true) noexcept { need_output_final_state = _need; }

	inline const char *get_type() const { return type; }
	inline const char *get_name() const { return name.c_str(); }
	inline size_t get_interval_num() const { return interval_num; }
	
	inline Model &get_model() const noexcept { return *model; }
	inline Step &get_step() const noexcept { return *step; }
	inline ResultFile &get_res_file() const noexcept { return *res_file; }

public:
	// Initialize before each steps
	virtual int init_per_step() { return 0; }
	// Finalize after each steps
	virtual void finalize_per_step() {}
	// Output functions
	TimeHistoryFunc output_func;
	inline int output() { return (*output_func)(*this); }

protected: // used by step
	double interval_time;
	double next_time;
	LinkListPointer<TimeHistory> pointer_by_step;

private: // non-copyable
	TimeHistory(TimeHistory& other);
};

#endif