#ifndef __Step_h__
#define __Step_h__

#include <string>
#include "LinkList.hpp"

#include "Model.h"
#include "TimeHistory.h"

typedef int (*CalSubstepFunc)(void *_self);
int cal_substep_base(void *_self);

/*===========================================================
Class Step:
	Functions needs to be rewritten in children classes:
	1. init()
	2. solve_substep()
	3. finalize()
 ===========================================================*/
class Step
{
protected:
	std::string name;
	const char *type;

	// model
	Model *model;

	// whether this is the first calculaton step
	bool is_first_step;

	// substep index at start of this step
	size_t prev_substep_num;
	// start time of this step
	double start_time;

	// time length of this step
	double step_time;

	// substep from the start of this step
	size_t substep_index;
	// time from the start of this step
	double current_time;
	double current_time_plus_tol;

	// time increment
	double dtime; // time increment
	double time_tol_ratio;
	double time_tol; // = dt * time_tol_ratio

public:
	Step(const char *_name, const char *_type = "Step",
		CalSubstepFunc _cal_substep_func = &cal_substep_base);
	~Step();

	inline void set_name(const char *_name) noexcept { name = _name; }
	inline void set_step_time(double _time) noexcept { step_time = _time; }
	inline void set_dtime(double _dtime, double t_tol_r = 1.0e-3) noexcept
	{
		dtime = _dtime;
		time_tol_ratio = t_tol_r;
		time_tol = dtime * t_tol_r;
	}
	inline void set_model(Model &md) noexcept { model = &md; }
	// continuate from prev step
	void set_prev_step(Step &prev_step)
	{
		model = prev_step.model;
		is_first_step = false;
		prev_substep_num = prev_step.get_total_substep_index();
		start_time = prev_step.get_total_time();
	}

	inline const char *get_type() const { return type; }
	inline const char *get_name() const noexcept { return name.c_str(); }
	inline Model &get_model() const noexcept { return *model; }

	// time length of this step
	inline double get_step_time() const noexcept { return step_time; }
	// time from the start of this step
	inline double get_current_time() const noexcept { return current_time; }
	// total time from the start of the whole simulation
	inline double get_total_time() const noexcept { return start_time + current_time; }
	// number of substep from the start of this step
	inline size_t get_substep_index() const noexcept { return substep_index; }
	// total number of substep from the start of the whole simulation
	inline size_t get_total_substep_index() const noexcept { return prev_substep_num + substep_index; }
	// size of time increment
	inline double get_dtime() const noexcept { return dtime; }

protected:
	// initialization before calculation
	virtual int init_calculation() { return 0; }
	// finalization after calculation
	virtual int finalize_calculation() { return 0; }
	// calculation of each substep
	CalSubstepFunc cal_substep_func;
	inline int cal_substep() { return (*cal_substep_func)(this); }

public: // solve this step
	virtual int solve();

	// =============== Time History ================
protected:
	typedef LinkList<TimeHistory, offsetof(TimeHistory, pointer_by_step)> TimeHistoryList;
	TimeHistoryList time_history_list;
	int init_time_history();
	int output_time_history();
	int finalize_time_history();

protected: // helper data and functions (called after init "next_time")
	double next_output_time;
	// get output with the latest "next_time"
	TimeHistory *get_latest_time_history(TimeHistoryList &th_list);
	// insert time_history in ascending order
	void insert_time_history_in_ascending_order(TimeHistoryList &th_list, TimeHistory *th);

public:
	void add_time_history(TimeHistory &th)
	{
		time_history_list.append(th);
		th.model = model;
		th.step = this;
	}
	inline void clear_time_history() { time_history_list.reset(); }

private: // non-copyable
	Step(Step &other);
};

// utilities functions
inline double get_sign(double var)
{
	if (var > 0.0)
		return 1.0;
	else if (var < 0.0)
		return -1.0;
	return 0.0;
}

#endif