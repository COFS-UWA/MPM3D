#include "Tests_pcp.h"

#include "Step.h"
#include "TimeHistory.h"
#include "TimeHistory_ConsoleProgressBar.h"

#include "test_simulations.h"

#include <Windows.h> // for Sleep(ms) func

int solve_substep_test(void *_self);
class Step_test : public Step
{
public:
	friend int solve_substep_test(void *_self);
	Step_test() : Step("test", "Step_test", &solve_substep_test) {}
};
int solve_substep_test(void *_self)
{
	//Step_test &self = *reinterpret_cast<Step_test *>(_self);
	//std::cout << "Substep " << self.substep_index
	//		  << " time: " << self.current_time
	//		  << " dt: " << self.dtime << " ";
	//
	//for (TimeHistory * pth = self.time_history_list.first();
	//	self.time_history_list.is_not_end(pth);
	//	pth = self.time_history_list.next(pth))
	//{
	//	std::cout << pth->get_name() << " ";
	//}
	//std::cout << "\n";

	//Sleep(500);
	return 0;
}

int time_history_output_func_test(TimeHistory &_self);
class TimeHistory_test : public TimeHistory
{
	static size_t cur_id;
	size_t id;
public:
	friend int time_history_output_func_test(TimeHistory &_self);
	TimeHistory_test(const char *_name) :
		TimeHistory(_name, "TimeHistory_test", &time_history_output_func_test), id(cur_id++) {}
};
size_t TimeHistory_test::cur_id = 0;
int time_history_output_func_test(TimeHistory &_self)
{
	TimeHistory_test &self = static_cast<TimeHistory_test &>(_self);
	std::cout << "TimeHistory: " << self.name << " time: " << self.get_step().get_current_time() << "\n"; return 0;
}

void test_solve_functions()
{
	TimeHistory_test th1("th1");
	th1.set_interval_num(5);
	th1.set_output_init_state();

	TimeHistory_test th2("th2");
	th2.set_interval_num(2);
	
	TimeHistory_test th3("th3");
	th3.set_interval_num(3);
	th3.set_output_init_state();

	TimeHistory_ConsoleProgressBar cpb;

	Step_test step;
	step.set_step_time(1.0);
	step.set_dtime(0.1);

	//step.add_time_history(th1);
	//step.add_time_history(th2);
	//step.add_time_history(th3);
	step.add_time_history(cpb);

	step.solve();
}
