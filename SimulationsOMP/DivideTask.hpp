#ifndef __Divide_Task_hpp__
#define __Divide_Task_hpp__

#include <assert.h>
#include <mutex>

#include "tbb/task_arena.h"
#include "tbb/task.h"

/* ========================================
Assumptions:
	class Work
	{
		void operator() (size_t task_id);
	};
  =========================================*/

std::mutex cout_lock2;

namespace DivideTaskUtils
{
	template <class Work>
	class WorkTask : public tbb::task
	{
	protected:
		const size_t tk_id;
		Work& work;
	public:
		WorkTask(size_t id, Work& wk) : tk_id(id), work(wk) {}
		~WorkTask() {}
		tbb::task* execute() override { work(tk_id); return nullptr; }
	};
	
	// (tk_id1 - tk_id0) < div_num
	template <class Work, size_t div_num = 2>
	class DivideTask2 : public tbb::task
	{
	protected:
		const size_t tk_id0, tk_id1;
		Work& work;
	public:
		DivideTask2(size_t start_id, size_t end_id, Work &wk) :
			tk_id0(start_id), tk_id1(end_id), work(wk) { assert(div_num > 1); }
		~DivideTask2() {}
		tbb::task* execute() override
		{
			const size_t tk_num = tk_id1 - tk_id0;
			if (tk_num > 1)
			{
				tbb::empty_task& c = *new(allocate_continuation()) tbb::empty_task;
				c.set_ref_count(tk_num - 1);
				for (size_t div_id = tk_id0 + 1; div_id < tk_id1; ++div_id)
					c.spawn(*new(c.allocate_child()) WorkTask<Work>(div_id, work));
			}
			work(tk_id0);
			return nullptr;
		}
	};
}

template <class Work, size_t div_num>
class DivideTask : public tbb::task
{
protected:
	size_t tk_id0, tk_id1;
	Work& work;
public:
	DivideTask(size_t start_id, size_t end_id, Work& wk) :
		tk_id0(start_id), tk_id1(end_id), work(wk)
	{ assert(div_num > 1); }
	~DivideTask() {}
	tbb::task* execute() override
	{
		size_t div_id, blk_start_id, blk_end_id;
		const size_t tk_num = tk_id1 - tk_id0;
		if (tk_num > (div_num * div_num))
		{
			tbb::empty_task &c = *new(allocate_continuation()) tbb::empty_task;
			c.set_ref_count(div_num);
			blk_start_id = tk_id0;
			for (div_id = 1; div_id < div_num; ++div_id)
			{
				blk_end_id = tk_id0 + div_id * tk_num / div_num;
				c.spawn(*new(c.allocate_child()) DivideTask<Work, div_num>(blk_start_id, blk_end_id, work));
				blk_start_id = blk_end_id;
			}
			recycle_as_child_of(c);
			tk_id0 = blk_start_id;
			return this;
		}
		else if (tk_num > div_num)
		{
			tbb::empty_task& c = *new(allocate_continuation()) tbb::empty_task;
			const size_t blk_num = (tk_num + div_num - 1) / div_num;
			c.set_ref_count(blk_num);
			blk_start_id = tk_id0;
			for (div_id = 1; div_id < blk_num; ++div_id)
			{
				blk_end_id = blk_start_id + div_num;
				c.spawn(*new(c.allocate_child())
					DivideTaskUtils::DivideTask2<Work, div_num>(
						blk_start_id,
						blk_end_id,
						work));
				blk_start_id = blk_end_id;
			}
			recycle_as_child_of(c);
			tk_id0 = blk_start_id;
			return this;
		}
		else if (tk_num > 1)
		{
			tbb::empty_task& c = *new(allocate_continuation()) tbb::empty_task;
			c.set_ref_count(tk_num - 1);
			for (div_id = tk_id0 + 1; div_id < tk_id1; ++div_id)
				c.spawn(*new(c.allocate_child()) DivideTaskUtils::WorkTask<Work>(div_id, work));
		}
		if (tk_num) // tk_num == 1
			work(tk_id0);
		return nullptr;
	}
};

template <class Work>
class DivideTask<Work, 2> : public tbb::task
{
protected:
	size_t tk_id0, tk_id1;
	Work &work;
public:
	DivideTask(size_t start_id, size_t end_id, Work &wk) :
		tk_id0(start_id), tk_id1(end_id), work(wk) {}
	~DivideTask() {}
	tbb::task* execute() override
	{
		const size_t tk_num = tk_id1 - tk_id0;
		if (tk_num > 2)
		{
			tbb::empty_task &c = *new(allocate_continuation()) tbb::empty_task;
			c.set_ref_count(2);
			const size_t mid_tk_id = tk_id0 + tk_num / 2;
			c.spawn(*new(c.allocate_child())
				DivideTask<Work, 2>(tk_id0, mid_tk_id, work));
			recycle_as_child_of(c);
			tk_id0 = mid_tk_id;
			return this;
		}
		else if (tk_num == 2)
		{
			tbb::empty_task &c = *new(allocate_continuation()) tbb::empty_task;
			c.set_ref_count(1);
			c.spawn(*new(c.allocate_child()) DivideTaskUtils::WorkTask<Work>(tk_id0 + 1, work));
		}
		if (tk_num) // tk_num == 1
			work(tk_id0);
		return nullptr;
	}
};

#endif