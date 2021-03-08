#ifndef __Merge_Task_hpp__
#define __Merge_Task_hpp__

#include <assert.h>
#include "tbb/task_arena.h"
#include "tbb/task.h"

/* ========================================
Assumptions:
	class Work
	{
		void operator() (size_t task_id, Result &res);
	};
  =========================================*/

namespace MergeTaskUtils
{
	template <class Result, size_t id>
	class MergeResult
	{
	public:
		inline static void merge(Result& res0, const Result* res)
		{
			res0 += res[id];
			MergeResult<Result, id - 1>::merge(res0, res);
		}
	};

	template <class Result>
	class MergeResult<Result, 0>
	{
	public:
		inline static void merge(Result& res0, const Result* res)
		{ res0 += res[0]; }
	};

	template <class Result, size_t div_num>
	class MergeResultTask : public tbb::task
	{
	public:
		Result res[div_num - 1];
		Result& parent_res;
		MergeResultTask(Result& pres) : parent_res(pres)
		{ assert(div_num > 2); }
		tbb::task* execute() override
		{
			MergeResult<Result, div_num - 2>::merge(parent_res, res);
			return nullptr;
		}
	};

	template <class Result>
	class MergeResultTask<Result, 2> : public tbb::task
	{
	public:
		Result& parent_res;
		Result res;
		MergeResultTask(Result& pres) : parent_res(pres) {}
		tbb::task* execute() override
		{
			parent_res += res;
			return nullptr;
		}
	};

	template <class Result, size_t div_num>
	class MergeResultTask2 : public tbb::task
	{
	public:
		Result res[div_num - 1];
		Result& parent_res;
		const size_t res_num_min_1;
		MergeResultTask2(Result& pres, size_t re_num) :
			parent_res(pres), res_num_min_1(re_num - 1)
		{
			assert(res_num_min_1);
			assert(res_num_min_1 < div_num);
		}
		tbb::task* execute() override
		{
			for (size_t i = 0; i < res_num_min_1; ++i)
				parent_res += res[i];
			return nullptr;
		}
	};

	template <class Work, class Result>
	class WorkTask : public tbb::task
	{
	protected:
		const size_t tk_id;
		size_t padding;
		Work& work;
		Result& res;
	public:
		WorkTask(size_t id, Work& wk, Result& re) :
			tk_id(id), work(wk), res(re) {}
		tbb::task* execute() override
		{ work(tk_id, res); return nullptr; }
	};

	template <class Work, class Result, size_t div_num, size_t div_id>
	class SpawnChildWorker
	{
	public:
		using ConTask = MergeResultTask<Result, div_num>;
		inline static void spawn(ConTask& c, Work& work, size_t start_id)
		{
			assert(div_id < div_num);
			c.spawn(*new(c.allocate_child())
				WorkTask<Work, Result>(start_id,
					work, c.res[div_num - div_id - 1]));
			SpawnChildWorker<Work, Result, div_num, div_id - 1>::spawn(c, work, start_id + 1);
		}
	};

	template <class Work, class Result, size_t div_num>
	class SpawnChildWorker<Work, Result, div_num, 1>
	{
	public:
		using ConTask = MergeResultTask<Result, div_num>;
		inline static void spawn(ConTask& c, Work& work, size_t start_id)
		{
			assert(1 < div_num);
			c.spawn(*new(c.allocate_child())
				WorkTask<Work, Result>(start_id,
					work, c.res[div_num - 2]));
		}
	};

	template <class Work, class Result, size_t div_num = 2>
	class MergeTask2 : public tbb::task
	{
	protected:
		using ConTask = MergeTaskUtils::MergeResultTask<Result, div_num>;

		const size_t tk_id0;
		size_t padding;
		Work& work;
		Result& res;

	public:
		MergeTask2(size_t start_id,
			Work& wk, Result& re) :
			tk_id0(start_id),
			work(wk), res(re)
		{
			assert(div_num > 1);
		}
		tbb::task* execute() override
		{
			ConTask& c = *new(allocate_continuation()) ConTask(res);
			c.set_ref_count(div_num);
			SpawnChildWorker<Work, Result, div_num, div_num - 1>::spawn(c, work, tk_id0);
			new (this) WorkTask<Work, Result>(tk_id0 + div_num - 1, work, res);
			recycle_as_child_of(c);
			return this;
		}
	};

	template <class Work, class Result, size_t div_num, size_t div_id>
	class SpawnChildren
	{
	public:
		using ConTask = MergeResultTask<Result, div_num>;
		static void spawn(ConTask& c, Work& work,
			size_t blk_start_id, size_t task_id0, size_t task_num);
	};

	template <class Work, class Result, size_t div_num>
	class SpawnChildren<Work, Result, div_num, 1>
	{
	public:
		using ConTask = MergeResultTask<Result, div_num>;
		static void spawn(ConTask& c, Work& work,
			size_t blk_start_id, size_t task_id0, size_t task_num);
	};
}

template <class Work, class Result, size_t div_num>
class MergeTask : public tbb::task
{
protected:
	using ConTask = MergeTaskUtils::MergeResultTask<Result, div_num>;
	using ConTask2 = MergeTaskUtils::MergeResultTask2<Result, div_num>;

	size_t tk_id0, tk_id1;
	Work& work;
	Result& res;

public:
	MergeTask(
		size_t start_id,
		size_t end_id,
		Work& wk, Result& re) :
		tk_id0(start_id),
		tk_id1(end_id),
		work(wk), res(re)
	{ assert(div_num > 2); }
	tbb::task* execute() override
	{
		size_t div_id;
		const size_t tk_num = tk_id1 - tk_id0;
		if (tk_num > (div_num * div_num))
		{
			ConTask& c = *new(allocate_continuation()) ConTask(res);
			c.set_ref_count(div_num);
			MergeTaskUtils::SpawnChildren<Work, Result, div_num, div_num - 1>::spawn(c, work, tk_id0, tk_id0, tk_num);
			recycle_as_child_of(c);
			tk_id0 = tk_id0 + (div_num - 1) * tk_num / div_num;
			return this;
		}
		else if (tk_num > div_num)
		{
			const size_t blk_num = (tk_num + div_num - 1) / div_num;
			ConTask2& c = *new(allocate_continuation()) ConTask2(res, blk_num);
			c.set_ref_count(blk_num);
			for (div_id = 0; div_id < blk_num - 1; ++div_id)
			{
				c.spawn(*new(c.allocate_child())
					MergeTaskUtils::MergeTask2<Work, Result, div_num>(
						tk_id0, work, c.res[div_id]));
				tk_id0 += div_num;
			}
			recycle_as_child_of(c);
			return this;
		}
		else if (tk_num == div_num)
		{
			ConTask& c = *new(allocate_continuation()) ConTask(res);
			c.set_ref_count(div_num);
			MergeTaskUtils::SpawnChildWorker<Work, Result, div_num, div_num - 1>::spawn(c, work, tk_id0);
			new (this) MergeTaskUtils::WorkTask<Work, Result>(tk_id0 + div_num - 1, work, res);
			recycle_as_child_of(c);
			return this;
		}
		else if (tk_num > 1)
		{
			ConTask2& c = *new(allocate_continuation()) ConTask2(res, tk_num);
			c.set_ref_count(tk_num);
			for (div_id = 0; div_id < tk_num - 1; ++div_id)
				c.spawn(*new(c.allocate_child())
					MergeTaskUtils::WorkTask<Work, Result>(
						tk_id0 + div_id, work, c.res[div_id]));
			new (this) MergeTaskUtils::WorkTask<Work, Result>(tk_id1 - 1, work, res);
			recycle_as_child_of(c);
			return this;
		}
		else if (tk_num == 1)
			work(tk_id1 - 1, res);
		return nullptr;
	}
};

namespace MergeTaskUtils
{
	template <class Work, class Result, size_t div_num, size_t div_id>
	inline void SpawnChildren<Work, Result, div_num, div_id>::spawn(
		ConTask& c, Work& work, size_t blk_start_id, size_t task_id0, size_t task_num)
	{
		assert(div_id < div_num);
		const size_t blk_end_id = task_id0 + (div_num - div_id) * task_num / div_num;
		c.spawn(*new(c.allocate_child())
			MergeTask<Work, Result, div_num>(
				blk_start_id, blk_end_id, work, c.res[div_num - div_id - 1]));
		SpawnChildren<Work, Result, div_num, div_id - 1>::spawn(
			c, work, blk_end_id, task_id0, task_num);
	}

	template <class Work, class Result, size_t div_num>
	inline void SpawnChildren<Work, Result, div_num, 1>::spawn(
		ConTask& c, Work& work, size_t blk_start_id, size_t task_id0, size_t task_num)
	{
		assert(1 < div_num);
		const size_t blk_end_id = task_id0 + (div_num - 1) * task_num / div_num;
		c.spawn(*new(c.allocate_child())
			MergeTask<Work, Result, div_num>(
				blk_start_id, blk_end_id, work, c.res[div_num - 2]));
	}
}

template <class Work, class Result>
class MergeTask<Work, Result, 2> : public tbb::task
{
protected:
	using ConTask = MergeTaskUtils::MergeResultTask<Result, 2>;
	size_t tk_id0, tk_id1;
	Work& work;
	Result& res;
public:
	MergeTask(size_t start_id, size_t end_id, Work& wk, Result& re) :
		tk_id0(start_id), tk_id1(end_id), work(wk), res(re) {}
	tbb::task* execute() override
	{
		const size_t tk_num = tk_id1 - tk_id0;
		if (tk_num > 2)
		{
			ConTask& c = *new(allocate_continuation()) ConTask(res);
			c.set_ref_count(2);
			const size_t mid_tk_id = tk_id0 + tk_num / 2;
			c.spawn(*new(c.allocate_child())
				MergeTask<Work, size_t, 2>(
					tk_id0, mid_tk_id, work, c.res));
			recycle_as_child_of(c);
			tk_id0 = mid_tk_id;
			return this;
		}
		else if (tk_num == 2)
		{
			ConTask& c = *new(allocate_continuation()) ConTask(res);
			c.set_ref_count(2);
			c.spawn(*new(c.allocate_child())
				MergeTaskUtils::WorkTask<Work, size_t>(
					tk_id0, work, c.res));
			new (this) MergeTaskUtils::WorkTask<Work, Result>(tk_id0 + 1, work, res);
			recycle_as_child_of(c);
			return this;
		}
		else if (tk_num == 1)
			work(tk_id0, res);
		return nullptr;
	}
};

#endif