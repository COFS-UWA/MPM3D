#ifndef __Parallel_Reduce_Task_hpp__
#define __Parallel_Reduce_Task_hpp__

#define TBB_SUPPRESS_DEPRECATED_MESSAGES 1
#include <tbb/tbb.h>

namespace ParaUtil
{
	template <class Result>
	class ParallelReduceContin2Task : public tbb::task
	{
	public:
		ParallelReduceContin2Task(Result& _res) : result0(_res) {}
		tbb::task* execute() override { result0.join(result1); return nullptr; }
		Result& result0, result1;
	};

	template <class Result>
	class ParallelReduceContin3Task : public tbb::task
	{
	public:
		ParallelReduceContin3Task(Result& _res) : result0(_res) {}
		tbb::task* execute() override
		{
			result0.join(result1); result0.join(result2); return nullptr;
		}
		Result& result0, result1, result2;
	};

	template <class Result>
	class ParallelReduceContin4Task : public tbb::task
	{
	public:
		ParallelReduceContin4Task(Result& _res) : result0(_res) {}
		tbb::task* execute() override
		{
			result0.join(result1);
			result0.join(result2);
			result0.join(result3);
			return nullptr;
		}
		Result& result0, result1, result2, result3;
	};

	template <class Operator, class Result>
	class ParallelReduceTask : public tbb::task
	{
	public:
		ParallelReduceTask(const Operator& _work, Result& _res, size_t _head, size_t _end) :
			work(_work), result(_res), head_id(_head), end_id(_end) {}

		inline static size_t blk_low_id(size_t blk_id, size_t blk_num, size_t data_num)
		{
			return blk_id * data_num / blk_num;
		}

		tbb::task* execute() override
		{
			constexpr size_t divide_num = 4;
			const size_t num = end_id - head_id;

			if (num == 1)
				return work(*this, head_id, result);

			if (num < (divide_num + 1))
			{
				end_id = head_id + 1;
				if (num == 2)
				{
					ParallelReduceContin2Task<Result>& continuation
						= *new (tbb::task::allocate_continuation())
						ParallelReduceContin2Task<Result>(result);
					continuation.set_ref_count(2);
					recycle_as_child_of(continuation);
					tbb::task::spawn(*new (continuation.allocate_child())
						ParallelReduceTask<Operator, Result>(work, continuation.result1, head_id + 1, head_id + 2));
					return this;
				}
				else if (num == 3)
				{
					ParallelReduceContin3Task<Result>& continuation
						= *new (tbb::task::allocate_continuation())
						ParallelReduceContin3Task<Result>(result);
					continuation.set_ref_count(3);
					recycle_as_child_of(continuation);
					tbb::task::spawn(*new (continuation.allocate_child())
						ParallelReduceTask<Operator, Result>(work, continuation.result2, head_id + 2, head_id + 3));
					tbb::task::spawn(*new (continuation.allocate_child())
						ParallelReduceTask<Operator, Result>(work, continuation.result1, head_id + 1, head_id + 2));
					return this;
				}
				ParallelReduceContin4Task<Result>& continuation
					= *new (tbb::task::allocate_continuation())
					ParallelReduceContin4Task<Result>(result);
				continuation.set_ref_count(4);
				recycle_as_child_of(continuation);
				tbb::task::spawn(*new (continuation.allocate_child())
					ParallelReduceTask<Operator, Result>(work, continuation.result3, head_id + 3, head_id + 4));
				tbb::task::spawn(*new (continuation.allocate_child())
					ParallelReduceTask<Operator, Result>(work, continuation.result2, head_id + 2, head_id + 3));
				tbb::task::spawn(*new (continuation.allocate_child())
					ParallelReduceTask<Operator, Result>(work, continuation.result1, head_id + 1, head_id + 2));
				return this;
			}

			ParallelReduceContin4Task<Result>& continuation
				= *new (tbb::task::allocate_continuation())
				ParallelReduceContin4Task<Result>(result);
			continuation.set_ref_count(4);
			recycle_as_child_of(continuation);
			const size_t head_id1 = head_id + blk_low_id(1, divide_num, num);
			const size_t head_id2 = head_id + blk_low_id(2, divide_num, num);
			const size_t head_id3 = head_id + blk_low_id(3, divide_num, num);
			tbb::task::spawn(*new (continuation.allocate_child()) ParallelReduceTask<Operator, Result>(work, continuation.result3, head_id3, end_id));
			tbb::task::spawn(*new (continuation.allocate_child()) ParallelReduceTask<Operator, Result>(work, continuation.result2, head_id2, head_id3));
			tbb::task::spawn(*new (continuation.allocate_child()) ParallelReduceTask<Operator, Result>(work, continuation.result1, head_id1, head_id2));
			end_id = head_id1;
			return this;
		}

		const Operator& work;
		Result& result;
		size_t head_id, end_id;
	};

	template <class Operator, class Result>
	void parallel_reduce(
		Operator& _work,
		Result& _result,
		size_t task_num)
	{
		if (task_num)
			tbb::task::spawn_root_and_wait(
				*new (tbb::task::allocate_root())
				ParallelReduceTask<Operator, Result>(_work, _result, 0, task_num));
	}

	template <class Operator, class Result>
	void parallel_for_within_task(
		tbb::task& parent,
		Operator& _work,
		Result& _result,
		size_t task_num)
	{
		if (task_num)
			parent.spawn(*new (parent.allocate_child())
				ParallelReduceTask<Operator, Result>(_work, _result, 0, task_num));
	}
}

#endif