#pragma once

#define TBB_SUPPRESS_DEPRECATED_MESSAGES 1
#include <tbb/tbb.h>

namespace ParallelUtils
{
	template <class Operator>
	class ParallelForTask : public tbb::task
	{
	public:
		ParallelForTask(const Operator& _work, size_t _head, size_t _end) :
			work(_work), head_id(_head), end_id(_end) {}

		tbb::task* execute() override
		{
			constexpr size_t divide_num = 4;
			const size_t num = end_id - head_id;

			if (num == 1)
				return work(*this, head_id);

			tbb::empty_task& continuation = *new (tbb::task::allocate_continuation()) tbb::empty_task;

			if (num < (divide_num + 1))
			{
				recycle_as_child_of(continuation);
				continuation.set_ref_count(num);
				end_id = head_id + 1;
				for (size_t i = 1; i < num; ++i)
					tbb::task::spawn(*new (continuation.allocate_child()) ParallelForTask<Operator>(work, head_id + i, head_id + i + 1));
				return this;
			}

			continuation.set_ref_count(4);
			recycle_as_child_of(continuation);
			const size_t head_id1 = head_id + blk_low_id(1, divide_num, num);
			const size_t head_id2 = head_id + blk_low_id(2, divide_num, num);
			const size_t head_id3 = head_id + blk_low_id(3, divide_num, num);
			tbb::task::spawn(*new (continuation.allocate_child()) ParallelForTask<Operator>(work, head_id1, head_id2));
			tbb::task::spawn(*new (continuation.allocate_child()) ParallelForTask<Operator>(work, head_id2, head_id3));
			tbb::task::spawn(*new (continuation.allocate_child()) ParallelForTask<Operator>(work, head_id3, end_id));
			end_id = head_id1;
			return this;
		}

	protected:
		const Operator& work;
		size_t head_id, end_id;

		inline static size_t blk_low_id(size_t blk_id, size_t blk_num, size_t data_num)
		{
			return blk_id * data_num / blk_num;
		}
	};

	template <class Operator>
	void parallel_for(Operator& _work, size_t task_num)
	{
		if (task_num)
			tbb::task::spawn_root_and_wait(
				*new (tbb::task::allocate_root())
				ParallelForTask<Operator>(_work, 0, task_num));
	}

	template <class Operator>
	void parallel_for_within_task(tbb::task& parent, Operator& _work, size_t task_num)
	{
		if (task_num)
			parent.spawn(*new (parent.allocate_child())
				ParallelForTask<Operator>(_work, 0, task_num));
	}
}