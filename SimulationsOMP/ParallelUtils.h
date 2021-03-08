#ifndef __Parallel_Utils_h__
#define __Parallel_Utils_h__

namespace ParallelUtils
{
	inline size_t block_low(size_t block_id, size_t block_num, size_t data_num) noexcept
	{ return block_id * data_num / block_num; }

	template <size_t min_data_num_per_task, size_t max_task_num_per_thread>
	inline size_t cal_task_num(size_t data_num, size_t thread_num) noexcept
	{
		const size_t task_num0 = (data_num + min_data_num_per_task - 1) / min_data_num_per_task;
		const size_t task_num1 = max_task_num_per_thread * thread_num;
		return task_num0 < task_num1 ? task_num0 : task_num1;
	}
}

#endif