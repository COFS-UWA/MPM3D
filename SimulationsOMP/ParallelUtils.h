#ifndef __Parallel_Utils_h__
#define __Parallel_Utils_h__

namespace ParallelUtils
{
	inline size_t block_low(size_t block_id, size_t block_num, size_t data_num) noexcept
	{ return block_id * data_num / block_num; }

    template <size_t min_data_num_per_task,
        size_t suggested_task_num_per_thread = 4>
        inline size_t cal_task_num(size_t thread_num, size_t data_num)
    {
        size_t task_num = thread_num * suggested_task_num_per_thread;
        if (data_num < task_num * min_data_num_per_task)
        {
            task_num = data_num / min_data_num_per_task;
            if (task_num == 0)
                task_num = 1;
        }
        return task_num;
    }
}

#endif