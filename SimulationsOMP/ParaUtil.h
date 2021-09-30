#ifndef __Para_Util__
#define __Para_Util__

namespace ParaUtil
{
    template <size_t min_data_num_per_task, size_t suggested_task_num_per_thread>
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

    template <size_t radix_bit_num>
    inline size_t cal_digit_num(size_t data)
    {
        size_t digit_num = 0;
        while (data)
        {
            digit_num += radix_bit_num;
            data >>= radix_bit_num;
        }
        return digit_num;
    }
}

#endif