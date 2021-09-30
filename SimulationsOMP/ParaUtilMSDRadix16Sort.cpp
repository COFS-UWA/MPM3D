#include "SimulationsOMP_pcp.h"

#include "ParaUtil.h"
#include "ParaUtilSerialSort.h"
#include "ParaUtilMSDRadix16Sort.h"

#define Data_Digit(num, disp, mask) (((num) >> (disp)) & (mask))
#define Data_Digit_16(num, disp) Data_Digit(num, disp, 0xF)

namespace ParaUtil
{
namespace Internal
{
	static constexpr uint64_t radix_num = 1 << msd_radix_16_bit_num;
	static constexpr uint64_t radix_mask = radix_num - 1;
	static constexpr uint64_t max_radix = radix_mask;

    inline void prefix_scan_16_serial(
        size_t keys[],
        size_t data_num,
        size_t bit_offset,
        size_t sum_bin[radix_num])
    {
        size_t i;

        memset(sum_bin, 0, radix_num * sizeof(size_t));

        for (i = 0; i < data_num; ++i)
            ++sum_bin[Data_Digit_16(keys[i], bit_offset)];

        for (i = 1; i < (radix_num - 1); ++i)
            sum_bin[i] += sum_bin[i - 1];

        for (i = radix_num - 1; i != 0; --i)
            sum_bin[i] = sum_bin[i - 1];
        sum_bin[0] = 0;
    };
    
    inline void swap_reorder_16_serial(
        size_t keys[],
        size_t values[],
        size_t data_num,
        size_t bit_offset,
        size_t sum_bin[radix_num])
    {
        const size_t end_bin[radix_num] = {
            sum_bin[1], sum_bin[2], sum_bin[3], sum_bin[4],
            sum_bin[5], sum_bin[6], sum_bin[7], sum_bin[8],
            sum_bin[9], sum_bin[10], sum_bin[11], sum_bin[12],
            sum_bin[13], sum_bin[14], sum_bin[15], data_num
        };

        for (uint64_t radix = 0; radix < radix_num; ++radix)
        {
            size_t& data_pos = sum_bin[radix];
            const size_t end_pos = end_bin[radix];

            while (data_pos < end_pos)
            {
                uint64_t tmp_radix = Data_Digit_16(keys[data_pos], bit_offset);
                if (tmp_radix != radix) // skip sorted part
                {
                    // fill misplaced item by swap reordering
                    do
                    {
                        size_t& new_data_pos = sum_bin[tmp_radix];
                        size_t new_radix;
                        while ((new_radix = Data_Digit_16(keys[new_data_pos], bit_offset)) == tmp_radix)
                            ++new_data_pos;
                        // swap data
                        size_t tmp = keys[data_pos];
                        keys[data_pos] = keys[new_data_pos];
                        keys[new_data_pos] = tmp;
                        tmp = values[data_pos];
                        values[data_pos] = values[new_data_pos];
                        values[new_data_pos] = tmp;
                        //
                        ++new_data_pos;
                        tmp_radix = new_radix;
                    } while (tmp_radix != radix);
                }
                ++data_pos;
            }
        }
    };

    tbb::task* MSDRadixSort16Task::execute()
    {
        if (!is_continuation)
        {
            if (data_num < min_data_num_for_msd_radix_sort_16 ||
                key_bit_num < (msd_radix_16_bit_num + 1))
            {
                serial_sort(keys, values, data_num,
                    key_bit_num, tmp_keys, tmp_values);
                return nullptr;
            }

            key_bit_num -= msd_radix_16_bit_num;
            size_t sum_bin[radix_num];
            prefix_scan_16_serial(keys, data_num, key_bit_num, sum_bin);
            swap_reorder_16_serial(keys, values, data_num, key_bit_num, sum_bin);

            is_continuation = true;
            recycle_as_continuation();

            size_t num = sum_bin[0] > 1 ? 1 : 0;
            for (size_t radix = 1; radix < radix_num; ++radix)
                if ((sum_bin[radix] - sum_bin[radix - 1]) > 1)
                    ++num;
            set_ref_count(num);

            tbb::task* child_task1 = nullptr;
            if (sum_bin[0] > 1)
                child_task1 = new (allocate_child())
                MSDRadixSort16Task(keys, values, sum_bin[0],
                    key_bit_num, tmp_keys, tmp_values);

            for (size_t radix = 1; radix < radix_num; ++radix)
                if ((num = sum_bin[radix] - sum_bin[radix - 1]) > 1)
                {
                    const size_t offset = sum_bin[radix - 1];
                    spawn(*new (allocate_child())
                        MSDRadixSort16Task(
                            keys + offset,
                            values + offset,
                            num, key_bit_num,
                            tmp_keys + offset,
                            tmp_values + offset));
                }
            return child_task1;
        }
        return nullptr;
    }

    void msd_radix_sort_16(
        size_t keys[],
        size_t values[],
        size_t data_num,
        size_t max_key,
        size_t tmp_keys[],
        size_t tmp_values[],
        int thread_num)
    {
        if (data_num < 2)
            return;

        tbb::task_scheduler_init init(thread_num);

        size_t key_bit_num = cal_digit_num<4>(max_key);
        if (key_bit_num)
        {
            MSDRadixSort16Task& task = *new(tbb::task::allocate_root())
                MSDRadixSort16Task(keys, values, data_num, key_bit_num, tmp_keys, tmp_values);
            tbb::task::spawn_root_and_wait(task);
        }
    }
}
}