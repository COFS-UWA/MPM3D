#include "SimulationsOMP_pcp.h"

#include <cstring>

#include "ParaUtilSerialSort.h"

#define Data_Digit(num, disp, mask) (((num) >> (disp)) & (mask))
#define Data_Digit_256(num, disp) Data_Digit(num, disp, 0xFF)

#define Swap_Address(a, b) \
		(a) = (size_t *)((size_t(a)) ^ (size_t(b))); \
		(b) = (size_t *)((size_t(a)) ^ (size_t(b))); \
		(a) = (size_t *)((size_t(a)) ^ (size_t(b)))

namespace ParaUtil
{
    namespace Internal
    {
        void insertion_sort(
            size_t keys[],
            size_t values[],
            size_t data_num)
        {
            for (size_t i = 0; i < data_num - 1; ++i)
            {
                size_t min_key = keys[i];
                size_t min_key_id = i;
                size_t j;
                for (j = i + 1; j < data_num; ++j)
                {
                    if (min_key > keys[j])
                    {
                        min_key = keys[j];
                        min_key_id = j;
                    }
                }
                keys[min_key_id] = keys[i];
                keys[i] = min_key;
                size_t tmp_value = values[min_key_id];
                values[min_key_id] = values[i];
                values[i] = tmp_value;
            }
        }

        void in_place_lsd_sort(
            size_t keys[],
            size_t values[],
            size_t data_num,
            size_t key_bit_offset,
            size_t radix_bit_num,
            size_t sum_bin[],
            size_t end_bin[])
        {
            const size_t radix_num = size_t(1) << radix_bit_num;
            const size_t radix_mask = radix_num - 1;
            size_t i;

            // prefix scan to form sum bin and end bin
            memset(sum_bin, 0, radix_num * sizeof(size_t));

            for (i = 0; i < data_num; ++i)
                ++sum_bin[Data_Digit(keys[i], key_bit_offset, radix_mask)];

            end_bin[0] = sum_bin[0];
            sum_bin[0] = 0;
            for (i = 1; i < radix_mask; ++i)
            {
                end_bin[i] = end_bin[i - 1] + sum_bin[i];
                sum_bin[i] = end_bin[i - 1];
            }
            end_bin[radix_mask] = data_num;
            sum_bin[radix_mask] = end_bin[radix_mask - 1];

            // swap reorder
            for (size_t radix = 0; radix < radix_num; ++radix)
            {
                size_t& data_pos = sum_bin[radix];
                const size_t end_pos = end_bin[radix];

                while (data_pos < end_pos)
                {
                    size_t tmp_radix = Data_Digit(keys[data_pos], key_bit_offset, radix_mask);
                    if (tmp_radix != radix) // skip sorted part
                    {
                        // fill misplaced item by swap reordering
                        do
                        {
                            size_t& new_data_pos = sum_bin[tmp_radix];
                            size_t new_radix;
                            while ((new_radix = Data_Digit(keys[new_data_pos], key_bit_offset, radix_mask)) == tmp_radix)
                                ++new_data_pos;
                            // swap
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
        }

        void serial_sort(
            size_t keys[],
            size_t values[],
            size_t data_num,
            size_t key_bit_num,
            size_t tmp_keys[],
            size_t tmp_values[])
        {
            if (data_num < max_data_num_for_inseration_sort)
            {
                insertion_sort(keys, values, data_num);
                return;
            }

            constexpr size_t radix_bit_num = 8;
            constexpr size_t radix_num = size_t(1) << radix_bit_num;
            constexpr size_t radix_mask = radix_num - 1;
            size_t sum_bin[radix_num], end_bin[radix_num];
            size_t sort_num = (key_bit_num + radix_bit_num - 1) / radix_bit_num;
            if (sort_num & 1)
            {
                --sort_num;
                key_bit_num -= sort_num * radix_bit_num;
                in_place_lsd_sort(keys, values, data_num, 0, key_bit_num, sum_bin, end_bin);
            }
            else
            {
                key_bit_num = 0;
            }

            size_t* key0 = keys, * key1 = tmp_keys;
            size_t* value0 = values, * value1 = tmp_values;
            size_t i;
            for (; sort_num-- > 0;)
            {
                memset(sum_bin, 0, radix_num * sizeof(size_t));

                for (i = 0; i < data_num; ++i)
                    ++sum_bin[Data_Digit_256(key0[i], key_bit_num)];

                for (i = 1; i < radix_num; ++i)
                    sum_bin[i] += sum_bin[i - 1];

                for (i = data_num; i-- > 0;)
                {
                    const size_t pos_id = --sum_bin[Data_Digit_256(key0[i], key_bit_num)];
                    key1[pos_id] = key0[i];
                    value1[pos_id] = value0[i];
                }

                key_bit_num += radix_bit_num;
                Swap_Address(key0, key1);
                Swap_Address(value0, value1);
            }
        }
    }
}