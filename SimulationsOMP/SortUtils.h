#ifndef __Sort_Utils_h__
#define __Sort_Utils_h__

#include "tbb/task.h"

namespace SortUtils
{
	namespace Internal
	{
		constexpr size_t radix_bucket_num = 0x100; // 256
	}

	struct SortBin
	{
		size_t count_bin[Internal::radix_bucket_num];
		size_t sum_bin[Internal::radix_bucket_num];
	};

	namespace Internal
	{
		inline void count_sort(
			size_t* const out_key,
			size_t* const out_val,
			const size_t* const in_key,
			const size_t* const in_val,
			const size_t data_num,
			SortBin& bin,
			const unsigned char digit_pos)
		{
#define __Cal_Key_Digit__(data, digit_pos) (((data) >> ((digit_pos) * 8)) & (0xFF))
			size_t* const c_bin = bin.count_bin;
			size_t* const s_bin = bin.sum_bin;
			memset(c_bin, 0, sizeof(size_t) * radix_bucket_num);
			size_t i;
			for (i = 0; i < data_num; ++i)
				++c_bin[__Cal_Key_Digit__(in_key[i], digit_pos)];
			s_bin[0] = c_bin[0];
			for (i = 1; i < radix_bucket_num; ++i)
				s_bin[i] = s_bin[i - 1] + c_bin[i];
			for (i = data_num; i--;)
			{
				const size_t pos = --s_bin[__Cal_Key_Digit__(in_key[i], digit_pos)];
				out_key[pos] = in_key[i];
				out_val[pos] = in_val[i];
			}
#undef __Cal_Key_Digit__
		}

		class CountSortTask : public tbb::task
		{
		protected:
			size_t*const out_key;
			size_t*const out_val;
			const size_t*const in_key;
			const size_t*const in_val;
			const size_t data_num;
			SortBin& bin;
			const unsigned char digit_pos;

		public:
			CountSortTask(
				size_t *const _out_key,
				size_t *const _out_val,
				const size_t *_in_key,
				const size_t *_in_val,
				const size_t _data_num,
				SortBin& _bin,
				const unsigned char _digit_pos);
			~CountSortTask();
			tbb::task* execute() override;
		};

		inline void serial_sort(
			size_t *out_key,
			size_t *out_val,
			const size_t* in_key,
			const size_t *in_val,
			const size_t data_num,
			const unsigned char digit_pos,
			SortBin& bin,
			size_t *const mid_key,
			size_t *const mid_val)
		{
			constexpr size_t insertion_sort_max_data_num = 64;
			if (data_num > insertion_sort_max_data_num) // radix sort
			{
				size_t *key0 = const_cast<size_t *>(in_key);
				size_t *val0 = const_cast<size_t *>(in_val);
				size_t *key1 = mid_key;
				size_t *val1 = mid_val;
				for (size_t d_id = 0; d_id < digit_pos; ++d_id)
				{
#define __Swap_Pointer__(a, b) \
	(*(size_t*)&(a)) = size_t(a) ^ size_t(b); \
	(*(size_t*)&(b)) = size_t(a) ^ size_t(b); \
	(*(size_t*)&(a)) = size_t(a) ^ size_t(b)
					count_sort(
						key1, val1,
						key0, val0,
						data_num,
						bin, d_id);
					__Swap_Pointer__(key0, key1);
					__Swap_Pointer__(val0, val1);
#undef __Swap_Pointer__
				}
				count_sort(
					out_key,
					out_val,
					key0, val0,
					data_num,
					bin,
					digit_pos);
			}
			else if (data_num > 2) // insertion sort
			{
				const size_t data_num_min_1 = data_num - 1;
				for (size_t i = 0; i < data_num_min_1; ++i)
				{
					size_t min_key = in_key[i];
					size_t min_key_id = i;
					for (size_t j = i + 1; j < data_num; ++j)
					{
						if (min_key > in_key[j])
						{
							min_key = in_key[j];
							min_key_id = j;
						}
					}
					out_key[i] = min_key;
					out_val[i] = in_val[min_key_id];
					const_cast<size_t *>(in_key)[min_key_id] = in_key[i];
					const_cast<size_t *>(in_val)[min_key_id] = in_val[i];
				}
				out_key[data_num_min_1] = in_key[data_num_min_1];
				out_val[data_num_min_1] = in_val[data_num_min_1];
			}
			else if (data_num == 2)
			{
				if (in_key[0] < in_key[1])
				{
					out_key[0] = in_key[0];
					out_key[1] = in_key[1];
					out_val[0] = in_val[0];
					out_val[1] = in_val[1];
				}
				else
				{
					out_key[1] = in_key[0];
					out_key[0] = in_key[1];
					out_val[1] = in_val[0];
					out_val[0] = in_val[1];
				}
			}
			else if (data_num == 1)
			{
				out_key[0] = in_key[0];
				out_val[0] = in_val[0];
			}
		}
	}
}

#endif