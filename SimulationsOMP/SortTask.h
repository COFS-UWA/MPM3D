#ifndef __Radix_Sort_Task_h__
#define __Radix_Sort_Task_h__

#include "CacheAlignedMem.h"
#include "SortUtils.h"

namespace SortUtils
{
	inline size_t max_digit_num(size_t data_num) noexcept
	{
		size_t max_num = 1 << 8;
		size_t digit_num = 1;
		while (data_num > max_num)
		{
			max_num <<= 8;
			++digit_num;
		}
		return digit_num;
	}

	struct SortMem
	{
	public:
		size_t** key_arrays; // digit_num + 1
		size_t** val_arrays; // digit_num + 1
		size_t** tmp_key_arrays; // digit_num
		size_t** tmp_val_arrays; // digit_num
		size_t *in_keys; // key_arrays[0]
		size_t *in_vals; // val_arrays[0]
		size_t *out_keys; // key_arrays[digit_num]
		size_t *out_vals; // val_arrays[digit_num]
		SortBin *root_bin;
		SortBin **thread_bins; // thread_num * block_num
	protected:
		CacheAlignedMem mem;
		inline void clear() noexcept { mem.free(); }
	public:
		SortMem() {}
		~SortMem() { clear(); }
		void init(
			const size_t thread_num,
			const size_t data_num,
			const size_t digit_num,
			const size_t min_data_num_per_block)
		{
			clear();
			const size_t init_block_num
				= (data_num + min_data_num_per_block - 1)
				  / min_data_num_per_block;
			char *cur_mem = (char *)mem.alloc(
				+ sizeof(size_t *) * (digit_num * 2 + 1) * 2
				+ sizeof(SortBin*) * thread_num
				+ sizeof(size_t) * (digit_num * 2 + 1) * 2 * data_num
				+ sizeof(SortBin) * (1 + thread_num) * init_block_num
				+ Cache_Alignment * (1 + (digit_num * 2 + 1) * 2) + 1);
			key_arrays = (size_t **)cur_mem; // digit_num + 1
			cur_mem += sizeof(size_t *) * (digit_num + 1);
			val_arrays = (size_t **)cur_mem; // digit_num + 1
			cur_mem += sizeof(size_t *) * (digit_num + 1);
			tmp_key_arrays = (size_t **)cur_mem; // digit_num
			cur_mem += sizeof(size_t*) * digit_num;
			tmp_val_arrays = (size_t **)cur_mem; // digit_num
			cur_mem += sizeof(size_t*) * digit_num;
			thread_bins = (SortBin **)cur_mem;
			cur_mem = cache_aligned(cur_mem + sizeof(SortBin *) * thread_num);
			
			// key and value memory
			size_t i;
			for (i = 0; i < digit_num + 1; ++i)
			{
				key_arrays[i] = (size_t *)cur_mem;
				cur_mem = cache_aligned(cur_mem + sizeof(size_t) * data_num);
			}
			for (i = 0; i < digit_num + 1; ++i)
			{
				val_arrays[i] = (size_t*)cur_mem;
				cur_mem = cache_aligned(cur_mem + sizeof(size_t) * data_num);
			}
			for (i = 0; i < digit_num; ++i)
			{
				tmp_key_arrays[i] = (size_t*)cur_mem;
				cur_mem = cache_aligned(cur_mem + sizeof(size_t) * data_num);
			}
			for (i = 0; i < digit_num; ++i)
			{
				tmp_val_arrays[i] = (size_t*)cur_mem;
				cur_mem = cache_aligned(cur_mem + sizeof(size_t));
			}
			in_keys = key_arrays[0];
			in_vals = val_arrays[0];
			out_keys = key_arrays[digit_num];
			out_vals = val_arrays[digit_num];

			// sort bin memory
			root_bin = (SortBin *)cur_mem;
			cur_mem += sizeof(SortBin) * init_block_num;
			for (i = 0; i < thread_num; ++i)
			{
				thread_bins[i] = (SortBin *)cur_mem;
				cur_mem += sizeof(SortBin) * init_block_num;
			}
		}
	};
	
	namespace Internal
	{
		constexpr size_t parallel_divide_min_pcl_num_per_block = 8;
		constexpr size_t serial_sort_max_data_num = 0;
	}
	
	class SortTask : public tbb::task
	{
	protected:
		const SortMem &sort_mem;
		const size_t start_id, data_num;
		const size_t digit_pos;

	public:
		SortTask(
			const SortMem& sm,
			const size_t st_id,
			const size_t d_num,
			const size_t dgt_pos);
		~SortTask() {}
		tbb::task* execute() override;
	};
}

#endif