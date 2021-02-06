#ifndef __Radix_Sort_Task_h__
#define __Radix_Sort_Task_h__

#include "CacheAlignedMem.h"
#include "SortUtils.h"

namespace SortUtils
{
	struct SortMem
	{
		static constexpr size_t parallel_divide_min_data_num_per_block = 8;
		static constexpr size_t serial_sort_max_data_num = 0;
		
		size_t max_block_num;
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
	};
	
	class SortTask : public tbb::task
	{
	protected:
		const SortMem &sort_mem;
		const size_t start_id, data_num;
		const size_t digit_pos;
	public:
		SortTask(
			const SortMem& sm,
			size_t st_id,
			size_t d_num,
			size_t dgt_pos) :
			sort_mem(sm),
			start_id(st_id),
			data_num(d_num),
			digit_pos(dgt_pos) {}
		~SortTask() {}
		tbb::task* execute() override;
	};
}

#endif