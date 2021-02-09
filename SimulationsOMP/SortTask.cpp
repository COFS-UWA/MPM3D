#include "SimulationsOMP_pcp.h"

#include "tbb/task_arena.h"
#include "SortTask.h"

#define Block_Low(block_id, block_num, data_num) ((block_id) * (data_num) / (block_num))

namespace SortUtils
{
	tbb::task* SortTask::execute()
	{
		const size_t my_th_id = tbb::task_arena::current_thread_index();
		SortBin* const th_bins = sort_mem.thread_bins[my_th_id];
		const size_t* const in_key = sort_mem.key_arrays[digit_pos + 1] + start_id;
		const size_t* const in_val = sort_mem.key_arrays[digit_pos + 1] + start_id;
		if (data_num > SortMem::serial_sort_max_data_num)
		{
			size_t*const out_key = sort_mem.key_arrays[digit_pos] + start_id;
			size_t*const out_val = sort_mem.val_arrays[digit_pos] + start_id;
			size_t bin_id;
			if (data_num > SortMem::parallel_divide_min_data_num_per_block)
			{
				size_t* const tmp_key = sort_mem.tmp_key_arrays[digit_pos] + start_id;
				size_t* const tmp_val = sort_mem.tmp_val_arrays[digit_pos] + start_id;
				// block number
				size_t block_num = (data_num + SortMem::parallel_divide_min_data_num_per_block - 1)
									/ SortMem::parallel_divide_min_data_num_per_block;
				if (block_num > sort_mem.max_block_num)
					block_num = sort_mem.max_block_num;
				set_ref_count(block_num + 1);
				size_t blk_id, child_end_id, child_start_id = 0;
				for (blk_id = 1; blk_id < block_num; ++blk_id)
				{
					child_end_id = Block_Low(blk_id, block_num, data_num);
					spawn(*new(allocate_child())
						Internal::CountSortTask(
							tmp_key + child_start_id,
							tmp_val + child_start_id,
							in_key + child_start_id,
							in_val + child_start_id,
							child_end_id - child_start_id,
							digit_pos,
							th_bins[blk_id - 1]));
					child_start_id = child_end_id;
				}
				// last child
				spawn_and_wait_for_all(*new(allocate_child())
					Internal::CountSortTask(
						tmp_key + child_start_id,
						tmp_val + child_start_id,
						in_key + child_start_id,
						in_val + child_start_id,
						data_num - child_start_id,
						digit_pos,
						th_bins[block_num - 1]));

				if (digit_pos) // not the last digit
				{
					child_start_id = 0;
					set_ref_count(Internal::radix_bucket_num + 1);
					for (bin_id = 0; bin_id < Internal::radix_bucket_num; ++bin_id)
					{
						size_t blk_start_id = 0;
						for (blk_id = 0; blk_id < block_num; ++blk_id)
						{
							const SortBin& blk_bin = th_bins[blk_id];
							const size_t blk_data_num = blk_bin.count_bin[bin_id];
							const size_t in_offset_id = Block_Low(blk_id, block_num, data_num) + blk_bin.sum_bin[bin_id];
							memcpy(out_key + child_start_id + blk_start_id,
								   tmp_key + in_offset_id,
								   blk_data_num * sizeof(size_t));
							memcpy(out_val + child_start_id + blk_start_id,
								   tmp_val + in_offset_id,
								   blk_data_num * sizeof(size_t));
							blk_start_id += blk_data_num;
						}
						if (blk_start_id)
						{
							spawn(*new(allocate_child())
								SortTask(sort_mem,
									start_id + child_start_id,
									blk_start_id,
									digit_pos - 1));
							child_start_id += blk_start_id;
						}
						else
							decrement_ref_count();
					}
					wait_for_all();
				}
				else // already last digit, merge the result and done
				{
					child_start_id = 0;
					for (bin_id = 0; bin_id < Internal::radix_bucket_num; ++bin_id)
					{
						for (blk_id = 0; blk_id < block_num; ++blk_id)
						{
							SortBin& blk_bin = th_bins[blk_id];
							const size_t blk_data_num = blk_bin.count_bin[bin_id];
							const size_t in_offset_id = Block_Low(blk_id, block_num, data_num) + blk_bin.sum_bin[bin_id];
							memcpy(out_key + child_start_id,
								   tmp_key + in_offset_id,
								   blk_data_num * sizeof(size_t));
							memcpy(out_val + child_start_id,
								   tmp_val + in_offset_id,
								   blk_data_num * sizeof(size_t));
							child_start_id += blk_data_num;
						}
					}
				}
			}
			else // divide serially and sort parallely
			{
				SortBin& th_bin = th_bins[0];
				Internal::count_sort(
					out_key, out_val,
					in_key, in_key,
					data_num,
					digit_pos,
					th_bin);
				if (digit_pos) // not the last digit
				{
					size_t* const c_bin = th_bin.count_bin;
					size_t* const s_bin = th_bin.sum_bin;
					set_ref_count(Internal::radix_bucket_num + 1);
					for (bin_id = 0; bin_id < Internal::radix_bucket_num; ++bin_id)
					{
						if (c_bin[bin_id])
						{
							spawn(*new(allocate_child())
								SortTask(sort_mem,
									start_id + s_bin[bin_id],
									c_bin[bin_id],
									digit_pos - 1));
						}
						else
							decrement_ref_count();
					}
					wait_for_all();
				}
			}
		}
		else // sort serially
		{
			Internal::serial_sort(
				sort_mem.out_keys + start_id,
				sort_mem.out_vals + start_id,
				in_key, in_val,
				data_num,
				digit_pos,
				th_bins[0],
				sort_mem.tmp_key_arrays[0],
				sort_mem.tmp_val_arrays[0]);
		}
		return nullptr;
	}
}
