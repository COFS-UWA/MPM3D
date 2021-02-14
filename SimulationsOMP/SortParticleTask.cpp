#include "SimulationsOMP_pcp.h"

#include <assert.h>
#include "SortParticleTask.h"

#define Block_Low(block_id, block_num, data_num) ((block_id) * (data_num) / (block_num))

namespace SortUtils
{
	void SortParticleMem::init(
		void* shared_mem,
		size_t thread_num,
		size_t pcl_num,
		size_t ori_pcl_num)
	{
		clear();
		valid_pcl_num = pcl_num;
		pcl_digit_num = max_digit_num(ori_pcl_num);
		max_block_num = thread_num * max_block_num_div_thread_num;
		size_t block_num = (pcl_num + parallel_divide_min_pcl_num_per_block - 1)
						/ parallel_divide_min_pcl_num_per_block;
		if (block_num > max_block_num)
			block_num = max_block_num;
		char* cur_mem = (char*)self_mem.alloc(
			  sizeof(size_t *) * (pcl_digit_num * 4 + 2)
			+ sizeof(SortBin *) * thread_num
			+ sizeof(size_t) * (pcl_num * 2 + 2)
			+ Cache_Alignment * 2);
		//
		pcl_in_elem_arrays = (size_t **)cur_mem; // pcl_digit_num + 1
		cur_mem += sizeof(size_t *) * (pcl_digit_num + 1);
		cur_to_prev_pcl_arrays = (size_t **)cur_mem; // pcl_digit_num + 1
		cur_mem += sizeof(size_t *) * (pcl_digit_num + 1);
		tmp_pcl_in_elem_arrays = (size_t **)cur_mem; // pcl_digit_num
		cur_mem += sizeof(size_t *) * pcl_digit_num;
		tmp_cur_to_prev_pcl_arrays = (size_t **)cur_mem; // pcl_digit_num
		cur_mem += sizeof(size_t *) * pcl_digit_num;
		thread_bins = (SortBin **)cur_mem;
		cur_mem = cache_aligned(cur_mem + sizeof(SortBin *) * thread_num);
		//
		pcl_in_elem = ((size_t*)cur_mem) + 1;
		cur_mem = cache_aligned(cur_mem + sizeof(size_t) * (pcl_num + 2));
		cur_to_prev_pcl = (size_t*)cur_mem;
		pcl_in_elem_arrays[0] = pcl_in_elem;
		cur_to_prev_pcl_arrays[0] = cur_to_prev_pcl;
		pcl_in_elem[-1] = SIZE_MAX;
		pcl_in_elem[pcl_num] = SIZE_MAX;

		cur_mem = (char*)cache_aligned(shared_mem);
		// key and value memory
		size_t i;
		for (i = 0; i < pcl_digit_num; ++i)
		{
			pcl_in_elem_arrays[i + 1] = (size_t*)cur_mem;
			cur_mem = cache_aligned(cur_mem + sizeof(size_t) * pcl_num);
		}
		for (i = 0; i < pcl_digit_num; ++i)
		{
			cur_to_prev_pcl_arrays[i + 1] = (size_t*)cur_mem;
			cur_mem = cache_aligned(cur_mem + sizeof(size_t) * pcl_num);
		}
		for (i = 0; i < pcl_digit_num; ++i)
		{
			tmp_pcl_in_elem_arrays[i] = (size_t*)cur_mem;
			cur_mem = cache_aligned(cur_mem + sizeof(size_t) * pcl_num);
		}
		for (i = 0; i < pcl_digit_num; ++i)
		{
			tmp_cur_to_prev_pcl_arrays[i] = (size_t*)cur_mem;
			cur_mem = cache_aligned(cur_mem + sizeof(size_t) * pcl_num);
		}
		ori_pcl_in_elem = pcl_in_elem_arrays[pcl_digit_num];
		ori_cur_to_prev_pcl = cur_to_prev_pcl_arrays[pcl_digit_num];

		// sort bin memory
		root_bin = (SortBin*)cur_mem;
		const size_t sort_bin_size = sizeof(SortBin) * block_num;
		for (i = 0; i < thread_num; ++i)
		{
			cur_mem += sort_bin_size;
			thread_bins[i] = (SortBin*)cur_mem;
		}
	}
	
	tbb::task* SortParticleTask::execute()
	{
		SortBin *const rt_bins = sort_mem.root_bin;
		const size_t* const in_pcl_in_elem = sort_mem.ori_pcl_in_elem;
		const size_t* const in_cur_to_prev_pcl = sort_mem.ori_cur_to_prev_pcl;
		const size_t pcl_digit_pos = sort_mem.pcl_digit_num - 1;
		const size_t pcl_num = sort_mem.valid_pcl_num;
		if (pcl_num > SortParticleMem::serial_sort_max_pcl_num)
		{
			size_t* const out_pcl_in_elem = sort_mem.pcl_in_elem_arrays[pcl_digit_pos];
			size_t* const out_cur_to_prev_pcl = sort_mem.cur_to_prev_pcl_arrays[pcl_digit_pos];
			size_t bin_id;
			if (pcl_num > SortParticleMem::parallel_divide_min_pcl_num_per_block)
			{
				size_t* const tmp_pcl_in_elem = sort_mem.tmp_pcl_in_elem_arrays[pcl_digit_pos];
				size_t* const tmp_cur_to_prev_pcl = sort_mem.tmp_cur_to_prev_pcl_arrays[pcl_digit_pos];
				size_t blk_num = (pcl_num + SortParticleMem::parallel_divide_min_pcl_num_per_block - 1)
									/ SortParticleMem::parallel_divide_min_pcl_num_per_block;
				if (blk_num > sort_mem.max_block_num)
					blk_num = sort_mem.max_block_num;
				size_t blk_id, end_id, start_id = 0;
				set_ref_count(blk_num + 1);
				for (blk_id = 1; blk_id < blk_num; ++blk_id)
				{
					end_id = Block_Low(blk_id, blk_num, pcl_num);
					spawn(*new(allocate_child())
						Internal::CountSortPclTask(
							tmp_pcl_in_elem + start_id,
							tmp_cur_to_prev_pcl + start_id,
							in_pcl_in_elem + start_id,
							in_cur_to_prev_pcl + start_id,
							end_id - start_id,
							pcl_digit_pos,
							rt_bins[blk_id - 1]));
					start_id = end_id;
				}
				// last child
				spawn_and_wait_for_all(*new(allocate_child())
					Internal::CountSortTask(
						tmp_pcl_in_elem + start_id,
						tmp_cur_to_prev_pcl + start_id,
						in_pcl_in_elem + start_id,
						in_cur_to_prev_pcl + start_id,
						pcl_num - start_id,
						pcl_digit_pos,
						rt_bins[blk_num - 1]));

				if (pcl_digit_pos) // not the last digit
				{
					set_ref_count(Internal::radix_bucket_num + 1);
					start_id = 0;
					for (bin_id = 0; bin_id < Internal::radix_bucket_num; ++bin_id)
					{
						size_t blk_start_id = 0;
						for (blk_id = 0; blk_id < blk_num; ++blk_id)
						{
							const SortBin& blk_bin = rt_bins[blk_id];
							const size_t blk_data_num = blk_bin.count_bin[bin_id];
							const size_t in_offset_id = Block_Low(blk_id, blk_num, pcl_num) + blk_bin.sum_bin[bin_id];
							memcpy(out_pcl_in_elem + start_id + blk_start_id,
								   tmp_pcl_in_elem + in_offset_id,
								   blk_data_num * sizeof(size_t));
							memcpy(out_cur_to_prev_pcl + start_id + blk_start_id,
								   tmp_cur_to_prev_pcl + in_offset_id,
								   blk_data_num * sizeof(size_t));
							blk_start_id += blk_data_num;
						}
						if (blk_start_id)
						{
							spawn(*new(allocate_child())
								SortTask(sort_mem._sort_mem,
									start_id,
									blk_start_id,
									pcl_digit_pos - 1));
							start_id += blk_start_id;
						}
						else
							decrement_ref_count();
					}
					sort_mem.valid_pcl_num = start_id;
					wait_for_all();
				}
				else // already last digit, merge the result and done
				{
					start_id = 0;
					for (bin_id = 0; bin_id < Internal::radix_bucket_num; ++bin_id)
						for (blk_id = 0; blk_id < blk_num; ++blk_id)
						{
							SortBin& blk_bin = rt_bins[blk_id];
							const size_t blk_data_num = blk_bin.count_bin[bin_id];
							const size_t in_offset_id = Block_Low(blk_id, blk_num, pcl_num) + blk_bin.sum_bin[bin_id];
							memcpy(out_pcl_in_elem + start_id,
								   tmp_pcl_in_elem + in_offset_id,
								   blk_data_num * sizeof(size_t));
							memcpy(out_cur_to_prev_pcl + start_id,
								   tmp_cur_to_prev_pcl + in_offset_id,
								   blk_data_num * sizeof(size_t));
							start_id += blk_data_num;
						}
					sort_mem.valid_pcl_num = start_id;
				}
			}
			else // divide serially and sort parallely
			{
				SortBin& rt_bin = rt_bins[0];
				sort_mem.valid_pcl_num
					= Internal::count_sort_pcl(
					out_pcl_in_elem,
					out_cur_to_prev_pcl,
					in_pcl_in_elem,
					in_cur_to_prev_pcl,
					pcl_num,
					pcl_digit_pos,
					rt_bin);
				if (pcl_digit_pos) // not the last digit
				{
					size_t* const c_bin = rt_bin.count_bin;
					size_t* const s_bin = rt_bin.sum_bin;
					set_ref_count(Internal::radix_bucket_num + 1);
					for (bin_id = 0; bin_id < Internal::radix_bucket_num; ++bin_id)
						if (c_bin[bin_id])
						{
							spawn(*new(allocate_child())
								SortTask(
									sort_mem._sort_mem,
									s_bin[bin_id],
									c_bin[bin_id],
									pcl_digit_pos - 1));
						}
						else
							decrement_ref_count();
					wait_for_all();
				}
			}
		}
		else // sort serially
		{
			sort_mem.valid_pcl_num = 0;
			for (size_t p_id = 0; p_id < pcl_num; ++p_id)
			{
				if (in_pcl_in_elem[p_id] != SIZE_MAX)
					++sort_mem.valid_pcl_num;
			}
			Internal::serial_sort(
				sort_mem.pcl_in_elem,
				sort_mem.cur_to_prev_pcl,
				in_pcl_in_elem,
				in_cur_to_prev_pcl,
				pcl_num,
				pcl_digit_pos,
				rt_bins[0],
				sort_mem.tmp_pcl_in_elem_arrays[0],
				sort_mem.tmp_cur_to_prev_pcl_arrays[0]);
		}
		sort_mem.pcl_in_elem[sort_mem.valid_pcl_num] = SIZE_MAX;
		return nullptr;
	}
}
