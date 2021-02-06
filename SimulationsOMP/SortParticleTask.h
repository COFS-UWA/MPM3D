#ifndef __Sort_Particle_Task_h__
#define __Sort_Particle_Task_h__

#include <assert.h>
#include "SortTask.h"

namespace SortUtils
{
	struct SortParticleMem
	{
		static constexpr size_t max_block_num_div_thread_num = 2;
		static constexpr size_t parallel_divide_min_pcl_num_per_block = 2;
		static constexpr size_t serial_sort_max_pcl_num = 1;
		
		union
		{
			struct
			{
				size_t max_block_num;
				size_t **pcl_in_elem_arrays; // pcl_digit_num + 1
				size_t **cur_to_prev_pcl_arrays; // pcl_digit_num + 1
				size_t **tmp_pcl_in_elem_arrays; // pcl_digit_num
				size_t **tmp_cur_to_prev_pcl_arrays; // pcl_digit_num
				size_t *ori_pcl_in_elem; // key_arrays[pcl_digit_num]
				size_t *ori_cur_to_prev_pcl; // val_arrays[pcl_digit_num]
				size_t *pcl_in_elem; // key_arrays[0]
				size_t *cur_to_prev_pcl; // val_arrays[0]
				SortBin *root_bin;
				SortBin **thread_bins; // thread_num * block_num
			};
			SortMem _sort_mem;
		};

		SortParticleMem() {}
		~SortParticleMem() { clear(); }
		// shared_mem used for :
		//  key_arrays
		//  val_arrays
		//  tmp_key_arrays
		//  tmp_val_arrays
		//  ori_pcl_in_elem
		//  ori_cur_to_prev_pcl
		//  root_bin
		//  thread_bins
		// self_mem used for:
		//  key_arrays *[]
		//  val_arrays *[]
		//  tmp_key_arrays *[]
		//  tmp_val_arrays *[]
		//  pcl_in_elem
		//  cur_to_prev_pcl
		inline size_t get_shared_memory_size(
			size_t thread_num, size_t pcl_num,
			size_t pcl_digit_num) noexcept
		{
			max_block_num = thread_num * SortParticleMem::max_block_num_div_thread_num;
			size_t block_num = (pcl_num + SortParticleMem::parallel_divide_min_pcl_num_per_block - 1)
				/ SortParticleMem::parallel_divide_min_pcl_num_per_block;
			if (block_num > max_block_num)
				block_num = max_block_num;
			return Cache_Alignment * pcl_digit_num * 2
				+ sizeof(size_t) * pcl_num * pcl_digit_num * 2
				+ sizeof(SortBin) * (1 + thread_num) * block_num;
		}
		void init(void* shared_mem, size_t thread_num,
				  size_t pcl_num, size_t pcl_digit_num);
	protected:
		CacheAlignedMem self_mem;
		inline void clear() noexcept { self_mem.free(); }
	};

	class SortParticleTask : public tbb::task
	{
	protected:
		const SortParticleMem& sort_mem;
		const size_t pcl_digit_num;
		size_t &valid_pcl_num;
	public:
		SortParticleTask(
			SortParticleMem &sm,
			size_t p_dgt_num,
			size_t &vp_num) :
			sort_mem(sm),
			pcl_digit_num(p_dgt_num),
			valid_pcl_num(vp_num)
		{ assert(vp_num > 0); }
		~SortParticleTask() {}
		tbb::task* execute() override;
	};
}

#endif