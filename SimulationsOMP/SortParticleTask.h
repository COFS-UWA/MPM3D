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

		size_t pcl_digit_num;
		size_t valid_pcl_num;

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
			size_t thread_num,
			size_t pcl_num,
			size_t ori_pcl_num) noexcept
		{
			const size_t _pcl_digit_num = max_digit_num(ori_pcl_num);
			max_block_num = thread_num * SortParticleMem::max_block_num_div_thread_num;
			size_t block_num = (pcl_num + SortParticleMem::parallel_divide_min_pcl_num_per_block - 1)
				/ SortParticleMem::parallel_divide_min_pcl_num_per_block;
			if (block_num > max_block_num)
				block_num = max_block_num;
			return sizeof(size_t) * pcl_num * _pcl_digit_num * 2
				+ sizeof(SortBin) * (1 + thread_num) * block_num
				+ Cache_Alignment * _pcl_digit_num * 2;
		}
		void init(void* shared_mem, size_t thread_num,
				  size_t pcl_num, size_t ori_pcl_num);
	protected:
		CacheAlignedMem self_mem;
		inline void clear() noexcept { self_mem.free(); }
	};

	namespace Internal
	{
		inline size_t count_sort_pcl(
			size_t* const out_pcl_in_elem,
			size_t* const out_cur_to_prev_pcl,
			const size_t* const in_pcl_in_elem,
			const size_t* const in_cur_to_prev_pcl,
			const size_t pcl_num,
			const size_t pcl_digit_pos,
			SortBin& bin)
		{
#define __Cal_Key_Digit__(data, digit_pos) (((data) >> ((digit_pos) * 8)) & (0xFF))
			size_t* const c_bin = bin.count_bin;
			size_t* const s_bin = bin.sum_bin;
			memset(c_bin, 0, sizeof(size_t) * radix_bucket_num);
			size_t i, vp_num = 0;
			for (i = 0; i < pcl_num; ++i)
			{
				if (in_pcl_in_elem[i] != SIZE_MAX)
				{
					++vp_num;
					++c_bin[__Cal_Key_Digit__(in_pcl_in_elem[i], pcl_digit_pos)];
				}
			}
			s_bin[0] = c_bin[0];
			for (i = 1; i < radix_bucket_num; ++i)
				s_bin[i] = s_bin[i - 1] + c_bin[i];
			for (i = pcl_num; i--;)
			{
				if (in_pcl_in_elem[i] != SIZE_MAX)
				{
					const size_t pos = --s_bin[__Cal_Key_Digit__(in_pcl_in_elem[i], pcl_digit_pos)];
					out_pcl_in_elem[pos] = in_pcl_in_elem[i];
					out_cur_to_prev_pcl[pos] = in_cur_to_prev_pcl[i];
				}
			}
			return vp_num;
#undef __Cal_Key_Digit__
		}
	
		class CountSortPclTask : public tbb::task
		{
		protected:
			size_t* const out_pcl_in_elem;
			size_t* const out_cur_to_prev_pcl;
			const size_t* const in_pcl_in_elem;
			const size_t* const in_cur_to_prev_pcl;
			const size_t pcl_num;
			const size_t pcl_digit_pos;
			SortBin& bin;
		public:
			CountSortPclTask(
				size_t *_out_pcl_in_elem,
				size_t *_out_cur_to_prev_pcl,
				const size_t *_in_pcl_in_elem,
				const size_t *_in_cur_to_prev_pcl,
				const size_t p_num,
				const size_t p_dgt_pos,
				SortBin& _bin) :
				out_pcl_in_elem(_out_pcl_in_elem),
				out_cur_to_prev_pcl(_out_cur_to_prev_pcl),
				in_pcl_in_elem(_in_pcl_in_elem),
				in_cur_to_prev_pcl(_in_cur_to_prev_pcl),
				pcl_num(p_num),
				pcl_digit_pos(p_dgt_pos),
				bin(_bin) {}
			~CountSortPclTask() {}
			tbb::task* execute() override
			{
				count_sort_pcl(
					out_pcl_in_elem, out_cur_to_prev_pcl,
					in_pcl_in_elem, in_cur_to_prev_pcl,
					pcl_num, pcl_digit_pos, bin);
				return nullptr;
			}
		};
	}

	class SortParticleTask : public tbb::task
	{
	protected:
		SortParticleMem& sort_mem;
	public:
		SortParticleTask(SortParticleMem &sm) : sort_mem(sm) {}
		~SortParticleTask() {}
		tbb::task* execute() override;
	};
}

#endif