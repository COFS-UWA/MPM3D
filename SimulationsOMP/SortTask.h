#ifndef __Radix_Sort_Task_h__
#define __Radix_Sort_Task_h__

#include "CacheAlignedMem.h"
#include "SortUtils.h"

namespace SortUtils
{
	namespace Internal
	{
		constexpr size_t cache_line_size = 0x40; // 64
		constexpr size_t block_div_thread_num_ratio = 2;
	}

	struct SortMem
	{
	public:
		// sorting array
		union
		{
			struct
			{
				size_t *in_key; // key0
				size_t *key1;
				size_t *key2;
				size_t *key3;
				size_t *key4;
				size_t *key5;
				size_t* key6;
				size_t* key7;
				size_t *out_key; // key8
				size_t *mid_key; // key9
				size_t *in_val; // val0
				size_t *val1;
				size_t *val2;
				size_t *val3;
				size_t* val4;
				size_t *val5;
				size_t* val6;
				size_t* val7;
				size_t *out_val; // val8
				size_t *mid_val; // val9
			};
			struct
			{
				size_t* keys[10];
				size_t* vals[10];
			};
		};
		// threadwise sort bin
		size_t init_block_num;
		SortBin **thread_bins;
		
		SortMem() {}
		~SortMem() { clear(); }
		void init(size_t data_num, size_t thread_num)
		{
			clear();
			char *cur_mem = (char *)mem.alloc(
				  sizeof(size_t) * data_num * 20
				+ sizeof(SortBin *) * thread_num
				+ sizeof(SortBin) * thread_num * init_block_num
				+ Internal::cache_line_size * 22);
			// key mem
			in_key = (size_t *)cur_mem;
			cur_mem = cache_aligned(cur_mem + sizeof(size_t) * data_num);
			key1 = (size_t*)cur_mem;
			cur_mem = cache_aligned(cur_mem + sizeof(size_t) * data_num);
			key2 = (size_t*)cur_mem;
			cur_mem = cache_aligned(cur_mem + sizeof(size_t) * data_num);
			key3 = (size_t*)cur_mem;
			cur_mem = cache_aligned(cur_mem + sizeof(size_t) * data_num);
			key4 = (size_t*)cur_mem;
			cur_mem = cache_aligned(cur_mem + sizeof(size_t) * data_num);
			key5 = (size_t*)cur_mem;
			cur_mem = cache_aligned(cur_mem + sizeof(size_t) * data_num);
			key6 = (size_t*)cur_mem;
			cur_mem = cache_aligned(cur_mem + sizeof(size_t) * data_num);
			key7 = (size_t*)cur_mem;
			cur_mem = cache_aligned(cur_mem + sizeof(size_t) * data_num);
			out_key = (size_t*)cur_mem;
			cur_mem = cache_aligned(cur_mem + sizeof(size_t) * data_num);
			mid_key = (size_t*)cur_mem;
			cur_mem = cache_aligned(cur_mem + sizeof(size_t) * data_num);
			// val mem
			in_val = (size_t *)cur_mem;
			cur_mem = cache_aligned(cur_mem + sizeof(size_t) * data_num);
			val1 = (size_t*)cur_mem;
			cur_mem = cache_aligned(cur_mem + sizeof(size_t) * data_num);
			val1 = (size_t*)cur_mem;
			cur_mem = cache_aligned(cur_mem + sizeof(size_t) * data_num);
			val2 = (size_t*)cur_mem;
			cur_mem = cache_aligned(cur_mem + sizeof(size_t) * data_num);
			val3 = (size_t*)cur_mem;
			cur_mem = cache_aligned(cur_mem + sizeof(size_t) * data_num);
			val4 = (size_t*)cur_mem;
			cur_mem = cache_aligned(cur_mem + sizeof(size_t) * data_num);
			val5 = (size_t*)cur_mem;
			cur_mem = cache_aligned(cur_mem + sizeof(size_t) * data_num);
			val6 = (size_t*)cur_mem;
			cur_mem = cache_aligned(cur_mem + sizeof(size_t) * data_num);
			val7 = (size_t*)cur_mem;
			cur_mem = cache_aligned(cur_mem + sizeof(size_t) * data_num);
			out_val = (size_t*)cur_mem;
			cur_mem = cache_aligned(cur_mem + sizeof(size_t) * data_num);
			mid_val = (size_t*)cur_mem;
			cur_mem = cache_aligned(cur_mem + sizeof(size_t) * data_num);
			// thread bin
			init_block_num = thread_num * Internal::block_div_thread_num_ratio;
			thread_bins = (SortBin **)cur_mem;
			cur_mem = cache_aligned(cur_mem + sizeof(SortBin *) * data_num);
			for (size_t th_id = 0; th_id < thread_num; ++th_id)
			{
				thread_bins[th_id] = (SortBin *)cur_mem;
				cur_mem += sizeof(SortBin) * init_block_num;
			}
		}
	protected:
		CacheAlignedMem mem;
		inline void clear() { mem.free(); }
	};
	
	namespace Internal
	{
		class SortContinuation : public tbb::task
		{
		public:
			SortContinuation() {}
			~SortContinuation() {}
			tbb::task* execute() override { return nullptr; }
		};
	}

	namespace Internal
	{
		constexpr size_t serial_divide_max_data_num = 2;
		constexpr size_t min_pcl_num_per_block = 2;
		constexpr size_t serial_sort_max_data_num = 0;
	}

	class SortTask : public tbb::task
	{
	protected:
		const SortMem &sort_mem;
		const size_t start_id, data_num;
		const unsigned char digit_pos;

	public:
		SortTask(
			const SortMem& _sort_mem,
			const size_t _start_id,
			const size_t _data_num,
			const unsigned char _digit_pos);
		~SortTask();
		tbb::task* execute() override;
	};
}

#endif