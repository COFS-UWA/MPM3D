#ifndef __MSD_Radix_Sort_Utils_h__
#define __MSD_Radix_Sort_Utils_h__

#include "tbb/task.h"
#include "CacheAlignedMem.h"

namespace MSDRadixSortUtils
{
	constexpr size_t cache_line_size = 0x40;
	
	constexpr size_t radix_digit_num = 8;
	constexpr size_t radix_bin_num = 1 << radix_digit_num;

	constexpr size_t max_block_num_per_thread = 2;
	constexpr size_t min_data_num_per_block = 100000;
	constexpr size_t serial_sort_max_data_num = 10000;
	constexpr size_t insertion_sort_max_data_num = 10;

	inline constexpr size_t max_digit_num(size_t max_data) noexcept
	{
		size_t digit_num = 0;
		while (max_data)
		{
			++digit_num;
			max_data >>= radix_digit_num;
		}
		return digit_num;
	}
	
	inline size_t digit(size_t data, size_t digit_pos) noexcept
	{
		constexpr size_t radix_digit_mask = ~(SIZE_MAX << radix_digit_num);
		return (data >> (digit_pos * radix_digit_num)) & radix_digit_mask;
	}
	
	inline size_t block_low(size_t block_id, size_t block_num, size_t data_num) noexcept
	{ return block_id * data_num / block_num; }

	struct RadixBin { size_t bin[radix_bin_num]; };

	class RadixBinBlockMem
	{
	protected:
		size_t bin_num_per_block;
		size_t block_num_per_page;
		char* pages_mem;
		RadixBin* free_blocks;
		char padding[cache_line_size
			- sizeof(bin_num_per_block)
			- sizeof(block_num_per_page)
			- sizeof(pages_mem)
			- sizeof(free_blocks)];

		void alloc_block_mem();

	public:
		RadixBinBlockMem();
		RadixBinBlockMem(size_t bin_num_per_blk,
			size_t blk_num_per_page = 1);
		~RadixBinBlockMem() { clear(); }

		void clear() noexcept;
		void reset_bin_num_per_block(size_t bin_num_per_blk) noexcept;
		void set_block_num_per_page(size_t blk_num_per_page) noexcept;

		inline RadixBin* alloc()
		{
			if (free_blocks == nullptr)
				alloc_block_mem();
			assert(free_blocks);
			RadixBin* block_tmp = free_blocks;
			free_blocks = *(RadixBin**)free_blocks;
			return block_tmp;
		}

		inline void free(RadixBin* bins) noexcept
		{
			*(RadixBin**)bins = free_blocks;
			free_blocks = bins;
		}
	};

	class RadixBinBlockMemArray
	{
	protected:
		CacheAlignedMem mem;
		size_t thread_num;
	public:
		RadixBinBlockMemArray() {}
		~RadixBinBlockMemArray() { clear(); }
		void init(size_t _thread_num, size_t block_per_page);
		void clear();
		inline RadixBinBlockMem* get_array() noexcept
		{ return (RadixBinBlockMem *)mem.aligned_address(); }
		inline size_t get_max_block_num() const noexcept
		{ return thread_num * max_block_num_per_thread; }
	};

	struct SortMem
	{
	public:
		union
		{
			size_t* radix_keys[2];
			struct { size_t* radix_keys0, *radix_keys1; };
			struct { size_t* _padding_keys, *res_keys; };
		};
		size_t* ori_keys;
		union
		{
			size_t* radix_vals[2];
			struct { size_t* radix_vals0, * radix_vals1; };
			struct { size_t* _padding_vals, * res_vals; };
		};
		size_t* ori_vals;
		size_t digit_num;

		size_t max_block_num;
		RadixBinBlockMem* thread_radix_bin_block; // thread_num

		void set_digit_num(size_t max_data);
		void set_radix_bin_block(RadixBinBlockMemArray &rbbs);
	};
}

#endif