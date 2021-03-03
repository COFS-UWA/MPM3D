#include "SimulationsOMP_pcp.h"

#include <assert.h>
#include "MSDRadixSortUtils.h"

namespace MSDRadixSortUtils
{
	RadixBinBlockMem::RadixBinBlockMem() :
		bin_num_per_block(1),
		block_num_per_page(1),
		pages_mem(nullptr),
		free_blocks(nullptr) {}

	RadixBinBlockMem::RadixBinBlockMem(
		size_t bin_num_per_blk,
		size_t blk_num_per_page) :
		bin_num_per_block(bin_num_per_blk),
		block_num_per_page(blk_num_per_page),
		pages_mem(nullptr),
		free_blocks(nullptr)
	{
		assert(bin_num_per_blk);
		assert(blk_num_per_page);
		alloc_block_mem();
	}

	void RadixBinBlockMem::clear() noexcept
	{
		while (pages_mem)
		{
			char* const page_tmp = pages_mem;
			pages_mem = *(char**)pages_mem;
			delete[] page_tmp;
		}
		free_blocks = nullptr;
	}

	void RadixBinBlockMem::reset_bin_num_per_block(
		size_t bin_num_per_blk) noexcept
	{
		assert(bin_num_per_blk);
		bin_num_per_block = bin_num_per_blk;
		clear();
	}

	void RadixBinBlockMem::set_block_num_per_page(
		size_t blk_num_per_page) noexcept
	{
		assert(blk_num_per_page);
		block_num_per_page = blk_num_per_page;
	}

	void RadixBinBlockMem::alloc_block_mem()
	{
#define Cache_Alignment_Padding(address) \
	((cache_line_size - (size_t(address) & (cache_line_size - 1))) & (cache_line_size - 1))
#define Cache_Aligned_Address(address) \
	((char *)(size_t(address) + Cache_Alignment_Padding(size_t(address))))

		char* const new_page = new char[cache_line_size + sizeof(char*)
			+ sizeof(RadixBin) * bin_num_per_block * block_num_per_page];
		// add to mem page list
		*(char**)new_page = pages_mem;
		pages_mem = new_page;
		// add to block list
		RadixBin* new_block = (RadixBin*)(Cache_Aligned_Address(new_page + sizeof(char*)));
		for (size_t blk_id = 0; blk_id < block_num_per_page; ++blk_id)
		{
			*(RadixBin**)new_block = free_blocks;
			free_blocks = new_block;
			new_block += bin_num_per_block;
		}
	}

	void RadixBinBlockMemArray::clear()
	{
		if (mem.raw_address())
		{
			RadixBinBlockMem *const thread_radix_bin_block
				= (RadixBinBlockMem *)mem.aligned_address();
			for (size_t th_id = 0; th_id < thread_num; ++th_id)
				thread_radix_bin_block[th_id].~RadixBinBlockMem();
			mem.free();
		}
	}

	void RadixBinBlockMemArray::init(
		size_t _thread_num,
		size_t block_per_page)
	{
		clear();
		if (_thread_num)
		{
			thread_num = _thread_num;
			const size_t max_block_num = thread_num * max_block_num_per_thread;
			RadixBinBlockMem* const thread_radix_bin_block
				= mem.alloc<RadixBinBlockMem>(thread_num);
			for (size_t th_id = 0; th_id < thread_num; ++th_id)
				new (thread_radix_bin_block + th_id)
					RadixBinBlockMem(max_block_num, block_per_page);
		}
	}

	void SortMem::set_digit_num(size_t max_data)
	{
		digit_num = max_digit_num(max_data);
		ori_keys = radix_keys[(digit_num & 1) ^ 1];
		ori_vals = radix_vals[(digit_num & 1) ^ 1];
	}

	void SortMem::set_radix_bin_block(RadixBinBlockMemArray& rbbs)
	{
		max_block_num = rbbs.get_max_block_num();
		thread_radix_bin_block = rbbs.get_array();
	}
}
