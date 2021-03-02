#include "SimulationsOMP_pcp.h"

#include "tbb/task_arena.h"

#include "DivideTask.hpp"
#include "MSDRadixSortTask.h"

namespace MSDRadixSortUtils
{
	void ScanAndFormBin::operator() (size_t blk_id)
	{
		RadixBin& bin = radix_bin_block[blk_id];
		memset(&bin, 0, sizeof(RadixBin));
		size_t i;
		const size_t data_id1 = block_low(blk_id + 1, block_num, data_num);
		for (i = block_low(blk_id, block_num, data_num); i < data_id1; ++i)
			++bin.bin[digit(keys[i], digit_pos)];
		for (i = 1; i < radix_bin_num; ++i)
			bin.bin[i] += bin.bin[i - 1];
	}

	void MoveAccToBin::operator() (size_t blk_id)
	{
		RadixBin bin;
		size_t i, o_blk_id;
		const RadixBin& my_bin = radix_bin_block[blk_id];
		for (i = 0; i < radix_bin_num; ++i)
			bin.bin[i] = my_bin.bin[i];
		for (o_blk_id = 0; o_blk_id < blk_id; ++o_blk_id)
		{
			const RadixBin& o_sbin = radix_bin_block[o_blk_id];
			for (i = 0; i < radix_bin_num; ++i)
				bin.bin[i] += o_sbin.bin[i];
		}
		for (o_blk_id = blk_id + 1; o_blk_id < block_num; ++o_blk_id)
		{
			const RadixBin& o_sbin = radix_bin_block[o_blk_id];
			for (i = 1; i < radix_bin_num; ++i)
				bin.bin[i] += o_sbin.bin[i - 1];
		}
		const size_t data_id0 = block_low(blk_id, block_num, data_num);
		for (i = block_low(blk_id + 1, block_num, data_num); i-- != data_id0;)
		{
			const size_t pos_id = --bin.bin[digit(in_keys[i], digit_pos)];
			out_keys[pos_id] = in_keys[i];
			out_vals[pos_id] = in_vals[i];
		}
	}

	void MoveAccToBin::operator() (size_t blk_id, RadixBin& bin)
	{
		size_t i, o_blk_id;
		const RadixBin& my_bin = radix_bin_block[blk_id];
		for (i = 0; i < radix_bin_num; ++i)
			bin.bin[i] = my_bin.bin[i];
		for (o_blk_id = 0; o_blk_id < blk_id; ++o_blk_id)
		{
			const RadixBin& o_sbin = radix_bin_block[o_blk_id];
			for (i = 0; i < radix_bin_num; ++i)
				bin.bin[i] += o_sbin.bin[i];
		}
		for (o_blk_id = blk_id + 1; o_blk_id < block_num; ++o_blk_id)
		{
			const RadixBin& o_sbin = radix_bin_block[o_blk_id];
			for (i = 1; i < radix_bin_num; ++i)
				bin.bin[i] += o_sbin.bin[i - 1];
		}
		const size_t data_id0 = block_low(blk_id, block_num, data_num);
		for (i = block_low(blk_id + 1, block_num, data_num); i-- != data_id0;)
		{
			const size_t pos_id = --bin.bin[digit(in_keys[i], digit_pos)];
			out_keys[pos_id] = in_keys[i];
			out_vals[pos_id] = in_vals[i];
		}
	}
}

tbb::task* MSDRadixSortTask::execute()
{
	using MSDRadixSortUtils::digit;
	using MSDRadixSortUtils::radix_bin_num;

	RadixBin bin;
	size_t i;
	const size_t radix_id = digit_pos & 1;
	size_t *in_keys = sort_mem.radix_keys[radix_id] + start_pos;
	size_t *in_vals = sort_mem.radix_vals[radix_id] + start_pos;
	size_t* out_keys, *out_vals;
	if (data_num > MSDRadixSortUtils::serial_sort_max_data_num)
	{
		out_keys = sort_mem.radix_keys[radix_id ^ 1] + start_pos;
		out_vals = sort_mem.radix_vals[radix_id ^ 1] + start_pos;
		if (data_num > MSDRadixSortUtils::min_data_num_per_block * 2)
		{
			size_t block_num = (data_num + MSDRadixSortUtils::min_data_num_per_block - 1)
							/ MSDRadixSortUtils::min_data_num_per_block;
			if (block_num > sort_mem.max_block_num)
				block_num = sort_mem.max_block_num;
			const size_t my_th_id = tbb::task_arena::current_thread_index();
			RadixBinBlockMem& th_radix_bin_block = sort_mem.thread_radix_bin_block[my_th_id];
			RadixBin *const radix_bin_block = th_radix_bin_block.alloc();
			// scan data
			ScanAndFormBin scan_and_fill_bin(
				block_num, data_num, digit_pos,
				radix_bin_block, in_keys);
			set_ref_count(2);
			spawn(*new(allocate_child())
				DivideTask<ScanAndFormBin, 2>(
					1, block_num,
					scan_and_fill_bin));
			scan_and_fill_bin(0);
			wait_for_all();
			// move data
			MoveAccToBin move_acc_to_bin(
				block_num, data_num, digit_pos,
				radix_bin_block,
				in_keys, in_vals,
				out_keys, out_vals);
			set_ref_count(2);
			spawn(*new(allocate_child())
				DivideTask<MoveAccToBin, 2>(
					1, block_num,
					move_acc_to_bin));
			move_acc_to_bin(0, bin);
			wait_for_all();
			th_radix_bin_block.free(radix_bin_block);
		}
		else // divide serially and sort parallely
		{
			memset(&bin, 0, sizeof(bin));
			for (i = 0; i < data_num; ++i)
				++bin.bin[digit(in_keys[i], digit_pos)];
			for (i = 1; i < radix_bin_num; ++i)
				bin.bin[i] += bin.bin[i - 1];
			for (i = data_num; i--;)
			{
				const size_t pos_id = --bin.bin[digit(in_keys[i], digit_pos)];
				out_keys[pos_id] = in_keys[i];
				out_vals[pos_id] = in_vals[i];
			}
		}

		if (digit_pos) // not the last digit
		{
			tbb::empty_task& c = *new(allocate_continuation()) tbb::empty_task;
			c.set_ref_count(4);
			c.spawn(*new(c.allocate_child()) ChildSpawner<3>(bin.bin, bin.bin[0x40], *this));
			c.spawn(*new(c.allocate_child()) ChildSpawner<3>(bin.bin + 0x40, bin.bin[0x40*2], *this));
			c.spawn(*new(c.allocate_child()) ChildSpawner<3>(bin.bin + 0x40*2, bin.bin[0x40*3], *this));
			return new(c.allocate_child()) ChildSpawner<3>(bin.bin + 0x40*3, data_num, *this);
		}
		return nullptr;
	}

	// sort serially
	if (data_num > MSDRadixSortUtils::insertion_sort_max_data_num) // radix sort
	{
		out_keys = sort_mem.radix_keys[radix_id ^ 1] + start_pos;
		out_vals = sort_mem.radix_vals[radix_id ^ 1] + start_pos;
		for (size_t d_id = 0; d_id < digit_pos + 1; ++d_id)
		{
			memset(&bin, 0, sizeof(bin));
			for (i = 0; i < data_num; ++i)
				++bin.bin[digit(in_keys[i], d_id)];
			for (i = 1; i < radix_bin_num; ++i)
				bin.bin[i] += bin.bin[i - 1];
			for (i = data_num; i--;)
			{
				const size_t pos_id = --bin.bin[digit(in_keys[i], d_id)];
				out_keys[pos_id] = in_keys[i];
				out_vals[pos_id] = in_vals[i];
			}
#define swap_pointer(pt1, pt2) \
			(pt1) = (size_t *)(size_t(pt1) ^ size_t(pt2)); \
			(pt2) = (size_t *)(size_t(pt1) ^ size_t(pt2)); \
			(pt1) = (size_t *)(size_t(pt1) ^ size_t(pt2))
			swap_pointer(in_keys, out_keys);
			swap_pointer(in_vals, out_vals);
		}
	}
	else if (data_num > 2) // insertion sort
	{
		out_keys = sort_mem.res_keys + start_pos;
		out_vals = sort_mem.res_vals + start_pos;
		const size_t data_num_min_1 = data_num - 1;
		for (size_t i = 0; i < data_num_min_1; ++i)
		{
			size_t min_key = in_keys[i];
			size_t min_val = in_vals[i];
			size_t min_key_id = i;
			for (size_t j = i + 1; j < data_num; ++j)
			{
				if (min_key > in_keys[j])
				{
					min_key = in_keys[j];
					min_val = in_vals[j];
					min_key_id = j;
				}
			}
			const_cast<size_t*>(in_keys)[min_key_id] = in_keys[i];
			out_keys[i] = min_key;
			const_cast<size_t*>(in_vals)[min_key_id] = in_vals[i];
			out_vals[i] = min_val;
		}
		out_keys[data_num_min_1] = in_keys[data_num_min_1];
		out_vals[data_num_min_1] = in_vals[data_num_min_1];
	}
	else if (data_num == 2)
	{
		out_keys = sort_mem.res_keys + start_pos;
		out_vals = sort_mem.res_vals + start_pos;
		if (in_keys[0] < in_keys[1])
		{
			out_keys[0] = in_keys[0];
			out_keys[1] = in_keys[1];
			out_vals[0] = in_vals[0];
			out_vals[1] = in_vals[1];
		}
		else
		{
			size_t tmp = in_keys[1];
			out_keys[1] = in_keys[0];
			out_keys[0] = tmp;
			tmp = in_vals[1];
			out_vals[1] = in_vals[0];
			out_vals[0] = tmp;
		}
	}
	else if (data_num == 1)
	{
		sort_mem.res_keys[start_pos] = in_keys[0];
		sort_mem.res_vals[start_pos] = in_vals[0];
	}
	return nullptr;
}
