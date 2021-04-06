#include "SimulationsOMP_pcp.h"

#include <assert.h>
#include <iostream>
#include "tbb/task_arena.h"

#include "DivideTask.hpp"
#include "SortTriMeshNodeTask.h"

void SortTriMeshNodeMem::init(
	size_t elem_num,
	size_t node_num,
	RadixBinBlockMemArray& rbbs)
{
	set_radix_bin_block(rbbs);
	const size_t three_elem_num = elem_num * 3;
	char *cur_mem = data_mem.alloc<char>(
		  max_block_num * sizeof(ValidElemBlock)
		+ elem_num * 2 * sizeof(size_t)
		+ (three_elem_num * 4 + 4) * sizeof(size_t)
		+ MSDRadixSortUtils::cache_line_size * 5);
	valid_elem_blocks = (ValidElemBlock *)cur_mem;
	cur_mem += max_block_num * sizeof(ValidElemBlock);
	ori_elems = (size_t *)cur_mem;
	cur_mem = cache_aligned(cur_mem + elem_num * sizeof(size_t));
	res_elems = (size_t *)cur_mem;
	cur_mem = cache_aligned(cur_mem + elem_num * sizeof(size_t));
	radix_keys0 = (size_t *)cur_mem + 1;
	cur_mem = cache_aligned(cur_mem + (three_elem_num + 2) * sizeof(size_t));
	radix_keys1 = (size_t*)cur_mem + 1;
	cur_mem = cache_aligned(cur_mem + (three_elem_num + 2) * sizeof(size_t));
	radix_vals0 = (size_t*)cur_mem;
	cur_mem = cache_aligned(cur_mem + three_elem_num * sizeof(size_t));
	radix_vals1 = (size_t*)cur_mem;
	radix_keys0[-1] = SIZE_MAX;
	radix_keys0[three_elem_num] = SIZE_MAX;
	radix_keys1[-1] = SIZE_MAX;
	radix_keys1[three_elem_num] = SIZE_MAX;
	set_digit_num(node_num);
}

void SortTriMeshNodeTask::ScanPcl::operator() (size_t blk_id) const
{
	using MSDRadixSortUtils::block_low;
	using MSDRadixSortUtils::digit;

	size_t e_id;
	size_t p_id0 = block_low(blk_id, block_num, pcl_num);
	e_id = pcl_in_elems[p_id0];
	while (p_id0 != SIZE_MAX && e_id == pcl_in_elems[--p_id0]);
	++p_id0;
	assert(p_id0 <= pcl_num);
	size_t p_id1 = block_low(blk_id + 1, block_num, pcl_num);
	e_id = pcl_in_elems[p_id1];
	while (p_id1 != SIZE_MAX && e_id == pcl_in_elems[--p_id1]);
	++p_id1;
	assert(p_id1 <= pcl_num);

	RadixBin& bin = radix_bin_block[blk_id];
	ValidElemBlock& ve_blk = valid_elem_blocks[blk_id];
	memset(&bin, 0, sizeof(RadixBin));
	e_id = pcl_in_elems[p_id0];
	size_t* const o_elems = ori_elems + e_id;
	ve_blk.ori_offset = e_id;
	ve_blk.num = 0;
	for (size_t p_id = p_id0; p_id < p_id1; ++p_id)
	{
		if (e_id != pcl_in_elems[p_id + 1])
		{
			o_elems[ve_blk.num++] = e_id;
			const ElemNodeIndex& eni = elem_node_ids[e_id];
			++bin.bin[digit(eni.n1, node_digit_pos)];
			++bin.bin[digit(eni.n2, node_digit_pos)];
			++bin.bin[digit(eni.n3, node_digit_pos)];
			e_id = pcl_in_elems[p_id + 1];
		}
	}
}

void SortTriMeshNodeTask::FormElemAndNodeArray::operator() (size_t blk_id) const
{
	using MSDRadixSortUtils::digit;

	size_t i, o_blk_id;
	RadixBin bin;
	const RadixBin& my_bin = radix_bin_block[blk_id];
	for (i = 0; i < MSDRadixSortUtils::radix_bin_num; ++i)
		bin.bin[i] = my_bin.bin[i];
	for (o_blk_id = 0; o_blk_id < blk_id; ++o_blk_id)
	{
		const RadixBin& o_sbin = radix_bin_block[o_blk_id];
		for (i = 0; i < MSDRadixSortUtils::radix_bin_num; ++i)
			bin.bin[i] += o_sbin.bin[i];
	}
	for (o_blk_id = blk_id + 1; o_blk_id < block_num; ++o_blk_id)
	{
		const RadixBin& o_sbin = radix_bin_block[o_blk_id];
		for (i = 1; i < MSDRadixSortUtils::radix_bin_num; ++i)
			bin.bin[i] += o_sbin.bin[i - 1];
	}

	ValidElemBlock& ve_blk = valid_elem_blocks[blk_id];
	const size_t* const o_elems = ori_elems + ve_blk.ori_offset;
	size_t ve_res_offset = ve_blk.res_offset;
	size_t pos_id;
	for (size_t ve_id = ve_blk.num; ve_id--;)
	{
		const size_t e_id = o_elems[ve_id];
		res_elems[--ve_res_offset] = e_id;
		const ElemNodeIndex& eni = elem_node_ids[e_id];
		pos_id = --bin.bin[digit(eni.n3, node_digit_pos)];
		out_node_has_elem[pos_id] = eni.n3;
		out_node_elem_pair[pos_id] = 3 * ve_id + 2;
		pos_id = --bin.bin[digit(eni.n2, node_digit_pos)];
		out_node_has_elem[pos_id] = eni.n2;
		out_node_elem_pair[pos_id] = 3 * ve_id + 1;
		pos_id = --bin.bin[digit(eni.n1, node_digit_pos)];
		out_node_has_elem[pos_id] = eni.n1;
		out_node_elem_pair[pos_id] = 3 * ve_id;
	}
}

tbb::task* SortTriMeshNodeTask::execute()
{
	using MSDRadixSortUtils::digit;

	SortTriMeshNodeMem& snm = static_cast<SortTriMeshNodeMem &>(sort_mem);
	const size_t radix_id = digit_pos & 1;
	res_elems = snm.res_elems;
	MSDRadixSortUtils::RadixBin bin;
	size_t i, j, e_id, pos_id;
	// scan in parallel
	if (pcl_num > serial_sort_max_data_num)
	{
		out_node_has_elem = snm.radix_keys[radix_id ^ 1];
		out_node_elem_pair = snm.radix_vals[radix_id ^ 1];
		if (pcl_num > min_pcl_num_per_block * 2)
		{
			block_num = (pcl_num + min_pcl_num_per_block - 1)
						/ min_pcl_num_per_block;
			if (block_num > snm.max_block_num)
				block_num = snm.max_block_num;
			ori_elems = snm.ori_elems;
			const size_t my_th_id = tbb::task_arena::current_thread_index();
			RadixBinBlockMem& th_radix_bin_block = snm.thread_radix_bin_block[my_th_id];
			radix_bin_block = th_radix_bin_block.alloc();
			valid_elem_blocks = snm.valid_elem_blocks;

			// scan data
			set_ref_count(2);
			spawn(*new(allocate_child())
				DivideTask<ScanPcl, 2>(0, block_num - 1, scan_pcl));
			scan_pcl(block_num - 1);
			wait_for_all();

			// summarize output
			valid_elem_blocks[0].res_offset = valid_elem_blocks[0].num;
			for (i = 1; i < block_num; ++i)
				valid_elem_blocks[i].res_offset = valid_elem_blocks[i - 1].res_offset + valid_elem_blocks[i].num;
			valid_elem_num = valid_elem_blocks[block_num - 1].res_offset;
			*const_cast<size_t *>(&data_num) = 3 * valid_elem_blocks[block_num - 1].res_offset;

			// move data
			set_ref_count(2);
			spawn(*new(allocate_child())
				DivideTask<FormElemAndNodeArray, 2>(
					0, block_num - 1, form_elem_and_node_array));
			form_elem_and_node_array(block_num - 1);
			wait_for_all();

			th_radix_bin_block.free(radix_bin_block);
		}
		else // scan serially
		{
			memset(&bin, 0, sizeof(bin));
			*const_cast<size_t *>(&data_num) = 0;
			e_id = pcl_in_elems[0];
			for (i = 0; i < pcl_num; ++i)
			{
				if (e_id != pcl_in_elems[i + 1])
				{
					res_elems[(*const_cast<size_t *>(&data_num))++] = e_id;
					const ElemNodeIndex& eni = elem_node_ids[e_id];
					++bin.bin[digit(eni.n1, node_digit_pos)];
					++bin.bin[digit(eni.n2, node_digit_pos)];
					++bin.bin[digit(eni.n3, node_digit_pos)];
					e_id = pcl_in_elems[i + 1];
				}
			}
			for (i = 1; i < MSDRadixSortUtils::radix_bin_num; ++i)
				bin.bin[i] += bin.bin[i - 1];
			for (i = data_num; i--;)
			{
				e_id = res_elems[i];
				const ElemNodeIndex& eni = elem_node_ids[e_id];
				pos_id = --bin.bin[digit(eni.n3, node_digit_pos)];
				out_node_has_elem[pos_id] = eni.n3;
				out_node_elem_pair[pos_id] = 3 * e_id + 2;
				pos_id = --bin.bin[digit(eni.n2, node_digit_pos)];
				out_node_has_elem[pos_id] = eni.n2;
				out_node_elem_pair[pos_id] = 3 * e_id + 1;
				pos_id = --bin.bin[digit(eni.n1, node_digit_pos)];
				out_node_has_elem[pos_id] = eni.n1;
				out_node_elem_pair[pos_id] = 3 * e_id;
			}
			valid_elem_num = data_num;
			*const_cast<size_t *>(&data_num) *= 3;

			if (digit_pos)
			{
				tbb::empty_task& c =
					*new(allocate_continuation()) tbb::empty_task;
				c.set_ref_count(4);
				c.spawn(*new(c.allocate_child())
					ChildSpawner<3>(
						bin.bin,
						bin.bin[0x40],
						start_pos,
						digit_pos - 1,
						sort_mem));
				c.spawn(*new(c.allocate_child())
					ChildSpawner<3>(
						bin.bin + 0x40,
						bin.bin[0x40 * 2],
						start_pos,
						digit_pos - 1,
						sort_mem));
				c.spawn(*new(c.allocate_child())
					ChildSpawner<3>(
						bin.bin + 0x40 * 2,
						bin.bin[0x40 * 3],
						start_pos,
						digit_pos - 1,
						sort_mem));
				return new(c.allocate_child())
					ChildSpawner<3>(
						bin.bin + 0x40 * 3,
						data_num,
						start_pos,
						digit_pos - 1,
						sort_mem);
			}
			return nullptr;
		}
	}
	
	// sort serially
	size_t *in_node_has_elem = snm.radix_keys[radix_id];
	size_t *in_node_elem_pair = snm.radix_vals[radix_id];
	*const_cast<size_t *>(&data_num) = 0;
	e_id = pcl_in_elems[0];
	for (i = 0; i < pcl_num; ++i)
	{
		if (e_id != pcl_in_elems[i + 1])
		{
			res_elems[data_num] = e_id;
			const ElemNodeIndex& eni = elem_node_ids[e_id];
			in_node_has_elem[data_num * 3] = eni.n1;
			in_node_has_elem[data_num * 3 + 1] = eni.n2;
			in_node_has_elem[data_num * 3 + 2] = eni.n3;
			in_node_elem_pair[data_num * 3] = data_num * 3;
			in_node_elem_pair[data_num * 3 + 1] = data_num * 3 + 1;
			in_node_elem_pair[data_num * 3 + 2] = data_num * 3 + 2;
			++(*const_cast<size_t *>(&data_num));
			e_id = pcl_in_elems[i + 1];
		}
	}
	valid_elem_num = data_num;
	*const_cast<size_t *>(&data_num) *= 3;
	if (data_num > insertion_sort_max_data_num) // radix sort
	{
		out_node_has_elem = snm.radix_keys[radix_id ^ 1];
		out_node_elem_pair = snm.radix_vals[radix_id ^ 1];
		for (size_t d_id = 0; d_id < digit_pos + 1; ++d_id)
		{
			memset(&bin, 0, sizeof(bin));
			for (i = 0; i < data_num; ++i)
				++bin.bin[digit(in_node_has_elem[i], d_id)];
			for (i = 1; i < MSDRadixSortUtils::radix_bin_num; ++i)
				bin.bin[i] += bin.bin[i - 1];
			for (i = data_num; i--;)
			{
				pos_id = --bin.bin[digit(in_node_has_elem[i], d_id)];
				out_node_has_elem[pos_id] = in_node_has_elem[i];
				out_node_elem_pair[pos_id] = in_node_elem_pair[i];
			}
#define swap_pointer(pt1, pt2) \
			(pt1) = (size_t *)(size_t(pt1) ^ size_t(pt2)); \
			(pt2) = (size_t *)(size_t(pt1) ^ size_t(pt2)); \
			(pt1) = (size_t *)(size_t(pt1) ^ size_t(pt2))
			swap_pointer(in_node_has_elem, out_node_has_elem);
			swap_pointer(in_node_elem_pair, out_node_elem_pair);
		}
	}
	else if (data_num > 2) // insertion sort
	{
		out_node_has_elem = sort_mem.res_keys;
		out_node_elem_pair = sort_mem.res_vals;
		const size_t data_num_min_1 = data_num - 1;
		for (i = 0; i < data_num_min_1; ++i)
		{
			size_t min_key = in_node_has_elem[i];
			size_t min_val = in_node_elem_pair[i];
			size_t min_key_id = i;
			for (j = i + 1; j < data_num; ++j)
			{
				if (min_key > in_node_has_elem[j])
				{
					min_key = in_node_has_elem[j];
					min_val = in_node_elem_pair[j];
					min_key_id = j;
				}
			}
			const_cast<size_t*>(in_node_has_elem)[min_key_id] = in_node_has_elem[i];
			out_node_has_elem[i] = min_key;
			const_cast<size_t*>(in_node_elem_pair)[min_key_id] = in_node_elem_pair[i];
			out_node_elem_pair[i] = min_val;
		}
		out_node_has_elem[data_num_min_1] = in_node_has_elem[data_num_min_1];
		out_node_elem_pair[data_num_min_1] = in_node_elem_pair[data_num_min_1];
	}
	else if (data_num == 2)
	{
		out_node_has_elem = sort_mem.res_keys;
		out_node_elem_pair = sort_mem.res_vals;
		if (in_node_has_elem[0] < in_node_has_elem[1])
		{
			out_node_has_elem[0] = in_node_has_elem[0];
			out_node_has_elem[1] = in_node_has_elem[1];
			out_node_elem_pair[0] = in_node_elem_pair[0];
			out_node_elem_pair[1] = in_node_elem_pair[1];
		}
		else
		{
			size_t tmp = in_node_has_elem[1];
			out_node_has_elem[1] = in_node_has_elem[0];
			out_node_has_elem[0] = tmp;
			tmp = in_node_elem_pair[1];
			out_node_elem_pair[1] = in_node_elem_pair[0];
			out_node_elem_pair[0] = tmp;
		}
	}
	else if (data_num == 1)
	{
		sort_mem.res_keys[0] = in_node_has_elem[0];
		sort_mem.res_vals[0] = in_node_elem_pair[0];
	}
	return nullptr;
}
