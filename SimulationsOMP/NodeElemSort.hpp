#ifndef __Node_Elem_Sort_hpp__
#define __Node_Elem_Sort_hpp__

#include <iostream>

#include "ParaUtil.h"
#include "ParallelReduceTask.hpp"
#include "ParallelForTask.hpp"
#include "ParaUtilSerialSort.h"
#include "ParaUtilMSDRadix16Sort.h"

#include "NodeElemSortInternal.hpp"

#define Block_Low(blk_id, blk_num, data_num) ((blk_id) * (data_num) / (blk_num))
#define Data_Digit(num, disp, mask) (((num) >> (disp)) & (mask))

namespace ParaUtil
{
namespace Internal
{
	template <size_t node_num_per_elem>
	inline void form_node_array(
		size_t node_ids[],
		size_t node_elem_offs[],
		size_t elem_id,
		const size_t elem_node_topologies[])
	{
		size_t ne_id = elem_id * node_num_per_elem;
		for (size_t n_id = 0; n_id < node_num_per_elem; ++n_id)
		{
			node_ids[n_id] = elem_node_topologies[ne_id];
			node_elem_offs[n_id] = ne_id++;
		}
	}

	template <>
	inline void form_node_array<4>(
		size_t node_ids[],
		size_t node_elem_offs[],
		size_t elem_id,
		const size_t elem_node_topologies[])
	{
		size_t ne_id = elem_id * 4;
		node_ids[0] = elem_node_topologies[ne_id];
		node_elem_offs[0] = ne_id;
		++ne_id;
		node_ids[1] = elem_node_topologies[ne_id];
		node_elem_offs[1] = ne_id;
		++ne_id;
		node_ids[2] = elem_node_topologies[ne_id];
		node_elem_offs[2] = ne_id;
		++ne_id;
		node_ids[3] = elem_node_topologies[ne_id];
		node_elem_offs[3] = ne_id;
	}

	template <size_t node_num_per_elem>
	tbb::task* ScanElem<node_num_per_elem>::operator()(
		tbb::task& parent,
		size_t task_id,
		ScanElemRes& result) const
	{
		// scan pcl to form node element array
		size_t e_id;
		// p_id 0
		size_t p_id0 = Block_Low(task_id, task_num, pcl_num);
		e_id = pcl_in_elems[p_id0];
		while (p_id0 != SIZE_MAX && e_id == pcl_in_elems[--p_id0]);
		++p_id0;
		// p_id1
		size_t p_id1 = Block_Low(task_id + 1, task_num, pcl_num);
		e_id = pcl_in_elems[p_id1];
		while (p_id1 != SIZE_MAX && e_id == pcl_in_elems[--p_id1]);
		++p_id1;			
		
		e_id = pcl_in_elems[p_id0];
		size_t &e_off = *node_elem_sort.task_off_sum_bin(task_id);
		e_off = e_id;
		size_t* const elem_ids = in_elem_ids + e_id;
		size_t* const node_ids = in_node_ids + e_id * node_num_per_elem;
		size_t* const node_elem_offs = in_node_elem_offs + e_id * node_num_per_elem;
		size_t elem_num = 0;
		for (size_t p_id = p_id0; p_id < p_id1; ++p_id)
		{
			if (e_id != pcl_in_elems[p_id + 1])
			{
				elem_ids[elem_num] = e_id;
				form_node_array<node_num_per_elem>(
					node_ids + elem_num * node_num_per_elem,
					node_elem_offs + elem_num * node_num_per_elem,
					e_id, elem_node_topologies);
				++elem_num;
				e_id = pcl_in_elems[p_id + 1];
			}
		}
		size_t& e_num = *((&e_off) + 1);
		e_num = elem_num;

		// sort my array
		size_t* comb_bin = node_elem_sort.task_combine_bin(task_id);
		size_t* sum_bin = (&e_num) + 1;
		ParaUtil::Internal::in_place_lsd_sort(
			node_ids, node_elem_offs,
			elem_num * node_num_per_elem,
			key_bit_offset, radix_bit_num,
			comb_bin, sum_bin);

		// form bins for merging
		for (size_t radix_id = radix_num; radix_id-- > 1;)
			comb_bin[radix_id] = sum_bin[radix_id] - sum_bin[radix_id - 1];

		// form result
		result.elem_num = elem_num;
		result.radix_num = radix_num;
		result.count_bin = comb_bin;

		return nullptr;
	}

	template <size_t node_num_per_elem>
	tbb::task* MergeToSortNE<node_num_per_elem>::operator() (
		tbb::task& parent,
		size_t radix_id) const
	{
		const size_t o_ne_offset = comb_bin[radix_id];
		size_t* const okeys = out_node_ids + o_ne_offset;
		size_t* const ovalues = out_node_elem_offs + o_ne_offset;
		size_t one_id = 0;
		if (radix_id == 0)
		{
			const size_t* const ielems = in_elem_ids;
			size_t oe_id = 0;
			for (size_t task_id = 0; task_id < task_num; ++task_id)
			{
				const size_t &i_e_offset = *node_elem_sort.task_off_sum_bin(task_id);
				const size_t& i_e_num = *((&i_e_offset) + 1);
				memcpy(out_elem_ids + oe_id, in_elem_ids + i_e_offset, i_e_num * sizeof(size_t));
				oe_id += i_e_num;
				const size_t* const sum_bin = (&i_e_num) + 1;
				const size_t i_ne_offset = i_e_offset * node_num_per_elem;
				const size_t* const ikeys = in_node_ids + i_ne_offset;
				const size_t* const ivalues = in_node_elem_offs + i_ne_offset;
				const size_t ne_num = sum_bin[0];
				memcpy(okeys + one_id, ikeys, ne_num * sizeof(size_t));
				memcpy(ovalues + one_id, ivalues, ne_num * sizeof(size_t));
				one_id += ne_num;
			}
		}
		else
		{
			for (size_t task_id = 0; task_id < task_num; ++task_id)
			{
				const size_t& i_e_offset = *node_elem_sort.task_off_sum_bin(task_id);
				const size_t* const sum_bin = (&i_e_offset) + 2;
				const size_t i_ne_offset = i_e_offset * node_num_per_elem + sum_bin[radix_id - 1];
				const size_t* const ikeys = in_node_ids + i_ne_offset;
				const size_t* const ivalues = in_node_elem_offs + i_ne_offset;
				const size_t ne_num = sum_bin[radix_id] - sum_bin[radix_id - 1];
				memcpy(okeys + one_id, ikeys, ne_num * sizeof(size_t));
				memcpy(ovalues + one_id, ivalues, ne_num * sizeof(size_t));
				one_id += ne_num;
			}
		}

		// sort in parallel
		if (one_id && key_bit_offset)
		{
			tbb::empty_task& contin = *new (parent.allocate_continuation()) tbb::empty_task;
			contin.set_ref_count(1);
			return new (contin.allocate_child())
				MSDRadixSort16Task(okeys, ovalues, one_id, key_bit_offset,
					tmp_node_ids + o_ne_offset,
					tmp_node_elem_offs + o_ne_offset);
		}
		return nullptr;
	}
}

	template <size_t node_num_per_elem>
	NodeElemSort<node_num_per_elem>::NodeElemSort() :
		thread_num(0),
		max_task_num(0),
		radix_num(0),
		radix_bit_num(0),
		key_bit_offset(0),
		in_elem_ids(nullptr),
		out_elem_ids(nullptr),
		sum_bins(nullptr),
		comb_bins(nullptr),
		scan_elem(*this),
		merge_to_sort(*this) {}

	template <size_t node_num_per_elem>
	void NodeElemSort<node_num_per_elem>::init(
		size_t max_elem_num,
		size_t max_node_num,
		size_t max_pcl_num,
		size_t th_num,
		const size_t p_in_es[],
		const size_t e_topos[])
	{
		thread_num = th_num;
		max_task_num = ParaUtil::cal_task_num<
			Internal::ne_sort_min_pcl_num_per_task,
			Internal::ne_sort_task_num_per_thread>(th_num, max_pcl_num);

		radix_bit_num = ParaUtil::cal_digit_num<1>(th_num * Internal::ne_sort_radix_num_per_thread);
		const size_t data_bit_num = ParaUtil::cal_digit_num<1>(max_node_num);
		if (radix_bit_num < data_bit_num)
		{
			key_bit_offset = data_bit_num - radix_bit_num;
		}
		else
		{
			key_bit_offset = 0;
			radix_bit_num = data_bit_num;
		}
		radix_num = size_t(1) << radix_bit_num;

		const size_t max_node_elem_num = max_elem_num * node_num_per_elem;
		in_elem_ids = node_elem_pair_mem.alloc<size_t>(
			  max_elem_num * (2 + node_num_per_elem * 6) + 2);
		out_elem_ids = in_elem_ids + max_elem_num;
		in_pairs.node_ids = out_elem_ids + max_elem_num;
		in_pairs.node_elem_offs = in_pairs.node_ids + max_node_elem_num;
		out_pairs.node_ids = in_pairs.node_elem_offs + max_node_elem_num + 1;
		out_pairs.node_elem_offs = out_pairs.node_ids + max_node_elem_num + 1;
		tmp_pairs.node_ids = out_pairs.node_elem_offs + max_node_elem_num;
		tmp_pairs.node_elem_offs = tmp_pairs.node_ids + max_node_elem_num;
		out_pairs.node_ids[-1] = SIZE_MAX;
		out_pairs.node_ids[max_elem_num * node_num_per_elem] = SIZE_MAX;

		sum_bins = sum_bin_mem.alloc<size_t>(max_task_num * (radix_num + 2));
		comb_bins = comb_bin_mem.alloc<size_t>(max_task_num * radix_num);

		scan_elem.init(p_in_es, e_topos);
		merge_to_sort.init();
	}

	template <size_t node_num_per_elem>
	void NodeElemSort<node_num_per_elem>::free()
	{
		node_elem_pair_mem.free();
		sum_bin_mem.free();
		comb_bin_mem.free();
		in_elem_ids = nullptr;
		out_elem_ids = nullptr;
		in_pairs.node_ids = nullptr;
		in_pairs.node_elem_offs = nullptr;
		out_pairs.node_ids = nullptr;
		out_pairs.node_elem_offs = nullptr;
		tmp_pairs.node_ids = nullptr;
		tmp_pairs.node_elem_offs = nullptr;
		sum_bins = nullptr;
		comb_bins = nullptr;

		thread_num = 0;
		max_task_num = 0;
		radix_num = 0;
		radix_bit_num = 0;
		key_bit_offset = 0;
	}

	template <size_t node_num_per_elem>
	void NodeElemSort<node_num_per_elem>::sort(size_t pcl_num)
	{
		const size_t task_num = ParaUtil::cal_task_num<
			Internal::ne_sort_min_pcl_num_per_task,
			Internal::ne_sort_task_num_per_thread>(
				thread_num, pcl_num);

		scan_elem.update(pcl_num, task_num);
		ParaUtil::parallel_reduce<Internal::ScanElem<node_num_per_elem>, Internal::ScanElemRes>(
			scan_elem, scan_elem_res, task_num);

		const size_t ne_num = scan_elem_res.elem_num * node_num_per_elem;
		comb_bins[radix_num - 1] = ne_num - comb_bins[radix_num - 1];
		for (size_t radix_id = radix_num - 1; radix_id-- > 1;)
			comb_bins[radix_id] = comb_bins[radix_id + 1] - comb_bins[radix_id];
		comb_bins[0] = 0;

		merge_to_sort.update(task_num);
		out_pairs.node_ids[ne_num] = SIZE_MAX;
		ParaUtil::parallel_for<Internal::MergeToSortNE<node_num_per_elem>>(merge_to_sort, radix_num);
	}
}

#endif