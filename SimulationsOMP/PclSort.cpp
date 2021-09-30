#include "SimulationsOMP_pcp.h"

#include "ParaUtil.h"
#include "ParallelReduceTask.hpp"
#include "ParallelForTask.hpp"
#include "ParaUtilSerialSort.h"
#include "ParaUtilMSDRadix16Sort.h"

#include "PclSort.h"

#define Block_Low(blk_id, blk_num, data_num) ((blk_id) * (data_num) / (blk_num))
#define Data_Digit(num, disp, mask) (((num) >> (disp)) & (mask))

namespace ParaUtil
{
namespace Internal
{
	tbb::task* ScanPcl::operator()(
		tbb::task& parent,
		size_t task_id,
		ScanPclRes &result) const
	{
		// sort my array
		size_t* comb_bin = pcl_sort.task_combine_bin(task_id);
		size_t* sum_bin = pcl_sort.task_sum_bin(task_id);
		const size_t data_id0 = Block_Low(task_id, task_num, pcl_num);
		ParaUtil::Internal::in_place_lsd_sort(
			in_pcl_in_elems + data_id0,
			in_prev_pcl_ids + data_id0,
			Block_Low(task_id + 1, task_num, pcl_num) - data_id0,
			key_bit_offset, radix_bit_num,
			comb_bin, sum_bin);
		
		// form bins for merging
		for (size_t radix_id = radix_num; radix_id-- > 1;)
			comb_bin[radix_id] = sum_bin[radix_id] - sum_bin[radix_id - 1];

		// form result
		result.radix_num = radix_num;
		result.count_bin = comb_bin;

		return nullptr;
	}
	
	tbb::task* MergeToSortPcl::operator() (tbb::task& parent, size_t radix_id) const
	{
		const size_t o_pcl_offset = comb_bin[radix_id];
		size_t *const okeys = out_pcl_in_elems + o_pcl_offset;
		size_t* const ovalues = out_prev_pcl_ids + o_pcl_offset;
		size_t opcl_id = 0;
		if (radix_id == 0)
		{
			for (size_t task_id = 0; task_id < task_num; ++task_id)
			{
				const size_t* const sum_bin = pcl_sort.task_sum_bin(task_id);
				const size_t p_num = sum_bin[0];
				const size_t ipcl_id = Block_Low(task_id, task_num, pcl_num);
				const size_t* const ikeys = in_pcl_in_elems + ipcl_id;
				const size_t* const ivalues = in_prev_pcl_ids + ipcl_id;
				memcpy(okeys + opcl_id, ikeys, p_num * sizeof(size_t));
				memcpy(ovalues + opcl_id, ivalues, p_num * sizeof(size_t));
				opcl_id += p_num;
			}
		}
		else
		{
			for (size_t task_id = 0; task_id < task_num; ++task_id)
			{
				const size_t* const sum_bin = pcl_sort.task_sum_bin(task_id);
				const size_t p_num = sum_bin[radix_id] - sum_bin[radix_id - 1];
				const size_t ipcl_id = Block_Low(task_id, task_num, pcl_num) + sum_bin[radix_id - 1];
				const size_t* const ikeys = in_pcl_in_elems + ipcl_id;
				const size_t* const ivalues = in_prev_pcl_ids + ipcl_id;
				memcpy(okeys + opcl_id, ikeys, p_num * sizeof(size_t));
				memcpy(ovalues + opcl_id, ivalues, p_num * sizeof(size_t));
				opcl_id += p_num;
			}
		}

		// sort in parallel
		if (opcl_id && key_bit_offset)
		{
			tbb::empty_task& contin = *new (parent.allocate_continuation()) tbb::empty_task;
			contin.set_ref_count(1);
			return new (contin.allocate_child())
				MSDRadixSort16Task(okeys, ovalues, opcl_id, key_bit_offset,
					tmp_pcl_in_elems + o_pcl_offset,
					tmp_prev_pcl_ids + o_pcl_offset);
		}
		return nullptr;
	}
}

	PclSort::PclSort() : 
		sum_bins(nullptr),
		comb_bins(nullptr),
		thread_num(0),
		max_task_num(0),
		radix_num(0),
		radix_bit_num(0),
		key_bit_offset(0),
		scan_pcl(*this),
		merge_to_sort(*this) {}

	void PclSort::init(size_t max_pcl_num, size_t th_num)
	{
		thread_num = th_num;
		max_task_num = ParaUtil::cal_task_num<
			Internal::pcl_sort_min_pcl_num_per_task,
			Internal::pcl_sort_task_num_per_thread>(
				th_num, max_pcl_num);

		radix_bit_num = ParaUtil::cal_digit_num<1>(th_num * Internal::pcl_sort_radix_num_per_thread);
		const size_t data_bit_num = ParaUtil::cal_digit_num<1>(max_pcl_num);
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
		
		in_pairs.pcl_in_elems = pcl_pair_mem.alloc<size_t>(max_pcl_num * 6 + 2);
		in_pairs.prev_pcl_ids = in_pairs.pcl_in_elems + max_pcl_num;
		out_pairs.pcl_in_elems = in_pairs.prev_pcl_ids + max_pcl_num + 1;
		out_pairs.prev_pcl_ids = out_pairs.pcl_in_elems + max_pcl_num + 1;
		tmp_pairs.pcl_in_elems = out_pairs.prev_pcl_ids + max_pcl_num;
		tmp_pairs.prev_pcl_ids = tmp_pairs.pcl_in_elems + max_pcl_num;
		out_pairs.pcl_in_elems[-1] = SIZE_MAX;
		out_pairs.pcl_in_elems[max_pcl_num] = SIZE_MAX;

		sum_bins = sum_bin_mem.alloc<size_t>(max_task_num * radix_num);
		comb_bins = comb_bin_mem.alloc<size_t>(max_task_num * radix_num);

		scan_pcl.init();
		merge_to_sort.init();
	}

	void PclSort::free()
	{
		pcl_pair_mem.free();
		sum_bin_mem.free();
		comb_bin_mem.free();
		in_pairs.pcl_in_elems = nullptr;
		in_pairs.prev_pcl_ids = nullptr;
		out_pairs.pcl_in_elems = nullptr;
		out_pairs.prev_pcl_ids = nullptr;
		tmp_pairs.pcl_in_elems = nullptr;
		tmp_pairs.prev_pcl_ids = nullptr;
		sum_bins = nullptr;
		comb_bins = nullptr;
		
		thread_num = 0;
		max_task_num = 0;
		radix_num = 0;
		radix_bit_num = 0;
		key_bit_offset = 0;
		pcl_num = 0;
		task_num = 0;
	}

	void PclSort::sort(size_t _pcl_num)
	{
		pcl_num = _pcl_num;
		task_num = ParaUtil::cal_task_num<
			Internal::pcl_sort_min_pcl_num_per_task,
			Internal::pcl_sort_task_num_per_thread>(
				thread_num, pcl_num);

		scan_pcl.update();
		ParaUtil::parallel_reduce<Internal::ScanPcl, Internal::ScanPclRes>(
			scan_pcl, scan_pcl_res, task_num);

		comb_bins[radix_num - 1] = pcl_num - comb_bins[radix_num - 1];
		for (size_t radix_id = radix_num - 1; radix_id-- > 1;)
			comb_bins[radix_id] = comb_bins[radix_id + 1] - comb_bins[radix_id];
		comb_bins[0] = 0;

		merge_to_sort.update();
		ParaUtil::parallel_for<Internal::MergeToSortPcl>(merge_to_sort, radix_num);
		out_pairs.pcl_in_elems[pcl_num] = SIZE_MAX;
	}
}