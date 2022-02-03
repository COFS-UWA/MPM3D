#ifndef __Pcl_Sort_h__
#define __Pcl_Sort_h__

#define TBB_SUPPRESS_DEPRECATED_MESSAGES 1
#include <tbb/tbb.h>

#include "DataMem.h"

namespace ParaUtil
{
	class PclSort;

namespace Internal
{
	constexpr size_t pcl_sort_min_pcl_num_per_task = 10;
	constexpr size_t pcl_sort_task_num_per_thread = 20;
	constexpr size_t pcl_sort_radix_num_per_thread = 16;

	struct PclPairs
	{
		size_t *pcl_in_elems; // key
		size_t *prev_pcl_ids; // value
		PclPairs() : pcl_in_elems(nullptr), prev_pcl_ids(nullptr) {}
		inline void swap(size_t i, size_t j) const noexcept
		{
			size_t tmp;
			tmp = pcl_in_elems[i];
			pcl_in_elems[i] = pcl_in_elems[j];
			pcl_in_elems[j] = tmp;
			tmp = prev_pcl_ids[i];
			prev_pcl_ids[i] = prev_pcl_ids[j];
			prev_pcl_ids[j] = tmp;
			return;
		}
	};

	struct ScanPclRes
	{
		size_t radix_num;
		union { size_t* count_bin, * sum_bin; };
		void join(const ScanPclRes& other)
		{
			size_t* const o_cbin = other.count_bin;
			for (size_t radix_id = 0; radix_id < radix_num; ++radix_id)
				count_bin[radix_id] += o_cbin[radix_id];
		}
	};

	class ScanPcl
	{
	public:
		ScanPcl(PclSort& _ps) : pcl_sort(_ps) {}
		void init();
		void update();
		tbb::task* operator()(tbb::task& parent, size_t task_id, ScanPclRes& result) const;

	protected:
		PclSort& pcl_sort;
		size_t* in_pcl_in_elems, * in_prev_pcl_ids;
		size_t key_bit_offset, radix_bit_num, radix_num;
		size_t pcl_num, task_num;
	};
	
	class MergeToSortPcl
	{
	public:
		MergeToSortPcl(PclSort& _ps) : pcl_sort(_ps) {}
		void init();
		void update();
		tbb::task* operator() (tbb::task& parent, size_t radix_id) const;

	protected:
		PclSort& pcl_sort;
		const size_t* comb_bin;
		const size_t* in_pcl_in_elems, * in_prev_pcl_ids;
		size_t* out_pcl_in_elems, * out_prev_pcl_ids;
		size_t* tmp_pcl_in_elems, * tmp_prev_pcl_ids;
		size_t key_bit_offset, radix_bit_num, radix_num;
		size_t pcl_num, task_num;
	};
}

class PclSort
{
	friend class Internal::ScanPcl;
	friend class Internal::MergeToSortPcl;

public:
	PclSort();
	~PclSort() { free(); }

	void init(size_t max_pcl_num, size_t thread_num);
	void free();

	inline size_t* in_pcl_in_elems() noexcept { return in_pairs.pcl_in_elems; }
	inline size_t* in_prev_pcl_ids() noexcept { return in_pairs.prev_pcl_ids; }
	inline const size_t* out_pcl_in_elems() const noexcept { return out_pairs.pcl_in_elems; }
	inline const size_t* out_prev_pcl_ids() const noexcept { return out_pairs.prev_pcl_ids; }

	void sort(size_t pcl_num);

protected: // vars updated during cal
	size_t pcl_num, task_num;

protected:
	typedef Internal::PclPairs PclPairs;
	
	size_t thread_num, max_task_num;
	size_t radix_num, radix_bit_num, key_bit_offset;

	Util::DataMem pcl_pair_mem;
	PclPairs in_pairs, out_pairs, tmp_pairs;

	Util::DataMem sum_bin_mem, comb_bin_mem;
	size_t *sum_bins, *comb_bins; // max_task_num * radix_num

	Internal::ScanPclRes scan_pcl_res;
	Internal::ScanPcl scan_pcl;
	Internal::MergeToSortPcl merge_to_sort;
	
	inline size_t* task_sum_bin(size_t tsk_id) noexcept
	{ return sum_bins + (tsk_id << radix_bit_num); }
	inline size_t* task_combine_bin(size_t tsk_id) noexcept
	{ return comb_bins + (tsk_id << radix_bit_num); }
};

namespace Internal
{
	inline void ScanPcl::init()
	{
		in_pcl_in_elems = pcl_sort.in_pcl_in_elems();
		in_prev_pcl_ids = pcl_sort.in_prev_pcl_ids();
		radix_num = pcl_sort.radix_num;
		radix_bit_num = pcl_sort.radix_bit_num;
		key_bit_offset = pcl_sort.key_bit_offset;
	}

	inline void ScanPcl::update()
	{
		pcl_num = pcl_sort.pcl_num;
		task_num = pcl_sort.task_num;
	}

	inline void MergeToSortPcl::init()
	{
		comb_bin = pcl_sort.comb_bins;
		in_pcl_in_elems = pcl_sort.in_pcl_in_elems();
		in_prev_pcl_ids = pcl_sort.in_prev_pcl_ids();
		out_pcl_in_elems = const_cast<size_t*>(pcl_sort.out_pcl_in_elems());
		out_prev_pcl_ids = const_cast<size_t*>(pcl_sort.out_prev_pcl_ids());
		tmp_pcl_in_elems = pcl_sort.tmp_pairs.pcl_in_elems;
		tmp_prev_pcl_ids = pcl_sort.tmp_pairs.prev_pcl_ids;
		radix_num = pcl_sort.radix_num;
		radix_bit_num = pcl_sort.radix_bit_num;
		key_bit_offset = pcl_sort.key_bit_offset;
	}

	inline void MergeToSortPcl::update()
	{
		pcl_num = pcl_sort.pcl_num;
		task_num = pcl_sort.task_num;
	}
}
}

#endif