#ifndef __Node_Elem_Sort_Internal_hpp__
#define __Node_Elem_Sort_Internal_hpp__

#define TBB_SUPPRESS_DEPRECATED_MESSAGES 1
#include <tbb/tbb.h>

#include "DataMem.h"

namespace ParaUtil
{
	template <size_t node_num_per_elem> class NodeElemSort;

	namespace Internal
	{
		constexpr size_t ne_sort_min_pcl_num_per_task = 10;
		constexpr size_t ne_sort_task_num_per_thread = 20;
		constexpr size_t ne_sort_radix_num_per_thread = 16;

		struct NodeElemPairs
		{
			size_t* node_ids; // key
			size_t* node_elem_offs; // value
			NodeElemPairs() : node_ids(nullptr), node_elem_offs(nullptr) {}
			inline void swap(size_t i, size_t j) const noexcept
			{
				size_t tmp;
				tmp = node_ids[i];
				node_ids[i] = node_ids[j];
				node_ids[j] = tmp;
				tmp = node_elem_offs[i];
				node_elem_offs[i] = node_elem_offs[j];
				node_elem_offs[j] = tmp;
				return;
			}
		};

		struct ScanElemRes
		{
			size_t elem_num;
			size_t radix_num;
			size_t* count_bin;
			void join(const ScanElemRes& other)
			{
				elem_num += other.elem_num;
				size_t* const o_cbin = other.count_bin;
				for (size_t radix_id = 0; radix_id < radix_num; ++radix_id)
					count_bin[radix_id] += o_cbin[radix_id];
			}
		};

		template <size_t node_num_per_elem>
		class ScanElem
		{
		public:
			ScanElem(NodeElemSort<node_num_per_elem>& _nes) : node_elem_sort(_nes) {}
			void init(const size_t p_in_es[], const size_t e_topos[]);
			void update(size_t p_num, size_t tsk_num);
			tbb::task* operator()(tbb::task& parent, size_t task_id, ScanElemRes& result) const;

		protected:
			NodeElemSort<node_num_per_elem> &node_elem_sort;
			const size_t *pcl_in_elems;
			const size_t *elem_node_topologies;
			size_t* in_elem_ids;
			size_t* in_node_ids, * in_node_elem_offs;
			size_t key_bit_offset, radix_bit_num, radix_num;
			size_t pcl_num, task_num;
		};

		template <size_t node_num_per_elem>
		class MergeToSortNE
		{
		public:
			MergeToSortNE(NodeElemSort<node_num_per_elem> & _nes) : node_elem_sort(_nes) {}
			void init();
			void update(size_t e_num);
			tbb::task* operator() (tbb::task& parent, size_t radix_id) const;

		protected:
			NodeElemSort<node_num_per_elem> &node_elem_sort;
			const size_t *comb_bin;
			const size_t *in_elem_ids;
			size_t* out_elem_ids;
			const size_t *in_node_ids, *in_node_elem_offs;
			size_t *out_node_ids, *out_node_elem_offs;
			size_t *tmp_node_ids, *tmp_node_elem_offs;
			size_t key_bit_offset, radix_bit_num, radix_num;
			size_t task_num;
		};
	}

	template <size_t node_num_per_elem>
	class NodeElemSort
	{
		template <size_t node_num_per_elem>
		friend class Internal::ScanElem;
		template <size_t node_num_per_elem>
		friend class Internal::MergeToSortNE;

	public:
		NodeElemSort();
		~NodeElemSort() { free(); }

		void init(size_t max_elem_num, size_t max_node_num,
				  size_t max_pcl_num, size_t th_num,
				  const size_t p_in_es[], const size_t e_topos[]);
		void free();

		inline const size_t* node_ids() const noexcept { return out_pairs.node_ids; }
		inline const size_t* node_elem_offs() const noexcept { return out_pairs.node_elem_offs; }
		inline const size_t* elem_ids() const noexcept { return out_elem_ids; }
		inline size_t elem_num() const noexcept { return scan_elem_res.elem_num; }

		void sort(size_t pcl_num);

	protected:
		typedef Internal::NodeElemPairs NodeElemPairs;

		size_t thread_num, max_task_num;
		size_t radix_num, radix_bit_num, key_bit_offset;

		Util::DataMem node_elem_pair_mem;
		size_t *in_elem_ids, *out_elem_ids;
		NodeElemPairs in_pairs, out_pairs, tmp_pairs;

		Util::DataMem sum_bin_mem, comb_bin_mem;
		size_t* sum_bins, * comb_bins; // max_task_num * radix_num

		Internal::ScanElemRes scan_elem_res;
		Internal::ScanElem<node_num_per_elem> scan_elem;
		Internal::MergeToSortNE<node_num_per_elem> merge_to_sort;

		inline size_t *task_off_sum_bin(size_t tsk_id) noexcept
		{ return sum_bins + (tsk_id << radix_bit_num) + tsk_id * 2; }
		inline size_t* task_combine_bin(size_t tsk_id) noexcept
		{ return comb_bins + (tsk_id << radix_bit_num); }
	};

	namespace Internal
	{
		template <size_t node_num_per_elem>
		inline void ScanElem<node_num_per_elem>::init(
			const size_t p_in_es[],
			const size_t e_topos[])
		{
			pcl_in_elems = p_in_es;
			elem_node_topologies = e_topos;
			in_elem_ids = node_elem_sort.in_elem_ids;
			in_node_ids = node_elem_sort.in_pairs.node_ids;
			in_node_elem_offs = node_elem_sort.in_pairs.node_elem_offs;
			radix_num = node_elem_sort.radix_num;
			radix_bit_num = node_elem_sort.radix_bit_num;
			key_bit_offset = node_elem_sort.key_bit_offset;
		}

		template <size_t node_num_per_elem>
		inline void ScanElem<node_num_per_elem>::update(
			size_t p_num, size_t tsk_num)
		{
			pcl_num = p_num;
			task_num = tsk_num;
		}

		template <size_t node_num_per_elem>
		inline void MergeToSortNE<node_num_per_elem>::init()
		{
			comb_bin = node_elem_sort.comb_bins;
			in_elem_ids = node_elem_sort.in_elem_ids;
			out_elem_ids = node_elem_sort.out_elem_ids;
			in_node_ids = node_elem_sort.in_pairs.node_ids;
			in_node_elem_offs = node_elem_sort.in_pairs.node_elem_offs;
			out_node_ids = node_elem_sort.out_pairs.node_ids;
			out_node_elem_offs = node_elem_sort.out_pairs.node_elem_offs;
			tmp_node_ids = node_elem_sort.tmp_pairs.node_ids;
			tmp_node_elem_offs = node_elem_sort.tmp_pairs.node_elem_offs;
			radix_num = node_elem_sort.radix_num;
			radix_bit_num = node_elem_sort.radix_bit_num;
			key_bit_offset = node_elem_sort.key_bit_offset;
		}

		template <size_t node_num_per_elem>
		inline void MergeToSortNE<node_num_per_elem>::update(size_t tsk_num)
		{
			task_num = tsk_num;
		}
	}
}

#endif