#ifndef __Sort_Tri_Mesh_Node_Task_h__
#define __Sort_Tri_Mesh_Node_Task_h__

#include <assert.h>
#include "SortTask.h"

namespace SortUtils
{
	struct SortTriMeshNodeMem
	{
		static constexpr size_t max_block_num_div_thread_num = 2;
		static constexpr size_t parallel_divide_min_pcl_num_per_block = 2;
		static constexpr size_t serial_sort_max_pcl_num = 1;

		union
		{
			struct
			{
				size_t max_block_num;
				size_t**node_has_elem_arrays; // (elem_num * 3) * node_digit_num + 1
				size_t**node_elem_pair_arrays; // (elem_num * 3) node_digit_num + 1
				size_t** tmp_node_has_elem_arrays; // (elem_num * 3) * node_digit_num
				size_t** tmp_node_elem_pair_arrays; // (elem_num * 3) * node_digit_num
				size_t* ori_node_has_elem; // key_arrays[node_digit_num]
				size_t* ori_node_elem_pair; // val_arrays[node_digit_num]
				size_t* node_has_elem; // key_arrays[0]
				size_t* node_elem_pair; // val_arrays[0]
				SortBin* root_bin;
				SortBin** thread_bins; // thread_num * block_num
			};
			SortMem _sort_mem;
		};

		size_t* tmp_valid_elems;
		size_t *valid_elems;

		SortTriMeshNodeMem() {}
		~SortTriMeshNodeMem() { clear(); }
		// shared_mem used for :
		//  key_arrays
		//  val_arrays
		//  tmp_key_arrays
		//  tmp_val_arrays
		//  ori_pcl_in_elem
		//  ori_cur_to_prev_pcl
		//  root_bin
		//  thread_bins
		// tmp_valid_elems
		// self_mem used for:
		//  key_arrays *[]
		//  val_arrays *[]
		//  tmp_key_arrays *[]
		//  tmp_val_arrays *[]
		//  pcl_in_elem
		//  cur_to_prev_pcl
		//  valid_elems
		inline size_t get_shared_memory_size(
			size_t thread_num, size_t elem_num,
			size_t node_digit_num) noexcept
		{
			max_block_num = thread_num * max_block_num_div_thread_num;
			const size_t three_elem_num = elem_num * 3;
			size_t block_num = (three_elem_num + parallel_divide_min_pcl_num_per_block - 1)
							  / parallel_divide_min_pcl_num_per_block;
			if (block_num > max_block_num)
				block_num = max_block_num;
			return Cache_Alignment * (node_digit_num * 2 + 1)
				+ sizeof(size_t) * three_elem_num * node_digit_num * 2
				+ sizeof(size_t) * elem_num
				+ sizeof(SortBin) * (1 + thread_num) * block_num;
		}
		void init(void* shared_mem, size_t thread_num,
				  size_t elem_num, size_t node_digit_num);
	protected:
		CacheAlignedMem self_mem;
		inline void clear() noexcept { self_mem.free(); }
	};

	namespace Internal
	{
		class CopyMemTask : public tbb::task
		{
		protected:
			void* const dst;
			const void* const src;
			const size_t len;
		public:
			CopyMemTask(
				void *const _dst,
				const void *const _src,
				const size_t _len) : 
				dst(_dst), src(_src), len(_len) {}
			~CopyMemTask() {}
			tbb::task* execute() override
			{ memcpy(dst, src, len); }
		};
		
		struct ElemNodeIndex { size_t n1, n2, n3; };

#define __Cal_Key_Digit__(data, digit_pos) (((data) >> ((digit_pos) * 8)) & (0xFF))	
		inline size_t count_sort_tri_mesh_node(
#ifdef _DEBUG
			const size_t elem_num,
#endif
			const ElemNodeIndex* const elem_node_id,
			const size_t* const pcl_in_elem,
			const size_t pcl_num,
			const size_t node_digit_pos,
			SortBin &bin,
			size_t* const out_valid_elems,
			size_t* const tmp_node_has_elem,
			size_t* const tmp_node_elem_pair,
			size_t* const out_node_has_elem,
			size_t* const out_node_elem_pair)
		{
			size_t* const c_bin = bin.count_bin;
			size_t* const s_bin = bin.sum_bin;
			memset(c_bin, 0, sizeof(size_t) * radix_bucket_num);
			size_t e_id = pcl_in_elem[0];
#ifdef _DEBUG
			assert(e_id < elem_num || e_id == SIZE_MAX);
#endif
			size_t valid_elem_num = 0;
			for (size_t p_id = 0; p_id < pcl_num; ++p_id)
			{
				if (e_id != pcl_in_elem[p_id + 1])
				{
					const ElemNodeIndex& eni = elem_node_id[e_id];
					++c_bin[__Cal_Key_Digit__(eni.n1, node_digit_pos)];
					++c_bin[__Cal_Key_Digit__(eni.n2, node_digit_pos)];
					++c_bin[__Cal_Key_Digit__(eni.n3, node_digit_pos)];
					tmp_node_has_elem[valid_elem_num * 3] = eni.n1;
					tmp_node_has_elem[valid_elem_num * 3 + 1] = eni.n2;
					tmp_node_has_elem[valid_elem_num * 3 + 2] = eni.n3;
					tmp_node_elem_pair[valid_elem_num * 3] = e_id * 3;
					tmp_node_elem_pair[valid_elem_num * 3 + 1] = e_id * 3 + 1;
					tmp_node_elem_pair[valid_elem_num * 3 + 2] = e_id * 3 + 2;
					out_valid_elems[valid_elem_num++] = e_id;
					e_id = pcl_in_elem[p_id + 1];
#ifdef _DEBUG
					assert(e_id < elem_num || e_id == SIZE_MAX);
#endif
				}
			}
			// form count and sort bin
			s_bin[0] = c_bin[0];
			for (size_t i = 1; i < radix_bucket_num; ++i)
				s_bin[i] = s_bin[i - 1] + c_bin[i];
			// reorder memory
			size_t pos;
			for (e_id = valid_elem_num; e_id--;)
			{
				pos = --s_bin[__Cal_Key_Digit__(tmp_node_has_elem[e_id * 3], node_digit_pos)];
				out_node_has_elem[pos] = tmp_node_has_elem[e_id * 3];
				out_node_elem_pair[pos] = tmp_node_elem_pair[e_id * 3];
				pos = --s_bin[__Cal_Key_Digit__(tmp_node_has_elem[e_id * 3 + 1], node_digit_pos)];
				out_node_has_elem[pos] = tmp_node_has_elem[e_id * 3 + 1];
				out_node_elem_pair[pos] = tmp_node_elem_pair[e_id * 3 + 1];
				pos = --s_bin[__Cal_Key_Digit__(tmp_node_has_elem[e_id * 3 + 2], node_digit_pos)];
				out_node_has_elem[pos] = tmp_node_has_elem[e_id * 3 + 2];
				out_node_elem_pair[pos] = tmp_node_elem_pair[e_id * 3 + 2];
			}
			return valid_elem_num;
		}
#undef __Cal_Key_Digit__

		class CountSortTriMeshNodeTask : public tbb::task
		{
		protected:
#ifdef _DEBUG
			size_t pcl_num;
			size_t elem_num;
#endif
			const ElemNodeIndex* const elem_node_id;
			const size_t* const pcl_in_elem;
			size_t p_id0, p_id1;
			const size_t node_digit_pos;
			SortBin &bin;
			size_t* const tmp_valid_elems;
			size_t* const tmp_node_has_elem;
			size_t* const tmp_node_elem_pair;
			size_t* const out_node_has_elem;
			size_t* const out_node_elem_pair;
		public:
			CountSortTriMeshNodeTask(
#ifdef _DEBUG
				const size_t p_num,
				const size_t e_num,
#endif
				const ElemNodeIndex* en_ids,
				const size_t* p_in_e,
				size_t _p_id0,
				size_t _p_id1,
				const size_t _n_dgt_pos,
				SortBin& _bin,
				size_t* tmp_ves,
				size_t* tmp_n_has_e,
				size_t* tmp_n_e_pair,
				size_t* out_n_has_e,
				size_t* out_n_e_pair) :
#ifdef _DEBUG
				pcl_num(p_num),
				elem_num(e_num),
#endif
				elem_node_id(en_ids),
				pcl_in_elem(p_in_e),
				p_id0(_p_id0),
				p_id1(_p_id1),
				node_digit_pos(_n_dgt_pos),
				bin(_bin),
				tmp_valid_elems(tmp_ves),
				tmp_node_has_elem(tmp_n_has_e),
				tmp_node_elem_pair(tmp_n_e_pair),
				out_node_has_elem(out_n_has_e),
				out_node_elem_pair(out_n_e_pair) {}
			~CountSortTriMeshNodeTask() {}
			tbb::task* execute() override;
		};
	}

	// Assumption:
	// pcl_num > 0
	// pcl_in_elem was sorted in accending order
	class SortTriMeshNodeTask : public tbb::task
	{
	protected:
		static constexpr size_t parallel_divide_min_pcl_num_per_block = 2;
		static constexpr size_t serial_sort_min_pcl_num = 1;
#ifdef _DEBUG
		const size_t elem_num;
#endif
		const size_t *const pcl_in_elem;
		const size_t pcl_num;
		const Internal::ElemNodeIndex *const elem_node_ids;
		const size_t node_digit_num;
		SortTriMeshNodeMem& sort_mem;
		size_t& valid_elem_num;
	public:
		SortTriMeshNodeTask(
#ifdef _DEBUG
			const size_t e_num,
#endif
			const size_t *p_in_e,
			const size_t p_num,
			const void *en_ids,
			const size_t n_dgt_num,
			SortTriMeshNodeMem& sm,
			size_t &ve_num) :
#ifdef _DEBUG
			elem_num(e_num),
#endif
			pcl_in_elem(p_in_e),
			pcl_num(p_num),
			elem_node_ids((const Internal::ElemNodeIndex *)en_ids),
			node_digit_num(n_dgt_num),
			sort_mem(sm),
			valid_elem_num(ve_num) {}
		~SortTriMeshNodeTask() {}
		tbb::task* execute() override;
	};
}

#endif