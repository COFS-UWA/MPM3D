#ifndef __Sort_Tri_Mesh_Node_Task_hpp__
#define __Sort_Tri_Mesh_Node_Task_hpp__

#include <assert.h>
#include "SortTask.h"

namespace SortUtils
{
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
		
#define __Cal_Key_Digit__(data, digit_pos) (((data) >> ((digit_pos) * 8)) & (0xFF))
		
		template <typename ElemNodeIndex>
		inline void count_sort_tri_mesh_node(
#ifdef _DEBUG
			const size_t elem_num,
#endif
			const ElemNodeIndex* const elem_node_id,
			const size_t* const pcl_in_elem,
			const size_t pcl_num,
			const size_t node_digit_pos,
			SortBin& bin,
			size_t* const out_valid_elems,
			size_t* const tmp_node_has_elem,
			size_t* const tmp_node_elem_pair,
			size_t* const out_node_has_elem,
			size_t* const out_node_elem_pair)
		{
			size_t* const c_bin = bin.count_bin;
			size_t* const s_bin = bin.sum_bin;
			memset(c_bin, 0, sizeof(size_t) * radix_bucket_num);
			size_t e_id = pcl_in_elem[p_id0];
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
		}
		
		template <typename ElemNodeIndex>
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
				const ElemNodeIndex* const en_ids,
				const size_t* const p_in_e,
				size_t _p_id0,
				size_t _p_id1,
				const size_t _n_dgt_pos,
				SortBin& _bin,
				size_t* const tmp_ves,
				size_t* const tmp_n_has_e,
				size_t* const tmp_n_e_pair,
				size_t* const out_n_has_e,
				size_t* const out_n_e_pair) :
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
			tbb::task *execute() override
			{
				size_t e_id;
				e_id = pcl_in_elem[p_id0];
				while (p_id0 != SIZE_MAX && e_id == pcl_in_elem[--p_id0]);
				++p_id0;
				e_id = pcl_in_elem[p_id1];
				while (p_id1 != SIZE_MAX && e_id == pcl_in_elem[--p_id1]);
				++p_id1;
#ifdef _DEBUG
				assert(p_id0 <= pcl_num);
				assert(p_id1 <= pcl_num);
#endif
				size_t* const c_bin = bin.count_bin;
				size_t* const s_bin = bin.sum_bin;
				memset(c_bin, 0, sizeof(size_t) * radix_bucket_num);
				size_t e_id = pcl_in_elem[p_id0];
#ifdef _DEBUG
				assert(e_id < elem_num || e_id == SIZE_MAX);
#endif
				size_t* const my_valid_elems = valid_elems + e_id;
				size_t* const my_node_has_elem = tmp_node_has_elem + e_id * 3;
				size_t* const my_node_elem_pair = tmp_node_elem_pair + e_id * 3;
				s_bin[0] = e_id * 3;
				size_t valid_elem_num = 0;
				for (size_t p_id = p_id0; p_id < p_id1; ++p_id)
				{
					if (e_id != pcl_in_elem[p_id + 1])
					{
						const ElemNodeIndex& eni = elem_node_id[e_id];
						++c_bin[__Cal_Key_Digit__(eni.n1, node_digit_pos)];
						++c_bin[__Cal_Key_Digit__(eni.n2, node_digit_pos)];
						++c_bin[__Cal_Key_Digit__(eni.n3, node_digit_pos)];
						my_node_has_elem[valid_elem_num * 3] = eni.n1;
						my_node_has_elem[valid_elem_num * 3 + 1] = eni.n2;
						my_node_has_elem[valid_elem_num * 3 + 2] = eni.n3;
						my_node_elem_pair[valid_elem_num * 3] = e_id * 3;
						my_node_elem_pair[valid_elem_num * 3 + 1] = e_id * 3 + 1;
						my_node_elem_pair[valid_elem_num * 3 + 2] = e_id * 3 + 2;
						my_valid_elems[valid_elem_num++] = e_id;
						e_id = pcl_in_elem[p_id + 1];
#ifdef _DEBUG
						assert(e_id < elem_num || e_id == SIZE_MAX);
#endif
					}
				}
				// form count and sort bin
				s_bin[0] += c_bin[0];
				for (size_t i = 1; i < radix_bucket_num; ++i)
					s_bin[i] = s_bin[i - 1] + c_bin[i];
				// reorder memory
				size_t pos;
				for (e_id = valid_elem_num; e_id--;)
				{
					pos = --s_bin[__Cal_Key_Digit__(my_node_has_elem[e_id * 3], node_digit_pos)];
					out_node_has_elem[pos] = my_node_has_elem[e_id * 3];
					out_node_elem_pair[pos] = my_node_elem_pair[e_id * 3];
					pos = --s_bin[__Cal_Key_Digit__(my_node_has_elem[e_id * 3 + 1], node_digit_pos)];
					out_node_has_elem[pos] = my_node_has_elem[e_id * 3 + 1];
					out_node_elem_pair[pos] = my_node_elem_pair[e_id * 3 + 1];
					pos = --s_bin[__Cal_Key_Digit__(my_node_has_elem[e_id * 3 + 2], node_digit_pos)];
					out_node_has_elem[pos] = my_node_has_elem[e_id * 3 + 2];
					out_node_elem_pair[pos] = my_node_elem_pair[e_id * 3 + 2];
				}
				return nullptr;
			}
		};
	}

#undef __Cal_Key_Digit__

	// Assumption:
	// pcl_num > 0
	// pcl_in_elem was sorted in accending order
	template <typename ElemNodeIndex>
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
		const size_t node_digit_num; // for max_digit_pos
		SortMem& sort_mem;
		size_t *const tmp_valid_elems;
		size_t *const out_valid_elems;
	public:
		SortTriMeshNodeTask(
#ifdef _DEBUG
			const size_t e_num,
#endif
			const size_t *const p_in_e,
			const size_t p_num,
			const size_t n_dgt_num,
			SortMem& sm,
			size_t *const tmp_ves,
			size_t *const out_ves) :
#ifdef _DEBUG
			elem_num(e_num),
#endif
			pcl_in_elem(p_in_e),
			pcl_num(p_num),
			node_digit_num(n_dgt_num),
			sort_mem(sm),
			tmp_valid_elems(tmp_ves),
			out_valid_elems(out_ves) {}
		~SortTriMeshNodeTask() {}
		tbb::task *execute() override
		{
			const size_t node_digit_pos = node_digit_num - 1;
			SortBin* const root_bins = sort_mem.root_bin;
			size_t* const in_node_has_elem = sort_mem.key_arrays[node_digit_num];
			size_t* const in_node_elem_pair = sort_mem.val_arrays[node_digit_num];
			if (pcl_num > serial_sort_min_pcl_num)
			{
				size_t* const tmp_node_has_elem = sort_mem.tmp_key_arrays[node_digit_pos];
				size_t* const tmp_node_elem_pair = sort_mem.tmp_val_arrays[node_digit_pos];
				size_t* const out_node_has_elem = sort_mem.key_arrays[node_digit_pos];
				size_t* const out_node_elem_pair = sort_mem.val_arrays[node_digit_pos];
				size_t bin_id;
				if (pcl_num > parallel_divide_min_pcl_num_per_block)
				{
					// block number
					const size_t blk_num = (pcl_num + parallel_divide_min_pcl_num_per_block - 1)
											 / parallel_divide_min_pcl_num_per_block;
					set_ref_count(blk_num + 1);
					size_t blk_id, pcl_end_id, pcl_start_id = 0;
					for (blk_id = 1; blk_id < blk_num; ++blk_id)
					{
#define Block_Low(block_id, block_num, data_num) ((block_id) * (data_num) / (block_num))
						pcl_end_id = Block_Low(blk_id, blk_num, pcl_num);
#undef Block_Low
						spawn(*new(allocate_child())
							Internal::CountSortTriMeshNodeTask<ElemNodeIndex>(
#ifdef _DEBUG
								elem_num,
#endif
								elem_node_ids,
								pcl_in_elem,
								pcl_start_id,
								pcl_end_id,
								node_digit_pos,
								root_bins[blk_id - 1],
								tmp_valid_elems,
								in_node_has_elem,
								in_node_elem_pair,
								tmp_node_has_elem,
								tmp_node_elem_pair));
						pcl_start_id = pcl_end_id;
					}
					// last child
					spawn_and_wait_for_all(*new(allocate_child())
						Internal::CountSortTriMeshNodeTask<ElemNodeIndex>(
#ifdef _DEBUG
							elem_num,
#endif
							elem_node_ids,
							pcl_in_elem,
							pcl_start_id,
							pcl_num,
							node_digit_pos,
							root_bins[blk_num - 1],
							tmp_valid_elems,
							in_node_has_elem,
							in_node_elem_pair,
							tmp_node_has_elem,
							tmp_node_elem_pair));

					size_t start_id;
					if (digit_pos) // not the last digit
					{
						set_ref_count(blk_num + Internal::radix_bucket_num + 1);

						// collect valid element into out_valid_elems
						start_id = 0;
						for (blk_id = 1; blk_id < blk_num; ++blk_id)
						{
							const SortBin& bin = root_bins[blk_id - 1];
							const size_t* sbin = bin.sum_bin;
							const size_t valid_elem_num = (sbin[Internal::radix_bucket_num - 1]
								+ bin.count_bin[Internal::radix_bucket_num - 1] - sbin[0]) / 3;
							spawn(*new(allocate_child())
								Internal::CopyMemTask(
									out_valid_elems + start_id,
									tmp_valid_elems + sbin[0] / 3,
									sizeof(size_t) * valid_elem_num));
							start_id += valid_elem_num;
						}

						// sort nodes
						start_id = 0;
						for (bin_id = 0; bin_id < Internal::radix_bucket_num; ++bin_id)
						{
							size_t child_start_id = 0;
							for (blk_id = 0; blk_id < blk_num; ++blk_id)
							{
								SortBin& blk_bin = root_bins[blk_id];
								const size_t blk_data_off = blk_bin.sum_bin[bin_id];
								const size_t blk_data_num = blk_bin.count_bin[bin_id];
								memcpy(out_node_has_elem + start_id + child_start_id,
									tmp_node_has_elem + blk_data_off,
									blk_data_num * sizeof(size_t));
								memcpy(out_node_elem_pair + start_id + child_start_id,
									tmp_node_elem_pair + blk_data_off,
									blk_data_num * sizeof(size_t));
								child_start_id += blk_data_num;
							}
							if (child_start_id)
							{
								spawn(*new(allocate_child())
									SortTask(
										sort_mem,
										start_id,
										child_start_id,
										node_digit_pos - 1));
								start_id += child_start_id;
							}
							else
								decrement_ref_count();
						}
					}
					else // already the last digit
					{
						set_ref_count(blk_num + 1);

						// collect valid element into out_valid_elems
						start_id = 0;
						for (blk_id = 1; blk_id < blk_num; ++blk_id)
						{
							const SortBin& bin = root_bins[blk_id - 1];
							const size_t* sbin = bin.sum_bin;
							const size_t valid_elem_num = (sbin[Internal::radix_bucket_num - 1]
								+ bin.count_bin[Internal::radix_bucket_num - 1] - sbin[0]) / 3;
							spawn(*new(allocate_child())
								Internal::CopyMemTask(
									out_valid_elems + start_id,
									tmp_valid_elems + sbin[0] / 3,
									sizeof(size_t) * valid_elem_num));
							start_id += valid_elem_num;
						}

						// sort nodes
						start_id = 0;
						for (bin_id = 0; bin_id < Internal::radix_bucket_num; ++bin_id)
						{
							size_t child_start_id = 0;
							for (blk_id = 0; blk_id < blk_num; ++blk_id)
							{
								SortBin& blk_bin = th_bins[blk_id];
								const size_t blk_data_off = blk_bin.sum_bin[bin_id];
								const size_t blk_data_num = blk_bin.count_bin[bin_id];
								memcpy(out_node_has_elem + start_id + child_start_id,
									mid_node_has_elem + blk_data_off,
									blk_data_num * sizeof(size_t));
								memcpy(out_node_elem_pair + start_id + child_start_id,
									mid_node_elem_pair + blk_data_off,
									blk_data_num * sizeof(size_t));
								child_start_id += blk_data_num;
							}
						}
					}
				}
				else // divide serially and sort parallely
				{
					Internal::count_sort_tri_mesh_node<ElemNodeIndex>(
#ifdef _DEBUG
						elem_num,
#endif
						elem_node_index,
						pcl_in_elem,
						pcl_num,
						node_digit_pos,
						root_bins[0],
						out_valid_elems,
						tmp_node_has_elem,
						tmp_node_elem_pair,
						out_node_has_elem,
						out_node_elem_pair);
					if (node_digit_pos) // not the last digit
					{
						size_t* const c_bin = root_bins[0].count_bin;
						size_t* const s_bin = root_bins[0].sum_bin;
						set_ref_count(Internal::radix_bucket_num + 1);
						for (bin_id = 0; bin_id < Internal::radix_bucket_num; ++bin_id)
						{
							if (c_bin[bin_id])
							{
								spawn(*new(allocate_child())
									SortTask(
										sort_mem,
										s_bin[bin_id],
										c_bin[bin_id],
										node_digit_pos - 1));
							}
							else
								decrement_ref_count();
						}
					}
				}

				wait_for_all();
			}
			else // sort serially
			{
				size_t e_id = pcl_in_elem[p_id + 1];
#ifdef _DEBUG
				assert(e_id < elem_num || e_id == SIZE_MAX);
#endif
				size_t valid_elem_num = 0;
				for (size_t p_id = 0; p_id < pcl_num; ++p_id)
				{
					if (e_id != pcl_in_elem[p_id + 1])
					{
						const ElemNodeIndex& eni = elem_node_id[e_id];
						in_node_has_elem[valid_elem_num * 3] = eni.n1;
						in_node_has_elem[valid_elem_num * 3 + 1] = eni.n2;
						in_node_has_elem[valid_elem_num * 3 + 2] = eni.n3;
						in_node_elem_pair[valid_elem_num * 3] = e_id * 3;
						in_node_elem_pair[valid_elem_num * 3 + 1] = e_id * 3 + 1;
						in_node_elem_pair[valid_elem_num * 3 + 2] = e_id * 3 + 2;
						out_valid_elems[valid_elem_num++] = e_id;
						e_id = pcl_in_elem[p_id + 1];
#ifdef _DEBUG
						assert(e_id < elem_num || e_id == SIZE_MAX);
#endif
					}
				}
				Internal::serial_sort(
					sort_mem.out_keys,
					sort_mem.out_vals,
					in_node_has_elem,
					in_node_elem_pair,
					valid_elem_num * 3,
					node_digit_pos,
					root_bins[0],
					sort_mem.tmp_key_arrays[0],
					sort_mem.tmp_val_arrays[0]);
			}
			return nullptr;
		}
	};
}

#endif