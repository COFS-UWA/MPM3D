#include "SimulationsOMP_pcp.h"

#include "SortTriMeshNodeTask.h"

#define __Cal_Key_Digit__(data, digit_pos) (((data) >> ((digit_pos) * 8)) & (0xFF))	
#define Block_Low(block_id, block_num, data_num) ((block_id) * (data_num) / (block_num))

namespace SortUtils
{
	void SortTriMeshNodeMem::init(
		void* shared_mem,
		size_t thread_num,
		size_t elem_num,
		size_t node_digit_num)
	{
		clear();
		max_block_num = thread_num * max_block_num_div_thread_num;
		const size_t three_elem_num = elem_num * 3;
		size_t block_num = (three_elem_num + parallel_divide_min_pcl_num_per_block - 1)
							/ parallel_divide_min_pcl_num_per_block;
		if (block_num > max_block_num)
			block_num = max_block_num;
		char* cur_mem = (char*)self_mem.alloc(
				sizeof(size_t *) * (node_digit_num * 4 + 2)
			+ sizeof(SortBin*) * thread_num
			+ sizeof(size_t) * (three_elem_num * 2 + 2)
			+ sizeof(size_t) * elem_num
			+ Cache_Alignment * 3);
		//
		node_has_elem_arrays = (size_t **)cur_mem; // node_digit_num + 1 
		cur_mem += sizeof(size_t *) * (node_digit_num + 1);
		node_elem_pair_arrays = (size_t**)cur_mem; // node_digit_num + 1 
		cur_mem += sizeof(size_t*) * (node_digit_num + 1);
		tmp_node_has_elem_arrays = (size_t**)cur_mem; // node_digit_num
		cur_mem += sizeof(size_t*) * node_digit_num;
		tmp_node_elem_pair_arrays = (size_t**)cur_mem; // pcl_digit_num
		cur_mem += sizeof(size_t*) * node_digit_num;
		thread_bins = (SortBin**)cur_mem;
		cur_mem = cache_aligned(cur_mem + sizeof(SortBin*) * thread_num);
		//
		node_has_elem = ((size_t *)cur_mem) + 1;
		cur_mem = cache_aligned(cur_mem + sizeof(size_t) * (three_elem_num + 2));
		node_elem_pair = (size_t *)cur_mem;
		cur_mem = cache_aligned(cur_mem + sizeof(size_t) * elem_num);
		node_has_elem_arrays[0] = node_has_elem;
		node_elem_pair_arrays[0] = node_elem_pair;
		node_elem_pair[-1] = SIZE_MAX;
		node_elem_pair[three_elem_num] = SIZE_MAX;
		valid_elems = (size_t *)cur_mem;
		
		cur_mem = (char*)cache_aligned(shared_mem);
		// key and value memory
		size_t i;
		for (i = 0; i < node_digit_num; ++i)
		{
			node_has_elem_arrays[i + 1] = (size_t*)cur_mem;
			cur_mem = cache_aligned(cur_mem + sizeof(size_t) * three_elem_num);
		}
		for (i = 0; i < node_digit_num; ++i)
		{
			node_elem_pair_arrays[i + 1] = (size_t*)cur_mem;
			cur_mem = cache_aligned(cur_mem + sizeof(size_t) * three_elem_num);
		}
		for (i = 0; i < node_digit_num; ++i)
		{
			tmp_node_has_elem_arrays[i] = (size_t*)cur_mem;
			cur_mem = cache_aligned(cur_mem + sizeof(size_t) * three_elem_num);
		}
		for (i = 0; i < node_digit_num; ++i)
		{
			tmp_node_elem_pair_arrays[i] = (size_t*)cur_mem;
			cur_mem = cache_aligned(cur_mem + sizeof(size_t) * three_elem_num);
		}
		ori_node_has_elem = node_has_elem_arrays[node_digit_num];
		ori_node_elem_pair = node_elem_pair_arrays[node_digit_num];
		tmp_valid_elems = (size_t*)cur_mem;
		cur_mem = cache_aligned(cur_mem + sizeof(size_t) * elem_num);

		// sort bin memory
		root_bin = (SortBin*)cur_mem;
		const size_t sort_bin_size = sizeof(SortBin) * block_num;
		for (i = 0; i < thread_num; ++i)
		{
			cur_mem += sort_bin_size;
			thread_bins[i] = (SortBin*)cur_mem;
		}
	}

	namespace Internal
	{
		tbb::task* CountSortTriMeshNodeTask::execute()
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
			e_id = pcl_in_elem[p_id0];
#ifdef _DEBUG
			assert(e_id < elem_num || e_id == SIZE_MAX);
#endif
			size_t* const my_valid_elems = tmp_valid_elems + e_id;
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
	}

	tbb::task* SortTriMeshNodeTask::execute()
	{
		SortBin* const rt_bins = sort_mem.root_bin;
		size_t* const in_node_has_elem = sort_mem.ori_node_has_elem;
		size_t* const in_node_elem_pair = sort_mem.ori_node_elem_pair;
		size_t* const valid_elems = sort_mem.valid_elems;
		const size_t node_digit_pos = node_digit_num - 1;
		if (pcl_num > serial_sort_min_pcl_num)
		{
			size_t* const out_node_has_elem = sort_mem.node_has_elem_arrays[node_digit_pos];
			size_t* const out_node_elem_pair = sort_mem.node_elem_pair_arrays[node_digit_pos];
			size_t bin_id;
			if (pcl_num > parallel_divide_min_pcl_num_per_block)
			{
				size_t* const tmp_node_has_elem = sort_mem.tmp_node_has_elem_arrays[node_digit_pos];
				size_t* const tmp_node_elem_pair = sort_mem.tmp_node_elem_pair_arrays[node_digit_pos];
				size_t* const tmp_valid_elems = sort_mem.tmp_valid_elems;
				size_t blk_num = (pcl_num + SortTriMeshNodeMem::parallel_divide_min_pcl_num_per_block - 1)
								/ SortTriMeshNodeMem::parallel_divide_min_pcl_num_per_block;
				if (blk_num > sort_mem.max_block_num)
					blk_num = sort_mem.max_block_num;
				set_ref_count(blk_num + 1);
				size_t blk_id, pcl_end_id, pcl_start_id = 0;
				for (blk_id = 1; blk_id < blk_num; ++blk_id)
				{
					pcl_end_id = Block_Low(blk_id, blk_num, pcl_num);
					spawn(*new(allocate_child())
						Internal::CountSortTriMeshNodeTask(
#ifdef _DEBUG
							pcl_num,
							elem_num,
#endif
							elem_node_ids,
							pcl_in_elem,
							pcl_start_id,
							pcl_end_id,
							node_digit_pos,
							rt_bins[blk_id - 1],
							sort_mem.tmp_valid_elems,
							in_node_has_elem,
							in_node_elem_pair,
							tmp_node_has_elem,
							tmp_node_elem_pair));
					pcl_start_id = pcl_end_id;
				}
				// last child
				spawn_and_wait_for_all(*new(allocate_child())
					Internal::CountSortTriMeshNodeTask(
#ifdef _DEBUG
						pcl_num,
						elem_num,
#endif
						elem_node_ids,
						pcl_in_elem,
						pcl_start_id,
						pcl_num,
						node_digit_pos,
						rt_bins[blk_num - 1],
						sort_mem.tmp_valid_elems,
						in_node_has_elem,
						in_node_elem_pair,
						tmp_node_has_elem,
						tmp_node_elem_pair));

				size_t start_id;
				if (node_digit_pos) // not the last digit
				{
					set_ref_count(blk_num + Internal::radix_bucket_num + 1);

					// collect valid element into out_valid_elems
					start_id = 0;
					for (blk_id = 0; blk_id < blk_num; ++blk_id)
					{
						const SortBin& bin = rt_bins[blk_id];
						const size_t* sbin = bin.sum_bin;
						const size_t ve_num = (bin.count_bin[Internal::radix_bucket_num - 1]
							+ sbin[Internal::radix_bucket_num - 1] - sbin[0]) / 3;
						spawn(*new(allocate_child())
							Internal::CopyMemTask(
								valid_elems + start_id,
								tmp_valid_elems + sbin[0] / 3,
								sizeof(size_t) * ve_num));
						start_id += ve_num;
					}
					valid_elem_num = start_id;

					// sort nodes
					start_id = 0;
					for (bin_id = 0; bin_id < Internal::radix_bucket_num; ++bin_id)
					{
						size_t child_start_id = 0;
						for (blk_id = 0; blk_id < blk_num; ++blk_id)
						{
							SortBin& blk_bin = rt_bins[blk_id];
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
									sort_mem._sort_mem,
									start_id,
									child_start_id,
									node_digit_pos - 1));
							start_id += child_start_id;
						}
						else
							decrement_ref_count();
					}
					wait_for_all();
				}
				else // already the last digit
				{
					set_ref_count(blk_num + 1);

					// collect valid element into out_valid_elems
					start_id = 0;
					for (blk_id = 0; blk_id < blk_num; ++blk_id)
					{
						const SortBin& bin = rt_bins[blk_id];
						const size_t* sbin = bin.sum_bin;
						const size_t ve_num = (bin.count_bin[Internal::radix_bucket_num - 1]
							+ sbin[Internal::radix_bucket_num - 1] - sbin[0]) / 3;
						spawn(*new(allocate_child())
							Internal::CopyMemTask(
								valid_elems + start_id,
								tmp_valid_elems + sbin[0] / 3,
								sizeof(size_t) * ve_num));
						start_id += ve_num;
					}
					valid_elem_num = start_id;

					// sort nodes
					start_id = 0;
					for (bin_id = 0; bin_id < Internal::radix_bucket_num; ++bin_id)
					{
						size_t child_start_id = 0;
						for (blk_id = 0; blk_id < blk_num; ++blk_id)
						{
							SortBin& blk_bin = rt_bins[blk_id];
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
					}
					wait_for_all();
				}
			}
			else // divide serially and sort parallely
			{
				SortBin& rt_bin = rt_bins[0];
				Internal::count_sort_tri_mesh_node(
#ifdef _DEBUG
					elem_num,
#endif
					elem_node_ids,
					pcl_in_elem,
					pcl_num,
					node_digit_pos,
					rt_bin,
					valid_elems,
					in_node_has_elem,
					in_node_elem_pair,
					out_node_has_elem,
					out_node_elem_pair);
				size_t* const c_bin = rt_bin.count_bin;
				size_t* const s_bin = rt_bin.sum_bin;
				valid_elem_num = s_bin[Internal::radix_bucket_num - 1]
							+ c_bin[Internal::radix_bucket_num - 1];
				if (node_digit_pos) // not the last digit
				{
					set_ref_count(Internal::radix_bucket_num + 1);
					for (bin_id = 0; bin_id < Internal::radix_bucket_num; ++bin_id)
					{
						if (c_bin[bin_id])
						{
							spawn(*new(allocate_child())
								SortTask(
									sort_mem._sort_mem,
									s_bin[bin_id],
									c_bin[bin_id],
									node_digit_pos - 1));
						}
						else
							decrement_ref_count();
					}
					wait_for_all();
				}
			}
		}
		else // sort serially
		{
			size_t e_id = pcl_in_elem[0];
#ifdef _DEBUG
			assert(e_id < elem_num || e_id == SIZE_MAX);
#endif
			size_t ve_num = 0;
			for (size_t p_id = 0; p_id < pcl_num; ++p_id)
			{
				if (e_id != pcl_in_elem[p_id + 1])
				{
					const Internal::ElemNodeIndex& eni = elem_node_ids[e_id];
					in_node_has_elem[ve_num * 3] = eni.n1;
					in_node_has_elem[ve_num * 3 + 1] = eni.n2;
					in_node_has_elem[ve_num * 3 + 2] = eni.n3;
					in_node_elem_pair[ve_num * 3] = e_id * 3;
					in_node_elem_pair[ve_num * 3 + 1] = e_id * 3 + 1;
					in_node_elem_pair[ve_num * 3 + 2] = e_id * 3 + 2;
					valid_elems[ve_num++] = e_id;
					e_id = pcl_in_elem[p_id + 1];
#ifdef _DEBUG
					assert(e_id < elem_num || e_id == SIZE_MAX);
#endif
				}
			}
			Internal::serial_sort(
				sort_mem.node_has_elem,
				sort_mem.node_elem_pair,
				in_node_has_elem,
				in_node_elem_pair,
				ve_num * 3,
				node_digit_pos,
				rt_bins[0],
				sort_mem.tmp_node_has_elem_arrays[0],
				sort_mem.tmp_node_elem_pair_arrays[0]);
			valid_elem_num = ve_num;
		}
		sort_mem.node_has_elem[valid_elem_num * 3] = SIZE_MAX;
		return nullptr;
	}
}
