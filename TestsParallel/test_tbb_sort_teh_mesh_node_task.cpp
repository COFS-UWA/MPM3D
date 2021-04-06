#include "TestsParallel_pcp.h"

#include <assert.h>
#include <memory>
#include "tbb/task_scheduler_init.h"

#include "SortTehMeshNodeTask.h"

#include "test_parallel_utils.h"
#include "test_simulations_omp.h"

void test_sort_teh_mesh_node_task()
{
	constexpr size_t pcl_num = 100;
	constexpr size_t elem_num = 30;
	constexpr size_t node_num = 200000;
	constexpr size_t thread_num = 6;

	tbb::task_scheduler_init init(thread_num);
	MSDRadixSortUtils::RadixBinBlockMemArray radix_bin_blocks;
	radix_bin_blocks.init(thread_num, 2);

	struct ElemNodeId { size_t n1, n2, n3, n4; };
	for (size_t i = 0; i < 1; i++)
	{
		std::unique_ptr<ElemNodeId> elem_node_ids_ptr(new ElemNodeId[elem_num]);
		ElemNodeId* elem_node_ids = elem_node_ids_ptr.get();
		srand(i);
		for (size_t e_id = 0; e_id < elem_num; ++e_id)
		{
			elem_node_ids[e_id].n1 = rand() % node_num;
			elem_node_ids[e_id].n2 = rand() % node_num;
			elem_node_ids[e_id].n3 = rand() % node_num;
			elem_node_ids[e_id].n4 = rand() % node_num;
		}

		std::unique_ptr<size_t> pcl_in_elem_ptr(new size_t[pcl_num]);
		size_t* pcl_in_elem = pcl_in_elem_ptr.get();
		const size_t pcl_per_elem = (pcl_num + elem_num - 1) / elem_num;
		size_t e_id = 0, p_in_e_id = 0;
		for (size_t i = 0; i < pcl_num; ++i, ++p_in_e_id)
		{
			pcl_in_elem[i] = e_id;
			if (p_in_e_id >= pcl_per_elem)
			{
				++e_id;
				p_in_e_id = 0;
			}
		}

		std::cout << "pcl_in_elem:\n";
		print_array(pcl_in_elem, pcl_num);

		SortTehMeshNodeMem node_sort_mem;
		node_sort_mem.init(elem_num, node_num, radix_bin_blocks);

		size_t valid_elem;
		tbb::task::spawn_root_and_wait(
			*new(tbb::task::allocate_root())
			SortTehMeshNodeTask(
				node_sort_mem,
				pcl_num,
				pcl_in_elem,
				elem_node_ids,
				valid_elem));
		std::cout << i << ", valid elem num: " << valid_elem << "\n\n";

		std::cout << "valid elem:\n";
		print_array(node_sort_mem.res_elems, valid_elem);
		std::cout << "nodes has elem:\n";
		print_array(node_sort_mem.res_keys, valid_elem * 4);
		std::cout << "node_elem_pair:\n";
		print_array(node_sort_mem.res_vals, valid_elem * 4);

		for (size_t i = 1; i < valid_elem * 4; ++i)
			assert(node_sort_mem.res_keys[i - 1] <= node_sort_mem.res_keys[i]);
	}
}
