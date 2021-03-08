#include "TestsParallel_pcp.h"

#include <assert.h>
#include <memory>
#include "tbb/task_scheduler_init.h"

#include "SortParticleTask.h"
#include "SortTriMeshNodeTask.h"
#include "TriangleMesh.h"

#include "test_parallel_utils.h"
#include "test_simulations_omp.h"

void test_sort_tri_mesh_node_task()
{
	constexpr size_t pcl_num = 30000;
	//constexpr size_t pcl_num = 100;
	//constexpr size_t pcl_num = 50;
	//constexpr size_t pcl_num = 10;
	//constexpr size_t pcl_num = 2;
	constexpr size_t thread_num = 6;

	tbb::task_scheduler_init init(thread_num);

	MSDRadixSortUtils::RadixBinBlockMemArray radix_bin_blocks;
	radix_bin_blocks.init(thread_num, 2);

	TriangleMesh tri_mesh;
	tri_mesh.load_mesh_from_hdf5("../../Asset/rect_mesh.h5");
	const size_t elem_num = tri_mesh.get_elem_num();
	std::cout << "elem: " << elem_num << ", node: " << tri_mesh.get_node_num() << "\n";
	auto* elems = tri_mesh.get_elems();
	struct ElemNodeId { size_t n1, n2, n3; };
	std::unique_ptr<ElemNodeId> elem_node_ids_ptr(new ElemNodeId[elem_num]);
	ElemNodeId* elem_node_ids = elem_node_ids_ptr.get();
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		elem_node_ids[e_id].n1 = elems[e_id].n1;
		elem_node_ids[e_id].n2 = elems[e_id].n2;
		elem_node_ids[e_id].n3 = elems[e_id].n3;
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
	//print_array(pcl_in_elem, pcl_num);

	SortTriMeshNodeMem node_sort_mem;
	node_sort_mem.init(
		tri_mesh.get_elem_num(),
		tri_mesh.get_node_num(),
		radix_bin_blocks);

	size_t valid_elem;
	tbb::task::spawn_root_and_wait(
		*new(tbb::task::allocate_root())
			SortTriMeshNodeTask(
				node_sort_mem,
				pcl_num,
				pcl_in_elem,
				elem_node_ids,
				valid_elem));
	std::cout << "valid elem num: " << valid_elem << "\n\n";

	std::cout << "valid elem:\n";
	//print_array(node_sort_mem.res_elems, valid_elem);
	//std::cout << "ori nodes has elem:\n";
	//print_array(node_sort_mem.ori_keys, valid_elem * 3);
	//std::cout << "ori node_elem_pair:\n";
	//print_array(node_sort_mem.ori_vals, valid_elem * 3);
	std::cout << "nodes has elem:\n";
	//print_array(node_sort_mem.res_keys, valid_elem * 3);
	std::cout << "node_elem_pair:\n";
	//print_array(node_sort_mem.res_vals, valid_elem * 3);

	//for (size_t i = 0; i < valid_elem * 3; ++i)
	//{
	//	size_t res_val = node_sort_mem.res_vals[i];
	//	size_t res_key = node_sort_mem.res_keys[i];
	//	size_t ori_key = node_sort_mem.ori_keys[res_val];
	//	assert(node_sort_mem.res_keys[i]
	//		== node_sort_mem.ori_keys[node_sort_mem.res_vals[i]]);
	//}

	for (size_t i = 1; i < valid_elem * 3; ++i)
		assert(node_sort_mem.res_keys[i - 1] <= node_sort_mem.res_keys[i]);
}

void test_sort_tri_mesh_node_task2()
{
	constexpr size_t pcl_num = 300000;
	constexpr size_t elem_num = 150000;
	constexpr size_t node_num = 200000;
	constexpr size_t thread_num = 6;

	tbb::task_scheduler_init init(thread_num);

	MSDRadixSortUtils::RadixBinBlockMemArray radix_bin_blocks;
	radix_bin_blocks.init(thread_num, 2);

	struct ElemNodeId { size_t n1, n2, n3; };
	for (size_t i = 0; i < 1000; i++)
	{
		std::unique_ptr<ElemNodeId> elem_node_ids_ptr(new ElemNodeId[elem_num]);
		ElemNodeId* elem_node_ids = elem_node_ids_ptr.get();
		srand(i);
		for (size_t e_id = 0; e_id < elem_num; ++e_id)
		{
			elem_node_ids[e_id].n1 = rand() % node_num;
			elem_node_ids[e_id].n2 = rand() % node_num;
			elem_node_ids[e_id].n3 = rand() % node_num;
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

		//std::cout << "pcl_in_elem:\n";
		//print_array(pcl_in_elem, pcl_num);

		SortTriMeshNodeMem node_sort_mem;
		node_sort_mem.init(elem_num, node_num, radix_bin_blocks);

		size_t valid_elem;
		tbb::task::spawn_root_and_wait(
			*new(tbb::task::allocate_root())
			SortTriMeshNodeTask(
				node_sort_mem,
				pcl_num,
				pcl_in_elem,
				elem_node_ids,
				valid_elem));
		std::cout << i << ", valid elem num: " << valid_elem << "\n\n";

		//std::cout << "valid elem:\n";
		//print_array(node_sort_mem.res_elems, valid_elem);
		//std::cout << "nodes has elem:\n";
		//print_array(node_sort_mem.res_keys, valid_elem * 3);
		//std::cout << "node_elem_pair:\n";
		//print_array(node_sort_mem.res_vals, valid_elem * 3);

		//for (size_t i = 0; i < valid_elem * 3; ++i)
		//{
		//	size_t res_val = node_sort_mem.res_vals[i];
		//	size_t res_key = node_sort_mem.res_keys[i];
		//	size_t ori_key = node_sort_mem.ori_keys[res_val];
		//	assert(node_sort_mem.res_keys[i]
		//		== node_sort_mem.ori_keys[node_sort_mem.res_vals[i]]);
		//}

		for (size_t i = 1; i < valid_elem * 3; ++i)
			assert(node_sort_mem.res_keys[i - 1] <= node_sort_mem.res_keys[i]);
	}
}
