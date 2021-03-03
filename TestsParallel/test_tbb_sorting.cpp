#include "TestsParallel_pcp.h"

#include <assert.h>
#include <memory>
#include "tbb/task_scheduler_init.h"

#include "SortParticleTask.h"
#include "SortTriMeshNodeTask.h"
#include "TriangleMesh.h"

#include "test_simulations_omp.h"

void print_array(const size_t *data, size_t num)
{
	constexpr size_t data_per_line = 20;
	for (size_t i = 0; i < num; ++i)
	{
		std::cout << data[i] << ", ";
		if (i % data_per_line == (data_per_line - 1))
			std::cout << "\n";
	}
	if ((num - 1) % data_per_line != (data_per_line - 1))
		std::cout << "\n";
	std::cout << "\n";
}

void test_sort_pcl_task()
{
	constexpr size_t pcl_num = 255;
	//constexpr size_t pcl_num = 100;
	//constexpr size_t pcl_num = 20;
	//constexpr size_t pcl_num = 10;
	//constexpr size_t pcl_num = 2;
	//constexpr size_t thread_num = 1;
	constexpr size_t thread_num = 6;

	MSDRadixSortUtils::RadixBinBlockMemArray radix_bin_blocks;
	SortParticleMem pcl_sort_mem;
	size_t *ori_keys, *ori_vals, *res_keys, *res_vals;

	tbb::task_scheduler_init init(thread_num);
	radix_bin_blocks.init(thread_num, 2);
	pcl_sort_mem.init(pcl_num, radix_bin_blocks);

	ori_keys = pcl_sort_mem.ori_keys;
	ori_vals = pcl_sort_mem.ori_vals;
	for (size_t i = 0; i < pcl_num; ++i)
	{
		ori_keys[i] = pcl_num - i + 1;
		ori_vals[i] = i;
	}
	ori_keys[0] = SIZE_MAX;
	ori_keys[1] = SIZE_MAX;
	print_array(ori_keys, pcl_num);

	tbb::task::spawn_root_and_wait(
		*new(tbb::task::allocate_root())
		SortParticleTask(pcl_sort_mem, pcl_num));

	res_keys = pcl_sort_mem.res_keys;
	res_vals = pcl_sort_mem.res_vals;
	print_array(res_keys, pcl_num);
	//print_array(res_vals, pcl_num);

	for (size_t i = 1; i < pcl_num; ++i)
		assert(res_keys[i - 1] <= res_keys[i]);

	pcl_sort_mem.update_key_and_val();
	ori_keys = pcl_sort_mem.ori_keys;
	ori_vals = pcl_sort_mem.ori_vals;
	// res will be uncorrect because (ori_key[i]) > pcl_num
	// but just for testing
	for (size_t i = 0; i < pcl_num; ++i)
	{
		ori_keys[i] = pcl_num - i + 2;
		ori_vals[i] = i;
	}
	std::cout << "new ori_keys:\n";
	print_array(ori_keys, pcl_num);
	std::cout << "old res_keys:\n";
	print_array(res_keys, pcl_num); // old res key

	tbb::task::spawn_root_and_wait(
		*new(tbb::task::allocate_root())
		SortParticleTask(pcl_sort_mem, pcl_num));
	
	res_keys = pcl_sort_mem.res_keys;
	res_vals = pcl_sort_mem.res_vals;
	std::cout << "new res_keys:\n";
	print_array(res_keys, pcl_num);
}

void test_sort_tri_mesh_node_task()
{
	//constexpr size_t pcl_num = 300;
	//constexpr size_t pcl_num = 100;
	//constexpr size_t pcl_num = 50;
	//constexpr size_t pcl_num = 10;
	constexpr size_t pcl_num = 2;
	constexpr size_t thread_num = 6;

	tbb::task_scheduler_init init(thread_num);

	MSDRadixSortUtils::RadixBinBlockMemArray radix_bin_blocks;
	radix_bin_blocks.init(thread_num, 2);

	TriangleMesh tri_mesh;
	tri_mesh.load_mesh_from_hdf5("../../Asset/rect_mesh.h5");
	const size_t elem_num = tri_mesh.get_elem_num();
	std::cout << "elem: " << elem_num << ", node: "
		<< tri_mesh.get_node_num() << "\n";
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
	size_t *pcl_in_elem = pcl_in_elem_ptr.get();
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
	print_array(node_sort_mem.res_elems, valid_elem);
	//std::cout << "ori nodes has elem:\n";
	//print_array(node_sort_mem.ori_keys, valid_elem * 3);
	//std::cout << "ori node_elem_pair:\n";
	//print_array(node_sort_mem.ori_vals, valid_elem * 3);
	std::cout << "nodes has elem:\n";
	print_array(node_sort_mem.res_keys, valid_elem * 3);
	std::cout << "node_elem_pair:\n";
	print_array(node_sort_mem.res_vals, valid_elem * 3);

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
