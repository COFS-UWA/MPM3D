#include "TestsParallel_pcp.h"

#include <assert.h>
#include <memory>
#include "tbb/task_scheduler_init.h"

#include "SortParticleTask.h"
#include "SortTriMeshNodeTask.h"
#include "TriangleMesh.h"

#include "test_parallel_utils.h"
#include "test_simulations_omp.h"

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
	pcl_sort_mem.init(pcl_num, pcl_num, radix_bin_blocks);

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

void test_sort_pcl_task2()
{
	//constexpr size_t thread_num = 1;
	constexpr size_t thread_num = 6;
	constexpr size_t pcl_num = 65536;
	constexpr size_t max_elem_num = 5500;

	constexpr size_t dgt_num0 = MSDRadixSortUtils::max_digit_num(0);
	constexpr size_t dgt_num1 = MSDRadixSortUtils::max_digit_num(255);
	constexpr size_t dgt_num2 = MSDRadixSortUtils::max_digit_num(256);
	constexpr size_t dgt_num3 = MSDRadixSortUtils::max_digit_num(65535);
	constexpr size_t dgt_num4 = MSDRadixSortUtils::max_digit_num(65536);

	tbb::task_scheduler_init init(thread_num);
	MSDRadixSortUtils::RadixBinBlockMemArray radix_bin_blocks;
	radix_bin_blocks.init(thread_num, 2);

	for (size_t i = 0; i < 1000; i++)
	{
		SortParticleMem pcl_sort_mem;
		pcl_sort_mem.init(pcl_num, max_elem_num, radix_bin_blocks);

		srand(i);
		size_t* ori_keys = pcl_sort_mem.ori_keys;
		size_t* ori_vals = pcl_sort_mem.ori_vals;
		for (size_t i = 0; i < pcl_num; ++i)
		{
			ori_keys[i] = rand() % max_elem_num;
			ori_vals[i] = i;
		}
		ori_keys[0] = SIZE_MAX;
		ori_keys[1] = SIZE_MAX;
		//print_array(ori_keys, pcl_num);

		//tbb::task::spawn_root_and_wait(
		//	*new(tbb::task::allocate_root())
		//	SortParticleTask(pcl_sort_mem, pcl_num));
		tbb::task::spawn_root_and_wait(
			*new(tbb::task::allocate_root())
				MSDRadixSortTask(pcl_sort_mem, 0,
					pcl_num, pcl_sort_mem.digit_num - 1));

		size_t* res_keys = pcl_sort_mem.res_keys;
		size_t* res_vals = pcl_sort_mem.res_vals;
		//print_array(res_keys, pcl_num);
		//print_array(res_vals, pcl_num);

		for (size_t i = 1; i < pcl_num; ++i)
			assert(res_keys[i - 1] <= res_keys[i]);

		std::cout << i << "\n";
	}
}
