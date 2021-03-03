#include "SimulationsOMP_pcp.h"

#include <assert.h>
#include "SortParticleTask.h"

void SortParticleMem::init(
	size_t pcl_num,
	RadixBinBlockMemArray &rbbs)
{
	set_radix_bin_block(rbbs);
	char *cur_mem = data_mem.alloc<char>(
		  (pcl_num * 4 + 4) * sizeof(size_t)
		+ MSDRadixSortUtils::cache_line_size * 3);
	radix_keys0 = (size_t *)cur_mem + 1;
	cur_mem = cache_aligned(cur_mem + (pcl_num + 2) * sizeof(size_t));
	radix_keys1 = (size_t*)cur_mem + 1;
	cur_mem = cache_aligned(cur_mem + (pcl_num + 2) * sizeof(size_t));
	radix_vals0 = (size_t *)cur_mem;
	cur_mem = cache_aligned(cur_mem + pcl_num * sizeof(size_t));
	radix_vals1 = (size_t *)cur_mem;
	radix_keys0[-1] = SIZE_MAX;
	radix_keys0[pcl_num] = SIZE_MAX;
	radix_keys1[-1] = SIZE_MAX;
	radix_keys1[pcl_num] = SIZE_MAX;
	set_digit_num(pcl_num);
}
