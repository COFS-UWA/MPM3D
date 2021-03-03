#ifndef __Sort_Particle_Task_h__
#define __Sort_Particle_Task_h__

#include "CacheAlignedMem.h"
#include "MSDRadixSortTask.h"

struct SortParticleMem : public MSDRadixSortUtils::SortMem
{
	using RadixBinBlockMemArray = MSDRadixSortUtils::RadixBinBlockMemArray;

	inline void update_key_and_val() noexcept
	{
		if (ori_keys == res_keys)
		{
			ori_keys = radix_keys0;
			radix_keys0 = radix_keys1;
			// radix_keys1 == res_keys
			radix_keys1 = ori_keys;
			ori_vals = radix_vals0;
			radix_vals0 = radix_vals1;
			// radix_vals1 == res_vals
			radix_vals1 = ori_vals;
		}
	}
	const size_t *get_prev_res_keys() const 
	{ return ori_keys == res_keys ? radix_vals0 : res_keys; }

	void init(size_t pcl_num, size_t elem_num, RadixBinBlockMemArray& rbbs);

protected:
	CacheAlignedMem data_mem;
};

// Assumption:
//     pcl_id < pcl_num (maximum pcl_id == pcl_num - 1)
//     => pcl_num <= SIZE_MAX, max(pcl_id) <= SIZE_MAX - 1
class SortParticleTask : public MSDRadixSortTask
{
public:
	SortParticleTask(SortParticleMem& _spm, size_t pcl_num) :
		MSDRadixSortTask(_spm, 0, pcl_num, _spm.digit_num-1) {}
};

#endif