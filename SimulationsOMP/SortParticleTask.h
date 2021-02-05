#ifndef __Sort_Particle_Task_h__
#define __Sort_Particle_Task_h__

#include "SortTask.h"

namespace SortUtils
{
	class SortParticleTask : public tbb::task
	{
	protected:
		static constexpr size_t parallel_divide_min_pcl_num_per_block = 2;
		static constexpr size_t serial_sort_min_pcl_num = 1;
		const size_t pcl_num;
		const size_t pcl_digit_num;
		SortMem& sort_mem;
	public:
		SortParticleTask(
			SortMem &sm,
			size_t p_num,
			size_t p_dgt_num);
		~SortParticleTask() {}
		tbb::task* execute() override;
	};
}

#endif