#include "SimulationsOMP_pcp.h"

#include "SortParticleTask.h"

namespace SortUtils
{
	SortParticleTask::SortParticleTask(
		SortMem &sm,
		size_t p_num,
		size_t p_dgt_num) :
		sort_mem(sm),
		pcl_num(p_num),
		pcl_digit_num(p_dgt_num) {}

	tbb::task *SortParticleTask::execute()
	{

		return nullptr;
	}
}
