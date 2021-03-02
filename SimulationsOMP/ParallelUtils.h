#ifndef __Parallel_Utils_h__
#define __Parallel_Utils_h__

namespace ParallelUtils
{
	inline size_t block_low(size_t block_id, size_t block_num, size_t data_num) noexcept
	{
		return block_id * data_num / block_num;
	}
}

#endif