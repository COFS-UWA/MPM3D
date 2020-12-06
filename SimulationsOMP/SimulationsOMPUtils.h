#ifndef __Simulation_OMP_Utils_h__
#define __Simulation_OMP_Utils_h__

namespace SimulationsOMP
{
	template <typename Data>
	inline void swap_sort_acc(Data* datas, size_t num)
	{
		for (size_t i = 0; i < num; ++i)
			for (size_t j = num; --j > i;)
			{
				if (datas[j - 1] > datas[j])
				{
					const Data tmp = datas[j];
					datas[j] = datas[j - 1];
					datas[j - 1] = tmp;
				}
			}
	}
}

#endif