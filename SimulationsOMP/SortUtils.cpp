#include "SimulationsOMP_pcp.h"

#include "SortUtils.h"

namespace SortUtils
{
	namespace Internal
	{
		CountSortTask::CountSortTask(
			size_t* const _out_key,
			size_t* const _out_val,
			const size_t* _in_key,
			const size_t* _in_val,
			const size_t _data_num,
			const size_t _digit_pos,
			SortBin& _bin) :
			out_key(_out_key),
			out_val(_out_val),
			in_key(_in_key),
			in_val(_in_val),
			data_num(_data_num),
			digit_pos(_digit_pos),
			bin(_bin) {}

		CountSortTask::~CountSortTask() {}

		tbb::task *CountSortTask::execute()
		{
			count_sort(
				out_key,
				out_val,
				in_key,
				out_val,
				data_num,
				digit_pos,
				bin);
			return nullptr;
		}
	}
}