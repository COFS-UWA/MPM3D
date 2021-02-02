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
			SortBin& _bin,
			const unsigned char _digit_pos) :
			out_key(_out_key),
			out_val(_out_val),
			in_key(_in_key),
			in_val(_in_val),
			data_num(_data_num),
			bin(_bin),
			digit_pos(_digit_pos) {}

		CountSortTask::~CountSortTask() {}

		tbb::task *CountSortTask::execute()
		{
			count_sort(
				out_key,
				out_val,
				in_key,
				out_val,
				data_num,
				bin,
				digit_pos);
			return nullptr;
		}
	}
}