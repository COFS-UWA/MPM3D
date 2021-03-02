#ifndef __MSD_Radix_Sort_Task_h__
#define __MSD_Radix_Sort_Task_h__

#include "CacheAlignedMem.h"
#include "MSDRadixSortUtils.h"

namespace MSDRadixSortUtils
{
	class ScanAndFormBin
	{
	protected:
		const size_t block_num;
		const size_t data_num;
		const size_t digit_pos;
		RadixBin* const radix_bin_block;
		const size_t* const keys;

	public:
		ScanAndFormBin(
			size_t blk_num,
			size_t dat_num,
			size_t dgt_pos,
			RadixBin* bin_blk,
			const size_t* ks) :
			block_num(blk_num),
			data_num(dat_num),
			digit_pos(dgt_pos),
			radix_bin_block(bin_blk),
			keys(ks) {}
		~ScanAndFormBin() {}
		void operator() (size_t blk_id);
	};

	class MoveAccToBin
	{
	protected:
		const size_t block_num;
		const size_t data_num;
		const size_t digit_pos;
		const RadixBin* const radix_bin_block;
		const size_t *const in_keys;
		const size_t *const in_vals;
		size_t *const out_keys;
		size_t *const out_vals;

	public:
		MoveAccToBin(
			size_t blk_num,
			size_t dat_num,
			size_t dgt_pos,
			const RadixBin* bin_blk,
			const size_t* in_ks,
			const size_t *in_vs,
			size_t *out_ks,
			size_t *out_vs) :
			block_num(blk_num),
			data_num(dat_num),
			digit_pos(dgt_pos),
			radix_bin_block(bin_blk),
			in_keys(in_ks),
			in_vals(in_vs),
			out_keys(out_ks),
			out_vals(out_vs) {}
		~MoveAccToBin() {}
		void operator() (size_t blk_id);
		void operator() (size_t blk_id, RadixBin& bin);
	};

	template <size_t level_num> class ChildSpawner;
}

class MSDRadixSortTask : public tbb::task
{
protected:
	template <size_t level_num>
	friend class MSDRadixSortUtils::ChildSpawner;

	using RadixBin = MSDRadixSortUtils::RadixBin;
	using RadixBinBlockMem = MSDRadixSortUtils::RadixBinBlockMem;
	using ScanAndFormBin = MSDRadixSortUtils::ScanAndFormBin;
	using MoveAccToBin = MSDRadixSortUtils::MoveAccToBin;
	template <size_t level_num>
	using ChildSpawner = MSDRadixSortUtils::ChildSpawner<level_num>;

	MSDRadixSortUtils::SortMem& sort_mem;
	const size_t start_pos, data_num;
	const size_t digit_pos;

public:
	MSDRadixSortTask(
		MSDRadixSortUtils::SortMem& _sm,
		size_t st_pos,
		size_t dat_num,
		size_t dgt_pos) :
		sort_mem(_sm),
		start_pos(st_pos),
		data_num(dat_num),
		digit_pos(dgt_pos) {}
	~MSDRadixSortTask() {}
	tbb::task* execute() override;
};

namespace MSDRadixSortUtils
{
	// div_num can be 2, 4, 16
	template <size_t level_num>
	class ChildSpawner : public tbb::task
	{
	protected:
		// bin_size = pow(4, level_num)
		static constexpr size_t bin_size = 1 << (level_num * 2);
		size_t bin[bin_size];
		const size_t max_id;
		MSDRadixSortTask &radix_sort;
	public:
		ChildSpawner(
			const size_t* _bin,
			size_t m_id,
			MSDRadixSortTask &_rs) :
			max_id(m_id),
			radix_sort(_rs)
		{ memcpy(bin, _bin, bin_size * sizeof(size_t)); }
		~ChildSpawner() {}
		tbb::task* execute()
		{
			tbb::empty_task& c = *new(allocate_continuation()) tbb::empty_task;
			c.set_ref_count(4);
			c.spawn(*new(c.allocate_child()) ChildSpawner<level_num - 1>(bin, bin[bin_size / 4], radix_sort));
			c.spawn(*new(c.allocate_child()) ChildSpawner<level_num - 1>(bin + bin_size / 4, bin[bin_size * 2 / 4], radix_sort));
			c.spawn(*new(c.allocate_child()) ChildSpawner<level_num - 1>(bin + bin_size * 2 / 4, bin[bin_size * 3 / 4], radix_sort));
			return new(c.allocate_child()) ChildSpawner<level_num - 1>(bin + bin_size * 3 / 4, max_id, radix_sort);
		}
	};

	template <>
	class ChildSpawner<1> : public tbb::task
	{
	protected:
		size_t bin[4];
		const size_t max_id;
		MSDRadixSortTask& radix_sort;
	public:
		ChildSpawner(
			const size_t* _bin,
			size_t m_id,
			MSDRadixSortTask& _rs) :
			max_id(m_id),
			radix_sort(_rs)
		{
			memcpy(bin, _bin, 4 * sizeof(size_t));
		}
		~ChildSpawner() {}
		tbb::task* execute()
		{
			if (bin[0] == max_id)
				return nullptr;
			tbb::empty_task& c = *new(allocate_continuation()) tbb::empty_task;
			if (bin[0] != bin[1])
			{
				c.increment_ref_count();
				c.spawn(*new(c.allocate_child())
					MSDRadixSortTask(
						radix_sort.sort_mem,
						radix_sort.start_pos + bin[0],
						bin[1] - bin[0],
						radix_sort.digit_pos - 1));
			}
			if (bin[1] != bin[2])
			{
				c.increment_ref_count();
				c.spawn(*new(c.allocate_child())
					MSDRadixSortTask(
						radix_sort.sort_mem,
						radix_sort.start_pos + bin[1],
						bin[2] - bin[1],
						radix_sort.digit_pos - 1));
			}
			if (bin[2] != bin[3])
			{
				c.increment_ref_count();
				c.spawn(*new(c.allocate_child())
					MSDRadixSortTask(
						radix_sort.sort_mem,
						radix_sort.start_pos + bin[2],
						bin[3] - bin[2],
						radix_sort.digit_pos - 1));
			}
			if (bin[3] != max_id)
			{
				c.increment_ref_count();
				return new(c.allocate_child())
					MSDRadixSortTask(
						radix_sort.sort_mem,
						radix_sort.start_pos + bin[3],
						max_id - bin[3],
						radix_sort.digit_pos - 1);
			}
			return nullptr;
		}
	};
}

#endif