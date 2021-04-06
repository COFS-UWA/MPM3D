#ifndef __MSD_Radix_Sort_Task_h__
#define __MSD_Radix_Sort_Task_h__

#include "CacheAlignedMem.h"
#include "MSDRadixSortUtils.h"

namespace MSDRadixSortUtils
{
	class ScanAndFormBin
	{
	protected:
		size_t data_num;
		size_t digit_pos;
		size_t block_num;
		RadixBin *radix_bin_block;
		const size_t *in_keys;
	public:
		void operator() (size_t blk_id) const;
	};

	class MoveAccToBin
	{
	protected:
		size_t data_num;
		size_t digit_pos;
		size_t block_num;
		const RadixBin *radix_bin_block;
		const size_t *in_keys;
		const size_t *in_vals;
		size_t *out_keys;
		size_t *out_vals;
	public:
		void operator() (size_t blk_id) const;
		void operator() (RadixBin& bin) const;
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
	const size_t start_pos;
	
	union
	{
		struct
		{
			const size_t data_num;
			const size_t digit_pos;
			size_t block_num;
			RadixBin* radix_bin_block;
			const size_t* in_keys;
			const size_t* in_vals;
			size_t* out_keys;
			size_t* out_vals;
		};
		ScanAndFormBin scan_and_form_bin;
		MoveAccToBin move_acc_to_bin;
	};

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
		const size_t start_pos;
		const size_t child_digit_pos;
		SortMem& sort_mem;
	public:
		ChildSpawner(
			const size_t* _bin,
			size_t m_id,
			size_t st_pos,
			size_t cld_dgt_pos,
			SortMem& _sm) :
			max_id(m_id),
			start_pos(st_pos),
			child_digit_pos(cld_dgt_pos),
			sort_mem(_sm)
		{ memcpy(bin, _bin, bin_size * sizeof(size_t)); }
		tbb::task* execute()
		{
			tbb::empty_task& c =
				*new(allocate_continuation()) tbb::empty_task;
			c.set_ref_count(4);
			c.spawn(*new(c.allocate_child())
				ChildSpawner<level_num - 1>(
					bin,
					bin[bin_size / 4],
					start_pos,
					child_digit_pos,
					sort_mem));
			c.spawn(*new(c.allocate_child())
				ChildSpawner<level_num - 1>(
					bin + bin_size / 4,
					bin[bin_size * 2 / 4],
					start_pos,
					child_digit_pos,
					sort_mem));
			c.spawn(*new(c.allocate_child())
				ChildSpawner<level_num - 1>(
					bin + bin_size * 2 / 4,
					bin[bin_size * 3 / 4],
					start_pos,
					child_digit_pos,
					sort_mem));
			return new(c.allocate_child())
				ChildSpawner<level_num - 1>(
					bin + bin_size * 3 / 4,
					max_id,
					start_pos,
					child_digit_pos,
					sort_mem);
		}
	};

	template <>
	class ChildSpawner<1> : public tbb::task
	{
	protected:
		size_t bin[4];
		const size_t max_id;
		const size_t start_pos;
		const size_t child_digit_pos;
		SortMem& sort_mem;
	public:
		ChildSpawner(
			const size_t* _bin,
			size_t m_id,
			size_t st_pos,
			size_t cld_dgt_pos,
			SortMem& _sm) :
			max_id(m_id),
			start_pos(st_pos),
			child_digit_pos(cld_dgt_pos),
			sort_mem(_sm)
		{ memcpy(bin, _bin, 4 * sizeof(size_t)); }
		tbb::task* execute()
		{
			if (bin[0] == max_id)
				return nullptr;
			tbb::empty_task& c = *new(allocate_continuation()) tbb::empty_task;
			c.set_ref_count(4);
			if (bin[0] != bin[1])
				c.spawn(*new(c.allocate_child())
					MSDRadixSortTask(
						sort_mem,
						start_pos + bin[0],
						bin[1] - bin[0],
						child_digit_pos));
			else
				c.decrement_ref_count();
			if (bin[1] != bin[2])
				c.spawn(*new(c.allocate_child())
					MSDRadixSortTask(
						sort_mem,
						start_pos + bin[1],
						bin[2] - bin[1],
						child_digit_pos));
			else
				c.decrement_ref_count();
			if (bin[2] != bin[3])
				c.spawn(*new(c.allocate_child())
					MSDRadixSortTask(
						sort_mem,
						start_pos + bin[2],
						bin[3] - bin[2],
						child_digit_pos));
			else
				c.decrement_ref_count();
			if (bin[3] != max_id)
				return new(c.allocate_child())
					MSDRadixSortTask(
						sort_mem,
						start_pos + bin[3],
						max_id - bin[3],
						child_digit_pos);
			c.decrement_ref_count();
			return nullptr;
		}
	};
}

#endif