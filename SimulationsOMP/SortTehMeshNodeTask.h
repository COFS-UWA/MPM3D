#ifndef __Sort_Teh_Mesh_Node_Task_h__
#define __Sort_Teh_Mesh_Node_Task_h__

#include "CacheAlignedMem.h"
#include "MSDRadixSortTask.h"

struct SortTehMeshNodeMem : public MSDRadixSortUtils::SortMem
{
	using RadixBinBlockMemArray = MSDRadixSortUtils::RadixBinBlockMemArray;

	struct ValidElemBlock
	{
		size_t ori_offset;
		size_t res_offset;
		size_t num;
		char _padding[MSDRadixSortUtils::cache_line_size / 2
					- sizeof(ori_offset)
					- sizeof(res_offset)
					- sizeof(num)];
	};

	ValidElemBlock *valid_elem_blocks;
	size_t *ori_elems, *res_elems;

	void init(size_t elem_num, size_t node_num, RadixBinBlockMemArray& rbbs);

protected:
	CacheAlignedMem data_mem;
};

// Assumption:
//     pcl_num > 0
//     pcl_in_elem has been sorted in accending order
class SortTehMeshNodeTask : public MSDRadixSortTask
{
protected:
	static constexpr size_t min_pcl_num_per_block = 50000;
	static constexpr size_t serial_sort_max_data_num = 500;
	static constexpr size_t insertion_sort_max_data_num = 20;

	using ValidElemBlock = SortTehMeshNodeMem::ValidElemBlock;

	struct ElemNodeIndex { size_t n1, n2, n3, n4; };

	class ScanPcl
	{
	protected:
		size_t block_num;
		size_t node_digit_pos;
		size_t pcl_num;
		const size_t *pcl_in_elems;
		const ElemNodeIndex* elem_node_ids;
		size_t *ori_elems;
		RadixBin *radix_bin_block;
		ValidElemBlock *valid_elem_blocks;
	public:
		void operator() (size_t blk_id) const;
	};

	class FormElemAndNodeArray
	{
	protected:
		size_t block_num;
		size_t node_digit_pos;
		size_t _padding_pcl_num;
		const size_t* _padding_pcl_in_elems;
		const ElemNodeIndex* elem_node_ids;
		const size_t *ori_elems;
		const RadixBin *radix_bin_block;
		ValidElemBlock* valid_elem_blocks;
		size_t* res_elems;
		size_t* out_node_has_elem;
		size_t* out_node_elem_pair;
	public:
		void operator() (size_t blk_id) const;
	};

	union
	{
		struct
		{
			size_t block_num;
			const size_t node_digit_pos;
			const size_t pcl_num;
			const size_t *const pcl_in_elems;
			const ElemNodeIndex *const elem_node_ids;
			const size_t* ori_elems;
			RadixBin* radix_bin_block;
			ValidElemBlock* valid_elem_blocks;
			size_t* res_elems;
			size_t* out_node_has_elem;
			size_t* out_node_elem_pair;
		};
		ScanPcl scan_pcl;
		FormElemAndNodeArray form_elem_and_node_array;
	};
	size_t &valid_elem_num;

public:
	SortTehMeshNodeTask(
		SortTehMeshNodeMem& _snm,
		const size_t p_num,
		const size_t *p_in_e,
		const void *en_ids,
		size_t &ve_num) :
		MSDRadixSortTask(_snm, 0, 0, _snm.digit_num - 1),
		node_digit_pos(_snm.digit_num - 1),
		pcl_num(p_num),
		pcl_in_elems(p_in_e),
		elem_node_ids((const ElemNodeIndex *)en_ids),
		valid_elem_num(ve_num) {}
	tbb::task* execute() override;
};

#endif