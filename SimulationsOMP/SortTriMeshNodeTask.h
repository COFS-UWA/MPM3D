#ifndef __Sort_Tri_Mesh_Node_Task_h__
#define __Sort_Tri_Mesh_Node_Task_h__

#include "CacheAlignedMem.h"
#include "MSDRadixSortTask.h"

struct SortTriMeshNodeMem : public MSDRadixSortUtils::SortMem
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
// pcl_num > 0
// pcl_in_elem was sorted in accending order
class SortTriMeshNodeTask : public MSDRadixSortTask
{
protected:
	static constexpr size_t min_pcl_num_per_block = 2;
	static constexpr size_t serial_sort_max_data_num = 3;
	static constexpr size_t insertion_sort_max_data_num = 1;

	using ValidElemBlock = SortTriMeshNodeMem::ValidElemBlock;
	struct ElemNodeIndex { size_t n1, n2, n3; };
	class ScanPcl
	{
	protected:
		const size_t block_num;
		const size_t pcl_num;
		const size_t*const pcl_in_elems;
		const size_t node_digit_pos;
		size_t* const ori_elems;
		const ElemNodeIndex *const elem_node_ids;
		RadixBin *const radix_bin_block;
		ValidElemBlock *const valid_elem_blocks;
	public:
		ScanPcl(
			size_t blk_num,
			size_t p_num,
			const size_t* p_in_e,
			size_t n_dgt_pos,
			size_t* ori_e,
			const ElemNodeIndex* en_id,
			RadixBin *rd_bin_blk,
			ValidElemBlock *ve_blks) :
			block_num(blk_num),
			pcl_num(p_num),
			pcl_in_elems(p_in_e),
			node_digit_pos(n_dgt_pos),
			ori_elems(ori_e),
			elem_node_ids(en_id),
			radix_bin_block(rd_bin_blk),
			valid_elem_blocks(ve_blks) {}
		void operator() (size_t blk_id);
	};
	class FormElemAndNodeArray
	{
	protected:
		const size_t block_num;
		const size_t node_digit_pos;
		const size_t *const ori_elems;
		size_t *const res_elems;
		const ElemNodeIndex* const elem_node_ids;
		const RadixBin* const radix_bin_block;
		ValidElemBlock* const valid_elem_blocks;
		size_t* out_node_has_elem;
		size_t* out_node_elem_pair;
	public:
		FormElemAndNodeArray(
			size_t blk_num,
			size_t n_dgt_pos,
			ValidElemBlock* ve_blks,
			const size_t* ori_e,
			size_t *res_e,
			const ElemNodeIndex* en_id,
			const RadixBin* rd_bin_blk,
			size_t *out_n_has_e,
			size_t *out_n_e_pair) :
			block_num(blk_num),
			node_digit_pos(n_dgt_pos),
			valid_elem_blocks(ve_blks),
			ori_elems(ori_e),
			res_elems(res_e),
			elem_node_ids(en_id),
			radix_bin_block(rd_bin_blk),
			out_node_has_elem(out_n_has_e),
			out_node_elem_pair(out_n_e_pair) {}
		void operator() (size_t blk_id);
	};

#ifdef _DEBUG
	const size_t elem_num;
#endif
	const size_t pcl_num;
	const size_t* const pcl_in_elems;
	const ElemNodeIndex *const elem_node_ids;
	size_t &valid_elem_num;

public:
	SortTriMeshNodeTask(
		SortTriMeshNodeMem& _snm,
#ifdef _DEBUG
		const size_t e_num,
#endif
		const size_t p_num,
		const size_t *p_in_e,
		const void *en_ids,
		size_t &ve_num) :
		MSDRadixSortTask(_snm, 0, 0, _snm.digit_num - 1),
#ifdef _DEBUG
		elem_num(e_num),
#endif
		pcl_num(p_num),
		pcl_in_elems(p_in_e),
		elem_node_ids((const ElemNodeIndex *)en_ids),
		valid_elem_num(ve_num) {}
	~SortTriMeshNodeTask() {}
	tbb::task* execute() override;
};

#endif