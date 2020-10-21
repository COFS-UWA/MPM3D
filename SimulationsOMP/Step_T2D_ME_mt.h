#ifndef __Step_T2D_ME_mt_h__
#define __Step_T2D_ME_mt_h__

#include "CacheAlignedMem.h"
#include "Step_OMP.h"
#include "Model_T2D_ME_mt.h"

class Model_T2D_ME_mt;
class Step_T2D_ME_mt;
namespace Model_T2D_ME_mt_hdf5_utilities
{
	int load_me_mt_model_from_hdf5_file(Model_T2D_ME_mt& md, Step_T2D_ME_mt& step, const char* hdf5_name, const char* th_name, size_t frame_id);
}

int substep_func_omp_T2D_ME_mt(void* _self, uint32_t my_th_id,
	float dt, float cur_time, uint32_t substp_id);

class Step_T2D_ME_mt : public Step_OMP
{
	friend int Model_T2D_ME_mt_hdf5_utilities::load_me_mt_model_from_hdf5_file(Model_T2D_ME_mt& md, Step_T2D_ME_mt& step, const char* hdf5_name, const char* th_name, size_t frame_id);
protected:
	typedef Model_T2D_ME_mt::PclBodyForce PclBodyForce;
	typedef Model_T2D_ME_mt::PclTraction PclTraction;
	typedef Model_T2D_ME_mt::PclPos PclPos;
	typedef Model_T2D_ME_mt::PclSortedVarArray PclSortedVarArray;
	typedef Model_T2D_ME_mt::PclDisp PclDisp;
	typedef Model_T2D_ME_mt::PclV PclV;
	typedef Model_T2D_ME_mt::PclShapeFunc PclShapeFunc;
	typedef Model_T2D_ME_mt::PclStress PclStress;
	typedef Model_T2D_ME_mt::ElemNodeIndex ElemNodeIndex;
	typedef Model_T2D_ME_mt::ElemShapeFuncAB ElemShapeFuncAB;
	typedef Model_T2D_ME_mt::ElemShapeFuncC ElemShapeFuncC;
	typedef Model_T2D_ME_mt::ElemStrainInc ElemStrainInc;
	typedef Model_T2D_ME_mt::ElemStress ElemStress;
	typedef Model_T2D_ME_mt::ElemNodeVM ElemNodeVM;
	typedef Model_T2D_ME_mt::ElemNodeForce ElemNodeForce;
	typedef Model_T2D_ME_mt::NodeA NodeA;
	typedef Model_T2D_ME_mt::NodeV NodeV;
	typedef Model_T2D_ME_mt::NodeHasVBC NodeHasVBC;
	typedef Model_T2D_ME_mt::NodePos NodePos;

	uint32_t pcl_num;
	uint32_t elem_num;
	uint32_t node_num;

	float* pcl_m;
	PclBodyForce* pcl_bf;
	PclTraction* pcl_t;
	PclPos* pcl_pos;
	MatModel::MaterialModel** pcl_mat_model;

	uint32_t pcl_sorted_var_id;
	PclSortedVarArray pcl_sorted_var_array[2];

	ElemNodeIndex* elem_node_id;
	float* elem_area;
	ElemShapeFuncAB* elem_sf_ab;
	ElemShapeFuncC* elem_sf_c;

	float* elem_density;
	float* elem_pcl_m;
	float* elem_pcl_vol;
	ElemStrainInc* elem_de;
	ElemStress* elem_stress;
	float* elem_m_de_vol;

	ElemNodeVM* elem_node_vm;
	ElemNodeForce* elem_node_force;

	uint32_t* elem_id_array;
	uint32_t* node_elem_id_array;
	uint32_t* node_elem_list;
	NodeA *node_a;
	NodeV *node_v;
	NodeHasVBC* node_has_vbc;
	float *node_am; // node_num
	float *node_de_vol; // node_num

protected:
	// task division
	CacheAlignedMem task_range_mem;
	union PclRange
	{
		uint32_t id;
		char padding[Cache_Alignment];
	};
	PclRange* pcl_range;
	uint32_t *elem_range;
	uint32_t *node_range;
	uint32_t* node_elem_range;

	// radix sort
	CacheAlignedMem sort_var_mem;
	uint32_t* new_to_ori_pcl_map0;
	uint32_t* new_to_ori_pcl_map1;
	uint32_t* pcl_in_elem_array0;
	uint32_t* pcl_in_elem_array1;
	CacheAlignedMem elem_bin_mem;
	uint32_t *elem_count_bin;
	uint32_t* elem_sum_bin;

	uint32_t new_pcl_num;

public:
	int init_calculation() override;
	friend int substep_func_omp_T2D_ME_mt(void* _self,
		uint32_t my_th_id, float dt, float cur_time, uint32_t substp_id);
	int finalize_calculation() override;

public:
	typedef Model_T2D_ME_mt::PclBodyForce PclBodyForce;
	typedef Model_T2D_ME_mt::PclTraction PclTraction;
	typedef Model_T2D_ME_mt::PclPos PclPos;

	typedef Model_T2D_ME_mt::PclSortedVarArray PclSortedVarArray;
	typedef Model_T2D_ME_mt::PclDisp PclDisp;
	typedef Model_T2D_ME_mt::PclV PclV;
	typedef Model_T2D_ME_mt::PclShapeFunc PclShapeFunc;
	typedef Model_T2D_ME_mt::PclStress PclStress;

	typedef Model_T2D_ME_mt::ElemNodeIndex ElemNodeIndex;
	typedef Model_T2D_ME_mt::ElemShapeFuncAB ElemShapeFuncAB;
	typedef Model_T2D_ME_mt::ElemShapeFuncC ElemShapeFuncC;

	typedef Model_T2D_ME_mt::ElemStrainInc ElemStrainInc;
	typedef Model_T2D_ME_mt::ElemStress ElemStress;

	typedef Model_T2D_ME_mt::ElemNodeVM ElemNodeVM;
	typedef Model_T2D_ME_mt::ElemNodeForce ElemNodeForce;

	Step_T2D_ME_mt(const char* _name);
	~Step_T2D_ME_mt();

	inline uint32_t get_pcl_num() const noexcept { return pcl_num; }
	inline uint32_t get_sorted_var_id() const noexcept { return pcl_sorted_var_id; }
	inline uint32_t* get_new_to_ori_pcl_map() const noexcept { return new_to_ori_pcl_map0; }
};

#endif