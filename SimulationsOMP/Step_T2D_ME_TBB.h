#ifndef __Step_T2D_ME_TBB_h__
#define __Step_T2D_ME_TBB_h__

#include "CacheAlignedMem.h"
#include "Model_T2D_ME_mt.h"
#include "Step_TBB.h"

class Step_T2D_ME_TBB;
namespace Model_T2D_ME_mt_hdf5_utilities
{
	int load_me_mt_model_from_hdf5_file(Model_T2D_ME_mt& md, Step_T2D_ME_TBB& step, const char* hdf5_name, const char* th_name, size_t frame_id);
}

int cal_substep_func_T2D_ME_TBB(void* _self);

class Step_T2D_ME_TBB : public Step_TBB
{
	friend int Model_T2D_ME_mt_hdf5_utilities::load_me_mt_model_from_hdf5_file(Model_T2D_ME_mt& md, Step_T2D_ME_TBB& step, const char* hdf5_name, const char* th_name, size_t frame_id);
protected:
	typedef Model_T2D_ME_mt::SortedPclVarArrays SortedPclVarArrays;
	typedef Model_T2D_ME_mt::Force Force;
	typedef Model_T2D_ME_mt::Acceleration Acceleration;
	typedef Model_T2D_ME_mt::Velocity Velocity;
	typedef Model_T2D_ME_mt::Displacement Displacement;
	typedef Model_T2D_ME_mt::Position Position;
	typedef Model_T2D_ME_mt::ShapeFunc ShapeFunc;
	typedef Model_T2D_ME_mt::ShapeFuncAB ShapeFuncAB;
	typedef Model_T2D_ME_mt::ShapeFuncC ShapeFuncC;
	typedef Model_T2D_ME_mt::Stress Stress;
	typedef Model_T2D_ME_mt::Strain Strain;
	typedef Model_T2D_ME_mt::StrainInc StrainInc;
	typedef Model_T2D_ME_mt::ElemNodeIndex ElemNodeIndex;
	typedef Model_T2D_ME_mt::ElemNodeVM ElemNodeVM;
	typedef Model_T2D_ME_mt::NodeHasVBC NodeHasVBC;

	size_t elem_num;
	size_t node_num;

	double* pcl_m;
	Force *pcl_bf;
	Force *pcl_t;
	Position *pcl_pos;
	double* pcl_vol;
	MatModel::MaterialModel** pcl_mat_model;

	SortedPclVarArrays pcl_sorted_var_arrays[2];

	ElemNodeIndex* elem_node_id;
	ShapeFuncAB* elem_dN_ab;
	ShapeFuncC* elem_dN_c;
	double* elem_area;

	double* elem_pcl_m;
	double* elem_density;
	StrainInc* elem_de;
	double* elem_m_de_vol;

	ElemNodeVM* elem_node_vm;
	Force* elem_node_force;

	Acceleration *node_a;
	Velocity *node_v;
	NodeHasVBC* node_has_vbc;
	double* node_am;
	double* node_de_vol;

	//int apply_rigid_rect(
	//	size_t p_id0, size_t p_id1,
	//	size_t* pcl_in_elem,
	//	SortedPclVarArrays& cur_spva,
	//	Force2D& rc_cf,
	//	size_t substp_id,
	//	ThreadData& thd) noexcept;

public:
	int init_calculation() override;
	friend int cal_substep_func_T2D_ME_TBB(void* _self);
	int finalize_calculation() override;

public:
	Step_T2D_ME_TBB(const char* _name);
	~Step_T2D_ME_TBB();

	//inline size_t get_pcl_num() const noexcept { return prev_valid_pcl_num; }
	//inline size_t get_sorted_pcl_var_id() const noexcept { return thread_datas[0].sorted_pcl_var_id; }
	//inline size_t* get_pcl_in_elem() const noexcept { return pcl_in_elems[thread_datas[0].sorted_pcl_in_elem_id]; }
};

#endif