#ifndef __Step_T2D_CHM_mt_h__
#define __Step_T2D_CHM_mt_h__

#include "CacheAlignedMem.h"
#include "Step_OMP.h"
#include "Model_T2D_CHM_mt.h"

class Model_T2D_CHM_mt;
class Step_T2D_CHM_mt;
//namespace Model_T2D_CHM_mt_hdf5_utilities
//{
//	int load_me_mt_model_from_hdf5_file(Model_T2D_CHM_mt& md, Step_T2D_CHM_mt& step, const char* hdf5_name, const char* th_name, size_t frame_id);
//}

int substep_func_omp_T2D_CHM_mt(void* _self, size_t my_th_id,
	double dt, double cur_time, size_t substp_id);

class Step_T2D_CHM_mt : public Step_OMP
{
protected:
	typedef Model_T2D_CHM_mt::PclBodyForce PclBodyForce;
	typedef Model_T2D_CHM_mt::PclTraction PclTraction;
	typedef Model_T2D_CHM_mt::PclPos PclPos;
	typedef Model_T2D_CHM_mt::PclSortedVarArray PclSortedVarArray;
	typedef Model_T2D_CHM_mt::PclDisp PclDisp;
	typedef Model_T2D_CHM_mt::PclV PclV;
	typedef Model_T2D_CHM_mt::PclShapeFunc PclShapeFunc;
	typedef Model_T2D_CHM_mt::PclStress PclStress;
	typedef Model_T2D_CHM_mt::ElemNodeIndex ElemNodeIndex;
	typedef Model_T2D_CHM_mt::ElemShapeFuncAB ElemShapeFuncAB;
	typedef Model_T2D_CHM_mt::ElemShapeFuncC ElemShapeFuncC;
	typedef Model_T2D_CHM_mt::ElemStrainInc ElemStrainInc;
	typedef Model_T2D_CHM_mt::ElemStress ElemStress;
	typedef Model_T2D_CHM_mt::ElemNodeVM ElemNodeVM;
	typedef Model_T2D_CHM_mt::ElemNodeForce ElemNodeForce;
	typedef Model_T2D_CHM_mt::NodeA NodeA;
	typedef Model_T2D_CHM_mt::NodeV NodeV;
	typedef Model_T2D_CHM_mt::NodeHasVBC NodeHasVBC;
	typedef Model_T2D_CHM_mt::NodePos NodePos;

	size_t elem_num;
	size_t node_num;

	double* pcl_m;
	PclBodyForce* pcl_bf;
	PclTraction* pcl_t;
	PclPos* pcl_pos;
	double* pcl_vol;
	MatModel::MaterialModel** pcl_mat_model;

	PclSortedVarArray pcl_sorted_var_array[2];

	ElemNodeIndex* elem_node_id;
	double* elem_area;
	ElemShapeFuncAB* elem_sf_ab;
	ElemShapeFuncC* elem_sf_c;

	double* elem_density;
	double* elem_pcl_m;
	double* elem_pcl_vol;
	ElemStrainInc* elem_de;
	ElemStress* elem_stress;
	double* elem_m_de_vol;

	ElemNodeVM* elem_node_vm;
	ElemNodeForce* elem_node_force;

	size_t* elem_id_array;
	size_t* node_elem_id_array;
	size_t* node_elem_list;
	NodeA *node_a;
	NodeV *node_v;
	NodeHasVBC* node_has_vbc;
	double *node_am; // node_num
	double *node_de_vol; // node_num

protected:
	// task division
	CacheAlignedMem task_range_mem;
	union PclRange
	{
		size_t id;
		char padding[Cache_Alignment];
	};
	PclRange* pcl_range;
	size_t *elem_range;
	size_t *node_range;
	size_t* node_elem_range;

	// radix sort
	CacheAlignedMem elem_bin_mem;
	size_t* elem_count_bin;
	size_t* elem_sum_bin;
	
	CacheAlignedMem radix_sort_var_mem;
	size_t *new_to_prev_pcl_maps[2];
	size_t *pcl_in_elem_arrays[2];

	size_t pcl_sorted_var_id;
	size_t radix_sort_var_id;
	size_t pcl_num;
	size_t new_pcl_num;
	
	double K_cont;
	double rr_fx_cont, rr_fy_cont, rr_m_cont;
	struct ContPos { double x, y; };
	CacheAlignedMem contact_mem;
	size_t* contact_state;
	ContPos *contact_pos;

	int apply_rigid_circle_avg(size_t my_th_id, double dt,
		size_t* pcl_in_elem_array, PclSortedVarArray& cur_pscv,
		RigidCircleForce& rr_force);

public:
	int init_calculation() override;
	friend int substep_func_omp_T2D_CHM_mt(void* _self,
		size_t my_th_id, double dt, double cur_time, size_t substp_id);
	int finalize_calculation() override;

public:
	typedef Model_T2D_CHM_mt::PclBodyForce PclBodyForce;
	typedef Model_T2D_CHM_mt::PclTraction PclTraction;
	typedef Model_T2D_CHM_mt::PclPos PclPos;

	typedef Model_T2D_CHM_mt::PclSortedVarArray PclSortedVarArray;
	typedef Model_T2D_CHM_mt::PclDisp PclDisp;
	typedef Model_T2D_CHM_mt::PclV PclV;
	typedef Model_T2D_CHM_mt::PclShapeFunc PclShapeFunc;
	typedef Model_T2D_CHM_mt::PclStress PclStress;

	typedef Model_T2D_CHM_mt::ElemNodeIndex ElemNodeIndex;
	typedef Model_T2D_CHM_mt::ElemShapeFuncAB ElemShapeFuncAB;
	typedef Model_T2D_CHM_mt::ElemShapeFuncC ElemShapeFuncC;

	typedef Model_T2D_CHM_mt::ElemStrainInc ElemStrainInc;
	typedef Model_T2D_CHM_mt::ElemStress ElemStress;

	typedef Model_T2D_CHM_mt::ElemNodeVM ElemNodeVM;
	typedef Model_T2D_CHM_mt::ElemNodeForce ElemNodeForce;

	Step_T2D_CHM_mt(const char* _name);
	~Step_T2D_CHM_mt();

	inline size_t get_pcl_num() const noexcept { return pcl_num; }
	inline size_t get_pcl_sorted_var_id() const noexcept { return pcl_sorted_var_id; }
	inline const size_t *get_new_to_prev_pcl_map() const noexcept { return new_to_prev_pcl_maps[radix_sort_var_id]; }
	inline double get_rr_fx_contact() const noexcept { return rr_fx_cont; }
	inline double get_rr_fy_contact() const noexcept { return rr_fy_cont; }
	inline double get_rr_m_contact() const noexcept { return rr_m_cont; }

	//friend int Model_T2D_CHM_mt_hdf5_utilities::load_me_mt_model_from_hdf5_file(Model_T2D_CHM_mt& md, Step_T2D_CHM_mt& step, const char* hdf5_name, const char* th_name, size_t frame_id);
};

#endif