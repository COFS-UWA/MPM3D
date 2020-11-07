#ifndef __Step_R2D_ME_mt_h__
#define __Step_R2D_ME_mt_h__

#include "CacheAlignedMem.h"
#include "Step_OMP.h"
#include "Model_R2D_ME_mt.h"

class Model_R2D_ME_mt;
class Step_R2D_ME_mt;
namespace Model_R2D_ME_mt_hdf5_utilities
{
	int load_me_mt_model_from_hdf5_file(Model_R2D_ME_mt& md, Step_R2D_ME_mt& step, const char* hdf5_name, const char* th_name, size_t frame_id);
}

int substep_func_omp_R2D_ME_mt(void* _self, size_t my_th_id,
	double dt, double cur_time, size_t substp_id);

class Step_R2D_ME_mt : public Step_OMP
{
	friend int Model_R2D_ME_mt_hdf5_utilities::load_me_mt_model_from_hdf5_file(Model_R2D_ME_mt& md, Step_R2D_ME_mt& step, const char* hdf5_name, const char* th_name, size_t frame_id);
protected:
	typedef Model_R2D_ME_mt::ShapeFunc ShapeFunc;
	typedef Model_R2D_ME_mt::DShapeFunc DShapeFunc;
	typedef Model_R2D_ME_mt::Stress Stress;
	typedef Model_R2D_ME_mt::StrainInc StrainInc;
	typedef Model_R2D_ME_mt::Stress Stress;
	typedef Model_R2D_ME_mt::PclBodyForce PclBodyForce;
	typedef Model_R2D_ME_mt::PclTraction PclTraction;
	typedef Model_R2D_ME_mt::PclPos PclPos;
	typedef Model_R2D_ME_mt::PclDisp PclDisp;
	typedef Model_R2D_ME_mt::PclV PclV;
	typedef Model_R2D_ME_mt::PclSortedVarArray PclSortedVarArray;
	typedef Model_R2D_ME_mt::ElemNodeVM ElemNodeVM;
	typedef Model_R2D_ME_mt::ElemNodeForce ElemNodeForce;
	typedef Model_R2D_ME_mt::NodeA NodeA;
	typedef Model_R2D_ME_mt::NodeV NodeV;
	typedef Model_R2D_ME_mt::NodeHasVBC NodeHasVBC;
	typedef Model_R2D_ME_mt::ContPos ContPos;

	size_t actual_elem_x_num, actual_elem_num;
	size_t node_x_num, node_num;
	double mh_xl, mh_xu, mh_yl, mh_yu;
	double inv_elem_hx, inv_elem_hy, elem_area;
	double dxi_dx, deta_dy;
	DShapeFunc elem_dN;

	double K_cont;
	double rr_fx_cont, rr_fy_cont, rr_m_cont;

	size_t pcl_sorted_var_id;
	size_t radix_sort_var_id;
	size_t pcl_num;
	
	double* pcl_m;
	PclBodyForce* pcl_bf;
	PclTraction* pcl_t;
	PclPos* pcl_pos;
	double* pcl_vol;
	ShapeFunc* pcl_N;
	DShapeFunc* pcl_dN;
	MatModel::MaterialModel** pcl_mat_model;
	PclSortedVarArray pcl_sorted_var_array[2];

	double* elem_density;
	double* elem_pcl_m;
	size_t *elem_substp_id;

	ElemNodeVM* elem_node_vm;
	ElemNodeForce* elem_node_force;

	NodeA *node_a;
	NodeV *node_v;
	NodeHasVBC* node_has_vbc;

	// contact
	size_t* contact_substep_id;
	ContPos* contact_pos;

	// task division
	union PclRange
	{
		size_t id;
		char padding[Cache_Alignment];
	};
	PclRange* pcl_range;
	size_t *elem_range;
	size_t *node_range;

	// radix sort
	size_t* elem_count_bin;
	size_t* elem_sum_bin;
	size_t* pcl_in_elem_arrays[2];

	// memory
	CacheAlignedMem task_range_mem;
	CacheAlignedMem elem_bin_mem;
	CacheAlignedMem radix_sort_var_mem;

	int apply_rigid_rect_avg(size_t p_id0, size_t p_id1,
		size_t *pcl_in_elem, PclSortedVarArray& cur_pscv,
		RigidRectForce& rr_force);

	int init_calculation() override;
	friend int substep_func_omp_R2D_ME_mt(void* _self,
		size_t my_th_id, double dt, double cur_time, size_t substp_id);
	int finalize_calculation() override;

public:
	Step_R2D_ME_mt(const char* _name);
	~Step_R2D_ME_mt();

	inline size_t get_pcl_num() const noexcept { return pcl_num; }
	inline size_t get_pcl_sorted_var_id() const noexcept { return pcl_sorted_var_id; }
	inline double get_rr_fx_contact() const noexcept { return rr_fx_cont; }
	inline double get_rr_fy_contact() const noexcept { return rr_fy_cont; }
	inline double get_rr_m_contact() const noexcept { return rr_m_cont; }
};

#endif