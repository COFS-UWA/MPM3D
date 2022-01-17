#ifndef __Step_T3D_CHM_ud_mt_subiter_h__
#define __Step_T3D_CHM_ud_mt_subiter_h__

#include "CacheAlignedMem.h"
#include "Step_OMP.h"
#include "Model_T3D_CHM_mt.h"

class Step_T3D_CHM_ud_mt_subiter;
namespace Model_T3D_CHM_mt_hdf5_utilities
{
	int load_model_from_hdf5_file(Model_T3D_CHM_mt& md, Step_T3D_CHM_ud_mt_subiter& step,
		const char* hdf5_name, const char* th_name, size_t frame_id);
}

int substep_func_omp_T3D_CHM_ud_mt_subiter(void* _self, size_t my_th_id,
	double dt, double cur_time, size_t substp_id);

class Step_T3D_CHM_ud_mt_subiter : public Step_OMP
{
protected:
	typedef Model_T3D_CHM_mt::ShapeFunc ShapeFunc;
	typedef Model_T3D_CHM_mt::DShapeFuncABC DShapeFuncABC;
	typedef Model_T3D_CHM_mt::DShapeFuncD DShapeFuncD;
	typedef Model_T3D_CHM_mt::Force Force;
	typedef Model_T3D_CHM_mt::Position Position;
	typedef Model_T3D_CHM_mt::Displacement Displacement;
	typedef Model_T3D_CHM_mt::Velocity Velocity;
	typedef Model_T3D_CHM_mt::Acceleration Acceleration;
	typedef Model_T3D_CHM_mt::Stress Stress;
	typedef Model_T3D_CHM_mt::Strain Strain;
	typedef Model_T3D_CHM_mt::StrainInc StrainInc;
	typedef Model_T3D_CHM_mt::SortedPclVarArrays SortedPclVarArrays;
	typedef Model_T3D_CHM_mt::ElemNodeIndex ElemNodeIndex;
	typedef Model_T3D_CHM_mt::ElemNodeVM ElemNodeVM;
	typedef Model_T3D_CHM_mt::NodeHasVBC NodeHasVBC;
	typedef Model_T3D_CHM_mt::NodeVBCVec NodeVBCVec;

	double* pcl_m_s; // ori_pcl_num
	double* pcl_density_s; // ori_pcl_num
	double* pcl_vol_s; // ori_pcl_num
	Force* pcl_bf_s; // ori_pcl_num
	Force* pcl_bf_f; // ori_pcl_num
	Force* pcl_t; // ori_pcl_num
	Position* pcl_pos; // ori_pcl_num
	MatModel::MaterialModel** pcl_mat_model;
	size_t* pcl_mat_model_copy_offset; // ori_pcl_num

	double *pcl_vol; // ori_pcl_num
	Stress* pcl_prev_stress; // ori_pcl_num
	StrainInc *pcl_destrain; // ori_pcl_num
	StrainInc *pcl_dpstrain; // ori_pcl_num

	SortedPclVarArrays sorted_pcl_var_arrays[2];

	ElemNodeIndex* elem_node_id;
	DShapeFuncABC* elem_N_abc;
	DShapeFuncD* elem_N_d;
	double* elem_vol;

	double* elem_pcl_int_vol; // elem_num
	double* elem_pcl_m; // elem_num
	double* elem_pcl_n; // elem_num
	double* elem_density_f; // elem_num
	StrainInc* elem_de; // elem_num
	double* elem_p; // elem_num
	double* elem_m_de_vol; // elem_num
	double* elem_m_pde_vol; // elem_num

	// element-node data
	ElemNodeVM* elem_node_vm; // elem_num * 4
	Force* elem_node_f_ext; // elem_num * 4
	Force* elem_node_f_int; // elem_num * 4

	double* node_am; // node_num, mass for cal acceleration
	Acceleration* node_a; // node_num
	Velocity* node_v; // node_num
	Displacement* node_du; // node_num, nodal disp inc
	Velocity* node_vn; // node_num, previous nodal velocity
	Velocity* node_pv; // node_num, pseudo velocity
	Displacement* node_pdu; // node_num, pseudo displacement
	double* node_de_vol; // node_num, strain enhancement
	double* node_pde_vol; // node_num, strain enhancement
	NodeHasVBC* node_has_vbc; // node_num
	NodeVBCVec* node_vbc_vec; // node_num
	
	size_t elem_num, node_num;
	double Kf, u_cav, m_cav, f_cav_end;
	double u_cav_off, u_div_u_cav_lim;

#ifdef _DEBUG
	size_t prev_valid_pcl_num_tmp;
#endif
	size_t prev_valid_pcl_num, valid_pcl_num;
	size_t valid_elem_num;
	Force3D cf_tmp;

	size_t* elem_count_bin;
	size_t* elem_sum_bin;

	size_t* prev_pcl_ids[2]; // (pcl_num + 2) * 2
	size_t* pcl_in_elems[2]; // pcl_num * 2
	size_t* node_has_elems[2]; // (elem_num * 4 + 2) * 2
	size_t* node_elem_pairs[2]; // elem_num * 4 * 2
	size_t* valid_elem_ids; // elem_num

	union ThreadData
	{
		struct
		{
			size_t sorted_pcl_var_id;
			size_t sorted_pcl_in_elem_id;
			double max_pcl_vol;
		};
		char padding[Cache_Alignment * 2];
	};
	ThreadData* thread_datas;

	// tmp buf to store material model
	char *mat_model_copy;

	CacheAlignedMem thread_mem;
	CacheAlignedMem cal_mem;
	CacheAlignedMem mat_model_copy_mem;
	
	// rigid bodies
	RigidCylinder* prcy;
	RigidObjectByT3DMesh* prm;
	ContactModel3D* pcm_s, * pcm_f;
	
	int apply_rigid_cylinder(
		size_t p_id0, size_t p_id1,
		const size_t* pcl_in_elem,
		const SortedPclVarArrays& cur_spva,
		Force3D& rc_cf,
		size_t substp_id,
		ThreadData& thd) noexcept;

	int apply_t3d_rigid_mesh(
		size_t p_id0, size_t p_id1,
		const size_t* pcl_in_elem,
		const SortedPclVarArrays& cur_spva,
		Force3D& rc_cf,
		size_t substp_id,
		ThreadData& thd) noexcept;
	
	ParticleVariablesGetter pv_place_holder;
	
	// subiteration parameter
	size_t max_subiter_num;
	double mass_factor;
	double converge_e_kin_ratio;
	double pdt; // pseudo dt
	// subiteration state
	size_t cal_status, subiter_index;
	double cur_e_kin, prev_e_kin, max_e_kin;

	int Step_T3D_CHM_ud_mt_subiter::subiteration(
		size_t my_th_id,
		size_t p_id0, size_t p_id1,
		size_t ve_id0, size_t ve_id1,
		const size_t* my_valid_elem_ids, size_t my_valid_elem_num,
		SortedPclVarArrays& spva0,
		const size_t* pcl_in_elem0,
		const size_t* node_has_elem0, const size_t* node_elem_pair0);

public:
	Step_T3D_CHM_ud_mt_subiter(const char* _name);
	~Step_T3D_CHM_ud_mt_subiter();
	
	int init_calculation() override;
	friend int substep_func_omp_T3D_CHM_ud_mt_subiter(void* _self,
		size_t my_th_id, double dt, double cur_time, size_t substp_id);
	int finalize_calculation() override;

	void set_max_subiter_num(size_t stp_num) noexcept { max_subiter_num = stp_num; }
	void set_mass_factor(double m_fac) noexcept { mass_factor = m_fac; }
	void set_converge_e_kin_ratio(double con_ratio) noexcept { converge_e_kin_ratio = con_ratio; }
	void set_pdt(double _pdt) noexcept { pdt = _pdt; }

	inline size_t get_pcl_num() const noexcept { return prev_valid_pcl_num; }
	inline size_t get_sorted_pcl_var_id() const noexcept { return thread_datas[0].sorted_pcl_var_id; }
	inline size_t* get_prev_pcl_id() const noexcept { return prev_pcl_ids[thread_datas[0].sorted_pcl_in_elem_id]; }
	inline size_t* get_pcl_in_elem() const noexcept { return pcl_in_elems[thread_datas[0].sorted_pcl_in_elem_id]; }
	inline size_t get_subiter_index() const noexcept { return subiter_index; }
	inline double get_cur_e_kin() const noexcept { return cur_e_kin; }
	inline double get_max_e_kin() const noexcept { return max_e_kin; }

	friend int Model_T3D_CHM_mt_hdf5_utilities::load_model_from_hdf5_file(
		Model_T3D_CHM_mt &md, Step_T3D_CHM_ud_mt_subiter &step, const char *hdf5_name,
		const char* th_name, size_t frame_id);
};

#endif