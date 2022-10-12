#ifndef __Step_T2D_ME_mt_Geo_h__
#define __Step_T2D_ME_mt_Geo_h__

#include "CacheAlignedMem.h"
#include "Step_OMP.h"
#include "Model_T2D_ME_mt.h"
#include "RigidObject/Force2D.h"

class Model_T2D_ME_mt;
class Step_T2D_ME_mt_Geo;
namespace Model_T2D_ME_mt_hdf5_utilities
{
	int load_model_from_hdf5_file(Model_T2D_ME_mt& md, const char* hdf5_name);
	int output_pcl_data_to_hdf5_file(Model_T2D_ME_mt& md, Step_T2D_ME_mt_Geo& stp, ResultFile_hdf5& rf, hid_t grp_id);
}

int substep_func_omp_T2D_ME_mt_Geo(void* _self, size_t my_th_id,
	double dt, double cur_time, size_t substp_id);

class Step_T2D_ME_mt_Geo : public Step_OMP
{
protected:
	typedef Model_T2D_ME_mt::ShapeFunc ShapeFunc;
	typedef Model_T2D_ME_mt::ShapeFuncAB DShapeFuncAB;
	typedef Model_T2D_ME_mt::ShapeFuncC DShapeFuncC;
	typedef Model_T2D_ME_mt::Force Force;
	typedef Model_T2D_ME_mt::Position Position;
	typedef Model_T2D_ME_mt::Displacement Displacement;
	typedef Model_T2D_ME_mt::Velocity Velocity;
	typedef Model_T2D_ME_mt::Acceleration Acceleration;
	typedef Model_T2D_ME_mt::Stress Stress;
	typedef Model_T2D_ME_mt::Strain Strain;
	typedef Model_T2D_ME_mt::StrainInc StrainInc;
	typedef Model_T2D_ME_mt::SortedPclVarArrays SortedPclVarArrays;
	typedef Model_T2D_ME_mt::ElemNodeIndex ElemNodeIndex;
	typedef Model_T2D_ME_mt::ElemNodeVM ElemNodeVM;
	typedef Model_T2D_ME_mt::NodeHasVBC NodeHasVBC;

	double* pcl_m; // ori_pcl_num
	Force* pcl_bf; // ori_pcl_num
	Force* pcl_t; // ori_pcl_num
	Position* pcl_pos; // ori_pcl_num
	double* pcl_vol; // ori_pcl_num
	MatModel::MaterialModel** pcl_mat_model; // ori_pcl_num

	SortedPclVarArrays sorted_pcl_var_arrays[2];

	ElemNodeIndex *elem_node_id;
	DShapeFuncAB *elem_N_ab;
	DShapeFuncC *elem_N_c;
	double* elem_area;

	double* elem_pcl_m; // elem_num
	StrainInc *elem_de; // elem_num
	double* elem_m_de_vol; // elem_num

	// element-node data
	ElemNodeVM* elem_node_vm; // elem_num * 3
	Force* elem_node_force; // elem_num * 3

	Acceleration* node_a; // node_num
	Velocity* node_v; // node_num
	NodeHasVBC* node_has_vbc; // node_num
	double* node_am; // node_num
	double* node_de_vol; // node_num

	size_t elem_num, node_num;
	double k, miu, Kf;
	
	size_t valid_pcl_num;
	size_t valid_elem_num;

	size_t* elem_count_bin;
	size_t* elem_sum_bin;

	size_t* prev_pcl_ids[2]; // (pcl_num + 2) * 2
	size_t* pcl_in_elems[2]; // pcl_num * 2
	size_t* valid_elem_id; // elem_num
	size_t* node_has_elems[2]; // (elem_num * 4 + 2) * 2
	size_t* node_elem_pairs[2]; // elem_num * 4 * 2
	
	union ThreadData
	{
		struct
		{
			size_t p_id0, p_id1;
			size_t valid_elem_num;
			size_t *valid_elem_id;
			size_t ve_id0, ve_id1;
		};
		char padding[Cache_Alignment];
	};
	ThreadData* thread_datas;
	
	CacheAlignedMem thread_mem;
	CacheAlignedMem cal_mem;
	
	// 0 continue
	// 1 reset velocity
	// 2 exit
	unsigned char cal_status;
	bool init_f_ub_is_init;
	bool e_kin_max_is_init;
	double e_kin_prev;
	double f_ub_ratio_limit;
	double e_kin_ratio_limit;

	double init_f_ub;
	double f_ub;
	double f_ub_ratio;

	double e_kin_max;
	double e_kin;
	double e_kin_ratio;
	
public:
	int init_calculation() override;
	friend int substep_func_omp_T2D_ME_mt_Geo(void* _self,
		size_t my_th_id, double dt, double cur_time, size_t substp_id);
	int finalize_calculation() override;

public:
	Step_T2D_ME_mt_Geo(const char* _name);
	~Step_T2D_ME_mt_Geo();

	inline size_t get_pcl_num() const noexcept { return valid_pcl_num; }
	inline double get_f_ub() const noexcept { return f_ub; }
	inline double get_e_kin() const noexcept { return e_kin; }
	inline double get_f_ub_ratio() const noexcept { return f_ub_ratio; }
	inline double get_e_kin_ratio() const noexcept { return e_kin_ratio; }
	inline void set_f_ub_ratio_limit(double _lim) noexcept { f_ub_ratio_limit = _lim; }
	inline void set_e_kin_ratio_limit(double _lim) noexcept { e_kin_ratio_limit = _lim; }

	static void reorder(Model_T2D_ME_mt& md, size_t valid_pcl_num, const size_t* prev_pcl_ids);

	friend int Model_T2D_ME_mt_hdf5_utilities::output_pcl_data_to_hdf5_file(Model_T2D_ME_mt& md, Step_T2D_ME_mt_Geo& stp, ResultFile_hdf5& rf, hid_t grp_id);
};

#endif