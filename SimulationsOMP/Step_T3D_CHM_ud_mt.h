#ifndef __Step_T3D_CHM_ud_mt_h__
#define __Step_T3D_CHM_ud_mt_h__

#include "Step_T3D_CHM_mt.h"

class Step_T3D_CHM_ud_mt;
namespace Model_T3D_CHM_mt_hdf5_utilities
{
	int load_model_from_hdf5_file(Model_T3D_CHM_mt& md, Step_T3D_CHM_ud_mt& step,
		const char* hdf5_name, const char* th_name, size_t frame_id);
}

int substep_func_omp_T3D_CHM_ud_mt(void* _self, size_t my_th_id,
	double dt, double cur_time, size_t substp_id);

class Step_T3D_CHM_ud_mt : public Step_T3D_CHM_mt
{
public:
	int init_calculation() override;
	friend int substep_func_omp_T3D_CHM_ud_mt(void* _self,
		size_t my_th_id, double dt, double cur_time, size_t substp_id);
	int finalize_calculation() override;

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
	
public:
	Step_T3D_CHM_ud_mt(const char* _name);
	~Step_T3D_CHM_ud_mt();

	friend int Model_T3D_CHM_mt_hdf5_utilities::load_model_from_hdf5_file(
		Model_T3D_CHM_mt &md, Step_T3D_CHM_ud_mt &step, const char *hdf5_name,
		const char* th_name, size_t frame_id);
};

#endif