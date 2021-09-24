#ifndef __Model_T3D_CHM_mt_h__
#define __Model_T3D_CHM_mt_h__

#include "BCs.h"
#include "macro_utils.h"
#include "Model.h"
#include "MatModelContainer.h"
#include "ParticleGenerator3D.hpp"
#include "TetrahedronMesh.h"
#include "RigidObject/RigidCylinder.h"
#include "RigidObject/RigidObjectByT3DMesh.h"
#include "RigidObject/SmoothContact3D.h"
#include "RigidObject/RoughContact3D.h"
#include "RigidObject/FrictionalContact3D.h"
#include "RigidObject/StickyContact3D.h"

class Model_T3D_CHM_mt;
class Step_T3D_CHM_mt;
int substep_func_omp_T3D_CHM_mt(void* _self, size_t my_th_id,
	double dt, double cur_time, size_t substp_id);
int substep_func_omp_T3D_CHM_mt2(void* _self, size_t my_th_id, double dt, double cur_time, size_t substp_id);
class Step_T3D_CHM_ud_mt;
int substep_func_omp_T3D_CHM_ud_mt(void* _self,
	size_t my_th_id, double dt, double cur_time, size_t substp_id);
class Step_T3D_CHM_ud_mt_subiter;
int substep_func_omp_T3D_CHM_ud_mt_subiter(void* _self,
	size_t my_th_id, double dt, double cur_time, size_t substp_id);
class Step_T3D_CHM_mt_Geo;
int substep_func_omp_T3D_CHM_mt_Geo(void* _self,
	size_t my_th_id, double dt, double cur_time, size_t substp_id);
class Step_T3D_CHM_TBB;
int substep_func_T3D_CHM_TBB(void* _self);
namespace Step_T3D_CHM_Task { class CalData; }

class ResultFile_hdf5;
namespace Model_T3D_CHM_mt_hdf5_utilities
{
	struct ParticleData;
	int output_background_mesh_to_hdf5_file(Model_T3D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_background_mesh_from_hdf5_file(Model_T3D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_search_mesh_to_hdf5_file(Model_T3D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_search_mesh_from_hdf5_file(Model_T3D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_boundary_condition_to_hdf5_file(Model_T3D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_boundary_condition_from_hdf5_file(Model_T3D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_ori_pcl_data_to_hdf5_file(Model_T3D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_pcl_data_to_hdf5_file(Model_T3D_CHM_mt& md, Step_T3D_CHM_mt &stp, ResultFile_hdf5& rf, hid_t grp_id);
	int load_pcl_data_from_hdf5_file(Model_T3D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_material_model_to_hdf5_file(Model_T3D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_material_model_from_hdf5_file(Model_T3D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_rigid_cylinder_to_hdf5_file(Model_T3D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_rigid_cylinder_from_hdf5_file(Model_T3D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_t3d_rigid_mesh_to_hdf5_file(Model_T3D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_t3d_rigid_mesh_state_to_hdf5_file(Model_T3D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_t3d_rigid_mesh_from_hdf5_file(Model_T3D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
}

struct Model_T3D_CHM_mt : public Model,
	public MatModel::MatModelContainer
{
	friend class Step_T3D_CHM_mt;
	friend int substep_func_omp_T3D_CHM_mt(void* _self,
		size_t my_th_id, double dt, double cur_time, size_t substp_id);
	friend int substep_func_omp_T3D_CHM_mt2(void* _self,
		size_t my_th_id, double dt, double cur_time, size_t substp_id);
	
	friend class Step_T3D_CHM_ud_mt;
	friend int substep_func_omp_T3D_CHM_ud_mt(void* _self,
		size_t my_th_id, double dt, double cur_time, size_t substp_id);

	friend class Step_T3D_CHM_ud_mt_subiter;
	friend int substep_func_omp_T3D_CHM_ud_mt_subiter(void* _self,
		size_t my_th_id, double dt, double cur_time, size_t substp_id);

	friend class Step_T3D_CHM_mt_Geo;
	friend int substep_func_omp_T3D_CHM_mt_Geo(void* _self,
		size_t my_th_id, double dt, double cur_time, size_t substp_id);

	friend class Step_T3D_CHM_TBB;
	friend int substep_func_T3D_CHM_TBB(void* _self);
	friend class Step_T3D_CHM_Task::CalData;

public:
	struct ShapeFunc { double N1, N2, N3, N4; };
	struct DShapeFuncABC
	{
		union { double a1; double dN1_dx; };
		union { double b1; double dN1_dy; };
		union { double c1; double dN1_dz; };
		union { double a2; double dN2_dx; };
		union { double b2; double dN2_dy; };
		union { double c2; double dN2_dz; };
		union { double a3; double dN3_dx; };
		union { double b3; double dN3_dy; };
		union { double c3; double dN3_dz; };
		union { double a4; double dN4_dx; };
		union { double b4; double dN4_dy; };
		union { double c4; double dN4_dz; };
	};
	struct DShapeFuncD { double d1, d2, d3, d4; };
	
	union Force
	{
		struct { double fx, fy, fz; };
		Vector3D vec;
		Force() {}
	};

	union Acceleration
	{
		struct { double ax, ay, az; };
		struct { size_t iax, iay, iaz; };
	};
	union Velocity
	{
		struct { double vx, vy, vz; };
		struct { size_t ivx, ivy, ivz; };
	};
	union Position
	{
		struct { double x, y, z; };
		Point3D pt;
	};
	union Displacement
	{
		struct { double ux, uy, uz; };
		Vector3D vec;
	};

	union Stress
	{
		struct { double s11, s22, s33, s12, s23, s31; };
		double s[6];
	};
	union Strain
	{
		struct { double e11, e22, e33, e12, e23, e31; };
		double e[6];
	};
	union StrainInc
	{
		struct { double de11, de22, de33, de12, de23, de31; };
		double de[6];
	};
	
	struct ElemNodeIndex { size_t n1, n2, n3, n4; };
	struct ElemNodeVM { double vm, vmx, vmy, vmz; };

	struct NodeHasVBC { bool has_vx_bc, has_vy_bc, has_vz_bc; };
	struct NodeVBCVec { double x, y, z; };

	struct SortedPclVarArrays
	{
		size_t *pcl_index; // ori_pcl_num
		double* pcl_n;
		double *pcl_density_f; // ori_pcl_num
		Velocity *pcl_v_s; // ori_pcl_num
		Velocity* pcl_v_f; // ori_pcl_num
		Displacement *pcl_u_s; // ori_pcl_num
		Displacement* pcl_u_f; // ori_pcl_num
		Stress* pcl_stress; // ori_pcl_num
		double* pcl_p; // ori_pcl_num
		Strain* pcl_strain; // ori_pcl_num
		Strain* pcl_estrain; // ori_pcl_num
		Strain* pcl_pstrain; // ori_pcl_num
		ShapeFunc* pcl_N; // ori_pcl_num
	};

protected:
	size_t pcl_num;
	size_t elem_num;
	size_t node_num;
	
	double *pcl_m_s; // ori_pcl_num
	double* pcl_density_s; // ori_pcl_num
	double* pcl_vol_s; // ori_pcl_num
	Force *pcl_bf_s; // ori_pcl_num
	Force* pcl_bf_f; // ori_pcl_num
	Force *pcl_t; // ori_pcl_num
	Position *pcl_pos; // ori_pcl_num
	MatModel::MaterialModel **pcl_mat_model; // ori_pcl_num
	size_t *pcl_mat_model_copy_offset; // ori_pcl_num
	size_t pcl_mat_model_total_size;

	double* pcl_vol; // ori_pcl_num
	Stress* pcl_prev_stress; // ori_pcl_num
	StrainInc* pcl_destrain; // ori_pcl_num
	StrainInc* pcl_dpstrain; // ori_pcl_num

	SortedPclVarArrays sorted_pcl_var_arrays[2];

	// element data
	ElemNodeIndex *elem_node_id; // elem_num
	DShapeFuncABC* elem_N_abc; // elem_num
	DShapeFuncD* elem_N_d; // elem_num
	double* elem_vol; // elem_num

	// element calculation data
	double* elem_pcl_int_vol; // elem_num
	double *elem_density_f; // elem_num
	double *elem_pcl_n; // elem_num
	double *elem_pcl_m_s; // elem_num
	double *elem_pcl_m_f; // elem_num
	StrainInc *elem_de; // elem_num
	double *elem_p; // elem_num
	double *elem_n2_miu_div_k_vol; // elem_num
	Force *elem_seep_force; // elem_num
	double *elem_m_de_vol_s; // elem_num
	double *elem_m_de_vol_f; // elem_num

	// element-node data
	ElemNodeVM* elem_node_vm_s; // elem_num * 4
	ElemNodeVM* elem_node_vm_f; // elem_num * 4
	Force* elem_node_force_s; // elem_num * 4
	Force* elem_node_force_f; // elem_num * 4
	
	// node data
	Position* node_pos; // node_num
	double* node_am_s; // node_num
	double* node_am_f; // node_num
	Acceleration *node_a_s; // node_num
	Acceleration *node_a_f; // node_num
	Velocity *node_v_s; // node_num
	Velocity *node_v_f; // node_num
	Displacement *node_du_s; // node_num
	Displacement* node_du_f; // node_num
	Velocity *node_vn_s; // node_num
	Velocity* node_vn_f; // node_num
	Velocity *node_pv_s; // node_num
	Velocity* node_pv_f; // node_num
	Displacement* node_pdu_s; // node_num
	Displacement* node_pdu_f;// node_num
	double* node_de_vol_s; // node_num
	double* node_de_vol_f; // node_num

	NodeHasVBC *node_has_vbc_s; // node_num
	NodeHasVBC *node_has_vbc_f; // node_num
	NodeVBCVec* node_vbc_vec_s; // node_num
	NodeVBCVec* node_vbc_vec_f; // node_num

	double Kf, k, miu;

// ======== Non calculation data ========
	size_t ori_pcl_num;
	char *pcl_mem_raw;

	char *mesh_mem_raw;

	// boundary conditions
	size_t bfx_s_num, bfy_s_num, bfz_s_num;
	BodyForceAtPcl *bfx_ss, *bfy_ss, *bfz_ss;
	size_t bfx_f_num, bfy_f_num, bfz_f_num;
	BodyForceAtPcl* bfx_fs, *bfy_fs, *bfz_fs;
	size_t tx_num, ty_num, tz_num;
	TractionBCAtPcl* txs, *tys, *tzs;

public:
	Model_T3D_CHM_mt();
	~Model_T3D_CHM_mt();

	inline double get_miu() const noexcept { return miu; }
	inline double get_k() const noexcept { return k; }
	inline double get_Kf() const noexcept { return Kf; }
	inline void set_fluid_props(double _Kf, double _k, double _miu) noexcept
	{ Kf = _Kf; k = _k; miu = _miu; }
	inline size_t get_ori_pcl_num() const noexcept { return ori_pcl_num; }
	inline size_t get_pcl_num() const noexcept { return pcl_num; }
	inline const size_t* get_pcl_index0() const noexcept { return sorted_pcl_var_arrays[0].pcl_index; }
	inline const Position *get_pcl_pos() const noexcept { return pcl_pos; }
	inline const double* get_pcl_m_s() const noexcept { return pcl_m_s; }
	inline const double* get_pcl_density_s() const noexcept { return pcl_density_s; }
	inline const double* get_pcl_n0() const noexcept { return sorted_pcl_var_arrays[0].pcl_n; }
	inline const double *get_pcl_vol() const noexcept { return pcl_vol; }
	inline Stress* get_pcl_stress0() noexcept { return sorted_pcl_var_arrays[0].pcl_stress; }
	inline size_t get_node_num() const noexcept { return node_num; }
	inline const Position *get_node_pos() const noexcept { return node_pos; }
	inline size_t get_elem_num() const noexcept { return elem_num; }
	inline const ElemNodeIndex* get_elem_node_index() const noexcept { return elem_node_id; }
	inline MatModel::MaterialModel **get_mat_models() noexcept { return pcl_mat_model; }
	Cube get_mesh_bbox();

	inline const NodeHasVBC *get_node_has_vbc_s() const noexcept { return node_has_vbc_s; }
	inline const NodeHasVBC* get_node_has_vbc_f() const noexcept { return node_has_vbc_f; }
	inline const NodeVBCVec* get_node_vbc_vec_s() const noexcept { return node_vbc_vec_s; }
	inline const NodeVBCVec* get_node_vbc_vec_f() const noexcept { return node_vbc_vec_f; }

	inline double get_bg_grid_xl() const noexcept { return grid_xl; }
	inline double get_bg_grid_yl() const noexcept { return grid_yl; }
	inline double get_bg_grid_zl() const noexcept { return grid_zl; }
	inline double get_bg_grid_hx() const noexcept { return grid_hx; }
	inline double get_bg_grid_hy() const noexcept { return grid_hy; }
	inline double get_bg_grid_hz() const noexcept { return grid_hz; }
	inline size_t get_grid_x_num() const noexcept { return grid_x_num; }
	inline size_t get_grid_y_num() const noexcept { return grid_y_num; }
	inline size_t get_grid_z_num() const noexcept { return grid_z_num; }
	inline const size_t* get_grid_elem_list_id_array() const noexcept { return grid_elem_list_id_array; }
	inline const size_t* get_grid_elem_list() const noexcept { return grid_elem_list; }
	
	inline void set_Kf(double _Kf) noexcept { Kf = _Kf; }
	inline void set_k(double _k) noexcept { k = _k; }
	inline void set_miu(double _miu) noexcept { miu = _miu; }

	void clear_mesh();
	void alloc_mesh(size_t n_num, size_t e_num);
	void init_mesh(const TetrahedronMesh &mesh);
	
	void clear_search_grid();
	int init_search_grid(TetrahedronMesh &mesh);
	
	void alloc_pcls(size_t num);
	void alloc_pcls(size_t num, size_t ori_num);
	void clear_pcls();
	int init_pcls(size_t num, double n, double m_s,
		double density_s, double density_f,
		double _Kf, double _k, double _miu);
	int init_pcls(ParticleGenerator3D<TetrahedronMesh> &pg,
		double n, double density_s, double density_f,
		double _Kf, double _k, double _miu);

	INIT_BC_TEMPLATE(bfx_s, BodyForceAtPcl)
	INIT_BC_TEMPLATE(bfy_s, BodyForceAtPcl)
	INIT_BC_TEMPLATE(bfz_s, BodyForceAtPcl)
	INIT_BC_TEMPLATE(bfx_f, BodyForceAtPcl)
	INIT_BC_TEMPLATE(bfy_f, BodyForceAtPcl)
	INIT_BC_TEMPLATE(bfz_f, BodyForceAtPcl)
	INIT_BC_TEMPLATE(tx, TractionBCAtPcl)
	INIT_BC_TEMPLATE(ty, TractionBCAtPcl)
	INIT_BC_TEMPLATE(tz, TractionBCAtPcl)
	void init_bfx_ss(size_t bf_num, const size_t *bf_pcls, const double *bfs);
	void init_bfy_ss(size_t bf_num, const size_t *bf_pcls, const double *bfs);
	void init_bfz_ss(size_t bf_num, const size_t* bf_pcls, const double* bfs);
	void init_bfx_fs(size_t bf_num, const size_t* bf_pcls, const double* bfs);
	void init_bfy_fs(size_t bf_num, const size_t* bf_pcls, const double* bfs);
	void init_bfz_fs(size_t bf_num, const size_t* bf_pcls, const double* bfs);
	void init_txs(size_t t_num, const size_t *t_pcls, const double *ts);
	void init_tys(size_t t_num, const size_t *t_pcls, const double *ts);
	void init_tzs(size_t t_num, const size_t* t_pcls, const double* ts);

	void init_fixed_vx_s_bc(size_t vx_bc_num, const size_t *vx_bcs);
	void init_fixed_vy_s_bc(size_t vy_bc_num, const size_t* vy_bcs);
	void init_fixed_vz_s_bc(size_t vz_bc_num, const size_t* vz_bcs);
	void init_fixed_vx_f_bc(size_t vx_bc_num, const size_t* vx_bcs);
	void init_fixed_vy_f_bc(size_t vy_bc_num, const size_t* vy_bcs);
	void init_fixed_vz_f_bc(size_t vz_bc_num, const size_t* vz_bcs);
	void set_vbc_vec_s(size_t n_id, double vecx, double vecy, double vecz);
	void set_vbc_vec_f(size_t n_id, double vecx, double vecy, double vecz);

	inline void add_mat_model(size_t pcl_id,
		MatModel::MaterialModel &mat_model,
		size_t model_size)
	{
		pcl_mat_model[pcl_id] = &mat_model;
		pcl_mat_model_copy_offset[pcl_id] = pcl_mat_model_total_size;
		pcl_mat_model_total_size += model_size;
	}

protected:
	// background grid for mesh
	double grid_xl, grid_yl, grid_zl;
	double grid_xu, grid_yu, grid_zu;
	double grid_hx, grid_hy, grid_hz;
	size_t grid_x_num, grid_y_num, grid_z_num;
	size_t grid_xy_num;
	size_t* grid_elem_list_id_array;
	size_t* grid_elem_list;

public:
	inline bool is_in_element(
		double pcl_x,
		double pcl_y,
		double pcl_z,
		size_t elem_id,
		ShapeFunc& p_N
	) const noexcept
	{
		const DShapeFuncABC& e_dN_abc = elem_N_abc[elem_id];
		const DShapeFuncD& e_dN_d = elem_N_d[elem_id];
		p_N.N1 = e_dN_abc.a1 * pcl_x + e_dN_abc.b1 * pcl_y + e_dN_abc.c1 * pcl_z + e_dN_d.d1;
		p_N.N2 = e_dN_abc.a2 * pcl_x + e_dN_abc.b2 * pcl_y + e_dN_abc.c2 * pcl_z + e_dN_d.d2;
		p_N.N3 = e_dN_abc.a3 * pcl_x + e_dN_abc.b3 * pcl_y + e_dN_abc.c3 * pcl_z + e_dN_d.d3;
		p_N.N4 = e_dN_abc.a4 * pcl_x + e_dN_abc.b4 * pcl_y + e_dN_abc.c4 * pcl_z + e_dN_d.d4;
		return p_N.N1 >= 0.0 && p_N.N1 <= 1.0 && p_N.N2 >= 0.0 && p_N.N2 <= 1.0
			&& p_N.N3 >= 0.0 && p_N.N3 <= 1.0 && p_N.N4 >= 0.0 && p_N.N4 <= 1.0;
	}

	inline bool is_in_element_tol(
		double pcl_x,
		double pcl_y,
		double pcl_z,
		size_t elem_id,
		ShapeFunc& p_N
	) const noexcept
	{
		const DShapeFuncABC& e_N_abc = elem_N_abc[elem_id];
		const DShapeFuncD& e_N_d = elem_N_d[elem_id];
		p_N.N1 = e_N_abc.a1 * pcl_x + e_N_abc.b1 * pcl_y + e_N_abc.c1 * pcl_z + e_N_d.d1;
		p_N.N2 = e_N_abc.a2 * pcl_x + e_N_abc.b2 * pcl_y + e_N_abc.c2 * pcl_z + e_N_d.d2;
		p_N.N3 = e_N_abc.a3 * pcl_x + e_N_abc.b3 * pcl_y + e_N_abc.c3 * pcl_z + e_N_d.d3;
		p_N.N4 = e_N_abc.a4 * pcl_x + e_N_abc.b4 * pcl_y + e_N_abc.c4 * pcl_z + e_N_d.d4;
#define in_elem_N_tol (1.0e-8)
		if (p_N.N1 >= -in_elem_N_tol && p_N.N1 <= (1.0 + in_elem_N_tol) &&
			p_N.N2 >= -in_elem_N_tol && p_N.N2 <= (1.0 + in_elem_N_tol) &&
			p_N.N3 >= -in_elem_N_tol && p_N.N3 <= (1.0 + in_elem_N_tol) &&
			p_N.N4 >= -in_elem_N_tol && p_N.N4 <= (1.0 + in_elem_N_tol))
		{
			if (p_N.N1 < 0.0)
				p_N.N1 = 0.0;
			if (p_N.N1 > 1.0)
				p_N.N1 = 1.0;
			if (p_N.N2 < 0.0)
				p_N.N2 = 0.0;
			if (p_N.N2 > 1.0)
				p_N.N2 = 1.0;
			if (p_N.N3 < 0.0)
				p_N.N3 = 0.0;
			if (p_N.N3 > 1.0)
				p_N.N3 = 1.0;
			if (p_N.N4 < 0.0)
				p_N.N4 = 0.0;
			if (p_N.N4 > 1.0)
				p_N.N4 = 1.0;
			return true;
		}
		return false;
#undef in_elem_N_tol
	}

	inline size_t find_pcl_in_which_elem(
		double pcl_x,
		double pcl_y,
		double pcl_z,
		ShapeFunc& p_N
		) const noexcept
	{
		if (pcl_x < grid_xl || pcl_x > grid_xu ||
			pcl_y < grid_yl || pcl_y > grid_yu ||
			pcl_z < grid_zl || pcl_z > grid_zu)
			return SIZE_MAX;
		const size_t x_id = (pcl_x - grid_xl) / grid_hx;
		const size_t y_id = (pcl_y - grid_yl) / grid_hy;
		const size_t z_id = (pcl_z - grid_zl) / grid_hz;
		const size_t g_id = grid_xy_num * z_id + grid_x_num * y_id + x_id;
		const size_t elem_end_id = grid_elem_list[g_id + 1];
		for (size_t el_id = grid_elem_list[g_id]; el_id < elem_end_id; ++el_id)
		{
			const size_t elem_id = grid_elem_list_id_array[el_id];
			if (is_in_element(pcl_x, pcl_y, pcl_z, elem_id, p_N))
				return elem_id;
		}
		return SIZE_MAX;
	}

	inline size_t find_pcl_in_which_elem_tol(
		double pcl_x,
		double pcl_y,
		double pcl_z,
		ShapeFunc& p_N
		) const noexcept
	{
		if (pcl_x < grid_xl || pcl_x > grid_xu ||
			pcl_y < grid_yl || pcl_y > grid_yu ||
			pcl_z < grid_zl || pcl_z > grid_zu)
			return SIZE_MAX;
		const size_t x_id = (pcl_x - grid_xl) / grid_hx;
		const size_t y_id = (pcl_y - grid_yl) / grid_hy;
		const size_t z_id = (pcl_z - grid_zl) / grid_hz;
		const size_t g_id = grid_xy_num * z_id + grid_x_num * y_id + x_id;
		const size_t elem_end_id = grid_elem_list[g_id + 1];
		for (size_t el_id = grid_elem_list[g_id];
			el_id < elem_end_id; ++el_id)
		{
			const size_t elem_id = grid_elem_list_id_array[el_id];
			if (is_in_element_tol(pcl_x, pcl_y, pcl_z, elem_id, p_N))
				return elem_id;
		}
		return SIZE_MAX;
	}

protected:
	bool rigid_cylinder_is_valid;
	RigidCylinder rigid_cylinder;
	bool rigid_t3d_mesh_is_valid;
	RigidObjectByT3DMesh rigid_t3d_mesh;
	
	// ad hoc design for output
	double Ksn_cont, Kst_cont;
	double fric_ratio, shear_strength;
	double Kfn_cont, Kft_cont;
	SmoothContact3D smooth_contact_s, smooth_contact_f;
	RoughContact3D rough_contact_s, rough_contact_f;
	FrictionalContact3D fric_contact_s;
	StickyContact3D sticky_contact_s;
	ContactModel3D *pcm_s, *pcm_f;

	// solid
	size_t* contact_substep_ids_s; // ori_pcl_num
	ContactModel3D::Position* prev_contact_poses_s; // ori_pcl_num
	ContactModel3D::Force* prev_contact_tan_forces_s; // ori_pcl_num
	double * prev_contact_dists_s; // ori_pcl_num
	// fluid
	size_t* contact_substep_ids_f; // ori_pcl_num
	ContactModel3D::Position* prev_contact_poses_f; // ori_pcl_num
	ContactModel3D::Force* prev_contact_tan_forces_f; // ori_pcl_num
	double * prev_contact_dists_f; // ori_pcl_num

	char* contact_mem;
	void clear_contact_mem();
	void alloc_contact_mem(size_t num);
	
public:
	// cylinder
	inline bool has_rigid_cylinder() const noexcept { return rigid_cylinder_is_valid; }
	inline RigidCylinder &get_rigid_cylinder() { return  rigid_cylinder; }
	void clear_rigid_cylinder() noexcept { rigid_cylinder_is_valid = false; }
	inline void init_rigid_cylinder(
		double x, double y, double z,
		double h, double r)
	{
		rigid_cylinder_is_valid = true;
		rigid_cylinder.init(x, y, z, h, r);
	}
	inline void set_rigid_cylinder_velocity(double vx, double vy, double vz)
	{ rigid_cylinder.set_vbc(vx, vy, vz); }

	// t3d rigid mesh
	bool has_t3d_rigid_mesh() const noexcept { return rigid_t3d_mesh_is_valid; }
	RigidObjectByT3DMesh &get_t3d_rigid_mesh() noexcept { return rigid_t3d_mesh; }
	void clear_t3d_rigid_mesh() noexcept { rigid_t3d_mesh_is_valid = false; }
	inline void init_t3d_rigid_mesh(
		double _density, const char *filename,
		double dx, double dy, double dz,
		double dx_ang, double dy_ang, double dz_ang,
		double ghx, double ghy, double ghz)
	{
		rigid_t3d_mesh_is_valid = true;
		rigid_t3d_mesh.init(_density, filename,
			dx, dy, dz, dx_ang, dy_ang, dz_ang, ghx, ghy, ghz);
	}
	inline void set_t3d_rigid_mesh_velocity(double vx, double vy, double vz)
	{ rigid_t3d_mesh.set_translation_velocity_bc(vx, vy, vz); }

	// for contact model
	inline void set_contact_param(
		double _Ksn_cont, double _Kst_cont,
		double _fric_ratio, double _shear_strength,
		double _Kfn_cont, double _Kft_cont)
	{
		Ksn_cont = _Ksn_cont;
		Kst_cont = _Kst_cont;
		fric_ratio = _fric_ratio;
		shear_strength = _shear_strength;
		Kfn_cont = _Kfn_cont;
		Kft_cont = _Kft_cont;
		smooth_contact_s.set_Kn_cont(_Ksn_cont);
		rough_contact_s.set_K_cont(_Ksn_cont, _Kst_cont);
		fric_contact_s.set_K_cont(_Ksn_cont, _Kst_cont);
		fric_contact_s.set_friction_ratio(_fric_ratio);
		sticky_contact_s.set_K_cont(_Ksn_cont, _Kst_cont);
		sticky_contact_s.set_shear_strength(_shear_strength);
		smooth_contact_f.set_Kn_cont(_Kfn_cont);
		rough_contact_f.set_K_cont(_Kfn_cont, _Kft_cont);
	}
	inline void set_smooth_contact_between_spcl_and_rect() noexcept { pcm_s = &smooth_contact_s; }
	inline void set_rough_contact_between_spcl_and_rect() noexcept { pcm_s = &rough_contact_s; }
	inline void set_frictional_contact_between_spcl_and_rect() noexcept { pcm_s = &fric_contact_s; }
	inline void set_sticky_contact_between_spcl_and_rect() noexcept { pcm_s = &sticky_contact_s; }
	inline void set_smooth_contact_between_fpcl_and_rect() noexcept { pcm_f = &smooth_contact_f; }
	inline void set_rough_contact_between_fpcl_and_rect() noexcept { pcm_f = &rough_contact_f; }
	inline ContactModel3D *get_contact_model_s() noexcept { return pcm_s; }
	inline ContactModel3D* get_contact_model_f() noexcept { return pcm_f; }

protected:
	friend struct Model_T3D_CHM_mt_hdf5_utilities::ParticleData;
	friend int Model_T3D_CHM_mt_hdf5_utilities::output_background_mesh_to_hdf5_file(Model_T3D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_CHM_mt_hdf5_utilities::load_background_mesh_from_hdf5_file(Model_T3D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_CHM_mt_hdf5_utilities::output_search_mesh_to_hdf5_file(Model_T3D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_CHM_mt_hdf5_utilities::load_search_mesh_from_hdf5_file(Model_T3D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_CHM_mt_hdf5_utilities::output_boundary_condition_to_hdf5_file(Model_T3D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_CHM_mt_hdf5_utilities::load_boundary_condition_from_hdf5_file(Model_T3D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_CHM_mt_hdf5_utilities::output_ori_pcl_data_to_hdf5_file(Model_T3D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_CHM_mt_hdf5_utilities::output_pcl_data_to_hdf5_file(Model_T3D_CHM_mt& md, Step_T3D_CHM_mt& stp, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_CHM_mt_hdf5_utilities::load_pcl_data_from_hdf5_file(Model_T3D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_CHM_mt_hdf5_utilities::output_material_model_to_hdf5_file(Model_T3D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_CHM_mt_hdf5_utilities::load_material_model_from_hdf5_file(Model_T3D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_CHM_mt_hdf5_utilities::output_rigid_cylinder_to_hdf5_file(Model_T3D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_CHM_mt_hdf5_utilities::load_rigid_cylinder_from_hdf5_file(Model_T3D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_CHM_mt_hdf5_utilities::output_t3d_rigid_mesh_to_hdf5_file(Model_T3D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_CHM_mt_hdf5_utilities::output_t3d_rigid_mesh_state_to_hdf5_file(Model_T3D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_CHM_mt_hdf5_utilities::load_t3d_rigid_mesh_from_hdf5_file(Model_T3D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
};

#endif