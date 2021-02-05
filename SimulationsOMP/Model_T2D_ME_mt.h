#ifndef __Model_T2D_ME_mt_h__
#define __Model_T2D_ME_mt_h__

#include <stdint.h>

#include "BCs.h"
#include "macro_utils.h"
#include "Model.h"
#include "MatModelContainer.h"
#include "ParticleGenerator2D.hpp"
#include "TriangleMesh.h"
#include "RigidBody/RigidRect.h"
#include "RigidObject/Force2D.h"

class Model_T2D_ME_mt;
class Step_T2D_ME_mt;
int substep_func_omp_T2D_ME_mt(void* _self, size_t my_th_id,
	double dt, double cur_time, size_t substp_id);
class Step_T2D_ME_TBB;
int cal_substep_func_T2D_ME_TBB(void* _self);

namespace Step_T2D_ME_Task
{
	class TaskData;
	class MapPclToBgMeshTask;
}

class ResultFile_hdf5;
namespace Model_T2D_ME_mt_hdf5_utilities
{
	struct ParticleData;
	int output_background_mesh_to_hdf5_file(Model_T2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_background_mesh_from_hdf5_file(Model_T2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_boundary_condition_to_hdf5_file(Model_T2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_boundary_condition_from_hdf5_file(Model_T2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_pcl_data_to_hdf5_file(Model_T2D_ME_mt& md, Step_T2D_ME_mt &stp, ResultFile_hdf5& rf, hid_t grp_id);
	int load_pcl_data_from_hdf5_file(Model_T2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_material_model_to_hdf5_file(Model_T2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_material_model_from_hdf5_file(Model_T2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_rigid_rect_to_hdf5_file(Model_T2D_ME_mt& md, Step_T2D_ME_mt& stp, ResultFile_hdf5& rf, hid_t grp_id);
	int output_rigid_rect_to_hdf5_file(Model_T2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_rigid_rect_from_hdf5_file(Model_T2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
}

struct Model_T2D_ME_mt : public Model,
	public MatModel::MatModelContainer
{
	friend class Step_T2D_ME_mt;
	friend int substep_func_omp_T2D_ME_mt(void* _self,
		size_t my_th_id, double dt, double cur_time, size_t substp_id);
	
	friend class Step_T2D_ME_TBB;
	friend int cal_substep_func_T2D_ME_TBB(void* _self);
	friend class Step_T2D_ME_Task::TaskData;
	friend class Step_T2D_ME_Task::MapPclToBgMeshTask;
	
public:
	struct ShapeFunc { double N1, N2, N3; };
	struct ShapeFuncAB
	{
		union { double a1; double dN1_dx; };
		union { double b1; double dN1_dy; };
		union { double a2; double dN2_dx; };
		union { double b2; double dN2_dy; };
		union { double a3; double dN3_dx; };
		union { double b3; double dN3_dy; };
	};
	struct ShapeFuncC { double c1, c2, c3; };

	union Force
	{
		struct { double fx, fy; };
		Vector2D vec;
		Force() {}
	};

	union Acceleration
	{
		struct { double ax, ay; };
		struct { size_t iax, iay; };
	};
	union Velocity
	{
		struct { double vx, vy; };
		struct { size_t ivx, ivy; };
	};
	union Position
	{
		struct { double x, y; };
		Point2D pt;
	};
	union Displacement
	{
		struct { double ux, uy; };
		Vector2D vec;
	};

	struct Stress { double s11, s22, s12; };
	struct Strain { double e11, e22, e12; };
	struct StrainInc { double de11, de22, de12; };
	
	struct ElemNodeIndex { size_t n1, n2, n3; };
	struct ElemNodeVM { double vm, vmx, vmy; };

	struct NodeHasVBC { bool has_vx_bc, has_vy_bc; };

	struct SortedPclVarArrays
	{
		size_t* pcl_index; // ori_pcl_num
		double* pcl_density; // ori_pcl_num
		Velocity* pcl_v; // ori_pcl_num
		Displacement* pcl_disp; // ori_pcl_num
		ShapeFunc* pcl_N; // ori_pcl_num
		Stress* pcl_stress; // ori_pcl_num
		Strain* pcl_strain; // ori_pcl_num
		Strain* pcl_estrain; // ori_pcl_num
		Strain* pcl_pstrain; // ori_pcl_num
	};
	
protected:
	size_t pcl_num;
	
	double *pcl_m; // ori_pcl_num
	Force *pcl_bf; // ori_pcl_num
	Force *pcl_t; // ori_pcl_num
	Position *pcl_pos; // ori_pcl_num
	double* pcl_vol; // ori_pcl_num
	MatModel::MaterialModel **pcl_mat_model; // ori_pcl_num

	SortedPclVarArrays sorted_pcl_var_arrays[2];

	// element data
	size_t elem_num;
	size_t node_num;

	ElemNodeIndex *elem_node_id; // elem_num
	ShapeFuncAB* elem_dN_ab; // elem_num
	ShapeFuncC* elem_dN_c; // elem_num
	double* elem_area; // elem_num

	// element calculation data
	double* elem_pcl_m; // elem_num
	double *elem_density; // elem_num
	StrainInc *elem_de; // elem_num
	double *elem_m_de_vol; // elem_num

	// element-node data
	ElemNodeVM* elem_node_vm; // elem_num * 3
	Force* elem_node_force; // elem_num * 3

	// node data
	Position* node_pos; // node_num
	Acceleration *node_a; // node_num
	Velocity *node_v; // node_num
	NodeHasVBC *node_has_vbc; // node_num
	double *node_am; // node_num
	double *node_de_vol; // node_num

// ======== Non calculation data ========
	size_t ori_pcl_num;
	char *pcl_mem_raw;
	char *mesh_mem_raw;

	// boundary conditions
	size_t bfx_num, bfy_num;
	BodyForceAtPcl* bfxs, * bfys;
	size_t tx_num, ty_num;
	TractionBCAtPcl* txs, * tys;

public:
	Model_T2D_ME_mt();
	~Model_T2D_ME_mt();

	inline size_t get_ori_pcl_num() const noexcept { return ori_pcl_num; }
	inline size_t get_pcl_num() const noexcept { return pcl_num; }
	inline const Position *get_pcl_pos() const noexcept { return pcl_pos; }
	inline const double *get_pcl_m() const noexcept { return pcl_m; }
	inline const double* get_pcl_density0() const noexcept { return sorted_pcl_var_arrays[0].pcl_density; }
	inline size_t get_node_num() const noexcept { return node_num; }
	inline const Position *get_node_pos() const noexcept { return node_pos; }
	inline size_t get_elem_num() const noexcept { return elem_num; }
	inline const ElemNodeIndex* get_elem_node_index() const noexcept { return elem_node_id; }
	inline MatModel::MaterialModel **get_mat_models() noexcept { return pcl_mat_model; }
	Rect get_mesh_bbox();

	inline const NodeHasVBC *get_has_vbcs() const noexcept { return node_has_vbc; }

	inline double get_bg_grid_xl() const noexcept { return grid_xl; }
	inline double get_bg_grid_yl() const noexcept { return grid_yl; }
	inline double get_bg_grid_hx() const noexcept { return grid_hx; }
	inline double get_bg_grid_hy() const noexcept { return grid_hy; }
	inline size_t get_grid_x_num() const noexcept { return grid_x_num; }
	inline size_t get_grid_y_num() const noexcept { return grid_y_num; }
	inline const size_t* get_grid_elem_list_id_array() const noexcept { return grid_elem_list_id_array; }
	inline const size_t* get_grid_elem_list() const noexcept { return grid_elem_list; }
	
	void clear_mesh();
	void alloc_mesh(size_t n_num, size_t e_num);
	void init_mesh(const TriangleMesh& mesh);
	
	void clear_search_grid();
	int init_search_grid(TriangleMesh& mesh, double _hx, double _hy);
	
	void alloc_pcls(size_t num);
	void alloc_pcls(size_t num, size_t ori_num);
	void clear_pcls();
	int init_pcls(size_t num, double m, double density);
	int init_pcls(ParticleGenerator2D<TriangleMesh>& pg, double density);

	INIT_BC_TEMPLATE(bfx, BodyForceAtPcl)
	INIT_BC_TEMPLATE(bfy, BodyForceAtPcl)
	INIT_BC_TEMPLATE(tx, TractionBCAtPcl)
	INIT_BC_TEMPLATE(ty, TractionBCAtPcl)
	void init_bfxs(size_t bf_num, const size_t *bf_pcls, const double *bfs);
	void init_bfys(size_t bf_num, const size_t *bf_pcls, const double *bfs);
	void init_txs(size_t t_num, const size_t *t_pcls, const double *ts);
	void init_tys(size_t t_num, const size_t *t_pcls, const double *ts);

	void init_fixed_vx_bc(size_t vx_bc_num, const size_t *vx_bcs);
	void init_fixed_vy_bc(size_t vy_bc_num, const size_t* vy_bcs);

protected:
	inline bool is_in_element(
		double pcl_x,
		double pcl_y,
		size_t elem_id,
		ShapeFunc &pcl_N
		) noexcept
	{
		const ShapeFuncAB& e_dN_ab = elem_dN_ab[elem_id];
		const ShapeFuncC& e_dN_c = elem_dN_c[elem_id];
		pcl_N.N1 = e_dN_ab.a1 * pcl_x + e_dN_ab.b1 * pcl_y + e_dN_c.c1;
		pcl_N.N2 = e_dN_ab.a2 * pcl_x + e_dN_ab.b2 * pcl_y + e_dN_c.c2;
		pcl_N.N3 = e_dN_ab.a3 * pcl_x + e_dN_ab.b3 * pcl_y + e_dN_c.c3;
		return pcl_N.N1 >= 0.0 && pcl_N.N1 <= 1.0
			&& pcl_N.N2 >= 0.0 && pcl_N.N2 <= 1.0
			&& pcl_N.N3 >= 0.0 && pcl_N.N3 <= 1.0;
	}
	
	inline bool is_in_element_tol(
		double pcl_x,
		double pcl_y,
		size_t elem_id,
		ShapeFunc& p_N
		) noexcept
	{
		const ShapeFuncAB& e_dN_ab = elem_dN_ab[elem_id];
		const ShapeFuncC& e_dN_c = elem_dN_c[elem_id];
		p_N.N1 = e_dN_ab.a1 * pcl_x + e_dN_ab.b1 * pcl_y + e_dN_c.c1;
		p_N.N2 = e_dN_ab.a2 * pcl_x + e_dN_ab.b2 * pcl_y + e_dN_c.c2;
		p_N.N3 = e_dN_ab.a3 * pcl_x + e_dN_ab.b3 * pcl_y + e_dN_c.c3;
#define in_elem_N_tol (1.0e-8)
		if (p_N.N1 >= -in_elem_N_tol && p_N.N1 <= (1.0 + in_elem_N_tol) &&
			p_N.N2 >= -in_elem_N_tol && p_N.N2 <= (1.0 + in_elem_N_tol) &&
			p_N.N3 >= -in_elem_N_tol && p_N.N3 <= (1.0 + in_elem_N_tol))
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
			return true;
		}
		return false;
#undef in_elem_N_tol
	}
	
	// background grid for mesh
	double grid_xl, grid_yl;
	double grid_xu, grid_yu;
	double grid_hx, grid_hy;
	size_t grid_x_num, grid_y_num;
	size_t* grid_elem_list_id_array;
	size_t* grid_elem_list;

public:
	inline size_t find_pcl_in_which_elem(
		double pcl_x,
		double pcl_y,
		ShapeFunc& pcl_N
		) noexcept
	{
		if (pcl_x < grid_xl || pcl_x > grid_xu ||
			pcl_y < grid_yl || pcl_y > grid_yu)
			return SIZE_MAX;
		const size_t x_id = (pcl_x - grid_xl) / grid_hx;
		const size_t y_id = (pcl_y - grid_yl) / grid_hy;
		const size_t g_id = grid_x_num * y_id + x_id;
		const size_t elem_end_id = grid_elem_list[g_id + 1];
		for (size_t el_id = grid_elem_list[g_id]; el_id < elem_end_id; ++el_id)
		{
			const size_t elem_id = grid_elem_list_id_array[el_id];
			if (is_in_element(pcl_x, pcl_y, elem_id, pcl_N))
				return elem_id;
		}
		return SIZE_MAX;
	}

	inline size_t find_pcl_in_which_elem_tol(
		double pcl_x,
		double pcl_y,
		ShapeFunc& pcl_N
		) noexcept
	{
		if (pcl_x < grid_xl || pcl_x > grid_xu ||
			pcl_y < grid_yl || pcl_y > grid_yu)
			return SIZE_MAX;
		const size_t x_id = (pcl_x - grid_xl) / grid_hx;
		const size_t y_id = (pcl_y - grid_yl) / grid_hy;
		const size_t g_id = grid_x_num * y_id + x_id;
		const size_t elem_end_id = grid_elem_list[g_id + 1];
		for (size_t el_id = grid_elem_list[g_id]; el_id < elem_end_id; ++el_id)
		{
			const size_t elem_id = grid_elem_list_id_array[el_id];
			if (is_in_element_tol(pcl_x, pcl_y, elem_id, pcl_N))
				return elem_id;
		}
		return SIZE_MAX;
	}

protected: // rigid object contact
	double K_cont;
	// rigid rect
	bool rigid_rect_is_valid;
	RigidRect rigid_rect;

public:
	inline bool has_rigid_rect() const noexcept { return rigid_rect_is_valid; }
	inline RigidRect& get_rigid_rect() { return rigid_rect; }
	inline void init_rigid_rect(double _K_cont,
		double x, double y, double hx, double hy, double density = 1.0)
	{
		rigid_rect_is_valid = true;
		K_cont = _K_cont;
		rigid_rect.init(x, y, hx, hy, density);
	}
	inline void set_rigid_rect_velocity(double vx, double vy, double v_ang)
	{
		rigid_rect.set_v_bc(vx, vy, v_ang);
	}

	friend class Model_T2D_ME_mt_hdf5_utilities::ParticleData;
	friend int Model_T2D_ME_mt_hdf5_utilities::output_background_mesh_to_hdf5_file(Model_T2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_ME_mt_hdf5_utilities::load_background_mesh_from_hdf5_file(Model_T2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_ME_mt_hdf5_utilities::output_boundary_condition_to_hdf5_file(Model_T2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_ME_mt_hdf5_utilities::load_boundary_condition_from_hdf5_file(Model_T2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_ME_mt_hdf5_utilities::output_pcl_data_to_hdf5_file(Model_T2D_ME_mt& md, Step_T2D_ME_mt& stp, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_ME_mt_hdf5_utilities::load_pcl_data_from_hdf5_file(Model_T2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_ME_mt_hdf5_utilities::output_material_model_to_hdf5_file(Model_T2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_ME_mt_hdf5_utilities::load_material_model_from_hdf5_file(Model_T2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_ME_mt_hdf5_utilities::output_rigid_rect_to_hdf5_file(Model_T2D_ME_mt& md, Step_T2D_ME_mt& stp, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_ME_mt_hdf5_utilities::output_rigid_rect_to_hdf5_file(Model_T2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_ME_mt_hdf5_utilities::load_rigid_rect_from_hdf5_file(Model_T2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
};

#endif