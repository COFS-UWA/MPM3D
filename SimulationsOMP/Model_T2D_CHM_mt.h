#ifndef __Model_T2D_CHM_mt_h__
#define __Model_T2D_CHM_mt_h__

#include "BCs.h"
#include "macro_utils.h"
#include "Model.h"
#include "MatModelContainer.h"
#include "ParticleGenerator2D.hpp"
#include "TriangleMesh.h"
#include "RigidBody/RigidCircle.h"

class Model_T2D_CHM_mt;
class Step_T2D_CHM_mt;
int substep_func_omp_T2D_CHM_mt(void* _self, size_t my_th_id,
	double dt, double cur_time, size_t substp_id);

class ResultFile_hdf5;
//namespace Model_T2D_CHM_mt_hdf5_utilities
//{
//	struct ParticleData;
//	int output_background_mesh_to_hdf5_file(Model_T2D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
//	int load_background_mesh_from_hdf5_file(Model_T2D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
//	int output_boundary_condition_to_hdf5_file(Model_T2D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
//	int load_boundary_condition_from_hdf5_file(Model_T2D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
//	int output_pcl_data_to_hdf5_file(Model_T2D_CHM_mt& md, Step_T2D_CHM_mt &stp, ResultFile_hdf5& rf, hid_t grp_id);
//	int load_pcl_data_from_hdf5_file(Model_T2D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
//	int output_material_model_to_hdf5_file(Model_T2D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
//	int load_material_model_from_hdf5_file(Model_T2D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
//	int output_rigid_rect_to_hdf5_file(Model_T2D_CHM_mt& md, Step_T2D_CHM_mt& stp, ResultFile_hdf5& rf, hid_t grp_id);
//	int output_rigid_rect_to_hdf5_file(Model_T2D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
//	int load_rigid_rect_from_hdf5_file(Model_T2D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
//}

struct Model_T2D_CHM_mt : public Model,
	public MatModel::MatModelContainer
{
	friend class Step_T2D_CHM_mt;
	friend int substep_func_omp_T2D_CHM_mt(void* _self,
		size_t my_th_id, double dt, double cur_time, size_t substp_id);

public:
	struct ShapeFunc { double N1, N2, N3; };
	struct DShapeFuncAB
	{
		union { double a1; double dN1_dx; };
		union { double b1; double dN1_dy; };
		union { double a2; double dN2_dx; };
		union { double b2; double dN2_dy; };
		union { double a3; double dN3_dx; };
		union { double b3; double dN3_dy; };
	};
	struct DShapeFuncC { double c1, c2, c3; };
	
	union Force
	{
		struct { double fx, fy; };
		Vector2D vec;
		inline Force() {}
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

	union Stress
	{
		struct { double s11, s22, s12; };
		double s[3];
	};
	union Strain
	{
		struct { double e11, e22, e12; };
		double e[3];
	};
	union StrainInc
	{
		struct { double de11, de22, de12; };
		double de[3];
	};

	struct ElemNodeIndex { size_t n1, n2, n3; };
	struct ElemNodeVM { double vm, vmx, vmy; };

	struct NodeHasVBC { bool has_vx_bc, has_vy_bc; };
	
	struct SortedPclVarArrays
	{
		size_t *pcl_index; // ori_pcl_num
		double* pcl_n;
		double *pcl_density_f; // ori_pcl_num
		Velocity *pcl_v; // ori_pcl_num
		Displacement *pcl_disp; // ori_pcl_num
		Stress* pcl_stress; // ori_pcl_num
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
	Force *pcl_bf; // ori_pcl_num
	Force *pcl_t; // ori_pcl_num
	Position *pcl_pos; // ori_pcl_num
	double *pcl_vol_s; // ori_pcl_num
	double *pcl_vol; // ori_pcl_num
	MatModel::MaterialModel **pcl_mat_model; // ori_pcl_num

	SortedPclVarArrays sorted_pcl_var_arrays[2];

	// element data
	ElemNodeIndex *elem_node_id; // elem_num
	//size_t* elem_id_array; // elem_num * 3 
	//size_t* node_elem_id_array; // elem_num * 3
	//size_t* node_elem_list;  // node_num
	DShapeFuncAB* elem_N_ab; // elem_num
	DShapeFuncC* elem_N_c; // elem_num
	double* elem_area; // elem_num
	// node data
	Position* node_pos; // node_num

	// element calculation data
	double *elem_density_f; // elem_num
	double *elem_pcl_n; // elem_num
	double *elem_pcl_m_s; // elem_num
	double *elem_pcl_m_f; // elem_num
	StrainInc *elem_de; // elem_num
	double *elem_p; // elem_num
	double *elem_m_de_vol_s; // elem_num
	double* elem_m_de_vol_f; // elem_num

	// element-node data
	ElemNodeVM* elem_node_vm_s; // elem_num * 3
	ElemNodeVM* elem_node_vm_f; // elem_num * 3
	Force* elem_node_force_s; // elem_num * 3
	Force* elem_node_force_f; // elem_num * 3
	
	Acceleration *node_a_s; // node_num
	Acceleration *node_a_f; // node_num
	Velocity *node_v_s; // node_num
	Velocity *node_v_f; // node_num
	NodeHasVBC *node_has_vbc_s;
	NodeHasVBC *node_has_vbc_f; // node_num
	double *node_am_s; // node_num
	double *node_am_f; // node_num
	double *node_de_vol_s; // node_num
	double * node_de_vol_f; // node_num

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
	Model_T2D_CHM_mt();
	~Model_T2D_CHM_mt();

	inline size_t get_ori_pcl_num() const noexcept { return ori_pcl_num; }
	inline size_t get_pcl_num() const noexcept { return pcl_num; }
	inline const Position *get_pcl_pos() const noexcept { return pcl_pos; }
	inline const double *get_pcl_m_s() const noexcept { return pcl_m_s; }
	inline size_t get_node_num() const noexcept { return node_num; }
	inline const Position *get_node_pos() const noexcept { return node_pos; }
	inline size_t get_elem_num() const noexcept { return elem_num; }
	inline const ElemNodeIndex* get_elem_node_index() const noexcept { return elem_node_id; }
	inline MatModel::MaterialModel **get_mat_models() noexcept { return pcl_mat_model; }
	Rect get_mesh_bbox();

	inline const NodeHasVBC *get_has_vbc_s() const noexcept { return node_has_vbc_s; }
	inline const NodeHasVBC* get_has_vbc_f() const noexcept { return node_has_vbc_f; }

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

	void init_fixed_vx_s_bc(size_t vx_bc_num, const size_t *vx_bcs);
	void init_fixed_vy_s_bc(size_t vy_bc_num, const size_t* vy_bcs);
	void init_fixed_vx_f_bc(size_t vx_bc_num, const size_t* vx_bcs);
	void init_fixed_vy_f_bc(size_t vy_bc_num, const size_t* vy_bcs);

protected:
	inline bool is_in_element(
		double pcl_x,
		double pcl_y,
		size_t elem_id,
		ShapeFunc &pcl_N
		) noexcept
	{
		DShapeFuncAB &e_sfab = elem_N_ab[elem_id];
		DShapeFuncC &e_sfc = elem_N_c[elem_id];
		pcl_N.N1 = e_sfab.a1 * pcl_x + e_sfab.b1 * pcl_y + e_sfc.c1;
		pcl_N.N2 = e_sfab.a2 * pcl_x + e_sfab.b2 * pcl_y + e_sfc.c2;
		pcl_N.N3 = e_sfab.a3 * pcl_x + e_sfab.b3 * pcl_y + e_sfc.c3;
		return pcl_N.N1 >= 0.0 && pcl_N.N1 <= 1.0
			&& pcl_N.N2 >= 0.0 && pcl_N.N2 <= 1.0
			&& pcl_N.N3 >= 0.0 && pcl_N.N3 <= 1.0;
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
		size_t x_id, y_id, elem_id;
		x_id = (pcl_x - grid_xl) / grid_hx;
		y_id = (pcl_y - grid_yl) / grid_hy;
		size_t g_id = grid_x_num * y_id + x_id;
		size_t elem_end_id = grid_elem_list[g_id + 1];
		for (size_t el_id = grid_elem_list[g_id];
			 el_id < elem_end_id; ++el_id)
		{
			elem_id = grid_elem_list_id_array[el_id];
			if (is_in_element(pcl_x, pcl_y, elem_id, pcl_N))
				return elem_id;
		}
		return SIZE_MAX;
	}

	//friend class Model_T2D_CHM_mt_hdf5_utilities::ParticleData;
	//friend int Model_T2D_CHM_mt_hdf5_utilities::output_background_mesh_to_hdf5_file(Model_T2D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	//friend int Model_T2D_CHM_mt_hdf5_utilities::load_background_mesh_from_hdf5_file(Model_T2D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	//friend int Model_T2D_CHM_mt_hdf5_utilities::output_boundary_condition_to_hdf5_file(Model_T2D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	//friend int Model_T2D_CHM_mt_hdf5_utilities::load_boundary_condition_from_hdf5_file(Model_T2D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	//friend int Model_T2D_CHM_mt_hdf5_utilities::output_pcl_data_to_hdf5_file(Model_T2D_CHM_mt& md, Step_T2D_CHM_mt& stp, ResultFile_hdf5& rf, hid_t grp_id);
	//friend int Model_T2D_CHM_mt_hdf5_utilities::load_pcl_data_from_hdf5_file(Model_T2D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	//friend int Model_T2D_CHM_mt_hdf5_utilities::output_material_model_to_hdf5_file(Model_T2D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	//friend int Model_T2D_CHM_mt_hdf5_utilities::load_material_model_from_hdf5_file(Model_T2D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	//friend int Model_T2D_CHM_mt_hdf5_utilities::output_rigid_rect_to_hdf5_file(Model_T2D_CHM_mt& md, Step_T2D_CHM_mt& stp, ResultFile_hdf5& rf, hid_t grp_id);
	//friend int Model_T2D_CHM_mt_hdf5_utilities::output_rigid_rect_to_hdf5_file(Model_T2D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	//friend int Model_T2D_CHM_mt_hdf5_utilities::load_rigid_rect_from_hdf5_file(Model_T2D_CHM_mt& md, ResultFile_hdf5& rf, hid_t grp_id);

protected: // rigid object contact
	double K_cont;
	// rigid rect
	bool rigid_circle_is_valid;
	RigidCircle rigid_circle;

	size_t* contact_state;
	Position* contact_pos;
	
public:
	inline bool has_rigid_circle() const noexcept { return rigid_circle_is_valid; }
	inline RigidCircle &get_rigid_circle() { return rigid_circle; }
	inline void init_rigid_circle(double _K_cont,
		double x, double y, double r, double density = 1.0)
	{
		rigid_circle_is_valid = true;
		K_cont = _K_cont;
		rigid_circle.init(x, y, r, density);
	}
	inline void set_rigid_circle_velocity(double vx, double vy, double v_ang)
	{ rigid_circle.set_v_bc(vx, vy, v_ang); }
};

#endif