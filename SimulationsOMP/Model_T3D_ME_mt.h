#ifndef __Model_T3D_ME_mt_h__
#define __Model_T3D_ME_mt_h__

#include <stdint.h>

#include "BCs.h"
#include "macro_utils.h"
#include "Model.h"
#include "MatModelContainer.h"
#include "ParticleGenerator3D.hpp"
#include "TetrahedronMesh.h"

class Model_T3D_ME_mt;
class Step_T3D_ME_mt;
int substep_func_omp_T3D_ME_mt(void* _self, size_t my_th_id,
	double dt, double cur_time, size_t substp_id);

class ResultFile_hdf5;
class Step_T3D_ME_mt;
namespace Model_T3D_ME_mt_hdf5_utilities
{
	struct ParticleData;
	int output_background_mesh_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_background_mesh_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_boundary_condition_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_boundary_condition_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_pcl_data_to_hdf5_file(Model_T3D_ME_mt& md, Step_T3D_ME_mt &stp, ResultFile_hdf5& rf, hid_t grp_id);
	int load_pcl_data_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_material_model_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_material_model_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_rigid_rect_to_hdf5_file(Model_T3D_ME_mt& md, Step_T3D_ME_mt& stp, ResultFile_hdf5& rf, hid_t grp_id);
	int output_rigid_rect_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_rigid_rect_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
}

struct Model_T3D_ME_mt : public Model,
	public MatModel::MatModelContainer
{
	friend class Step_T3D_ME_mt;
	friend int substep_func_omp_T3D_ME_mt(void* _self,
		size_t my_th_id, double dt, double cur_time, size_t substp_id);

public:
	struct PclBodyForce { double bfx, bfy, bfz; };
	struct PclTraction { double tx, ty, tz; };
	struct PclPos { double x, y, z; };

	struct PclDisp { double ux, uy, uz; };
	struct PclV { double vx, vy, vz; };
	struct PclShapeFunc { double N1, N2, N3, N4; };
	struct PclStress { double s11, s22, s33, s12, s23, s31; };

	struct ElemNodeIndex { size_t n1, n2, n3, n4; };
	struct ElemShapeFuncABC
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
	struct ElemShapeFuncD { double d1, d2, d3, d4; };

	union ElemStrainInc
	{
		struct { double de11, de22, de33, de12, de23, de31; };
		double de[6];
	};
	struct ElemStress { double s11, s22, s33, s12, s23, s31; };

	struct ElemNodeVM { double vm, vmx, vmy, vmz; };
	struct ElemNodeForce { double fx, fy, fz; };

	union NodeA
	{
		struct { double ax, ay, az; };
		struct { size_t ax_ui, ay_ui, az_ui; };
	};
	union NodeV
	{
		struct { double vx, vy, vz; };
		struct { size_t vx_ui, vy_ui, vz_ui; };
	};
	struct NodePos { double x, y, z; };
	struct NodeHasVBC { bool has_vx_bc, has_vy_bc, has_vz_bc; };

	struct PclSortedVarArray
	{
		size_t* pcl_index; // ori_pcl_num
		double* pcl_density; // ori_pcl_num
		PclDisp* pcl_disp; // ori_pcl_num
		PclV* pcl_v; // ori_pcl_num
		PclShapeFunc* pcl_N; // ori_pcl_num
		PclStress* pcl_stress; // ori_pcl_num
	};
	
protected:
	size_t pcl_num;
	size_t elem_num;
	size_t node_num;
	
	double *pcl_m; // ori_pcl_num
	PclBodyForce* pcl_bf; // ori_pcl_num
	PclTraction* pcl_t; // ori_pcl_num
	PclPos* pcl_pos; // ori_pcl_num
	double* pcl_vol; // ori_pcl_num
	MatModel::MaterialModel **pcl_mat_model; // ori_pcl_num

	PclSortedVarArray pcl_sorted_var_array[2];

	// element data
	ElemNodeIndex *elem_node_id; // elem_num
	double *elem_vol; // elem_num
	ElemShapeFuncABC* elem_sf_abc; // elem_num
	ElemShapeFuncD* elem_sf_d; // elem_num

	// element calculation data
	double *elem_density; // elem_num
	double *elem_pcl_m; // elem_num
	double *elem_pcl_vol; // elem_num
	ElemStrainInc *elem_de; // elem_num
	ElemStress *elem_stress; // elem_num
	double *elem_m_de_vol; // elem_num

	// element-node data
	ElemNodeVM* elem_node_vm; // elem_num * 4
	ElemNodeForce* elem_node_force; // elem_num * 4
	
	size_t *elem_id_array; // elem_num * 4 
	size_t *node_elem_id_array; // elem_num * 4
	size_t *node_elem_list;  // node_num

	NodeA *node_a; // node_num
	NodeV *node_v; // node_num
	NodeHasVBC *node_has_vbc; // node_num
	double *node_am; // node_num
	double *node_de_vol; // node_num

// ======== Non calculation data ========
	size_t ori_pcl_num;
	char *pcl_mem_raw;

	char *mesh_mem_raw;

	NodePos *node_pos;

	// boundary conditions
	size_t bfx_num, bfy_num, bfz_num;
	BodyForceAtPcl *bfxs, *bfys, *bfzs;
	size_t tx_num, ty_num, tz_num;
	TractionBCAtPcl *txs, *tys, *tzs;

public:
	Model_T3D_ME_mt();
	~Model_T3D_ME_mt();

	inline size_t get_ori_pcl_num() const noexcept { return ori_pcl_num; }
	inline size_t get_pcl_num() const noexcept { return pcl_num; }
	inline const PclPos* get_pcl_pos() const noexcept { return pcl_pos; }
	inline const double *get_pcl_m() const noexcept { return pcl_m; }
	inline const double* get_pcl_density0() const noexcept { return pcl_sorted_var_array[0].pcl_density; }
	inline size_t get_node_num() const noexcept { return node_num; }
	inline const NodePos* get_node_pos() const noexcept { return node_pos; }
	inline size_t get_elem_num() const noexcept { return elem_num; }
	inline const ElemNodeIndex* get_elem_node_index() const noexcept { return elem_node_id; }
	inline MatModel::MaterialModel **get_mat_models() noexcept { return pcl_mat_model; }
	Cube get_mesh_bbox();

	inline const NodeHasVBC *get_has_vbcs() const noexcept { return node_has_vbc; }

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
	
	void clear_mesh();
	void alloc_mesh(size_t n_num, size_t e_num);
	void init_mesh(const TetrahedronMesh &mesh);
	
	void clear_search_grid();
	int init_search_grid(TetrahedronMesh &mesh, double _hx, double _hy, double _hz);
	
	void alloc_pcls(size_t num);
	void alloc_pcls(size_t num, size_t ori_num);
	void clear_pcls();
	int init_pcls(size_t num, double m, double density);
	int init_pcls(ParticleGenerator3D<TetrahedronMesh>& pg, double density);

	INIT_BC_TEMPLATE(bfx, BodyForceAtPcl)
	INIT_BC_TEMPLATE(bfy, BodyForceAtPcl)
	INIT_BC_TEMPLATE(bfz, BodyForceAtPcl)
	INIT_BC_TEMPLATE(tx, TractionBCAtPcl)
	INIT_BC_TEMPLATE(ty, TractionBCAtPcl)
	INIT_BC_TEMPLATE(tz, TractionBCAtPcl)
	void init_bfxs(size_t bf_num, const size_t *bf_pcls, const double *bfs);
	void init_bfys(size_t bf_num, const size_t *bf_pcls, const double *bfs);
	void init_bfzs(size_t bf_num, const size_t *bf_pcls, const double *bfs);
	void init_txs(size_t t_num, const size_t *t_pcls, const double *ts);
	void init_tys(size_t t_num, const size_t *t_pcls, const double *ts);
	void init_tzs(size_t t_num, const size_t *t_pcls, const double *ts);

	void init_fixed_vx_bc(size_t vx_bc_num, const size_t *vx_bcs);
	void init_fixed_vy_bc(size_t vy_bc_num, const size_t *vy_bcs);
	void init_fixed_vz_bc(size_t vz_bc_num, const size_t *vz_bcs);

protected:
	inline bool is_in_element(
		double pcl_x,
		double pcl_y,
		double pcl_z,
		size_t elem_id,
		PclShapeFunc &pcl_sf
		) noexcept
	{
		ElemShapeFuncABC& e_sfabc = elem_sf_abc[elem_id];
		ElemShapeFuncD& e_sfd = elem_sf_d[elem_id];
		pcl_sf.N1 = e_sfabc.a1 * pcl_x + e_sfabc.b1 * pcl_y + e_sfabc.c1 * pcl_z + e_sfd.d1;
		pcl_sf.N2 = e_sfabc.a2 * pcl_x + e_sfabc.b2 * pcl_y + e_sfabc.c2 * pcl_z + e_sfd.d2;
		pcl_sf.N3 = e_sfabc.a3 * pcl_x + e_sfabc.b3 * pcl_y + e_sfabc.c3 * pcl_z + e_sfd.d3;
		pcl_sf.N4 = e_sfabc.a4 * pcl_x + e_sfabc.b4 * pcl_y + e_sfabc.c4 * pcl_z + e_sfd.d4;
		// for numerical accuracy
		return pcl_sf.N1 >= 0.0 && pcl_sf.N1 <= 1.0
			&& pcl_sf.N2 >= 0.0 && pcl_sf.N2 <= 1.0
			&& pcl_sf.N3 >= 0.0 && pcl_sf.N3 <= 1.0
			&& pcl_sf.N4 >= 0.0 && pcl_sf.N4 <= 1.0;
	}
	
	// background grid for mesh
	double grid_xl, grid_yl, grid_zl;
	double grid_xu, grid_yu, grid_zu;
	double grid_hx, grid_hy, grid_hz;
	size_t grid_x_num, grid_y_num, grid_z_num;
	size_t grid_xy_num;
	size_t *grid_elem_list_id_array;
	size_t *grid_elem_list;

public:
	inline size_t find_pcl_in_which_elem(
		double pcl_x,
		double pcl_y,
		double pcl_z,
		PclShapeFunc& pcl_sf
		) noexcept
	{
		if (pcl_x < grid_xl || pcl_x > grid_xu ||
			pcl_y < grid_yl || pcl_y > grid_yu ||
			pcl_z < grid_zl || pcl_z > grid_zu)
			return elem_num;
		size_t x_id, y_id, z_id, elem_id;
		x_id = (pcl_x - grid_xl) / grid_hx;
		y_id = (pcl_y - grid_yl) / grid_hy;
		z_id = (pcl_z - grid_zl) / grid_hz;
		size_t g_id = grid_xy_num * z_id + grid_x_num * y_id + x_id;
		size_t elem_end_id = grid_elem_list[g_id + 1];
		for (size_t el_id = grid_elem_list[g_id];
			 el_id < elem_end_id; ++el_id)
		{
			elem_id = grid_elem_list_id_array[el_id];
			if (is_in_element(pcl_x, pcl_y, pcl_z, elem_id, pcl_sf))
				return elem_id;
		}
		return elem_num;
	}

protected:
	double K_cont;

	friend class Model_T3D_ME_mt_hdf5_utilities::ParticleData;
	friend int Model_T3D_ME_mt_hdf5_utilities::output_background_mesh_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::load_background_mesh_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::output_boundary_condition_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::load_boundary_condition_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::output_pcl_data_to_hdf5_file(Model_T3D_ME_mt& md, Step_T3D_ME_mt& stp, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::load_pcl_data_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::output_material_model_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::load_material_model_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::output_rigid_rect_to_hdf5_file(Model_T3D_ME_mt& md, Step_T3D_ME_mt& stp, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::output_rigid_rect_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::load_rigid_rect_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
};

#endif