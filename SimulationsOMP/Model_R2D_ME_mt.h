#ifndef __Model_R2D_ME_mt_h__
#define __Model_R2D_ME_mt_h__

#include <stdint.h>
#include <hdf5.h>

#include "BCs.h"
#include "macro_utils.h"
#include "Model.h"
#include "MatModelContainer.h"
#include "ParticleGenerator2D.hpp"
#include "RigidBody/RigidRect.h"

class Model_R2D_ME_mt;
class Step_R2D_ME_mt;
int substep_func_omp_R2D_ME_mt(void* _self, size_t my_th_id,
	double dt, double cur_time, size_t substp_id);

class ResultFile_hdf5;
class Step_R2D_ME_mt;
namespace Model_R2D_ME_mt_hdf5_utilities
{
	struct ParticleData;
	int output_background_mesh_to_hdf5_file(Model_R2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_background_mesh_from_hdf5_file(Model_R2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_boundary_condition_to_hdf5_file(Model_R2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_boundary_condition_from_hdf5_file(Model_R2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_pcl_data_to_hdf5_file(Model_R2D_ME_mt& md, Step_R2D_ME_mt& stp, ResultFile_hdf5& rf, hid_t grp_id);
	int load_pcl_data_from_hdf5_file(Model_R2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_material_model_to_hdf5_file(Model_R2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_material_model_from_hdf5_file(Model_R2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_rigid_rect_to_hdf5_file(Model_R2D_ME_mt& md, Step_R2D_ME_mt& stp, ResultFile_hdf5& rf, hid_t grp_id);
	int output_rigid_rect_to_hdf5_file(Model_R2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_rigid_rect_from_hdf5_file(Model_R2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
}

struct Model_R2D_ME_mt : public Model,
	public MatModel::MatModelContainer
{
	friend class Step_R2D_ME_mt;
	friend int substep_func_omp_R2D_ME_mt(void* _self,
		size_t my_th_id, double dt, double cur_time, size_t substp_id);

public:
	struct ShapeFunc { double N1, N2, N3, N4; };
	struct DShapeFunc
	{
		double dN1_dx, dN1_dy;
		double dN2_dx, dN2_dy;
		double dN3_dx, dN3_dy;
		double dN4_dx, dN4_dy;
	};
	struct Stress { double s11, s22, s12; };
	struct Strain { double e11, e22, e12; };
	struct StrainInc { double de11, de22, de12; };

	struct PclBodyForce { double bfx, bfy; };
	struct PclTraction { double tx, ty; };
	struct PclPos { double x, y; };

	struct PclDisp { double ux, uy; };
	struct PclV { double vx, vy; };

	struct ElemNodeVM { double vm, vmx, vmy; };
	struct ElemNodeForce { double fx, fy; };

	union NodeA
	{
		struct { double ax, ay; };
		struct { size_t ax_ui, ay_ui; };
	};
	union NodeV
	{
		struct { double vx, vy; };
		struct { size_t vx_ui, vy_ui; };
	};
	struct NodeHasVBC { bool has_vx_bc, has_vy_bc; };

	struct PclSortedVarArray
	{
		size_t* pcl_index; // ori_pcl_num
		double* pcl_density; // ori_pcl_num
		PclDisp* pcl_disp; // ori_pcl_num
		PclV* pcl_v; // ori_pcl_num
		Stress* pcl_stress; // ori_pcl_num
		Strain* pcl_strain; // ori_pcl_num
		Strain* pcl_estrain; // ori_pcl_num
		Strain* pcl_pstrain; // ori_pcl_num
	};

protected:
	size_t pcl_num;
	double* pcl_m; // ori_pcl_num
	PclBodyForce* pcl_bf; // ori_pcl_num
	PclTraction* pcl_t; // ori_pcl_num
	PclPos* pcl_pos; // ori_pcl_num
	double* pcl_vol; // ori_pcl_num
	ShapeFunc* pcl_N; // ori_pcl_num
	DShapeFunc* pcl_dN; // ori_pcl_num
	MatModel::MaterialModel** pcl_mat_model; // ori_pcl_num
	PclSortedVarArray pcl_sorted_var_array[2];

	// bg mesh data
	size_t elem_x_num, elem_y_num, elem_num;
	size_t node_x_num, node_y_num, node_num;
	size_t acutal_elem_x_num;
	double mh_xl, mh_yl, mh_xu, mh_yu;
	double elem_hx, elem_hy;
	double inv_elem_hx, inv_elem_hy, elem_area;
	double dxi_dx, deta_dy;
	DShapeFunc elem_dN;

	double* elem_density; // elem_num
	double* elem_pcl_m; // elem_num
	size_t* elem_substp_id; // elem_num

	ElemNodeVM* elem_node_vm; // elem_num * 3
	ElemNodeForce* elem_node_force; // elem_num * 3

	NodeA* node_a; // node_num
	NodeV* node_v; // node_num
	NodeHasVBC* node_has_vbc; // node_num

// ======== Non calculation data ========
	size_t ori_pcl_num;
	char* pcl_mem_raw;

	char* mesh_mem_raw;
	
	// boundary conditions
	size_t bfx_num, bfy_num;
	BodyForceAtPcl *bfxs, *bfys;
	size_t tx_num, ty_num;
	TractionBCAtPcl *txs, *tys;

public:
	Model_R2D_ME_mt();
	~Model_R2D_ME_mt();

	inline size_t get_ori_pcl_num() const noexcept { return ori_pcl_num; }
	inline size_t get_pcl_num() const noexcept { return pcl_num; }
	inline const PclPos* get_pcl_pos() const noexcept { return pcl_pos; }
	inline const double* get_pcl_m() const noexcept { return pcl_m; }
	inline const double* get_pcl_density0() const noexcept { return pcl_sorted_var_array[0].pcl_density; }
	inline size_t get_node_num() const noexcept { return node_num; }
	inline size_t get_elem_num() const noexcept { return elem_num; }
	inline MatModel::MaterialModel** get_mat_models() noexcept { return pcl_mat_model; }
	Rect get_mesh_bbox() { return Rect(mh_xl, mh_xu, mh_yl, mh_yu); }
	inline const NodeHasVBC* get_has_vbcs() const noexcept { return node_has_vbc; }

	void clear_mesh();
	void alloc_mesh(size_t n_num, size_t e_num);
	int init_mesh(double _xl, double _yl,
				  double _xu, double _yu,
				  size_t _x_num, size_t _y_num);

	void alloc_pcls(size_t num);
	void alloc_pcls(size_t num, size_t ori_num);
	void clear_pcls();
	int init_pcls(size_t num, double m, double density);
	void init_pcls(double _xl, double _yl, double _xn, double _yn,
				   size_t _x_num, size_t _y_num, double density);

	INIT_BC_TEMPLATE(bfx, BodyForceAtPcl)
	INIT_BC_TEMPLATE(bfy, BodyForceAtPcl)
	INIT_BC_TEMPLATE(tx, TractionBCAtPcl)
	INIT_BC_TEMPLATE(ty, TractionBCAtPcl)
	void init_bfxs(size_t bf_num, const size_t* bf_pcls, const double* bfs);
	void init_bfys(size_t bf_num, const size_t* bf_pcls, const double* bfs);
	void init_txs(size_t t_num, const size_t* t_pcls, const double* ts);
	void init_tys(size_t t_num, const size_t* t_pcls, const double* ts);

	void init_fixed_vx_bc(size_t vx_bc_num, const size_t* vx_bcs);
	void init_fixed_vy_bc(size_t vy_bc_num, const size_t* vy_bcs);

protected: // rigid object contact
	double K_cont;
	// rigid rect
	bool rigid_rect_is_valid;
	RigidRect rigid_rect;

	struct ContPos
	{
		double x, y;
	};
	size_t* contact_substep_id;
	ContPos* contact_pos;
	CacheAlignedMem contact_mem;

public:
	inline bool has_rigid_rect() const noexcept { return rigid_rect_is_valid; }
	inline RigidRect& get_rigid_rect() { return rigid_rect; }
	inline void init_rigid_rect(
		double _K_cont,
		double x, double y,
		double hx, double hy,
		double density = 1.0
		)
	{
		rigid_rect_is_valid = true;
		K_cont = _K_cont;
		rigid_rect.init(x, y, hx, hy, density);
	}
	inline void set_rigid_rect_velocity(double vx, double vy, double v_ang)
	{
		rigid_rect.set_v_bc(vx, vy, v_ang);
	}

protected:
	friend class Model_R2D_ME_mt_hdf5_utilities::ParticleData;
	friend int Model_R2D_ME_mt_hdf5_utilities::output_background_mesh_to_hdf5_file(Model_R2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_R2D_ME_mt_hdf5_utilities::load_background_mesh_from_hdf5_file(Model_R2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_R2D_ME_mt_hdf5_utilities::output_boundary_condition_to_hdf5_file(Model_R2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_R2D_ME_mt_hdf5_utilities::load_boundary_condition_from_hdf5_file(Model_R2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_R2D_ME_mt_hdf5_utilities::output_pcl_data_to_hdf5_file(Model_R2D_ME_mt& md, Step_R2D_ME_mt& stp, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_R2D_ME_mt_hdf5_utilities::load_pcl_data_from_hdf5_file(Model_R2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_R2D_ME_mt_hdf5_utilities::output_material_model_to_hdf5_file(Model_R2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_R2D_ME_mt_hdf5_utilities::load_material_model_from_hdf5_file(Model_R2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_R2D_ME_mt_hdf5_utilities::output_rigid_rect_to_hdf5_file(Model_R2D_ME_mt& md, Step_R2D_ME_mt& stp, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_R2D_ME_mt_hdf5_utilities::output_rigid_rect_to_hdf5_file(Model_R2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_R2D_ME_mt_hdf5_utilities::load_rigid_rect_from_hdf5_file(Model_R2D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
};

#endif