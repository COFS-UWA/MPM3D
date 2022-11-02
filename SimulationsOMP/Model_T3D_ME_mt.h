#ifndef __Model_T3D_ME_mt_h__
#define __Model_T3D_ME_mt_h__

#include "BCs.h"
#include "macro_utils.h"
#include "Model.h"
#include "MatModelContainer.h"
#include "ParticleGenerator3D.hpp"
#include "TetrahedronMesh.h"
#include "RigidObject/RigidCylinder.h"
#include "RigidObject/RigidCone.h"
#include "RigidObject/RigidCube.h"
#include "RigidObject/RigidObjectByT3DMesh.h"
#include "RigidObject/SmoothContact3D.h"
#include "RigidObject/SmoothContact3D_Quad.h"
#include "RigidObject/RoughContact3D.h"
#include "RigidObject/FrictionalContact3D.h"
#include "RigidObject/StickyContact3D.h"
#include "RigidObject/NonSeperateContact3D.h"

class Model_T3D_ME_mt;
class Step_T3D_ME_mt;
int substep_func_omp_T3D_ME_mt(void* _self, size_t my_th_id,
	double dt, double cur_time, size_t substp_id);
class Step_T3D_ME_mt_Geo;
int substep_func_omp_T3D_ME_mt_Geo(void* _self,
	size_t my_th_id, double dt, double cur_time, size_t substp_id);

class Step_T3D_ME_TBB;
int cal_substep_func_T3D_ME_TBB(void* _self);
namespace Step_T3D_ME_Task { class CalData; }

class ResultFile_hdf5;
namespace Model_T3D_ME_mt_hdf5_utilities
{
	struct ParticleData;
	int output_background_mesh_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_background_mesh_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_boundary_condition_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_boundary_condition_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_pcl_data_to_hdf5_file(Model_T3D_ME_mt& md, Step_T3D_ME_mt& stp, ResultFile_hdf5& rf, hid_t grp_id);
	int load_pcl_data_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_material_model_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_material_model_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_rigid_cylinder_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_rigid_cylinder_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_rigid_cone_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_rigid_cone_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_rigid_cube_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_rigid_cube_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_t3d_rigid_mesh_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_t3d_rigid_mesh_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
}

class PclVar_T3D_ME_mt;

struct Model_T3D_ME_mt : public Model,
	public MatModel::MatModelContainer
{
	friend class Step_T3D_ME_mt;
	friend int substep_func_omp_T3D_ME_mt(void* _self,
		size_t my_th_id, double dt, double cur_time, size_t substp_id);
	
	friend class Step_T3D_ME_mt_Geo;
	friend int substep_func_omp_T3D_ME_mt_Geo(void* _self,
		size_t my_th_id, double dt, double cur_time, size_t substp_id);

	friend class Step_T3D_ME_TBB;
	friend int cal_substep_func_T3D_ME_TBB(void* _self);
	friend class Step_T3D_ME_Task::CalData;

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
		size_t* pcl_index; // ori_pcl_num
		double* pcl_density; // ori_pcl_num
		Velocity* pcl_v; // ori_pcl_num
		Displacement* pcl_disp; // ori_pcl_num
		Stress* pcl_stress; // ori_pcl_num
		Strain* pcl_strain; // ori_pcl_num
		Strain* pcl_estrain; // ori_pcl_num
		Strain* pcl_pstrain; // ori_pcl_num
		ShapeFunc* pcl_N; // ori_pcl_num
	};

protected:
	// pcl data
	size_t pcl_num;

	double* pcl_m; // ori_pcl_num
	Force* pcl_bf; // ori_pcl_num
	Force* pcl_t; // ori_pcl_num
	Position* pcl_pos; // ori_pcl_num
	double* pcl_vol; // ori_pcl_num
	MatModel::MaterialModel** pcl_mat_model; // ori_pcl_num

	SortedPclVarArrays sorted_pcl_var_arrays[2];

	// mesh data
	size_t elem_num;
	size_t node_num;

	// mesh topology
	ElemNodeIndex* elem_node_id; // elem_num
	size_t* elem_id_array; // elem_num * 4 
	size_t* node_elem_id_array; // elem_num * 4
	size_t* node_elem_list;  // node_num + 1

	// element data
	DShapeFuncABC* elem_dN_abc; // elem_num
	DShapeFuncD* elem_dN_d; // elem_num
	double* elem_vol; // elem_num
	
	// element calculation data
	double* elem_pcl_m; // elem_num
	double* elem_density; // elem_num
	StrainInc* elem_de; // elem_num
	double* elem_m_de_vol; // elem_num

	// element-node calculation data
	ElemNodeVM* elem_node_vm; // elem_num * 4
	Force* elem_node_force; // elem_num * 4

	// node calculation data
	Position* node_pos; // node_num
	Acceleration* node_a; // node_num
	Velocity* node_v; // node_num
	NodeHasVBC* node_has_vbc; // node_num
	NodeVBCVec *node_vbc_vec;
	double* node_am; // node_num
	double* node_de_vol; // node_num

// ======== Non calculation data ========
	size_t ori_pcl_num;
	char* pcl_mem_raw;

	char* mesh_mem_raw;

	// boundary conditions
	size_t bfx_num, bfy_num, bfz_num;
	BodyForceAtPcl* bfxs, * bfys, * bfzs;
	size_t tx_num, ty_num, tz_num;
	TractionBCAtPcl* txs, * tys, * tzs;

public:
	Model_T3D_ME_mt();
	~Model_T3D_ME_mt();

	inline size_t get_ori_pcl_num() const noexcept { return ori_pcl_num; }
	inline size_t get_pcl_num() const noexcept { return pcl_num; }
	inline const Position* get_pcl_pos() const noexcept { return pcl_pos; }
	inline const double* get_pcl_m() const noexcept { return pcl_m; }
	inline const double* get_pcl_density0() const noexcept { return sorted_pcl_var_arrays[0].pcl_density; }
	inline Stress* get_pcl_stress0() noexcept { return sorted_pcl_var_arrays[0].pcl_stress; }
	inline size_t get_node_num() const noexcept { return node_num; }
	inline const Position* get_node_pos() const noexcept { return node_pos; }
	inline size_t get_elem_num() const noexcept { return elem_num; }
	inline const ElemNodeIndex* get_elem_node_index() const noexcept { return elem_node_id; }
	inline MatModel::MaterialModel** get_mat_models() noexcept { return pcl_mat_model; }
	Cube get_mesh_bbox();

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
	inline const NodeHasVBC* get_node_has_vbc() const noexcept { return node_has_vbc; }
	inline const NodeVBCVec *get_node_vbc_vec() const noexcept { return node_vbc_vec; }

	void clear_mesh();
	void alloc_mesh(size_t n_num, size_t e_num);
	void init_mesh(const TetrahedronMesh& mesh);

	void clear_search_grid();
	int init_search_grid(TetrahedronMesh& mesh);

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
	void init_bfxs(size_t bf_num, const size_t* bf_pcls, const double* bfs);
	void init_bfys(size_t bf_num, const size_t* bf_pcls, const double* bfs);
	void init_bfzs(size_t bf_num, const size_t* bf_pcls, const double* bfs);
	void init_txs(size_t t_num, const size_t* t_pcls, const double* ts);
	void init_tys(size_t t_num, const size_t* t_pcls, const double* ts);
	void init_tzs(size_t t_num, const size_t* t_pcls, const double* ts);

	void init_fixed_vx_bc(size_t vx_bc_num, const size_t* vx_bcs);
	void init_fixed_vy_bc(size_t vy_bc_num, const size_t* vy_bcs);
	void init_fixed_vz_bc(size_t vz_bc_num, const size_t* vz_bcs);
	void set_vbc_vec(size_t n_id, double vecx, double vecy, double vecz);

	inline void add_mat_model(size_t pcl_id,
		MatModel::MaterialModel& mat_model,
		size_t model_size)
	{
		pcl_mat_model[pcl_id] = &mat_model;
		//pcl_mat_model_copy_offset[pcl_id] = pcl_mat_model_total_size;
		//pcl_mat_model_total_size += model_size;
	}

protected:
public:
	inline bool is_in_element(
		double pcl_x,
		double pcl_y,
		double pcl_z,
		size_t elem_id,
		ShapeFunc& p_N
		) noexcept
	{
		const DShapeFuncABC& e_dN_abc = elem_dN_abc[elem_id];
		const DShapeFuncD& e_dN_d = elem_dN_d[elem_id];
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
		) noexcept
	{
		const DShapeFuncABC& e_dN_abc = elem_dN_abc[elem_id];
		const DShapeFuncD& e_dN_d = elem_dN_d[elem_id];
		p_N.N1 = e_dN_abc.a1 * pcl_x + e_dN_abc.b1 * pcl_y + e_dN_abc.c1 * pcl_z + e_dN_d.d1;
		p_N.N2 = e_dN_abc.a2 * pcl_x + e_dN_abc.b2 * pcl_y + e_dN_abc.c2 * pcl_z + e_dN_d.d2;
		p_N.N3 = e_dN_abc.a3 * pcl_x + e_dN_abc.b3 * pcl_y + e_dN_abc.c3 * pcl_z + e_dN_d.d3;
		p_N.N4 = e_dN_abc.a4 * pcl_x + e_dN_abc.b4 * pcl_y + e_dN_abc.c4 * pcl_z + e_dN_d.d4;
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

	// background grid for mesh
	double grid_xl, grid_yl, grid_zl;
	double grid_xu, grid_yu, grid_zu;
	double grid_hx, grid_hy, grid_hz;
	size_t grid_x_num, grid_y_num, grid_z_num;
	size_t grid_xy_num;
	size_t* grid_elem_list_id_array;
	size_t* grid_elem_list;

	inline size_t find_pcl_in_which_elem(
		double pcl_x,
		double pcl_y,
		double pcl_z,
		ShapeFunc& p_N
	) noexcept
	{
		if (pcl_x < grid_xl || pcl_x > grid_xu ||
			pcl_y < grid_yl || pcl_y > grid_yu ||
			pcl_z < grid_zl || pcl_z > grid_zu)
			return SIZE_MAX;
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
	) noexcept
	{
		if (pcl_x < grid_xl || pcl_x > grid_xu ||
			pcl_y < grid_yl || pcl_y > grid_yu ||
			pcl_z < grid_zl || pcl_z > grid_zu)
			return SIZE_MAX;
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
			if (is_in_element_tol(pcl_x, pcl_y, pcl_z, elem_id, p_N))
				return elem_id;
		}
		return SIZE_MAX;
	}

protected:
	bool rigid_cylinder_is_valid;
	bool rigid_cone_is_valid;
	bool rigid_cube_is_valid;
	bool rigid_t3d_mesh_is_valid;
	RigidCylinder rigid_cylinder;
	RigidCube rigid_cube;
	RigidObjectByT3DMesh rigid_t3d_mesh;
	RigidCone rigid_cone;

	// ad hoc design for output
	double Kn_cont, Kt_cont;
	double fric_ratio, shear_strength;
	SmoothContact3D smooth_contact;
	SmoothContact3D_Quad smh_contact_quad;
	RoughContact3D rough_contact;
	FrictionalContact3D fric_contact;
	StickyContact3D sticky_contact;
	ContactModel3D* pcm;

	size_t *contact_substep_ids; // ori_pcl_num
	ContactModel3D::Position *prev_contact_poses; // ori_pcl_num
	ContactModel3D::Force *prev_contact_tan_forces; // ori_pcl_num
	double *prev_contact_dists; // ori_pcl_num
	char *contact_mem;
	void clear_contact_mem();
	void alloc_contact_mem(size_t num);

public:
	// cylinder
	bool has_rigid_cylinder() const noexcept { return rigid_cylinder_is_valid; }
	RigidCylinder& get_rigid_cylinder() noexcept { return rigid_cylinder; }
	void clear_rigid_cylinder() noexcept { rigid_cylinder_is_valid = false; }
	inline void init_rigid_cylinder(
		double x, double y, double z,
		double h, double r, double den = 1.0)
	{
		rigid_cylinder_is_valid = true;
		rigid_cylinder.init(x, y, z, h, r, den);
	}
	inline void fix_rigid_cylinder_vx() { rigid_cylinder.set_vx_bc(0.0); }
	inline void fix_rigid_cylinder_vy() { rigid_cylinder.set_vy_bc(0.0); }
	inline void set_rigid_cylinder_velocity(
		double vx, double vy, double vz)
	{ rigid_cylinder.set_vbc(vx, vy, vz); }
	inline void set_rigid_cylinder_init_velocity(
		double vx, double vy, double vz)
	{ rigid_cylinder.set_init_v(vx, vy, vz); }
	inline void set_rigid_cylinder_ext_force(
		double fx, double fy, double fz)
	{ rigid_cylinder.set_ext_force(fx, fy, fz, 0.0, 0.0, 0.0); }
	inline double get_rigid_cylinder_m() { return rigid_cylinder.get_m(); }
	// cone
	bool has_rigid_cone() const noexcept { return rigid_cone_is_valid; }
	RigidCone& get_rigid_cone() { return rigid_cone; }
	inline void init_rigid_cone(
		double x, double y, double z,
		double h_tip, double h_shaft, double r)
	{
		rigid_cone_is_valid = true;
		rigid_cone.init(x, y, z, r, h_tip, h_shaft);
	}
	inline void set_rigid_cone_velocity(
		double vx, double vy, double vz)
	{
		rigid_cone.set_vbc(vx, vy, vz);
	}

	// cube
	bool has_rigid_cube() const noexcept { return rigid_cube_is_valid; }
	RigidCube& get_rigid_cube() { return rigid_cube; }
	inline void init_rigid_cube(
		double x, double y, double z,
		double hx, double hy, double hz, double density)
	{
		rigid_cube_is_valid = true;
		rigid_cube.init(x, y, z, hx, hy, hz, density);
	}
	inline void set_rigid_cube_force(double fx, double fy, double fz)
	{
		rigid_cube.set_fx_ext(fx);
		rigid_cube.set_fy_ext(fy);
		rigid_cube.set_fz_ext(fz);
	}

	// t3d rigid mesh
	bool has_t3d_rigid_mesh() const noexcept { return rigid_t3d_mesh_is_valid; }
	RigidObjectByT3DMesh &get_t3d_rigid_mesh() noexcept { return rigid_t3d_mesh; }
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
	{ rigid_t3d_mesh.set_translation_velocity_bc(vx, vy, vz);	}
	inline void set_t3d_rigid_mesh_init_velocity(double vx, double vy, double vz)
	{ rigid_t3d_mesh.set_velocity(vx, vy, vz); }
	inline void set_t3d_rigid_mesh_ext_force(double fx, double fy, double fz)
	{ rigid_t3d_mesh.set_ext_force(fx, fy, fz); }

	// for contact model
	// for contact model
	inline void set_contact_param(
		double _Kn_cont,
		double _Kt_cont,
		double _fric_ratio,
		double _shear_strength)
	{
		Kn_cont = _Kn_cont;
		Kt_cont = _Kt_cont;
		fric_ratio = _fric_ratio;
		shear_strength = _shear_strength;
		smooth_contact.set_Kn_cont(_Kn_cont);
		smh_contact_quad.set_Kn_cont(_Kn_cont);
		fric_contact.set_K_cont(_Kn_cont, _Kt_cont);
		fric_contact.set_friction_ratio(_fric_ratio);
		sticky_contact.set_K_cont(_Kn_cont, _Kt_cont);
		sticky_contact.set_shear_strength(_shear_strength);
		rough_contact.set_K_cont(_Kn_cont, _Kt_cont);
	}
	inline void set_smooth_contact_between_pcl_and_rect() noexcept { pcm = &smooth_contact; }
	inline void set_smooth_quad_contact_between_pcl_and_rect() noexcept { pcm = &smh_contact_quad; }
	inline void set_rough_contact_between_pcl_and_rect() noexcept { pcm = &rough_contact; }
	inline void set_frictional_contact_between_pcl_and_rect() noexcept { pcm = &fric_contact; }
	inline void set_sticky_contact_between_pcl_and_rect() noexcept { pcm = &sticky_contact; }
	inline ContactModel3D *get_contact_model() noexcept { return pcm; }
	inline bool is_smooth_contact_between_pcl_and_rect() noexcept { return pcm == &smooth_contact; }
	inline bool is_smooth_quad_contact_between_pcl_and_rect() noexcept { return pcm == &smh_contact_quad; }
	inline bool is_frictional_contact_between_pcl_and_rect() noexcept { return pcm == &fric_contact; }
	inline bool is_sticky_contact_between_pcl_and_rect() noexcept { return pcm == &sticky_contact; }
	inline bool is_rough_contact_between_pcl_and_rect() noexcept { return pcm == &rough_contact; }

	inline void set_cylinder_vx_bc_ramp_up_time(double ramp_up_time) { rigid_cylinder.set_ramp_up_vx_bc(ramp_up_time); }
	inline void set_cylinder_vy_bc_ramp_up_time(double ramp_up_time) { rigid_cylinder.set_ramp_up_vy_bc(ramp_up_time); }
	inline void set_cylinder_vz_bc_ramp_up_time(double ramp_up_time) { rigid_cylinder.set_ramp_up_vz_bc(ramp_up_time); }
	
	friend class Model_T3D_ME_mt_hdf5_utilities::ParticleData;
	friend int Model_T3D_ME_mt_hdf5_utilities::output_background_mesh_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::load_background_mesh_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::output_boundary_condition_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::load_boundary_condition_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::output_pcl_data_to_hdf5_file(Model_T3D_ME_mt& md, Step_T3D_ME_mt& stp, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::load_pcl_data_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::output_material_model_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::load_material_model_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::output_rigid_cylinder_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::load_rigid_cylinder_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::output_rigid_cone_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::load_rigid_cone_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::output_rigid_cube_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::load_rigid_cube_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::output_t3d_rigid_mesh_to_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T3D_ME_mt_hdf5_utilities::load_t3d_rigid_mesh_from_hdf5_file(Model_T3D_ME_mt& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend class PclVar_T3D_ME_mt;
};

#include "ParticleVariablesGetter.h"

class PclVar_T3D_ME_mt : public ParticleVariablesGetter
{
public:
	Model_T3D_ME_mt* pmodel;
	Model_T3D_ME_mt::SortedPclVarArrays* cur_sorted_pcl_vars;
	size_t pcl_id;

	size_t get_index() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_index[pcl_id];
	}
	double get_m() const noexcept override
	{
		return pmodel->pcl_m[cur_sorted_pcl_vars->pcl_index[pcl_id]];
	}
	double get_bfx() const noexcept override
	{
		return pmodel->pcl_bf[cur_sorted_pcl_vars->pcl_index[pcl_id]].fx;
	}
	double get_bfy() const noexcept override
	{
		return pmodel->pcl_bf[cur_sorted_pcl_vars->pcl_index[pcl_id]].fy;
	}
	double get_bfz() const noexcept override
	{
		return pmodel->pcl_bf[cur_sorted_pcl_vars->pcl_index[pcl_id]].fz;
	}
	double get_tx() const noexcept override
	{
		return pmodel->pcl_t[cur_sorted_pcl_vars->pcl_index[pcl_id]].fx;
	}
	double get_ty() const noexcept override
	{
		return pmodel->pcl_t[cur_sorted_pcl_vars->pcl_index[pcl_id]].fy;
	}
	double get_tz() const noexcept override
	{
		return pmodel->pcl_t[cur_sorted_pcl_vars->pcl_index[pcl_id]].fz;
	}
	double get_x() const noexcept override
	{
		return pmodel->pcl_pos[cur_sorted_pcl_vars->pcl_index[pcl_id]].x;
	}
	double get_y() const noexcept override
	{
		return pmodel->pcl_pos[cur_sorted_pcl_vars->pcl_index[pcl_id]].y;
	}
	double get_z() const noexcept override
	{
		return pmodel->pcl_pos[cur_sorted_pcl_vars->pcl_index[pcl_id]].z;
	}
	double get_vol() const noexcept override
	{
		return pmodel->pcl_vol[cur_sorted_pcl_vars->pcl_index[pcl_id]];
	}
	double get_square_r() const noexcept override
	{
		double p_vol = pmodel->pcl_vol[cur_sorted_pcl_vars->pcl_index[pcl_id]];
		return 0.5 * pow(p_vol, 1.0 / 3.0);
	}
	double get_circle_r() const noexcept override
	{
		double p_vol = pmodel->pcl_vol[cur_sorted_pcl_vars->pcl_index[pcl_id]];
		return pow(0.75 / 3.14159265359 * p_vol, 1.0 / 3.0);
	}
	MatModel::MaterialModel* get_mat_model() const noexcept override
	{
		return pmodel->pcl_mat_model[cur_sorted_pcl_vars->pcl_index[pcl_id]];
	}
	double get_density() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_density[pcl_id];
	}
	double get_vx() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_v[pcl_id].vx;
	}
	double get_vy() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_v[pcl_id].vy;
	}
	double get_vz() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_v[pcl_id].vz;
	}
	double get_ux() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_disp[pcl_id].ux;
	}
	double get_uy() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_disp[pcl_id].uy;
	}
	double get_uz() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_disp[pcl_id].uz;
	}
	double get_s11() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_stress[pcl_id].s11;
	}
	double get_s22() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_stress[pcl_id].s22;
	}
	double get_s33() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_stress[pcl_id].s33;
	}
	double get_s12() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_stress[pcl_id].s12;
	}
	double get_s23() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_stress[pcl_id].s23;
	}
	double get_s31() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_stress[pcl_id].s31;
	}
	double get_e11() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_strain[pcl_id].e11;
	}
	double get_e22() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_strain[pcl_id].e22;
	}
	double get_e33() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_strain[pcl_id].e33;
	}
	double get_e12() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_strain[pcl_id].e12;
	}
	double get_e23() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_strain[pcl_id].e23;
	}
	double get_e31() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_strain[pcl_id].e31;
	}
	double get_ee11() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_estrain[pcl_id].e11;
	}
	double get_ee22() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_estrain[pcl_id].e22;
	}
	double get_ee33() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_estrain[pcl_id].e33;
	}
	double get_ee12() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_estrain[pcl_id].e12;
	}
	double get_ee23() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_estrain[pcl_id].e23;
	}
	double get_ee31() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_estrain[pcl_id].e31;
	}
	double get_pe11() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_pstrain[pcl_id].e11;
	}
	double get_pe22() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_pstrain[pcl_id].e22;
	}
	double get_pe33() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_pstrain[pcl_id].e33;
	}
	double get_pe12() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_pstrain[pcl_id].e12;
	}
	double get_pe23() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_pstrain[pcl_id].e23;
	}
	double get_pe31() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_pstrain[pcl_id].e31;
	}
	double get_N1() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_N[pcl_id].N1;
	}
	double get_N2() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_N[pcl_id].N2;
	}
	double get_N3() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_N[pcl_id].N3;
	}
	double get_N4() const noexcept override
	{
		return cur_sorted_pcl_vars->pcl_N[pcl_id].N4;
	}
};

#endif