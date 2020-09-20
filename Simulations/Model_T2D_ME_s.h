#ifndef __Model_T2D_ME_s_h__
#define __Model_T2D_ME_s_h__

#include "BCs.h"
#include "macro_utils.h"
#include "Model.h"
#include "TriangleMeshTemplate.hpp"
#include "MatModelContainer.h"
#include "SearchingGrid2D.hpp"
#include "ParticleGenerator2D.hpp"
#include "RigidBody/RigidCircle.h"
#include "RigidBody/RigidRect.h"

namespace Model_T2D_ME_s_Internal
{
struct Node
{
	size_t id;
	double x, y;

	bool has_mp;
	double m;
	double ax, ay, am;
	double vx, vy, vm;
	double dux, duy;
	double fx_ext, fy_ext;
	double fx_int, fy_int;

	// strain enhancement
	double de_vol_by_3;
	double se_pcl_vol;
};

struct Element;
struct Particle
{
	size_t id;
	double x, y;

	double ux, uy;
	double vx, vy;

	double m, density;

	double e11, e22, e12;
	double s11, s22, s12;

	// calculation variables
	double x_ori, y_ori;
	double vol;
	inline double get_vol() { return m / density; }
	
	Element* pe;
	double N1, N2, N3;
	
	Particle* next; // used by Element

	MatModel::MaterialModel* mm;
	inline void set_mat_model(MatModel::MaterialModel &_mm)
	{
		_mm.ext_data_pt = this;
		mm = &_mm;
	}
};

struct Element
{
	size_t id;
	size_t n1, n2, n3;
	double area;

	// for shape function calculation
	// Ni = ai * p.x + bi * p.y + coefi
	double a1, b1, coef1;
	double a2, b2, coef2;
	double a3, b3, coef3;
	// shape function at element centre
	// dN1_dx = a1, dN1_dy = b1
	// dN2_dx = a2, dN2_dy = b2
	// dN3_dx = a3, dN3_dy = b3
	double dN1_dx, dN1_dy;
	double dN2_dx, dN2_dy;
	double dN3_dx, dN3_dy;

	bool has_pcl;

	// particle list
	Particle* pcls;
	inline void add_pcl(Particle& pcl) noexcept
	{
		pcl.next = pcls;
		pcls = &pcl;
	}

	// m
	double pcl_m, pcl_density;

	// mixed integration
	double mi_pcl_vol, s11, s22, s12;

	// strain enhancement apprach
	double dde11, dde22, de12;
	double de_vol_by_3;
};

struct Edge { size_t n1, n2; };

typedef TriangleMeshTemplate<Node, Element, Edge> BgMesh;

}

class Model_T2D_ME_s;

class Step_T2D_ME_s_Geo;
int solve_substep_T2D_ME_s_Geo(void* _self);
int solve_substep_T2D_ME_s_Geo_avg(void* _self);
class Step_T2D_ME_s;
int solve_substep_T2D_ME_s(void* _self);
int solve_substep_T2D_ME_s_avg(void* _self);

class ResultFile_hdf5;
namespace Model_T2D_ME_s_hdf5_utilities
{
	int output_background_mesh_to_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_background_mesh_from_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_boundary_condition_to_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_boundary_condition_from_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_pcl_data_to_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_pcl_data_from_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_material_model_to_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_material_model_from_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_rigid_circle_to_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_rigid_circle_from_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_rigid_rect_to_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_rigid_rect_from_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);
}

struct Model_T2D_ME_s : public Model,
	public Model_T2D_ME_s_Internal::BgMesh,
	public MatModel::MatModelContainer
{
	friend class Step_T2D_ME_s_Geo;
	friend int solve_substep_T2D_ME_s_Geo(void* _self);
	friend int solve_substep_T2D_ME_s_Geo_avg(void* _self);
	friend class Step_T2D_ME_s;
	friend int solve_substep_T2D_ME_s(void *_self);
	friend int solve_substep_T2D_ME_s_avg(void* _self);

public:
	typedef Model_T2D_ME_s_Internal::BgMesh BgMesh;
	typedef Model_T2D_ME_s_Internal::Node Node;
	typedef Model_T2D_ME_s_Internal::Element Element;
	typedef Model_T2D_ME_s_Internal::Edge Edge;
	typedef Model_T2D_ME_s_Internal::Particle Particle;

protected:
	size_t pcl_num;
	Particle *pcls;
	
	// boundary conditions
	size_t bfx_num, bfy_num;
	BodyForceAtPcl* bfxs, * bfys;
	size_t tx_num, ty_num;
	TractionBCAtPcl* txs, * tys;
	size_t ax_num, ay_num;
	AccelerationBC* axs, * ays;
	size_t vx_num, vy_num;
	VelocityBC* vxs, * vys;

	// background grid to accelerate searching
	SearchingGrid2D<Model_T2D_ME_s> search_bg_grid;
	
public:
	Model_T2D_ME_s();
	~Model_T2D_ME_s();

	inline double get_bg_grid_xl() { return search_bg_grid.get_x_min(); }
	inline double get_bg_grid_xu() { return search_bg_grid.get_x_max(); }
	inline double get_bg_grid_yl() { return search_bg_grid.get_y_min(); }
	inline double get_bg_grid_yu() { return search_bg_grid.get_y_max(); }
	inline double get_bg_grid_hx() { return search_bg_grid.get_hx(); }
	inline double get_bg_grid_hy() { return search_bg_grid.get_hy(); }
	
	inline size_t get_pcl_num() { return pcl_num; }
	inline Particle* get_pcls() { return pcls; }

	int init_mesh(double* node_coords, size_t n_num,
				  size_t* elem_n_ids, size_t e_num);
	int load_mesh_from_hdf5(const char* file_name);

	int init_search_grid(double _hx, double _hy);

	int init_pcls(size_t num, double m, double density);
	int init_pcls(ParticleGenerator2D<Model_T2D_ME_s>& pg, double density);

	void alloc_pcls(size_t num);
	void clear_pcls();

	INIT_BC_TEMPLATE(bfx, BodyForceAtPcl)
	INIT_BC_TEMPLATE(bfy, BodyForceAtPcl)
	INIT_BC_TEMPLATE(tx, TractionBCAtPcl)
	INIT_BC_TEMPLATE(ty, TractionBCAtPcl)
	INIT_BC_TEMPLATE(ax, AccelerationBC)
	INIT_BC_TEMPLATE(ay, AccelerationBC)
	INIT_BC_TEMPLATE(vx, VelocityBC)
	INIT_BC_TEMPLATE(vy, VelocityBC)

protected: // helper functions
	void init_mesh_shape_funcs();

public:
	inline double cal_N1(Element& e, double x, double y)
	{ return e.a1 * x + e.b1 * y + e.coef1; }
	inline double cal_N2(Element& e, double x, double y)
	{ return e.a2 * x + e.b2 * y + e.coef2; }
	inline double cal_N3(Element& e, double x, double y)
	{ return e.a3 * x + e.b3 * y + e.coef3; }

	using BgMesh::is_in_triangle;

	inline bool is_in_triangle(Element &e, Particle &p)
	{
		double a = e.a1 * p.x + e.b1 * p.y + e.coef1;
		double b = e.a2 * p.x + e.b2 * p.y + e.coef2;
		//double c = e.a3 * p.x + e.b3 * p.y + e.coef3;
		double c = 1.0 - a - b;
		bool res = 0.0 <= a && a <= 1.0 &&
				   0.0 <= b && b <= 1.0 &&
				   0.0 <= c && c <= 1.0;
		if (res)
		{
			p.N1 = a > N_tol ? a : N_tol;
			p.N2 = b > N_tol ? b : N_tol;
			p.N3 = c > N_tol ? c : N_tol;
		}
		return res;
	}

	// search using background grid
	template <typename Point2D>
	inline Element* find_in_which_element(Point2D& pcl)
	{
		return search_bg_grid.find_in_which_element<Point2D>(pcl);
	}

	// brute force searching
	template <typename Point2D>
	inline Element* find_in_which_element_bf(Point2D &pcl)
	{
		for (size_t e_id = 0; e_id < elem_num; ++e_id)
		{
			Element& e = elems[e_id];
			if (is_in_triangle(e, pcl))
				return &e;
		}
		return nullptr;
	}

protected:
	double K_cont; // contact stiffness

	// rigid circle
	bool rigid_circle_is_init;
	RigidCircle rigid_circle;
	int apply_rigid_circle(double dt);

	// rigid rect
	bool rigid_rect_is_init;
	RigidRect rigid_rect;
	int apply_rigid_rect(double dt);

public:
	// interaction with rigid circle
	inline bool rigid_circle_is_valid() { return rigid_circle_is_init; }
	inline RigidCircle& get_rigid_circle() { return rigid_circle; }
	inline void init_rigid_circle(double _K_cont, double x, double y, double r, double density = 1.0)
	{
		rigid_circle_is_init = true;
		K_cont = _K_cont;
		rigid_circle.init(x, y, r, density);
	}
	inline void set_rigid_circle_velocity(double vx, double vy, double v_ang)
	{ rigid_circle.set_v_bc(vx, vy, v_ang); }

	// interaction with rigid rect
	inline bool rigid_rect_is_valid() { return rigid_rect_is_init; }
	inline RigidRect& get_rigid_rect() { return rigid_rect; }
	inline void init_rigid_rect(double _K_cont,
		double x, double y, double hx, double hy, double density = 1.0)
	{
		rigid_rect_is_init = true;
		K_cont = _K_cont;
		rigid_rect.init(x, y, hx, hy, density);
	}
	inline void set_rigid_rect_velocity(double vx, double vy, double v_ang)
	{ rigid_rect.set_v_bc(vx, vy, v_ang); }
	
protected:
	friend int Model_T2D_ME_s_hdf5_utilities::output_background_mesh_to_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_ME_s_hdf5_utilities::load_background_mesh_from_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_ME_s_hdf5_utilities::output_boundary_condition_to_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_ME_s_hdf5_utilities::load_boundary_condition_from_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_ME_s_hdf5_utilities::output_pcl_data_to_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_ME_s_hdf5_utilities::load_pcl_data_from_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_ME_s_hdf5_utilities::output_material_model_to_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_ME_s_hdf5_utilities::load_material_model_from_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_ME_s_hdf5_utilities::output_rigid_circle_to_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_ME_s_hdf5_utilities::load_rigid_circle_from_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_ME_s_hdf5_utilities::output_rigid_rect_to_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_ME_s_hdf5_utilities::load_rigid_rect_from_hdf5_file(Model_T2D_ME_s& md, ResultFile_hdf5& rf, hid_t grp_id);

public: // for debug
	void sum_vol_for_all_elements();
};

#endif