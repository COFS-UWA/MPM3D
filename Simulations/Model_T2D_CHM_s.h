#ifndef __Model_T2D_CHM_s_h__
#define __Model_T2D_CHM_s_h__

#include "BCs.h"
#include "macro_utils.h"
#include "Model.h"
#include "TriangleMeshTemplate.hpp"
#include "MatModelContainer.h"
#include "SearchingGrid2D.hpp"
#include "ParticleGenerator2D.hpp"
#include "RigidCircle.h"

namespace Model_T2D_CHM_s_Internal
{
struct Node
{
	size_t id;
	double x, y;

	bool has_mp;
	// solid phase
	double m_s;
	double ax_s, ay_s;
	double vx_s, vy_s;
	double dux_s, duy_s;
	double fx_ext_s, fy_ext_s;
	double fx_int_s, fy_int_s;
	// fluid phase
	double m_f;
	double ax_f, ay_f;
	double vx_f, vy_f;
	double dux_f, duy_f;
	double fx_ext_f, fy_ext_f;
	double fx_int_f, fy_int_f;
	// solid - fluid interaction
	double fx_drag, fy_drag;

	// for strain enhancement approach
	double pcl_vol, de_vol_s, de_vol_f;
};

struct Element;
struct Particle
{
	size_t id;
	double x, y;

	double ux_s, uy_s;
	double vx_s, vy_s;

	double ux_f, uy_f;
	double vx_f, vy_f;

	double n; // porosity
	double m_s; // mass of solid
	double density_s;
	double density_f;

	double e11, e22, e12;
	double s11, s22, s12;
	double p;

	// calculation variables
	double x_ori, y_ori;
	double vol_s, vol, m_f;
	
	Element* pe;
	double N1, N2, N3;

	Particle* next; // used by Element

	MatModel::MaterialModel* mm;
	inline void set_mat_model(MatModel::MaterialModel& _mm)
	{
		_mm.ext_data = this;
		mm = &_mm;
	}

	inline double get_vol() { return m_s / (density_s * (1.0 - n)); }
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

	// particle list
	Particle* pcls;
	inline void add_pcl(Particle& pcl) noexcept
	{
		pcl.next = pcls;
		pcls = &pcl;
	}

	// calculation variables
	double pcl_vol, pcl_n, s11, s22, s12, p;
	
	// strain enhancement apprach
	double dde11, dde22, de12;
	double de_vol_s, de_vol_f;
};

struct Edge { size_t n1, n2; };

typedef TriangleMeshTemplate<Node, Element, Edge> BgMesh;

}

class Step_T2D_CHM_s_Geo;
int solve_substep_T2D_CHM_s_Geo(void* _self);
class Step_T2D_CHM_s;
int solve_substep_T2D_CHM_s(void* _self);

struct Model_T2D_CHM_s : public Model,
	public Model_T2D_CHM_s_Internal::BgMesh,
	public MatModel::MatModelContainer
{
	friend class Step_T2D_CHM_s_Geo;
	friend int solve_substep_T2D_CHM_s_Geo(void* _self);
	friend class Step_T2D_CHM_s;
	friend int solve_substep_T2D_CHM_s(void *_self);

public:
	typedef Model_T2D_CHM_s_Internal::BgMesh BgMesh;
	typedef Model_T2D_CHM_s_Internal::Node Node;
	typedef Model_T2D_CHM_s_Internal::Element Element;
	typedef Model_T2D_CHM_s_Internal::Particle Particle;
	typedef Model_T2D_CHM_s_Internal::Edge Edge;

	size_t pcl_num;
	Particle *pcls;
	
	// boundary conditions
	size_t bfx_num, bfy_num;
	BodyForceAtPcl *bfxs, *bfys;
	size_t tx_num, ty_num;
	TractionBCAtPcl *txs, *tys;
	size_t asx_num, asy_num;
	AccelerationBC *asxs, *asys;
	size_t afx_num, afy_num;
	AccelerationBC *afxs, *afys;
	size_t vsx_num, vsy_num;
	VelocityBC* vsxs, * vsys;
	size_t vfx_num, vfy_num;
	VelocityBC *vfxs, *vfys;

	// water constitutive parameters
	double Kf;  // Bulk modulus of water
	double k;   // Permeability
	double miu; // Dynamic viscosity
	
	// background grid to accelerate searching
	SearchingGrid2D<Model_T2D_CHM_s> search_bg_grid;

	// rigid object
	RigidCircle rigid_circle;
	// contact stiffness
	double Ks_cont; // solid - circle
	double Kf_cont; // fluid - circle

public:
	Model_T2D_CHM_s();
	~Model_T2D_CHM_s();
	
	inline double get_bg_grid_hx() { return search_bg_grid.get_hx(); }
	inline double get_bg_grid_hy() { return search_bg_grid.get_hy(); }

	inline size_t get_pcl_num() { return pcl_num; }
	inline Particle* get_pcls() { return pcls; }

	int init_mesh(double *node_coords, size_t n_num,
				   size_t *elem_n_ids,  size_t e_num);
	int load_mesh_from_hdf5(const char* file_name);

	int init_search_grid(double _hx, double _hy);

	int init_pcls(size_t num, double n, double m_s,
		double density_s, double density_f,
		double _Kf, double _k, double _miu);
	int init_pcls(ParticleGenerator2D<Model_T2D_CHM_s>& pg,
		double n, double density_s, double density_f,
		double _Kf, double _k, double _miu);

	void alloc_pcls(size_t num);
	void clear_pcls();

	INIT_BC_TEMPLATE(bfx, BodyForceAtPcl)
	INIT_BC_TEMPLATE(bfy, BodyForceAtPcl)
	INIT_BC_TEMPLATE(tx, TractionBCAtPcl)
	INIT_BC_TEMPLATE(ty, TractionBCAtPcl)
	INIT_BC_TEMPLATE(asx, AccelerationBC)
	INIT_BC_TEMPLATE(asy, AccelerationBC)
	INIT_BC_TEMPLATE(afx, AccelerationBC)
	INIT_BC_TEMPLATE(afy, AccelerationBC)
	INIT_BC_TEMPLATE(vsx, VelocityBC)
	INIT_BC_TEMPLATE(vsy, VelocityBC)
	INIT_BC_TEMPLATE(vfx, VelocityBC)
	INIT_BC_TEMPLATE(vfy, VelocityBC)

protected: // helper functions
	void init_mesh_shape_funcs();

public:
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
	inline Element* find_in_which_element_bf(Point2D& pcl)
	{
		for (size_t e_id = 0; e_id < elem_num; ++e_id)
		{
			Element& e = elems[e_id];
			if (is_in_triangle(e, pcl))
				return &e;
		}
		return nullptr;
	}

public: // interaction with rigid circle
	inline RigidCircle& get_rigid_circle() { return rigid_circle; }
	inline void init_rigid_circle(double _r, double _x, double _y, double max_pcl_size)
	{
		rigid_circle.init(_r, _x, _y, max_pcl_size);
	}
	inline void set_rigid_circle_velocity(double _vx, double _vy, double _w)
	{
		rigid_circle.set_velocity(_vx, _vy, _w);
	}
	inline void set_contact_stiffness(double _Ks_cont, double _Kf_cont)
	{ Ks_cont = _Ks_cont; Kf_cont = _Kf_cont; }

protected:
	int apply_contact_force_to_bg_mesh(double dtime);

public: // for debug
	void sum_vol_for_each_elements();
};

#endif