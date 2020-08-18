#ifndef __Model_T2D_CHM_p_h__
#define __Model_T2D_CHM_p_h__

#include "BCs.h"
#include "macro_utils.h"
#include "Model.h"
#include "TriangleMeshTemplate.hpp"
#include "MatModelContainer.h"
#include "SearchingGrid2D.hpp"
#include "ParticleGenerator2D.hpp"
#include "RigidCircle.h"

namespace Model_T2D_CHM_p_Internal
{
	class Particle;

	struct NodeToElem
	{
		size_t e_id;
		unsigned char n_id;
	};

	struct Node
	{
		size_t id;
		double x, y;

		// elements that this node associated with
		size_t n2e_num;
		NodeToElem* n2es;
		
		bool has_mp;
		// solid phase
		double m_s;
		double ax_s, ay_s;
		double vx_s, vy_s;
		double vx_s_ori, vy_s_ori;
		double dux_s, duy_s;
		double fx_s, fy_s;
		// fluid phase
		double m_f;
		double ax_f, ay_f;
		double vx_f, vy_f;
		double vx_f_ori, vy_f_ori;
		double dux_f, duy_f;
		double fx_f, fy_f;

		// for strain enhancement approach
		double de_vol_s;
		double de_vol_f;
		double se_pcl_vol;
	};

	struct NodeVarAtElem
	{
		double m_s;
		double vx_s, vy_s;
		double fx_s, fy_s;
		double m_f;
		double vx_f, vy_f;
		double fx_f, fy_f;
		// strain enhancement
		double de_vol_s;
		double de_vol_f;
		double se_pcl_vol;
		inline void init()
		{
			m_s = 0.0;
			vx_s = 0.0;
			vy_s = 0.0;
			fx_s = 0.0;
			fy_s = 0.0;
			m_f = 0.0;
			vx_f = 0.0;
			vy_f = 0.0;
			fx_f = 0.0;
			fy_f = 0.0;
			de_vol_s = 0.0;
			de_vol_f = 0.0;
			se_pcl_vol = 0.0;
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

		// calculation varibles
		bool has_mp;
		Particle* pcls;
		inline Particle* first_pcl() { return pcls; }
		Particle* next_pcl(Particle* pcl);
		inline bool not_last_pcl(Particle* pcl) { return pcl != nullptr; }
		void add_pcl(Particle& pcl);

		// mixed integration
		double mi_pcl_vol, mi_pcl_n, s11, s22, s12, p;

		// strain enhancement apprach
		double dde11, dde22, de12;
		double de_vol_s, de_vol_f;

		NodeVarAtElem node_vars[3];
	};

	struct Edge { size_t n1, n2; };

	typedef TriangleMeshTemplate<Node, Element, Edge> BgMesh;

	struct Particle
	{
		size_t id;
		double x, y;
		double x_f, y_f; // position of fluid phase
		
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
		double x_f_ori, y_f_ori;
		double vol_s, vol, m_f;

		inline double get_vol() { return m_s / (density_s * (1.0 - n)); }

		Element* pe;
		double N1, N2, N3;

		Particle* next; // used by Element

		// material model
		MatModel::MaterialModel* mm;
		inline void set_mat_model(MatModel::MaterialModel& _mm)
		{
			_mm.ext_data_pt = this;
			mm = &_mm;
		}

		Particle* next_pcl_in_elem;

		bool has_fx_s_ext, has_fy_s_ext;
		double fx_s_ext, fy_s_ext;
		bool has_fx_f_ext, has_fy_f_ext;
		double fx_f_ext, fy_f_ext;
	};

	inline Particle* Element::next_pcl(Particle* pcl) { return pcl->next_pcl_in_elem; }

	inline void Element::add_pcl(Particle& pcl)
	{
		pcl.next = pcls;
		pcls = &pcl;
	}
}

class Model_T2D_CHM_p;
class ResultFile_hdf5;
namespace Model_T2D_CHM_p_hdf5_utilities
{
	int output_background_mesh_to_hdf5_file(Model_T2D_CHM_p& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_background_mesh_from_hdf5_file(Model_T2D_CHM_p& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_boundary_condition_to_hdf5_file(Model_T2D_CHM_p& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_boundary_condition_from_hdf5_file(Model_T2D_CHM_p& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_pcl_data_to_hdf5_file(Model_T2D_CHM_p& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_pcl_data_from_hdf5_file(Model_T2D_CHM_p& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_material_model_to_hdf5_file(Model_T2D_CHM_p& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_material_model_from_hdf5_file(Model_T2D_CHM_p& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_rigid_circle_to_hdf5_file(Model_T2D_CHM_p& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_rigid_circle_from_hdf5_file(Model_T2D_CHM_p& md, ResultFile_hdf5& rf, hid_t grp_id);
}

class Step_T2D_CHM_p;
int solve_substep_T2D_CHM_p(void* _self);
class Step_T2D_CHM_p_Geo;
int solve_substep_T2D_CHM_p_Geo(void* _self);

struct Model_T2D_CHM_p : public Model,
	public Model_T2D_CHM_p_Internal::BgMesh,
	public MatModel::MatModelContainer
{
	friend class Step_T2D_CHM_p;
	friend int solve_substep_T2D_CHM_p(void *_self);
	friend class Step_T2D_CHM_p_Geo;
	friend int solve_substep_T2D_CHM_p_Geo(void* _self);

public:
	typedef Model_T2D_CHM_p_Internal::NodeToElem NodeToElem;
	typedef Model_T2D_CHM_p_Internal::Node Node;
	typedef Model_T2D_CHM_p_Internal::NodeVarAtElem NodeVarAtElem;
	typedef Model_T2D_CHM_p_Internal::Element Element;
	typedef Model_T2D_CHM_p_Internal::Edge Edge;
	typedef Model_T2D_CHM_p_Internal::BgMesh BgMesh;
	typedef Model_T2D_CHM_p_Internal::Particle Particle;

protected:
	NodeToElem* node2elems;
	
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
	SearchingGrid2D<Model_T2D_CHM_p> search_bg_grid;

public:
	explicit Model_T2D_CHM_p();
	~Model_T2D_CHM_p();
	
	Model_T2D_CHM_p(const Model_T2D_CHM_p &other) = delete;
	Model_T2D_CHM_p& operator=(const Model_T2D_CHM_p& other) = delete;

	inline double get_bg_grid_xl() { return search_bg_grid.get_x_min(); }
	inline double get_bg_grid_xu() { return search_bg_grid.get_x_max(); }
	inline double get_bg_grid_yl() { return search_bg_grid.get_y_min(); }
	inline double get_bg_grid_yu() { return search_bg_grid.get_y_max(); }
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
	int init_pcls(ParticleGenerator2D<Model_T2D_CHM_p>& pg,
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
	void init_node_to_elem_info();

public:
	void alloc_n2e_info(size_t num);
	void clear_n2e_info();
	void init_bcs();

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
	{ return search_bg_grid.find_in_which_element<Point2D>(pcl); }

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

// ===================== rigid circle =====================
protected:
	bool rigid_circle_is_init;
	// contact stiffness
	double Ks_cont, Kf_cont;
	RigidCircle rigid_circle;

public: // interaction with rigid circle
	inline bool rigid_circle_is_valid() { return rigid_circle_is_init; }
	inline void init_rigid_circle(double _Ks_cont, double _Kf_cont,
		double x, double y, double r, double density = 1.0)
	{
		rigid_circle_is_init = true;
		Ks_cont = _Ks_cont;
		Kf_cont = _Kf_cont;
		rigid_circle.init(x, y, r, density);
	}
	inline void set_rigid_circle_velocity(double vx, double vy, double w)
	{ rigid_circle.set_v_bc(vx, vy, w); }
	inline RigidCircle& get_rigid_circle() { return rigid_circle; }
	
protected:
	friend int Model_T2D_CHM_p_hdf5_utilities::output_background_mesh_to_hdf5_file(Model_T2D_CHM_p& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_CHM_p_hdf5_utilities::load_background_mesh_from_hdf5_file(Model_T2D_CHM_p& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_CHM_p_hdf5_utilities::output_boundary_condition_to_hdf5_file(Model_T2D_CHM_p& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_CHM_p_hdf5_utilities::load_boundary_condition_from_hdf5_file(Model_T2D_CHM_p& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_CHM_p_hdf5_utilities::output_pcl_data_to_hdf5_file(Model_T2D_CHM_p& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_CHM_p_hdf5_utilities::load_pcl_data_from_hdf5_file(Model_T2D_CHM_p& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_CHM_p_hdf5_utilities::output_material_model_to_hdf5_file(Model_T2D_CHM_p& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_CHM_p_hdf5_utilities::load_material_model_from_hdf5_file(Model_T2D_CHM_p& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_CHM_p_hdf5_utilities::output_rigid_circle_to_hdf5_file(Model_T2D_CHM_p& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_CHM_p_hdf5_utilities::load_rigid_circle_from_hdf5_file(Model_T2D_CHM_p& md, ResultFile_hdf5& rf, hid_t grp_id);

public: // for debug
	void sum_vol_for_each_elements();
};

#endif