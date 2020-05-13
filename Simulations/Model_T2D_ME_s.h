#ifndef __Model_T2D_ME_s_H__
#define __Model_T2D_ME_s_H__

#include "BC.h"
#include "Model.h"
#include "ModelContainer.h"

#include "TriangleMesh.h"
#include "TriangleMeshToParticles.h"
#include "DispConRigidCircle.h"

#define SHAPE_FUNC_VALUE_CONTENT   \
	double N1, N2, N3;             \
	double dN1_dx, dN2_dx, dN3_dx; \
	double dN1_dy, dN2_dy, dN3_dy

#define N1(xi, eta) (xi)
#define N2(xi, eta) (eta)
#define N3(xi, eta) (1.0 - (xi) - (eta))
#define dN1_dxi   (1.0)
#define dN2_dxi   (0.0)
#define dN3_dxi  (-1.0)
#define dN1_deta  (0.0)
#define dN2_deta  (1.0)
#define dN3_deta (-1.0)
#define N_tol (1.0e-30)

class ConstitutiveModel;

int solve_substep_T2D_ME_s(void *_self);

struct Model_T2D_ME_s : public Model
{
	friend int solve_substep_T2D_ME_s(void *_self);

public: // Node, Element and Particle data structures
	struct Node
	{
		size_t id;
		double x, y;

		// material point
		bool has_mp;
		// solid phase
		double m;
		double ax, ay;
		double vx, vy;
		double dux, duy;
		double fx_ext, fy_ext;
		double fx_int, fy_int;

		// rigid body
		bool has_rb;
		double vx_rb, vy_rb, vol_rb;

		// for strain enhancement approach
		double vol, de_vol;
	};
	
	struct Element;
	struct Particle
	{
		size_t id;
		double x, y;

		double ux, uy;
		double vx, vy;

		double m, density, vol;

		// total strain
		double e11, e22, e12;
		// effective stress
		double s11, s22, s12;

		// calculation variables
		double x_ori, y_ori;
		Element *pe;

		// shape function value
		double N1, N2, N3;
		
		// Used by Element
		Particle *next;

		ConstitutiveModel *cm;
		inline void set_cm(ConstitutiveModel &_cm)
		{ 
			cm = &_cm;
			_cm.ext_data = this;
		}
	};

	struct Element
	{
		// index
		size_t id;
		//topology
		size_t n1, n2, n3;
		double area, area_2; // 2 * A
		
		// shape function of element centre
		double dN1_dx, dN1_dy;
		double dN2_dx, dN2_dy;
		double dN3_dx, dN3_dy;
		
		// calculation variables
		double vol, s11, s22, s12;
		
		// particles list
		Particle *pcls;
		inline void add_pcl(Particle &pcl) noexcept
		{
			pcl.next = pcls;
			pcls = &pcl;
		}

		// strain enhancement apprach
		double dde11, dde22, de12;
		double de_vol;
	};

public:
	size_t node_num;
	Node *nodes;
	size_t elem_num;
	Element *elems;

	size_t pcl_num;
	Particle *pcls;
	
	// boundary conditions
	size_t bfx_num, bfy_num;
	BodyForce *bfxs, *bfys;
	size_t tx_num, ty_num;
	TractionBC_MPM *txs, *tys;

	size_t ax_num, ay_num;
	AccelerationBC *axs, *ays;
	size_t vx_num, vy_num;
	VelocityBC *vxs, *vys;

	// Constitutive model
	double E;   // Elastic modulus
	double niu; // Poisson ratio
	
	double K_cont; // contact stiffness

	ModelContainer model_container;

public:
	Model_T2D_ME_s();
	~Model_T2D_ME_s();
	
	void init_mesh(double *node_coords, size_t n_num,
				   size_t *elem_n_ids,  size_t e_num);
	void clear_mesh(void);

	void init_pcls(size_t num, double m, double density, double _E, double _niu);
	void alloc_pcls(size_t num);
	void clear_pcls(void);

	void init_mesh(TriangleMesh &tri_mesh);
	void init_pcls(TriangleMeshToParticles &mh_2_pcl, double density, double E, double niu);

#define INIT_BC_TEMPLATE(name, type)    \
	void init_ ## name ## s(size_t num) \
	{                                   \
		if (name ## s)                  \
		{                               \
			if (name ## _num < num)		\
				delete[] name ## s;		\
			else                        \
			{                           \
				name ## _num = num;	    \
				return;                 \
			}                           \
		}                               \
		name ## s = new type ## [num];  \
		name ## _num = num;             \
	}                                   \
	void clear_ ## name ## s(void)      \
	{                                   \
		if (name ## s)                  \
		{                               \
			delete[] name ## s;         \
			name ## s = nullptr;        \
		}                               \
		name ## _num = 0;               \
	}

	INIT_BC_TEMPLATE(bfx, BodyForce)
	INIT_BC_TEMPLATE(bfy, BodyForce)
	INIT_BC_TEMPLATE(tx, TractionBC_MPM)
	INIT_BC_TEMPLATE(ty, TractionBC_MPM)
	INIT_BC_TEMPLATE(ax, AccelerationBC)
	INIT_BC_TEMPLATE(ay, AccelerationBC)
	INIT_BC_TEMPLATE(vx, VelocityBC)
	INIT_BC_TEMPLATE(vy, VelocityBC)

public:
	inline bool is_in_triangle(Element &elem, double x, double y)
	{
		Node &n1 = nodes[elem.n1];
		Node &n2 = nodes[elem.n2];
		Node &n3 = nodes[elem.n3];
		double a = (n2.x - n1.x) * (y - n1.y) - (x - n1.x) * (n2.y - n1.y);
		double b = (x - n1.x) * (n3.y - n1.y) - (n3.x - n1.x) * (y - n1.y);
		double c = elem.area_2 - a - b;
		return 0.0 <= a && a <= elem.area_2
			&& 0.0 <= b && b <= elem.area_2
			&& 0.0 <= c && c <= elem.area_2;
	}
	
public:
	inline bool is_in_triangle(Element &e, Particle &p)
	{
		Node &n1 = nodes[e.n1];
		Node &n2 = nodes[e.n2];
		Node &n3 = nodes[e.n3];
		double a = (n2.x - p.x) * (n3.y - p.y) - (n3.x - p.x) * (n2.y - p.y);
		double b = (n3.x - p.x) * (n1.y - p.y) - (n1.x - p.x) * (n3.y - p.y);
		double c = e.area_2 - a - b;
		bool res = 0.0 <= a && a <= e.area_2
				&& 0.0 <= b && b <= e.area_2
				&& 0.0 <= c && c <= e.area_2;
		if (res)
		{
			p.N1 = a / e.area_2;
			p.N2 = b / e.area_2;
			p.N3 = 1.0 - p.N1 - p.N2;
			if (p.N1 < N_tol)
				p.N1 = N_tol;
			if (p.N2 < N_tol)
				p.N2 = N_tol;
			if (p.N3 < N_tol)
				p.N3 = N_tol;
		}
		return res;
	}

	// find particle is in which element (to be accelerated)
	// calculate shape function
	// return true if particle is in mesh and false vise versa.
	inline bool init_pcl_cal_var(Particle &pcl)
	{
		for (size_t e_id = 0; e_id < elem_num; ++e_id)
		{
			Element &e = elems[e_id];
			if (is_in_triangle(e, pcl))
			{
				pcl.pe = &e;
				return true;
			}
		}
		pcl.pe = nullptr;
		return false;
	}

protected:
	DispConRigidCircle rigid_circle;
	typedef DispConRigidCircle::Particle RCParticle;
	inline bool is_in_triangle(Element &e, RCParticle &p)
	{
		Node &n1 = nodes[e.n1];
		Node &n2 = nodes[e.n2];
		Node &n3 = nodes[e.n3];
		double a = (n2.x - p.x) * (n3.y - p.y) - (n3.x - p.x) * (n2.y - p.y);
		double b = (n3.x - p.x) * (n1.y - p.y) - (n1.x - p.x) * (n3.y - p.y);
		double c = e.area_2 - a - b;
		bool res = 0.0 <= a && a <= e.area_2
				&& 0.0 <= b && b <= e.area_2
				&& 0.0 <= c && c <= e.area_2;
		if (res)
		{
			p.N1 = a / e.area_2;
			p.N2 = b / e.area_2;
			p.N3 = 1.0 - p.N1 - p.N2;
			if (p.N1 < N_tol)
				p.N1 = N_tol;
			if (p.N2 < N_tol)
				p.N2 = N_tol;
			if (p.N3 < N_tol)
				p.N3 = N_tol;
			// map rigid body velocity to nodes
			double n_vol;
			// node 1
			n_vol = p.vol * p.N1;
			n1.vol_rb += n_vol;
			n1.vx_rb  += n_vol * p.vx;
			n1.vy_rb  += n_vol * p.vy;
			n1.has_rb = true;
			// node 2
			n_vol = p.vol * p.N2;
			n2.vol_rb += n_vol;
			n2.vx_rb  += n_vol * p.vx;
			n2.vy_rb  += n_vol * p.vy;
			n2.has_rb = true;
			// node 3
			n_vol = p.vol * p.N3;
			n3.vol_rb += n_vol;
			n3.vx_rb  += n_vol * p.vx;
			n3.vy_rb  += n_vol * p.vy;
			n3.has_rb = true;
		}
		return res;
	}
	inline bool apply_pcl_to_mesh(RCParticle &pcl)
	{
		return find_in_which_element(pcl) ? true : false;	
	}
	// velocity control
	int apply_rigid_body_to_bg_mesh(double dtime);
	// force control (panelty function)
	int apply_contact_force_to_bg_mesh(double dtime);

public:
	inline void init_rigid_circle(double _r, double _x, double _y, double max_pcl_size)
	{
		rigid_circle.init(_r, _x, _y, max_pcl_size);
	}
	inline void set_rigid_circle_velocity(double _vx, double _vy, double _w)
	{
		rigid_circle.set_velocity(_vx, _vy, _w);
	}
	inline DispConRigidCircle &get_rigid_circle(void) noexcept { return rigid_circle; }
	inline void set_contact_stiffness(double _K_cont) noexcept { K_cont = _K_cont; }

	// =========================================================
	// ==================== background grid ====================
protected:
	struct PElement
	{
		Element *e;
		PElement *next;
	};

	struct Grid
	{
		size_t x_id, y_id;
		PElement *pelems;
	};
	
	double grid_x_min, grid_x_max, grid_hx;
	double grid_y_min, grid_y_max, grid_hy;
	size_t grid_x_num, grid_y_num, grid_num;
	Grid *bg_grids;

	MemoryUtilities::ItemBuffer<PElement> pe_buffer;
	inline void add_elem_to_grid(Grid &g, Element &e)
	{
		PElement *pe = pe_buffer.alloc();
		pe->e = &e;
		pe->next = g.pelems;
		g.pelems = pe;
	}

	void add_elem_to_bg_grid(Element &e);

public:
	template <typename Point>
	Element *find_in_which_element(Point &p);

	int init_bg_mesh(double hx, double hy);
	inline void clear_bg_mesh(void)
	{
		if (bg_grids)
		{
			delete[] bg_grids;
			bg_grids = nullptr;
		}
		grid_x_min = 0.0;
		grid_x_max = 0.0;
		grid_y_min = 0.0;
		grid_y_max = 0.0;
		grid_num = 0;
		grid_x_num = 0;
		grid_y_num = 0;
	}
};

template <typename Point>
inline Model_T2D_ME_s::Element *Model_T2D_ME_s::find_in_which_element(Point &p)
{
	if (p.x < grid_x_min || p.x > grid_x_max ||
		p.y < grid_y_min || p.y > grid_y_max)
		return nullptr;
	size_t x_id = size_t((p.x - grid_x_min) / grid_hx);
	size_t y_id = size_t((p.y - grid_y_min) / grid_hy);
	Grid &g = bg_grids[grid_x_num * y_id + x_id];
	PElement *pelem = g.pelems;
	Element *elem;
	while (pelem)
	{
		elem = pelem->e;
		if (is_in_triangle(*elem, p))
			return elem;
		pelem = pelem->next;
	}
	return nullptr;
}

#undef SHAPE_FUNC_VALUE_CONTENT
#undef N1
#undef N2
#undef N3
#undef dN1_dxi
#undef dN2_dxi
#undef dN3_dxi
#undef dN1_deta
#undef dN2_deta
#undef dN3_deta
#undef N_tol

#endif