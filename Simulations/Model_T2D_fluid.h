#ifndef __Model_T2D_fluid_s_H__
#define __Model_T2D_fluid_s_H__

#include "BC.h"
#include "Model.h"

#include "TriangleMesh.h"
#include "TriangleMeshToParticles.h"

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

int solve_substep_T2D_fluid(void *_self);

struct Model_T2D_fluid : public Model
{
	friend int solve_substep_T2D_fluid(void *_self);

public: // Node, Element and Particle data structures
	struct Node
	{
		size_t id;
		double x, y;

		bool has_mp;
		double m;
		double ax, ay;
		double vx, vy;
		double dux, duy;
		double fx_ext, fy_ext;
		double fx_int, fy_int;

		// strain enhancement
		double vol, de_vol, de_vol_dt;
	};

	struct Element;
	struct Particle
	{
		size_t id;
		double m, density;
		double x, y;
		double vx, vy;

		// effective stress
		double t11, t22, t12, p;

		// calculation variables
		double vol;
		double ux, uy;
		double x_ori, y_ori;
		// the element it locates
		Element *pe;
		// shape function value
		double N1, N2, N3;

		// Used by Element
		Particle *next;
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
		
		double t11, t22, t12;
		// volumetric strain
		double de_vol;

		// particles list
		Particle *pcls;
		inline void add_pcl(Particle &pcl) noexcept
		{
			pcl.next = pcls;
			pcls = &pcl;
		}

		// mixed integration
		double vol, p;
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

	size_t ax_num, ay_num;
	AccelerationBC *axs, *ays;
	size_t vx_num, vy_num;
	VelocityBC *vxs, *vys;

	// Constitutive model
	double miu, lambda;
	double Kf; // bulk modulus
	
public:
	Model_T2D_fluid();
	~Model_T2D_fluid();

	void init_mesh(double *node_coords, size_t n_num, size_t *elem_n_ids, size_t e_num);
	void clear_mesh(void);

	void init_pcls(size_t num, double m, double density,
				   double _miu, double _lambda, double _Kf);
	void alloc_pcls(size_t num);
	void clear_pcls(void);

	void init_mesh(TriangleMesh &tri_mesh);
	void init_pcls(TriangleMeshToParticles &mh_2_pcl, double density,
				   double _miu, double _lambda, double _Kf);

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
	INIT_BC_TEMPLATE(ax, AccelerationBC)
	INIT_BC_TEMPLATE(ay, AccelerationBC)
	INIT_BC_TEMPLATE(vx, VelocityBC)
	INIT_BC_TEMPLATE(vy, VelocityBC)

public:
	inline bool is_in_triangle(Element &e, double x, double y)
	{
		Node &n1 = nodes[e.n1];
		Node &n2 = nodes[e.n2];
		Node &n3 = nodes[e.n3];
		double a = (n2.x - n1.x) * (y - n1.y) - (x - n1.x) * (n2.y - n1.y);
		double b = (x - n1.x) * (n3.y - n1.y) - (n3.x - n1.x) * (y - n1.y);
		double c = e.area_2 - a - b;
		return 0.0 <= a && a <= e.area_2
			&& 0.0 <= b && b <= e.area_2
			&& 0.0 <= c && c <= e.area_2;
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
inline Model_T2D_fluid::Element *Model_T2D_fluid::find_in_which_element(Point &p)
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