#ifndef __Model_T3D_ME_s_h__
#define __Model_T3D_ME_s_h__

#include "BCs.h"
#include "Model.h"
#include "TetrahedronMeshTemplate.hpp"
#include "SearchingGrid3D.hpp"
#include "MatModelContainer.h"
#include "ParticleGenerator3D.hpp"

namespace Model_T3D_ME_s_Internal
{
	struct Node
	{
		size_t id;
		double x, y, z;

		bool has_mp;
		double m;
		double ax, ay, az;
		double vx, vy, vz;
		double dux, duy, duz;
		double fx_ext, fy_ext, fz_ext;
		double fx_int, fy_int, fz_int;

		// strain enhancement
		double pcl_vol, de_vol;
	};

	struct Element;
	struct Particle
	{
		size_t id;
		double x, y, z;
		double vx, vy, vz;
		double m, density, vol;

		// total strain
		double e11, e22, e33, e12, e23, e31;
		// effective stress
		double s11, s22, s33, s12, s23, s31;

		// calculation variables
		double x_ori, y_ori, z_ori;
		double ux, uy, uz;
		Element *pe;

		// shape function value
		double N1, N2, N3, N4;

		Particle *next; // Used by Element

		MatModel::MaterialModel *mm;
		inline void set_mat_model(MatModel::MaterialModel &_mm) noexcept
		{
			_mm.ext_data = this;
			mm = &_mm;
		}
	};

	struct Element
	{
		// index
		size_t id;
		//topology
		size_t n1, n2, n3, n4;
		// volume
		double vol;

		// for shape function calculation
		// Ni = ai * p.x + bi * p.y + ci * p.z - coefi
		double a1, b1, c1, coef1;
		double a2, b2, c2, coef2;
		double a3, b3, c3, coef3;
		double a4, b4, c4, coef4;
		// derivative of shape functions at element centre
		// dN1_dx = a1, dN1_dy = b1, dN1_dz = c1
		// dN2_dx = a2, dN2_dy = b2, dN2_dz = c2
		// dN3_dx = a3, dN3_dy = b3, dN3_dz = c3
		// dN4_dx = a4, dN4_dy = b4, dN4_dz = c4
		double dN1_dx, dN1_dy, dN1_dz;
		double dN2_dx, dN2_dy, dN2_dz;
		double dN3_dx, dN3_dy, dN3_dz;
		double dN4_dx, dN4_dy, dN4_dz;

		// particles in this element
		Particle *pcls;
		inline void add_pcl(Particle &pcl) noexcept
		{
			pcl.next = pcls;
			pcls = &pcl;
		}

		// mixed integration
		double s11, s22, s33, s12, s23, s31;
		double pcl_vol;

		// strain enhancement
		double dde11, dde22, dde33, de12, de23, de31;
		double de_vol;
	};

	struct Edge { size_t n1, n2; };

	typedef TetrahedronMeshTemplate<Node, Element, Edge> BgMesh;
};

class Step_T3D_ME_s;
int solve_substep_T3D_ME_s(void *_self);

class Model_T3D_ME_s : public Model, 
	public Model_T3D_ME_s_Internal::BgMesh,
	public MatModel::MatModelContainer
{
	friend Step_T3D_ME_s;
	friend int solve_substep_T3D_ME_s(void *_self);
protected:
	typedef Model_T3D_ME_s_Internal::BgMesh BgMesh;
public:
	typedef Model_T3D_ME_s_Internal::Node Node;
	typedef Model_T3D_ME_s_Internal::Element Element;
	typedef Model_T3D_ME_s_Internal::Particle Particle;
	typedef Model_T3D_ME_s_Internal::Edge Edge;

public:
	// material particles
	size_t pcl_num;
	Particle *pcls;

	SearchingGrid3D<Model_T3D_ME_s> search_bg_grid;

	// boundary conditions
	size_t bfx_num, bfy_num, bfz_num;
	BodyForceAtPcl *bfxs, *bfys, *bfzs;
	size_t tx_num, ty_num, tz_num;
	TractionBCAtPcl *txs, *tys, *tzs;
	size_t ax_num, ay_num, az_num;
	AccelerationBC *axs, *ays, *azs;
	size_t vx_num, vy_num, vz_num;
	VelocityBC *vxs, *vys, *vzs;

public: // for init model data
	Model_T3D_ME_s();
	~Model_T3D_ME_s();

	size_t get_pcl_num() { return pcl_num; }
	Particle *get_pcls() { return pcls; }

protected:
	void init_mesh_shape_funcs();

public:
	int init_mesh(double *node_coords, size_t n_num,
				   size_t *elem_n_ids, size_t e_num);
	int load_mesh_from_hdf5(const char *file_name);

	int init_search_grid(double _hx, double _hy, double _hz);

	void init_pcls(size_t num, double m, double density);
	void init_pcls(ParticleGenerator3D<Model_T3D_ME_s> &pg, double density);
	void clear_pcls();

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
	void clear_ ## name ## s()      \
	{                                   \
		if (name ## s)                  \
		{                               \
			delete[] name ## s;         \
			name ## s = nullptr;        \
		}                               \
		name ## _num = 0;               \
	}

	INIT_BC_TEMPLATE(bfx, BodyForceAtPcl)
	INIT_BC_TEMPLATE(bfy, BodyForceAtPcl)
	INIT_BC_TEMPLATE(bfz, BodyForceAtPcl)
	INIT_BC_TEMPLATE(tx, TractionBCAtPcl)
	INIT_BC_TEMPLATE(ty, TractionBCAtPcl)
	INIT_BC_TEMPLATE(tz, TractionBCAtPcl)
	INIT_BC_TEMPLATE(ax, AccelerationBC)
	INIT_BC_TEMPLATE(ay, AccelerationBC)
	INIT_BC_TEMPLATE(az, AccelerationBC)
	INIT_BC_TEMPLATE(vx, VelocityBC)
	INIT_BC_TEMPLATE(vy, VelocityBC)
	INIT_BC_TEMPLATE(vz, VelocityBC)

public: // for calculation
	using BgMesh::is_in_tetrahedron;

	inline bool is_in_tetrahedron(Element &e, Particle &p)
	{
		p.N1 = e.a1 * p.x + e.b1 * p.y + e.c1 * p.z - e.coef1;
		p.N2 = e.a2 * p.x + e.b2 * p.y + e.c2 * p.z - e.coef2;
		p.N3 = e.a3 * p.x + e.b3 * p.y + e.c3 * p.z - e.coef3;
		//p.N4 = 1.0 - p.N1 - p.N2 - p.N3;
		p.N4 = e.a4 * p.x + e.b4 * p.y + e.c4 * p.z - e.coef4;

		if (p.N1 < 0.0 || p.N1 > 1.0 || p.N2 < 0.0 || p.N2 > 1.0 ||
			p.N3 < 0.0 || p.N3 > 1.0 || p.N4 < 0.0 || p.N4 > 1.0)
			return false;

#define N_tol 1.0e-10
		if (p.N1 < N_tol)
			p.N1 = N_tol;
		if (p.N2 < N_tol)
			p.N2 = N_tol;
		if (p.N3 < N_tol)
			p.N3 = N_tol;
		if (p.N4 < N_tol)
			p.N4 = N_tol;
		return true;
#undef N_tol
	}

public:
	// brute force searching
	inline Element *find_in_which_element_bf(Particle &pcl)
	{
		for (size_t e_id = 0; e_id < elem_num; ++e_id)
		{
			Element &e = elems[e_id];
			if (is_in_tetrahedron(e, pcl))
				return &e;
		}
		return nullptr;
	}

	inline Element *find_in_which_element(Particle &pcl)
	{
		return search_bg_grid.find_in_which_element<Particle>(pcl);
	}
};

#undef INIT_BC_TEMPLATE

#endif