#ifndef __Model_T3D_ME_s_h__
#define __Model_T3D_ME_s_h__

#include "BCs.h"
#include "Model.h"
#include "TetrahedronMeshTemplate.hpp"
#include "MatModelContainer.h"
#include "ParticleGenerator3D.h"

int solve_substep_T3D_ME_s(void *_self);

namespace
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


};

struct Model_T3D_ME_s : public Model
{
	friend int solve_substep_T3D_ME_s(void *_self);
public:
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

public:
	// model geometry
	size_t node_num;
	Node *nodes;
	size_t elem_num;
	Element *elems;

	size_t pcl_num;
	Particle *pcls;

	size_t edge_num;
	Edge *edges;

	// boundary conditions
	size_t bfx_num, bfy_num, bfz_num;
	BodyForceAtPcl *bfxs, *bfys, *bfzs;
	size_t tx_num, ty_num, tz_num;
	TractionBCAtPcl *txs, *tys, *tzs;

	size_t ax_num, ay_num, az_num;
	AccelerationBC *axs, *ays, *azs;
	size_t vx_num, vy_num, vz_num;
	VelocityBC *vxs, *vys, *vzs;

	// material models
	MatModel::MatModelContainer model_container;

	size_t get_node_num() const noexcept { return node_num; }
	Node *get_nodes() const noexcept { return nodes; }
	size_t get_elem_num() const noexcept { return elem_num; }
	Element *get_elems() const noexcept { return elems; }
	size_t get_edge_num() const noexcept { return edge_num; }
	Edge *get_edges() const noexcept { return edges; }

public:
	Model_T3D_ME_s();
	~Model_T3D_ME_s();

	void init_mesh(double *node_coords, size_t n_num,
				   size_t *elem_n_ids, size_t e_num);
	void init_mesh(TetrahedronMesh &tri_mesh);
	void clear_mesh();

	void init_pcls(size_t num, double m, double density);
	void init_pcls(ParticleGenerator3D &pg, double density);
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

	bool is_in_tetrahedron(Element &e, Particle &p);

	// brute force searching
	inline Element *find_in_which_element(Particle &pcl)
	{
		for (size_t e_id = 0; e_id < elem_num; ++e_id)
		{
			Element &e = elems[e_id];
			if (is_in_tetrahedron(e, pcl))
				return &e;
		}
		return nullptr;
	}
};

#undef INIT_BC_TEMPLATE

#endif