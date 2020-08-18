#ifndef __Model_FEM_T3D_ME_s_h__
#define __Model_FEM_T3D_ME_s_h__

#include "macro_utils.h"
#include "BCs.h"
#include "Model.h"
#include "TetrahedronMeshTemplate.hpp"
#include "MatModelContainer.h"
#include "ParticleGenerator3D.hpp"

namespace Model_FEM_T3D_ME_s_Internal
{
	struct Node
	{
		size_t id;
		double x, y, z;
		double ux, uy, uz;
		double vx, vy, vz;

		double m;
		double ax, ay, az;
		double dux, duy, duz;
		double fx_ext, fy_ext, fz_ext;
		double fx_int, fy_int, fz_int;

		// strain enhancement
		double pcl_w, de_vol;
	};

	class Particle;
	// element face id:
	// face 0: n1, n2, n3
	// face 1: n1, n4, n2
	// face 2: n1, n3, n4
	// face 3: n2, n4, n3
	struct Element
	{
		// index
		size_t id;
		//topology
		size_t n1, n2, n3, n4;
		// volume
		double vol;

		double density;

		// particles in this element
		Particle* p1, * p2, * p3, * p4;

		// mixed integration
		double s11, s22, s33, s12, s23, s31;

		// strain enhancement
		double dde11, dde22, dde33, de12, de23, de31;
		double de_vol;

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
	};
	
	// gauss point
	struct Particle
	{
		size_t id;
		double x, y, z;
		inline double get_vol() { return w * pe->vol; }

		// total strain
		double e11, e22, e33, e12, e23, e31;
		// effective stress
		double s11, s22, s33, s12, s23, s31;

		// element
		Element *pe;
		// weight of gauss points
		double w;
		// shape function value
		double N1, N2, N3, N4;

		MatModel::MaterialModel* mm;
		inline void set_mat_model(MatModel::MaterialModel& _mm) noexcept
		{
			_mm.ext_data_pt = this;
			mm = &_mm;
		}
	};

	struct Edge { size_t n1, n2; };

	typedef TetrahedronMeshTemplate<Node, Element, Edge> BgMesh;
};

class Step_FEM_T3D_ME_s;
int solve_substep_FEM_T3D_ME_s(void* _self);

class Model_FEM_T3D_ME_s : public Model,
	public Model_FEM_T3D_ME_s_Internal::BgMesh,
	public MatModel::MatModelContainer
{
	friend Step_FEM_T3D_ME_s;
	friend int solve_substep_FEM_T3D_ME_s(void* _self);
public:
	typedef Model_FEM_T3D_ME_s_Internal::BgMesh BgMesh;
	typedef Model_FEM_T3D_ME_s_Internal::Node Node;
	typedef Model_FEM_T3D_ME_s_Internal::Element Element;
	typedef Model_FEM_T3D_ME_s_Internal::Particle Particle;
	typedef Model_FEM_T3D_ME_s_Internal::Edge Edge;

public:
	// material particles
	size_t pcl_num;
	Particle* pcls;

	// boundary conditions
	size_t bfx_num, bfy_num, bfz_num;
	BodyForceAtElem *bfxs, *bfys, *bfzs;
	size_t tx_num, ty_num, tz_num;
	TractionBCAtFace *txs, *tys, *tzs;
	size_t ax_num, ay_num, az_num;
	AccelerationBC* axs, * ays, * azs;
	size_t vx_num, vy_num, vz_num;
	VelocityBC* vxs, * vys, * vzs;

public:
	Model_FEM_T3D_ME_s();
	~Model_FEM_T3D_ME_s();

	inline size_t get_pcl_num() { return pcl_num; }
	inline Particle *get_pcls() { return pcls; }

	int init_mesh(double* node_coords, size_t n_num,
				  size_t* elem_n_ids, size_t e_num,
				  double density);
	int load_mesh_from_hdf5(const char* file_name, double density);

	void clear_pcls();

	INIT_BC_TEMPLATE(bfx, BodyForceAtElem)
	INIT_BC_TEMPLATE(bfy, BodyForceAtElem)
	INIT_BC_TEMPLATE(bfz, BodyForceAtElem)
	INIT_BC_TEMPLATE(tx, TractionBCAtFace)
	INIT_BC_TEMPLATE(ty, TractionBCAtFace)
	INIT_BC_TEMPLATE(tz, TractionBCAtFace)
	INIT_BC_TEMPLATE(ax, AccelerationBC)
	INIT_BC_TEMPLATE(ay, AccelerationBC)
	INIT_BC_TEMPLATE(az, AccelerationBC)
	INIT_BC_TEMPLATE(vx, VelocityBC)
	INIT_BC_TEMPLATE(vy, VelocityBC)
	INIT_BC_TEMPLATE(vz, VelocityBC)

protected: // helper functions
	void init_mesh(double density);
};

#endif