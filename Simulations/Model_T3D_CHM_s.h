#ifndef __Model_T3D_CHM_s_h__
#define __Model_T3D_CHM_s_h__

#include "macro_utils.h"
#include "BCs.h"
#include "Model.h"
#include "TetrahedronMeshTemplate.hpp"
#include "SearchingGrid3D.hpp"
#include "MatModelContainer.h"
#include "ParticleGenerator3D.hpp"

namespace Model_T3D_CHM_s_Internal
{
	struct Node
	{
		size_t id;
		double x, y, z;

		bool has_mp;
		// solid phase
		double m_s;
		double ax_s, ay_s, az_s;
		double vx_s, vy_s, vz_s;
		double dux_s, duy_s, duz_s;
		double fx_ext_s, fy_ext_s, fz_ext_s;
		double fx_int_s, fy_int_s, fz_int_s;
		// fluid phase
		double m_f;
		double ax_f, ay_f, az_f;
		double vx_f, vy_f, vz_f;
		double dux_f, duy_f, duz_f;
		double fx_ext_f, fy_ext_f, fz_ext_f;
		double fx_int_f, fy_int_f, fz_int_f;
		// solid - fluid interaction
		double fx_drag, fy_drag, fz_drag;

		// strain enhancement
		double pcl_vol, de_vol_s, de_vol_f;
	};

	struct Element;
	struct Particle
	{
		size_t id;
		double x, y, z;

		double n; // porosity
		double m_s; // solid mass 
		double density_s;
		double density_f;

		double vx_s, vy_s, vz_s;
		double vx_f, vy_f, vz_f;

		// total strain
		double e11, e22, e33, e12, e23, e31;
		// effective stress
		double s11, s22, s33, s12, s23, s31;
		// pore pressure
		double p;

		// calculation variables
		double vol_s; // solid volume
		double vol;
		double m_f;
		double x_ori, y_ori, z_ori;
		double ux_s, uy_s, uz_s;
		double ux_f, uy_f, uz_f;

		Element *pe;

		// shape function value
		double N1, N2, N3, N4;

		Particle *next; // Used by Element

		MatModel::MaterialModel *mm;
		inline void set_mat_model(MatModel::MaterialModel &_mm) noexcept
		{
			_mm.ext_data_pt = this;
			mm = &_mm;
		}

		inline double get_vol() const noexcept { return m_s / (density_s * (1.0 - n)); }
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
		double s11, s22, s33, s12, s23, s31, p;
		double pcl_vol, n;

		// strain enhancement
		double dde11, dde22, dde33, de12, de23, de31;
		double de_vol_s, de_vol_f;
	};

	struct Edge { size_t n1, n2; };

	typedef TetrahedronMeshTemplate<Node, Element, Edge> BgMesh;
};

class Step_T3D_CHM_s;
int solve_substep_T3D_CHM_s(void *_self);

struct Model_T3D_CHM_s : public Model,
	public Model_T3D_CHM_s_Internal::BgMesh,
	public MatModel::MatModelContainer
{
	friend Step_T3D_CHM_s;
	friend int solve_substep_T3D_CHM_s(void *_self);
public:
	typedef Model_T3D_CHM_s_Internal::BgMesh BgMesh;
	typedef Model_T3D_CHM_s_Internal::Node Node;
	typedef Model_T3D_CHM_s_Internal::Element Element;
	typedef Model_T3D_CHM_s_Internal::Particle Particle;
	typedef Model_T3D_CHM_s_Internal::Edge Edge;

public:
	// background grid to accelerate searching
	SearchingGrid3D<Model_T3D_CHM_s> search_bg_grid;

	// material particles
	size_t pcl_num;
	Particle *pcls;

	// boundary conditions
	size_t bfx_num, bfy_num, bfz_num;
	BodyForceAtPcl *bfxs, *bfys, *bfzs;
	size_t tx_num, ty_num, tz_num;
	TractionBCAtPcl *txs, *tys, *tzs;
	size_t asx_num, asy_num, asz_num;
	AccelerationBC *asxs, *asys, *aszs;
	size_t vsx_num, vsy_num, vsz_num;
	VelocityBC *vsxs, *vsys, *vszs;
	size_t afx_num, afy_num, afz_num;
	AccelerationBC *afxs, *afys, *afzs;
	size_t vfx_num, vfy_num, vfz_num;
	VelocityBC *vfxs, *vfys, *vfzs;

	double Kf;
	double k, miu;

public:
	Model_T3D_CHM_s();
	~Model_T3D_CHM_s();

	inline size_t get_pcl_num() { return pcl_num; }
	inline Particle *get_pcls() { return pcls; }

	int init_mesh(double *node_coords, size_t n_num,
				  size_t *elem_n_ids,  size_t e_num);
	int load_mesh_from_hdf5(const char *file_name);

	int init_search_grid(double _hx, double _hy, double _hz);

	int init_pcls(size_t num, double n, double m_s, double density_s, double density_f,
				   double _Kf, double _k, double _miu);
	int init_pcls(ParticleGenerator3D<Model_T3D_CHM_s> &pg,
				   double n, double density_s, double density_f,
				   double _Kf, double _k, double _miu);
	void alloc_pcls(size_t num);
	void clear_pcls();
	
	INIT_BC_TEMPLATE(bfx, BodyForceAtPcl)
	INIT_BC_TEMPLATE(bfy, BodyForceAtPcl)
	INIT_BC_TEMPLATE(bfz, BodyForceAtPcl)
	INIT_BC_TEMPLATE(tx, TractionBCAtPcl)
	INIT_BC_TEMPLATE(ty, TractionBCAtPcl)
	INIT_BC_TEMPLATE(tz, TractionBCAtPcl)
	INIT_BC_TEMPLATE(asx, AccelerationBC)
	INIT_BC_TEMPLATE(asy, AccelerationBC)
	INIT_BC_TEMPLATE(asz, AccelerationBC)
	INIT_BC_TEMPLATE(vsx, VelocityBC)
	INIT_BC_TEMPLATE(vsy, VelocityBC)
	INIT_BC_TEMPLATE(vsz, VelocityBC)
	INIT_BC_TEMPLATE(afx, AccelerationBC)
	INIT_BC_TEMPLATE(afy, AccelerationBC)
	INIT_BC_TEMPLATE(afz, AccelerationBC)
	INIT_BC_TEMPLATE(vfx, VelocityBC)
	INIT_BC_TEMPLATE(vfy, VelocityBC)
	INIT_BC_TEMPLATE(vfz, VelocityBC)

protected: // helper functions
	void init_mesh_shape_funcs();

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

		if (p.N1 < N_tol)
			p.N1 = N_tol;
		if (p.N2 < N_tol)
			p.N2 = N_tol;
		if (p.N3 < N_tol)
			p.N3 = N_tol;
		if (p.N4 < N_tol)
			p.N4 = N_tol;
		return true;
	}

public: // search with background grid
	inline Element *find_in_which_element(Particle &pcl)
	{
		return search_bg_grid.find_in_which_element<Particle>(pcl);
	}
};

#endif