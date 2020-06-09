#include "Simulations_pcp.h"

#include "ParticleGenerator3D.hpp"
#include "Model_FEM_T3D_ME_s.h"

Model_FEM_T3D_ME_s::Model_FEM_T3D_ME_s() :
	pcl_num(0), pcls(nullptr),
	bfx_num(0), bfxs(nullptr),
	bfy_num(0), bfys(nullptr),
	bfz_num(0), bfzs(nullptr),
	tx_num(0), txs(nullptr),
	ty_num(0), tys(nullptr),
	tz_num(0), tzs(nullptr),
	ax_num(0), axs(nullptr),
	ay_num(0), ays(nullptr),
	az_num(0), azs(nullptr),
	vx_num(0), vxs(nullptr),
	vy_num(0), vys(nullptr),
	vz_num(0), vzs(nullptr) {}

Model_FEM_T3D_ME_s::~Model_FEM_T3D_ME_s()
{
	clear_pcls();

	clear_bfxs();
	clear_bfys();
	clear_bfzs();
	clear_txs();
	clear_tys();
	clear_tzs();
	clear_axs();
	clear_ays();
	clear_azs();
	clear_vxs();
	clear_vys();
	clear_vzs();
}

namespace
{

inline void init_pcl(Model_FEM_T3D_ME_s::Particle &pcl)
{
	pcl.s11 = 0.0;
	pcl.s22 = 0.0;
	pcl.s33 = 0.0;
	pcl.s12 = 0.0;
	pcl.s23 = 0.0;
	pcl.s31 = 0.0;
	pcl.e11 = 0.0;
	pcl.e22 = 0.0;
	pcl.e33 = 0.0;
	pcl.e12 = 0.0;
	pcl.e23 = 0.0;
	pcl.e31 = 0.0;
	pcl.mm = nullptr;
}

}

void Model_FEM_T3D_ME_s::init_mesh(double density)
{
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node& n = nodes[n_id];
		n.vx = 0.0;
		n.vy = 0.0;
		n.vz = 0.0;
		n.ux = 0.0;
		n.uy = 0.0;
		n.uz = 0.0;
	}

	double inv_vol_6;
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		Element& e = elems[e_id];

		e.density = density;

		Node& n1 = nodes[e.n1];
		Node& n2 = nodes[e.n2];
		Node& n3 = nodes[e.n3];
		Node& n4 = nodes[e.n4];
		inv_vol_6 = 1.0 / (6.0 * e.vol);
		// N1
		e.a1 = ((n4.y - n2.y) * (n3.z - n2.z) - (n3.y - n2.y) * (n4.z - n2.z)) * inv_vol_6;
		e.b1 = ((n4.z - n2.z) * (n3.x - n2.x) - (n3.z - n2.z) * (n4.x - n2.x)) * inv_vol_6;
		e.c1 = ((n4.x - n2.x) * (n3.y - n2.y) - (n3.x - n2.x) * (n4.y - n2.y)) * inv_vol_6;
		e.coef1 = e.a1 * n2.x + e.b1 * n2.y + e.c1 * n2.z;
		// N2
		e.a2 = ((n4.y - n3.y) * (n1.z - n3.z) - (n1.y - n3.y) * (n4.z - n3.z)) * inv_vol_6;
		e.b2 = ((n4.z - n3.z) * (n1.x - n3.x) - (n1.z - n3.z) * (n4.x - n3.x)) * inv_vol_6;
		e.c2 = ((n4.x - n3.x) * (n1.y - n3.y) - (n1.x - n3.x) * (n4.y - n3.y)) * inv_vol_6;
		e.coef2 = e.a2 * n3.x + e.b2 * n3.y + e.c2 * n3.z;
		// N3
		e.a3 = ((n2.y - n4.y) * (n1.z - n4.z) - (n1.y - n4.y) * (n2.z - n4.z)) * inv_vol_6;
		e.b3 = ((n2.z - n4.z) * (n1.x - n4.x) - (n1.z - n4.z) * (n2.x - n4.x)) * inv_vol_6;
		e.c3 = ((n2.x - n4.x) * (n1.y - n4.y) - (n1.x - n4.x) * (n2.y - n4.y)) * inv_vol_6;
		e.coef3 = e.a3 * n4.x + e.b3 * n4.y + e.c3 * n4.z;
		// N4
		e.a4 = ((n2.y - n1.y) * (n3.z - n1.z) - (n3.y - n1.y) * (n2.z - n1.z)) * inv_vol_6;
		e.b4 = ((n2.z - n1.z) * (n3.x - n1.x) - (n3.z - n1.z) * (n2.x - n1.x)) * inv_vol_6;
		e.c4 = ((n2.x - n1.x) * (n3.y - n1.y) - (n3.x - n1.x) * (n2.y - n1.y)) * inv_vol_6;
		e.coef4 = e.a4 * n1.x + e.b4 * n1.y + e.c4 * n1.z;
		// derivatives
		e.dN1_dx = e.a1;
		e.dN1_dy = e.b1;
		e.dN1_dz = e.c1;
		e.dN2_dx = e.a2;
		e.dN2_dy = e.b2;
		e.dN2_dz = e.c2;
		e.dN3_dx = e.a3;
		e.dN3_dy = e.b3;
		e.dN3_dz = e.c3;
		e.dN4_dx = e.a4;
		e.dN4_dy = e.b4;
		e.dN4_dz = e.c4;
	}

	// gauss point
	clear_pcls();
	pcl_num = elem_num * 4;
	pcls = new Particle[pcl_num];
	Particle* ppcl = pcls;
	size_t pcl_id = 0;
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		Element& e = elems[e_id];
		Node& n1 = nodes[e.n1];
		Node& n2 = nodes[e.n2];
		Node& n3 = nodes[e.n3];
		Node& n4 = nodes[e.n4];
		Particle& p1 = ppcl[0];
		Particle& p2 = ppcl[1];
		Particle& p3 = ppcl[2];
		Particle& p4 = ppcl[3];

#define alpha 0.58541020
#define beta  0.13819660
		// gauss point 1
		e.p1 = &p1;
		p1.id = pcl_id++;
		p1.pe = &e;
		p1.w = 0.25;
		p1.N1 = alpha;
		p1.N2 = beta;
		p1.N3 = beta;
		p1.N4 = beta;
		p1.x = n1.x * p1.N1 + n2.x * p1.N2 + n3.x * p1.N3 + n4.x * p1.N4;
		p1.y = n1.y * p1.N1 + n2.y * p1.N2 + n3.y * p1.N3 + n4.y * p1.N4;
		p1.z = n1.z * p1.N1 + n2.z * p1.N2 + n3.z * p1.N3 + n4.z * p1.N4;
		init_pcl(p1);

		// gauss point 2
		e.p2 = &p2;
		p2.id = pcl_id++;
		p2.pe = &e;
		p2.w = 0.25;
		p2.N1 = beta;
		p2.N2 = alpha;
		p2.N3 = beta;
		p2.N4 = beta;
		p2.x = n1.x * p2.N1 + n2.x * p2.N2 + n3.x * p2.N3 + n4.x * p2.N4;
		p2.y = n1.y * p2.N1 + n2.y * p2.N2 + n3.y * p2.N3 + n4.y * p2.N4;
		p2.z = n1.z * p2.N1 + n2.z * p2.N2 + n3.z * p2.N3 + n4.z * p2.N4;
		init_pcl(p2);

		// gauss point 3
		e.p3 = &p3;
		p3.id = pcl_id++;
		p3.pe = &e;
		p3.w = 0.25;
		p3.N1 = beta;
		p3.N2 = beta;
		p3.N3 = alpha;
		p3.N4 = beta;
		p3.x = n1.x * p3.N1 + n2.x * p3.N2 + n3.x * p3.N3 + n4.x * p3.N4;
		p3.y = n1.y * p3.N1 + n2.y * p3.N2 + n3.y * p3.N3 + n4.y * p3.N4;
		p3.z = n1.z * p3.N1 + n2.z * p3.N2 + n3.z * p3.N3 + n4.z * p3.N4;
		init_pcl(p3);

		// gauss point 4
		e.p4 = &p4;
		p4.id = pcl_id++;
		p4.pe = &e;
		p4.w = 0.25;
		p4.N1 = beta;
		p4.N2 = beta;
		p4.N3 = beta;
		p4.N4 = alpha;
		p4.x = n1.x * p4.N1 + n2.x * p4.N2 + n3.x * p4.N3 + n4.x * p4.N4;
		p4.y = n1.y * p4.N1 + n2.y * p4.N2 + n3.y * p4.N3 + n4.y * p4.N4;
		p4.z = n1.z * p4.N1 + n2.z * p4.N2 + n3.z * p4.N3 + n4.z * p4.N4;
		init_pcl(p4);

		ppcl += 4;
	}
}

int Model_FEM_T3D_ME_s::init_mesh(
	double* node_coords, size_t n_num,
	size_t* elem_n_ids, size_t e_num,
	double density
)
{
	int res = BgMesh::init_mesh(node_coords, n_num, elem_n_ids, e_num);
	if (res < 0)
		return res;
	init_mesh(density);
	return 0;
}

int Model_FEM_T3D_ME_s::load_mesh_from_hdf5(
	const char* file_name,
	double density
	)
{
	int res = BgMesh::load_mesh_from_hdf5(file_name);
	if (res < 0)
		return res;
	init_mesh(density);
	return 0;
}

void Model_FEM_T3D_ME_s::clear_pcls()
{
	if (pcls)
	{
		delete[] pcls;
		pcls = nullptr;
	}
	pcl_num = 0;
}
