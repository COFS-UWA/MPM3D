#include "Simulations_pcp.h"

#include "Model_T3D_ME_s.h"

Model_T3D_ME_s::Model_T3D_ME_s() :
	node_num(0), nodes(nullptr),
	elem_num(0), elems(nullptr),
	edge_num(0), edges(nullptr),
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
	vz_num(0), vzs(nullptr)
{

}

Model_T3D_ME_s::~Model_T3D_ME_s()
{
	clear_mesh();
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

void Model_T3D_ME_s::init_mesh(
	double *node_coords, size_t n_num,
	size_t *elem_n_ids, size_t e_num
)
{
	clear_mesh();
	if (n_num == 0 || e_num == 0)
		return;
	
	node_num = n_num;
	nodes = new Node[n_num];
	double *ncrds = node_coords;
	for (size_t n_id = 0; n_id < n_num; ++n_id)
	{
		Node &n = nodes[n_id];
		n.id = n_id;
		n.x = *ncrds;
		++ncrds;
		n.y = *ncrds;
		++ncrds;
		n.z = *ncrds;
		++ncrds;
	}

	elem_num = e_num;
	elems = new Element[e_num];
	size_t *e_n_id = elem_n_ids;
	for (size_t e_id = 0; e_id < e_num; ++e_id)
	{
		Element &e = elems[e_id];
		e.id = e_id;
		e.n1 = *e_n_id;
		++e_n_id;
		e.n2 = *e_n_id;
		++e_n_id;
		e.n3 = *e_n_id;
		++e_n_id;
		e.n4 = *e_n_id;
		++e_n_id;
	}
}

void Model_T3D_ME_s::init_mesh(TetrahedronMesh &tri_mesh)
{
	clear_mesh();
	node_num = tri_mesh.get_node_num();
	elem_num = tri_mesh.get_elem_num();
	if (node_num == 0 || elem_num == 0)
		return;

	TetrahedronMesh::Node *mh_nodes = tri_mesh.get_nodes();
	nodes = new Node[node_num];
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node &n = nodes[n_id];
		TetrahedronMesh::Node &mh_n = mh_nodes[n_id];
		n.id = mh_n.id;
		n.x = mh_n.x;
		n.y = mh_n.y;
		n.z = mh_n.z;
	}

	double inv_vol_6;
	TetrahedronMesh::Element *mh_elems = tri_mesh.get_elems();
	elems = new Element[elem_num];
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		Element &e = elems[e_id];
		TetrahedronMesh::Element &mh_e = mh_elems[e_id];
		e.id = mh_e.id;
		e.n1 = mh_e.n1;
		e.n2 = mh_e.n2;
		e.n3 = mh_e.n3;
		e.n4 = mh_e.n4;
		e.vol = mh_e.vol;
		Node &n1 = nodes[e.n1];
		Node &n2 = nodes[e.n2];
		Node &n3 = nodes[e.n3];
		Node &n4 = nodes[e.n4];
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
}

void Model_T3D_ME_s::clear_mesh()
{
	if (nodes)
	{
		delete[] nodes;
		nodes = nullptr;
	}
	node_num = 0;

	if (elems)
	{
		delete[] elems;
		elems = nullptr;
	}
	elem_num = 0;

	if (edges)
	{
		delete[] edges;
		edges = nullptr;
	}
	edge_num = 0;
}

void Model_T3D_ME_s::init_pcls(size_t num, double m, double density)
{
	clear_pcls();
	pcls = new Particle[num];
	pcl_num = num;
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		Particle &pcl = pcls[pcl_id];
		pcl.id = pcl_id;
		pcl.m = m;
		pcl.density = density;
		pcl.vx = 0.0;
		pcl.vy = 0.0;
		pcl.e11 = 0.0;
		pcl.e22 = 0.0;
		pcl.e12 = 0.0;
		pcl.s11 = 0.0;
		pcl.s12 = 0.0;
		pcl.s22 = 0.0;
	}
}

void Model_T3D_ME_s::init_pcls(ParticleGenerator3D &pg, double density)
{
	clear_pcls();
	init_pcls(pg.get_num(), density, density);
	ParticleGenerator3D::Particle *p_iter = pg.first();
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		Particle &pcl = pcls[pcl_id];
		pcl.x = p_iter->x;
		pcl.y = p_iter->y;
		pcl.z = p_iter->z;
		pcl.m *= p_iter->vol;
		p_iter = pg.next(p_iter);
	}
}

void Model_T3D_ME_s::clear_pcls()
{
	if (pcls)
	{
		delete[] pcls;
		pcls = nullptr;
	}
	pcl_num = 0;
}

inline bool Model_T3D_ME_s::is_in_tetrahedron(Element &e, Particle &p)
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
