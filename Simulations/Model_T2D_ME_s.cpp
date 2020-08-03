#include "Simulations_pcp.h"

#include "Geometry.h"
#include "Model_T2D_ME_s.h"

Model_T2D_ME_s::Model_T2D_ME_s() :
	Model("Model_T2D_ME_s"),
	pcls(nullptr), pcl_num(0),
	bfx_num(0), bfxs(nullptr),
	bfy_num(0), bfys(nullptr),
	tx_num(0), txs(nullptr),
	ty_num(0), tys(nullptr),
	ax_num(0), axs(nullptr),
	ay_num(0), ays(nullptr),
	vx_num(0), vxs(nullptr),
	vy_num(0), vys(nullptr),
	rigid_circle_is_init(false) {}

Model_T2D_ME_s::~Model_T2D_ME_s()
{
	clear_pcls();

	clear_bfxs();
	clear_bfys();
	clear_txs();
	clear_tys();
	clear_axs();
	clear_ays();
	clear_vxs();
	clear_vys();
}

void Model_T2D_ME_s::init_mesh_shape_funcs()
{
	double area2;
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		Element& e = elems[e_id];
		Node& n1 = nodes[e.n1];
		Node& n2 = nodes[e.n2];
		Node& n3 = nodes[e.n3];
		area2 = n1.x * n2.y - n2.x * n1.y
			  + n2.x * n3.y - n3.x * n2.y
			  + n3.x * n1.y - n1.x * n3.y;
		e.a1 = (n2.y - n3.y) / area2;
		e.b1 = (n3.x - n2.x) / area2;
		e.coef1 = (n2.x * n3.y - n3.x * n2.y) / area2;
		e.a2 = (n3.y - n1.y) / area2;
		e.b2 = (n1.x - n3.x) / area2;
		e.coef2 = (n3.x * n1.y - n1.x * n3.y) / area2;
		e.a3 = (n1.y - n2.y) / area2;
		e.b3 = (n2.x - n1.x) / area2;
		e.coef3 = (n1.x * n2.y - n2.x * n1.y) / area2;
		// shape function derivatives
		e.dN1_dx = e.a1;
		e.dN1_dy = e.b1;
		e.dN2_dx = e.a2;
		e.dN2_dy = e.b2;
		e.dN3_dx = e.a3;
		e.dN3_dy = e.b3;
	}
}

int Model_T2D_ME_s::init_mesh(
	double *node_coords, size_t n_num,
	size_t *elem_n_ids,  size_t e_num
	)
{
	int res = BgMesh::init_mesh(node_coords, n_num, elem_n_ids, e_num);
	if (res < 0)
		return res;
	init_mesh_shape_funcs();
	return 0;
}

int Model_T2D_ME_s::load_mesh_from_hdf5(const char* file_name)
{
	int res = BgMesh::load_mesh_from_hdf5(file_name);
	if (res < 0)
		return res;
	init_mesh_shape_funcs();
	return 0;
}

int Model_T2D_ME_s::init_search_grid(double _hx, double _hy)
{
	return search_bg_grid.init(*this, _hx, _hy);
}

int Model_T2D_ME_s::init_pcls(size_t num, double m, double density)
{
	clear_pcls();

	if (num == 0)
		return -1;
	alloc_pcls(num);
	for (size_t pcl_id = 0; pcl_id < num; ++pcl_id)
	{
		Particle &pcl = pcls[pcl_id];
		pcl.id = pcl_id;
		pcl.ux = 0.0;
		pcl.uy = 0.0;
		pcl.vx = 0.0;
		pcl.vy = 0.0;
		pcl.m = m;
		pcl.density = density;
		pcl.e11 = 0.0;
		pcl.e22 = 0.0;
		pcl.e12 = 0.0;
		pcl.s11 = 0.0;
		pcl.s12 = 0.0;
		pcl.s22 = 0.0;
	}
	return 0;
}

int Model_T2D_ME_s::init_pcls(
	ParticleGenerator2D<Model_T2D_ME_s>& pg,
	double density
	)
{
	int res = init_pcls(pg.get_num(), density, density);
	if (res)
		return res;

	typedef ParticleGenerator2D<Model_T2D_ME_s>::Particle PgPcl;
	PgPcl* pg_pcl = pg.first();
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		Particle& pcl = pcls[pcl_id];
		pcl.x = pg_pcl->x;
		pcl.y = pg_pcl->y;
		pcl.m *= pg_pcl->area;
		pg_pcl = pg.next(pg_pcl);
	}

	return 0;
}


void Model_T2D_ME_s::alloc_pcls(size_t num)
{
	pcls = new Particle[num];
	pcl_num = num;
}

void Model_T2D_ME_s::clear_pcls()
{
	if (pcls)
	{
		delete[] pcls;
		pcls = nullptr;
	}
	pcl_num = 0;
}

// rigid circle - soil interaction
int Model_T2D_ME_s::apply_rigid_circle(double dt)
{
	rigid_circle.reset_rf();

	double dist, norm_x, norm_y;
	double f_cont, fx_cont, fy_cont;
	double nfx_cont, nfy_cont, ndax, nday;
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		Particle &pcl = pcls[p_id];
		if (pcl.pe && rigid_circle.detect_collision_with_point(
							pcl.x, pcl.y, pcl.vol, dist, norm_x, norm_y))
		{
			f_cont = K_cont * dist;
			fx_cont = f_cont * norm_x;
			fy_cont = f_cont * norm_y;
			// add reaction force to rigid object
			rigid_circle.add_rf(pcl.x, pcl.y, -fx_cont, -fy_cont);
			// adjust velocity at nodes
			Element &e = *pcl.pe;
			// node 1
			Node &n1 = nodes[e.n1];
			nfx_cont = pcl.N1 * fx_cont;
			nfy_cont = pcl.N1 * fy_cont;
			ndax = nfx_cont / n1.m;
			nday = nfy_cont / n1.m;
			n1.ax += ndax;
			n1.ay += nday;
			n1.vx += ndax * dt;
			n1.vy += nday * dt;
			// node 2
			Node &n2 = nodes[e.n2];
			nfx_cont = pcl.N2 * fx_cont;
			nfy_cont = pcl.N2 * fy_cont;
			ndax = nfx_cont / n2.m;
			nday = nfy_cont / n2.m;
			n2.ax += ndax;
			n2.ay += nday;
			n2.vx += ndax * dt;
			n2.vy += nday * dt;
			// node 3
			Node &n3 = nodes[e.n3];
			nfx_cont = pcl.N3 * fx_cont;
			nfy_cont = pcl.N3 * fy_cont;
			ndax = nfx_cont / n3.m;
			nday = nfy_cont / n3.m;
			n3.ax += ndax;
			n3.ay += nday;
			n3.vx += ndax * dt;
			n3.vy += nday * dt;
		}
	}

	rigid_circle.update_motion(dt);
	return 0;
}

#include <fstream>

void Model_T2D_ME_s::sum_vol_for_all_elements()
{
	// init elements
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		Element& e = elems[e_id];
		e.pcls = nullptr;
		e.mi_pcl_vol = 0.0;
	}

	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		Particle& pcl = pcls[p_id];
		if (!(pcl.pe = find_in_which_element(pcl)))
			continue;
		pcl.pe->add_pcl(pcl);
		pcl.vol = pcl.m / pcl.density;
		pcl.pe->mi_pcl_vol += pcl.vol;
	}

	std::fstream out_file;
	out_file.open("area.txt", std::ios::binary | std::ios::out);
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		Element& e = elems[e_id];
		if (e.pcls)
		{
			if (e.mi_pcl_vol > e.area)
			{
				out_file << "id: " << e.id << ", elem_a: " << e.area
					<< ", pcl_a:" << e.mi_pcl_vol
					<< " * +" << (e.mi_pcl_vol - e.area) / e.area * 100.0 << "% *\n";
			}
			else if (e.mi_pcl_vol < e.area)
			{
				out_file << "id: " << e.id << ", elem_a: " << e.area
					<< ", pcl_a:" << e.mi_pcl_vol
					<< " * -" << (e.area - e.mi_pcl_vol) / e.area * 100.0 << "% *\n";
			}
			else
			{
				out_file << "id: " << e.id << ", elem_a: " << e.area
					<< ", pcl_a:" << e.mi_pcl_vol << "\n";
			}
		}
	}
	out_file.close();
}
