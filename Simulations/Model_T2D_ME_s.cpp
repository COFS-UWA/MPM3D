#include "Simulations_pcp.h"

#include "Geometry.h"

#include "Model_T2D_ME_s.h"

Model_T2D_ME_s::Model_T2D_ME_s() :
	Model("Model_T2D_ME_s"),
	elems(nullptr), elem_num(0), nodes(nullptr), node_num(0),
	pcls(nullptr), pcl_num(0),
	bfx_num(0), bfy_num(0), bfxs(nullptr), bfys(nullptr),
	tx_num(0),  ty_num(0),  txs(nullptr),  tys(nullptr),
	ax_num(0), ay_num(0), axs(nullptr), ays(nullptr),
	vx_num(0), vy_num(0), vxs(nullptr), vys(nullptr),
	grid_x_min(0.0), grid_x_max(0.0), grid_y_min(0.0), grid_y_max(0.0),
	grid_x_num(0), grid_y_num(0), grid_num(0), bg_grids(nullptr) {}

Model_T2D_ME_s::~Model_T2D_ME_s()
{
	clear_mesh();
	clear_pcls();
	clear_bfxs();
	clear_bfys();
	clear_txs();
	clear_tys();
	clear_axs();
	clear_ays();
	clear_vxs();
	clear_vys();
	clear_bg_mesh();
}

void Model_T2D_ME_s::init_mesh(double *node_coords, size_t n_num,
								size_t *elem_n_ids,  size_t e_num)
{
	clear_mesh();
	if (n_num == 0 || e_num == 0)
		return;
	
	node_num = n_num;
	nodes = new Node[node_num];
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node &n = nodes[n_id];
		n.id = n_id;
		n.x = node_coords[2*n_id];
		n.y = node_coords[2*n_id + 1];
	}
	
	elem_num = e_num;
	elems = new Element[elem_num];
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		Element &e = elems[e_id];
		e.id = e_id;
		e.n1 = elem_n_ids[3*e_id];
		e.n2 = elem_n_ids[3*e_id + 1];
		e.n3 = elem_n_ids[3*e_id + 2];
		Node &n1 = nodes[e.n1];
		Node &n2 = nodes[e.n2];
		Node &n3 = nodes[e.n3];
		e.area_2 = (n2.x - n1.x) * (n3.y - n1.y) - (n3.x - n1.x) * (n2.y - n1.y);
		// Ensure that nodes index of element is in counter-clockwise sequence.
		if (e.area_2 < 0.0)
		{
			e.area_2 = -e.area_2;
			size_t n_tmp = e.n2;
			e.n2 = e.n3;
			e.n3 = n_tmp;
		}
		e.area = e.area_2 * 0.5;
		Node &_n2 = nodes[e.n2];
		Node &_n3 = nodes[e.n3];
		e.dN1_dx = (_n2.y - _n3.y) / e.area_2;
		e.dN1_dy = (_n3.x - _n2.x) / e.area_2;
		e.dN2_dx = (_n3.y -  n1.y) / e.area_2;
		e.dN2_dy = ( n1.x - _n3.x) / e.area_2;
		e.dN3_dx = ( n1.y - _n2.y) / e.area_2;
		e.dN3_dy = (_n2.x -  n1.x) / e.area_2;
	}
}

void Model_T2D_ME_s::clear_mesh(void)
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
}

void Model_T2D_ME_s::init_pcls(size_t num, double m, double density, double _E, double _niu)
{
	clear_pcls();
	pcl_num = num;
	pcls = new Particle[pcl_num];
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

	// constitutive model
	E = _E;
	niu = _niu;
	LinearElasticity *cms = model_container.add_LinearElasticity(pcl_num);
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		LinearElasticity &cm = cms[p_id];
		cm.set_param(_E, _niu);
		pcls[p_id].set_cm(cm);
	}
}

void Model_T2D_ME_s::alloc_pcls(size_t num)
{
	pcls = new Particle[num];
	pcl_num = num;
}

void Model_T2D_ME_s::clear_pcls(void)
{
	if (pcls)
	{
		delete[] pcls;
		pcls = nullptr;
	}
	pcl_num = 0;
}

void Model_T2D_ME_s::init_mesh(TriangleMesh &tri_mesh)
{
	clear_mesh();
	node_num = tri_mesh.get_node_num();
	elem_num = tri_mesh.get_elem_num();
	if (node_num == 0 || elem_num == 0)
		return;

	TriangleMesh::Node *ext_nodes = tri_mesh.get_nodes();
	nodes = new Node[node_num];
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		TriangleMesh::Node &ext_n = ext_nodes[n_id];
		Node &n = nodes[n_id];
		n.id = ext_n.id;
		n.x  = ext_n.x;
		n.y  = ext_n.y;
	}

	TriangleMesh::Element *ext_elems = tri_mesh.get_elems();
	elems = new Element[elem_num];
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		TriangleMesh::Element &ext_e = ext_elems[e_id];
		Element &e = elems[e_id];
		e.id = ext_e.id;
		e.n1 = ext_e.n1;
		e.n2 = ext_e.n2;
		e.n3 = ext_e.n3;
		e.area_2 = ext_e.area;
		e.area = e.area_2 * 0.5;
		Node &n1 = nodes[e.n1];
		Node &n2 = nodes[e.n2];
		Node &n3 = nodes[e.n3];
		e.dN1_dx = (n2.y - n3.y) / e.area_2;
		e.dN1_dy = (n3.x - n2.x) / e.area_2;
		e.dN2_dx = (n3.y - n1.y) / e.area_2;
		e.dN2_dy = (n1.x - n3.x) / e.area_2;
		e.dN3_dx = (n1.y - n2.y) / e.area_2;
		e.dN3_dy = (n2.x - n1.x) / e.area_2;
	}
}

void Model_T2D_ME_s::init_pcls(TriangleMeshToParticles &mh_2_pcl,
								double density, double E, double niu)
{
	pcl_num = mh_2_pcl.get_pcl_num();
	if (pcl_num == 0)
		return;
	init_pcls(pcl_num, density, density, E, niu);
	TriangleMeshToParticles::Particle *ext_ppcl = mh_2_pcl.first();
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		Particle &pcl = pcls[pcl_id];
		pcl.x = ext_ppcl->x;
		pcl.y = ext_ppcl->y;
		pcl.m *= ext_ppcl->vol;
		ext_ppcl = mh_2_pcl.next(ext_ppcl);
	}
}

int Model_T2D_ME_s::apply_rigid_body_to_bg_mesh(double dtime)
{
	// reset reaction force
	rigid_circle.rfx = 0.0;
	rigid_circle.rfy = 0.0;
	rigid_circle.rm = 0.0;
	// update state
	rigid_circle.update(dtime);

	for (size_t pcl_id = 0; pcl_id < rigid_circle.pcl_num; ++pcl_id)
	{
		RCParticle &pcl = rigid_circle.pcls[pcl_id];
		apply_pcl_to_mesh(pcl);
	}
	
	double fx, fy;
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node &n = nodes[n_id];
		if (n.has_rb && n.has_mp && n.vol_rb != 0.0)
		{
			n.vx_rb /= n.vol_rb;
			n.vy_rb /= n.vol_rb;
			// rough contact for both solid and fluid phase
			fx = -(n.vx_rb - n.vx) * n.m / dtime;
			fy = -(n.vy_rb - n.vy) * n.m / dtime;
			rigid_circle.add_reaction_force(n.x, n.y, fx, fy);
			// modify nodal velocity
			n.vx = n.vx_rb;
			n.vy = n.vy_rb;
		}
	}

	return 0;
}

int Model_T2D_ME_s::apply_contact_force_to_bg_mesh(double dtime)
{
	// reset reaction force
	rigid_circle.rfx = 0.0;
	rigid_circle.rfy = 0.0;
	rigid_circle.rm = 0.0;
	// update state
	rigid_circle.update(dtime);

	double dist, norm_x, norm_y;
	double f_cont, fx_cont, fy_cont;
	double nfx_cont, nfy_cont, ndax, nday;
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		Particle &pcl = pcls[p_id];
		if (pcl.pe && rigid_circle.is_in_circle(pcl.x, pcl.y))
		{
			dist = rigid_circle.get_overlay_and_norm(pcl.x, pcl.y, norm_x, norm_y);
			f_cont = K_cont * dist;
			fx_cont = f_cont * norm_x;
			fy_cont = f_cont * norm_y;
			// reaction force by the rigid object
			rigid_circle.add_reaction_force(pcl.x, pcl.y, fx_cont, fy_cont);
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
			n1.vx += ndax * dtime;
			n1.vy += nday * dtime;
			// node 2
			Node &n2 = nodes[e.n2];
			nfx_cont = pcl.N2 * fx_cont;
			nfy_cont = pcl.N2 * fy_cont;
			ndax = nfx_cont / n2.m;
			nday = nfy_cont / n2.m;
			n2.ax += ndax;
			n2.ay += nday;
			n2.vx += ndax * dtime;
			n2.vy += nday * dtime;
			// node 3
			Node &n3 = nodes[e.n3];
			nfx_cont = pcl.N3 * fx_cont;
			nfy_cont = pcl.N3 * fy_cont;
			ndax = nfx_cont / n3.m;
			nday = nfy_cont / n3.m;
			n3.ax += ndax;
			n3.ay += nday;
			n3.vx += ndax * dtime;
			n3.vy += nday * dtime;
		}
	}

	return 0;
}

int Model_T2D_ME_s::init_bg_mesh(double hx, double hy)
{
	if (elem_num == 0)
		return -1;

	clear_bg_mesh();
	pe_buffer.reset();

	grid_hx = hx;
	grid_hy = hy;
	pe_buffer.set_page_size(elem_num);

	// get bounding box
	grid_x_min = nodes[0].x;
	grid_x_max = grid_x_min;
	grid_y_min = nodes[0].y;
	grid_y_max = grid_y_min;
	for (size_t n_id = 1; n_id < node_num; ++n_id)
	{
		Node &n = nodes[n_id];
		if (grid_x_min > n.x)
			grid_x_min = n.x;
		if (grid_x_max < n.x)
			grid_x_max = n.x;
		if (grid_y_min > n.y)
			grid_y_min = n.y;
		if (grid_y_max < n.y)
			grid_y_max = n.y;
	}

	// init background grid
	double x_padding, y_padding;
	x_padding = (grid_x_max - grid_x_min) * 0.001;
	y_padding = (grid_y_max - grid_y_min) * 0.001;
	grid_x_min -= x_padding;
	grid_x_max += x_padding;
	grid_y_min -= y_padding;
	grid_y_max += y_padding;
	grid_x_num = size_t(ceil((grid_x_max - grid_x_min) / grid_hx));
	grid_y_num = size_t(ceil((grid_y_max - grid_y_min) / grid_hy));
	grid_num = grid_x_num * grid_y_num;
	grid_x_max = grid_x_min + double(grid_x_num) * grid_hx;
	grid_y_max = grid_y_min + double(grid_y_num) * grid_hy;
	bg_grids = new Grid[grid_num];
	size_t grid_id = 0;
	for (size_t row_id = 0; row_id < grid_x_num; ++row_id)
		for (size_t col_id = 0; col_id < grid_y_num; ++col_id)
		{
			Grid &g = bg_grids[grid_id];
			g.x_id = col_id;
			g.y_id = row_id;
			g.pelems = nullptr;
			++grid_id;
		}

	// add element to grids
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
		add_elem_to_bg_grid(elems[e_id]);

	return 0;
}

void Model_T2D_ME_s::add_elem_to_bg_grid(Element &e)
{
	// get bounding box of element
	double e_x_min, e_x_max, e_y_min, e_y_max;
	// node 1
	Node &n1 = nodes[e.n1];
	e_x_min = n1.x;
	e_x_max = e_x_min;
	e_y_min = n1.y;
	e_y_max = e_y_min;
	// node 2
	Node &n2 = nodes[e.n2];
	if (e_x_min > n2.x)
		e_x_min = n2.x;
	if (e_x_max < n2.x)
		e_x_max = n2.x;
	if (e_y_min > n2.y)
		e_y_min = n2.y;
	if (e_y_max < n2.y)
		e_y_max = n2.y;
	// node 3
	Node &n3 = nodes[e.n3];
	if (e_x_min > n3.x)
		e_x_min = n3.x;
	if (e_x_max < n3.x)
		e_x_max = n3.x;
	if (e_y_min > n3.y)
		e_y_min = n3.y;
	if (e_y_max < n3.y)
		e_y_max = n3.y;

	// test elem intersection with bg grid
	size_t min_x_id, max_x_id, min_y_id, max_y_id, grid_row_id;
	double grid_xl, grid_xu, grid_yl, grid_yu;
	min_x_id = size_t(floor((e_x_min - grid_x_min) / grid_hx));
	max_x_id = size_t( ceil((e_x_max - grid_x_min) / grid_hx)) + 1;
	min_y_id = size_t(floor((e_y_min - grid_y_min) / grid_hy));
	max_y_id = size_t( ceil((e_y_max - grid_y_min) / grid_hy)) + 1;
	for (size_t y_id = min_y_id; y_id < max_y_id; ++y_id)
	{
		grid_row_id = grid_x_num * y_id;
		grid_yl = grid_y_min + double(y_id) * grid_hy;
		grid_yu = grid_yl + grid_hy;
		for (size_t x_id = min_x_id; x_id < max_x_id; ++x_id)
		{
			Grid &g = bg_grids[grid_row_id + x_id];
			grid_xl = grid_x_min + double(x_id) * grid_hx;
			grid_xu = grid_xl + grid_hx;
			if (test_AABB_triangle_intersection(grid_xl, grid_xu, grid_yl, grid_yu,
												n1.x, n1.y, n2.x, n2.y, n3.x, n3.y))
				add_elem_to_grid(g, e);
		}
	}
}
