#include "SimulationsOMP_pcp.h"

#include <iostream>

#include "TriangleUtils.h"
#include "SearchingGrid2D.hpp"
#include "Model_T2D_ME_mt.h"

Model_T2D_ME_mt::Model_T2D_ME_mt() :
	ori_pcl_num(0), pcl_num(0),
	pcl_mem_raw(nullptr),
	node_num(0), elem_num(0),
	mesh_mem_raw(nullptr),
	bfx_num(0), bfxs(nullptr),
	bfy_num(0), bfys(nullptr),
	tx_num(0), txs(nullptr),
	ty_num(0), tys(nullptr),
	grid_x_num(0), grid_y_num(0),
	grid_elem_list(nullptr),
	grid_elem_list_id_array(nullptr),
	rigid_rect_is_valid(false),
	contact_mem(nullptr),
	//pcm(&smooth_contact)
	pcm(&rough_contact)
	//pcm(&fric_contact)
	//pcm(&sticky_contact)
	{}

Model_T2D_ME_mt::~Model_T2D_ME_mt()
{
	clear_mesh();
	clear_search_grid();
	clear_pcls();
	clear_bfxs();
	clear_bfys();
	clear_txs();
	clear_tys();
	clear_contact_mem();
}

Rect Model_T2D_ME_mt::get_mesh_bbox()
{
	if (!node_num)
		return Rect(0.0, 0.0, 0.0, 0.0);

	Rect res(node_pos[0].x, node_pos[0].x,
			 node_pos[0].y, node_pos[0].y);
	for (size_t n_id = 1; n_id < node_num; ++n_id)
	{
		Position &np = node_pos[n_id];
		if (res.xl > np.x)
			res.xl = np.x;
		if (res.xu < np.x)
			res.xu = np.x;
		if (res.yl > np.y)
			res.yl = np.y;
		if (res.yu < np.y)
			res.yu = np.y;
	}
	return res;
}

void Model_T2D_ME_mt::clear_mesh()
{
	if (mesh_mem_raw)
	{
		delete[] mesh_mem_raw;
		mesh_mem_raw = nullptr;
	}
}

void Model_T2D_ME_mt::alloc_mesh(
	size_t n_num,
	size_t e_num)
{
	node_num = n_num;
	elem_num = e_num;

	size_t mem_len = (sizeof(ElemNodeIndex)
		+ sizeof(ShapeFuncAB) + sizeof(ShapeFuncC)
		+ sizeof(double) * 4 + sizeof(StrainInc)
		+ sizeof(ElemNodeVM) * 3 + sizeof(Force) * 3) * e_num
		+ (sizeof(Position) + sizeof(Acceleration)
		+ sizeof(Velocity) + sizeof(NodeHasVBC)
		+ sizeof(double) * 2) * n_num;
	mesh_mem_raw = new char[mem_len];

	char* cur_mem = mesh_mem_raw;
	elem_node_id = (ElemNodeIndex*)cur_mem;
	cur_mem += sizeof(ElemNodeIndex) * elem_num;
	elem_dN_ab = (ShapeFuncAB*)cur_mem;
	cur_mem += sizeof(ShapeFuncAB) * elem_num;
	elem_dN_c = (ShapeFuncC*)cur_mem;
	cur_mem += sizeof(ShapeFuncC) * elem_num;
	elem_area = (double*)cur_mem;
	cur_mem += sizeof(double) * elem_num;

	elem_pcl_m = (double*)cur_mem;
	cur_mem += sizeof(double) * elem_num;
	elem_density = (double*)cur_mem;
	cur_mem += sizeof(double) * elem_num;
	elem_de = (StrainInc *)cur_mem;
	cur_mem += sizeof(StrainInc) * elem_num;
	elem_m_de_vol = (double*)cur_mem;
	cur_mem += sizeof(double) * elem_num;

	elem_node_vm = (ElemNodeVM*)cur_mem;
	cur_mem += sizeof(ElemNodeVM) * elem_num * 3;
	elem_node_force = (Force*)cur_mem;
	cur_mem += sizeof(Force) * elem_num * 3;

	node_pos = (Position*)cur_mem;
	cur_mem += sizeof(Position) * node_num;
	node_a = (Acceleration *)cur_mem;
	cur_mem += sizeof(Acceleration) * node_num;
	node_v = (Velocity *)cur_mem;
	cur_mem += sizeof(Velocity) * node_num;
	node_has_vbc = (NodeHasVBC *)cur_mem;
	cur_mem += sizeof(NodeHasVBC) * node_num;
	node_am = (double*)cur_mem;
	cur_mem += sizeof(double) * node_num;
	node_de_vol = (double*)cur_mem;
	cur_mem += sizeof(double) * node_num;
}

void Model_T2D_ME_mt::init_mesh(const TriangleMesh &mesh)
{
	clear_mesh();

	if (mesh.get_elem_num() == 0 ||
		mesh.get_node_num() == 0)
		return;
	
	alloc_mesh(mesh.get_node_num(), mesh.get_elem_num());

	// init node coordinates
	const TriangleMesh::Node* nodes = mesh.get_nodes();
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		const TriangleMesh::Node& n = nodes[n_id];
		Position &np = node_pos[n_id];
		np.x = n.x;
		np.y = n.y;
		NodeHasVBC& n_vbc = node_has_vbc[n_id];
		n_vbc.has_vx_bc = false;
		n_vbc.has_vy_bc = false;
	}

	// init elem connectivity, area and shape functions
	PointInTriangle pit;
	const TriangleMesh::Element *elems = mesh.get_elems();
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		const TriangleMesh::Element &e = elems[e_id];
		ElemNodeIndex& eni = elem_node_id[e_id];
		eni.n1 = e.n1;
		eni.n2 = e.n2;
		eni.n3 = e.n3;
		Position &n1_pos = node_pos[e.n1];
		Position &n2_pos = node_pos[e.n2];
		Position &n3_pos = node_pos[e.n3];
		elem_area[e_id] = cal_triangle_area(n1_pos, n2_pos, n3_pos);
		// shape functions
		pit.init_triangle(n1_pos, n2_pos, n3_pos, elem_area[e_id]);
		ShapeFuncAB &esfab = elem_dN_ab[e_id];
		esfab.dN1_dx = pit.dN1_dx();
		esfab.dN1_dy = pit.dN1_dy();
		esfab.dN2_dx = pit.dN2_dx();
		esfab.dN2_dy = pit.dN2_dy();
		esfab.dN3_dx = pit.dN3_dx();
		esfab.dN3_dy = pit.dN3_dy();
		ShapeFuncC& esfc = elem_dN_c[e_id];
		esfc.c1 = pit.get_coef1();
		esfc.c2 = pit.get_coef2();
		esfc.c3 = pit.get_coef3();
	}
}

void Model_T2D_ME_mt::clear_search_grid()
{
	if (grid_elem_list)
	{
		delete[] grid_elem_list;
		grid_elem_list = nullptr;
	}
	if (grid_elem_list_id_array)
	{
		delete[] grid_elem_list_id_array;
		grid_elem_list_id_array = nullptr;
	}
	grid_x_num = 0;
	grid_y_num = 0;
}

int Model_T2D_ME_mt::init_search_grid(
	TriangleMesh& mesh,
	double _hx,
	double _hy)
{
	clear_search_grid();
	
	SearchingGrid2D<TriangleMesh> search_grid;
	search_grid.init(mesh, _hx, _hy);

	grid_xl = search_grid.get_x_min();
	grid_yl = search_grid.get_y_min();
	grid_xu = search_grid.get_x_max();
	grid_yu = search_grid.get_y_max();
	grid_hx = search_grid.get_hx();
	grid_hy = search_grid.get_hy();
	grid_x_num = search_grid.get_x_num();
	grid_y_num = search_grid.get_y_num();

	const size_t num = grid_x_num * grid_y_num;
	typedef SearchingGrid2D<TriangleMesh>::Grid Grid;
	Grid *sg = search_grid.get_grids();
	grid_elem_list = new size_t[num + 1];
	grid_elem_list[0] = 0;
	size_t cur_id = 0;
	for (size_t g_id = 0; g_id < num; ++g_id)
	{
		Grid &grid = sg[g_id];
		for (auto pe = grid.pelems; pe; pe = pe->next)
			++cur_id;
		grid_elem_list[g_id+1] = cur_id;
	}

	grid_elem_list_id_array = new size_t[cur_id];
	cur_id = 0;
	for (size_t g_id = 0; g_id < num; ++g_id)
	{
		Grid& grid = sg[g_id];
		for (auto pe = grid.pelems; pe; pe = pe->next)
		{
			grid_elem_list_id_array[cur_id] = pe->e->id;
			++cur_id;
		}
	}

	return 0;
}

void Model_T2D_ME_mt::clear_pcls()
{
	if (pcl_mem_raw)
	{
		delete[] pcl_mem_raw;
		pcl_mem_raw = nullptr;
	}
	ori_pcl_num = 0;
	pcl_num = 0;
}

void Model_T2D_ME_mt::alloc_pcls(size_t num)
{
	clear_pcls();

	if (num == 0)
		return;

	ori_pcl_num = num;
	pcl_num = ori_pcl_num;
	size_t mem_len = (sizeof(double) + sizeof(Force) * 2
			 + sizeof(Position) + sizeof(double)
			 + sizeof(MatModel::MaterialModel *)
			+ (sizeof(size_t) + sizeof(double)
			 + sizeof(Velocity) + sizeof(Displacement)
			 + sizeof(ShapeFunc) + sizeof(Stress)
			 + sizeof(Strain) * 3) * 2
			) * num;
	pcl_mem_raw = new char[mem_len];

	char *cur_mem = pcl_mem_raw;
	pcl_m = (double *)cur_mem;
	cur_mem += sizeof(double) * num;
	pcl_bf = (Force *)(cur_mem);
	cur_mem += sizeof(Force) * num;
	pcl_t = (Force *)cur_mem;
	cur_mem += sizeof(Force) * num;
	pcl_pos = (Position *)cur_mem;
	cur_mem += sizeof(Position) * num;
	pcl_vol = (double *)cur_mem;
	cur_mem += sizeof(double) * num;

	SortedPclVarArrays& spva0 = sorted_pcl_var_arrays[0];
	spva0.pcl_index = (size_t *)cur_mem;
	cur_mem += sizeof(size_t) * num;
	spva0.pcl_density = (double *)cur_mem;
	cur_mem += sizeof(double) * num;
	spva0.pcl_v = (Velocity*)cur_mem;
	cur_mem += sizeof(Velocity) * num;
	spva0.pcl_disp = (Displacement *)cur_mem;
	cur_mem += sizeof(Displacement) * num;
	spva0.pcl_N = (ShapeFunc *)cur_mem;
	cur_mem += sizeof(ShapeFunc) * num;
	spva0.pcl_stress = (Stress*)cur_mem;
	cur_mem += sizeof(Stress) * num;
	spva0.pcl_strain = (Strain*)cur_mem;
	cur_mem += sizeof(Strain) * num;
	spva0.pcl_estrain = (Strain*)cur_mem;
	cur_mem += sizeof(Strain) * num;
	spva0.pcl_pstrain = (Strain*)cur_mem;
	cur_mem += sizeof(Strain) * num;

	SortedPclVarArrays &spva1 = sorted_pcl_var_arrays[1];
	spva1.pcl_index = (size_t *)cur_mem;
	cur_mem += sizeof(size_t) * num;
	spva1.pcl_density = (double *)cur_mem;
	cur_mem += sizeof(double) * num;
	spva1.pcl_v = (Velocity *)cur_mem;
	cur_mem += sizeof(Velocity) * num;
	spva1.pcl_disp = (Displacement *)cur_mem;
	cur_mem += sizeof(Displacement) * num;
	spva1.pcl_N = (ShapeFunc*)cur_mem;
	cur_mem += sizeof(ShapeFunc) * num;
	spva1.pcl_stress = (Stress*)cur_mem;
	cur_mem += sizeof(Stress) * num;
	spva1.pcl_strain = (Strain*)cur_mem;
	cur_mem += sizeof(Strain) * num;
	spva1.pcl_estrain = (Strain*)cur_mem;
	cur_mem += sizeof(Strain) * num;
	spva1.pcl_pstrain = (Strain*)cur_mem;
	cur_mem += sizeof(Strain) * num;

	pcl_mat_model = new MatModel::MaterialModel*[num];

	alloc_contact_mem(pcl_num);
}

void Model_T2D_ME_mt::alloc_pcls(
	size_t num,
	size_t ori_num)
{
	clear_pcls();

	if (num == 0 || ori_num == 0)
		return;
	
	alloc_pcls(ori_num);
	pcl_num = num;
}

int Model_T2D_ME_mt::init_pcls(size_t num, double m, double density)
{
	alloc_pcls(num);

	SortedPclVarArrays &spva0 = sorted_pcl_var_arrays[0];
	size_t p_id;
	for (p_id = 0; p_id < num; ++p_id)
	{
		pcl_m[p_id] = m;
		Force &p_bf = pcl_bf[p_id];
		p_bf.fx = 0.0;
		p_bf.fy = 0.0;
		Force &p_t = pcl_t[p_id];
		p_t.fx = 0.0;
		p_t.fy = 0.0;
		spva0.pcl_index[p_id] = p_id;
		spva0.pcl_density[p_id] = density;
		Displacement& p_disp = spva0.pcl_disp[p_id];
		p_disp.ux = 0.0;
		p_disp.uy = 0.0;
		Velocity& p_v = spva0.pcl_v[p_id];
		p_v.vx = 0.0;
		p_v.vy = 0.0;
		Stress& p_s = spva0.pcl_stress[p_id];
		p_s.s11 = 0.0;
		p_s.s22 = 0.0;
		p_s.s12 = 0.0;
		Strain& p_e = spva0.pcl_strain[p_id];
		p_e.e11 = 0.0;
		p_e.e22 = 0.0;
		p_e.e12 = 0.0;
		Strain& p_ee = spva0.pcl_estrain[p_id];
		p_ee.e11 = 0.0;
		p_ee.e22 = 0.0;
		p_ee.e12 = 0.0;
		Strain& p_pe = spva0.pcl_pstrain[p_id];
		p_pe.e11 = 0.0;
		p_pe.e22 = 0.0;
		p_pe.e12 = 0.0;
		pcl_mat_model[p_id] = nullptr;
	}
	
	if (bfxs && bfx_num)
	{
		for (size_t bf_id = 0; bf_id < bfx_num; ++bf_id)
		{
			BodyForceAtPcl& bf = bfxs[bf_id];
			p_id = bf.pcl_id;
			pcl_bf[p_id].fx += pcl_m[p_id] * bf.bf;
		}
		clear_bfxs();
	}

	if (bfys && bfy_num)
	{
		for (size_t bf_id = 0; bf_id < bfy_num; ++bf_id)
		{
			BodyForceAtPcl& bf = bfys[bf_id];
			p_id = bf.pcl_id;
			pcl_bf[p_id].fy += pcl_m[p_id] * bf.bf;
		}
		clear_bfys();
	}

	if (txs && tx_num)
	{
		for (size_t t_id = 0; t_id < tx_num; ++t_id)
		{
			TractionBCAtPcl& t = txs[t_id];
			p_id = t.pcl_id;
			pcl_t[p_id].fx += t.t;
		}
		clear_txs();
	}

	if (tys && ty_num)
	{
		for (size_t t_id = 0; t_id < ty_num; ++t_id)
		{
			TractionBCAtPcl& t = tys[t_id];
			p_id = t.pcl_id;
			pcl_t[p_id].fy += t.t;
		}
		clear_tys();
	}

	return 0;
}

int Model_T2D_ME_mt::init_pcls(
	ParticleGenerator2D<TriangleMesh>& pg,
	double density)
{
	int res = init_pcls(pg.get_num(), density, density);
	if (res)
		return res;

	typedef ParticleGenerator2D<TriangleMesh>::Particle PgPcl;
	PgPcl* pg_pcl = pg.first();
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		Position &p_p = pcl_pos[p_id];
		p_p.x = pg_pcl->x;
		p_p.y = pg_pcl->y;
		pcl_m[p_id] *= pg_pcl->area;
		pcl_vol[p_id] = pg_pcl->area;
		Force &bf = pcl_bf[p_id];
		bf.fx *= pcl_vol[p_id];
		bf.fy *= pcl_vol[p_id];
		pg_pcl = pg.next(pg_pcl);
	}

	return 0;
}

void Model_T2D_ME_mt::init_bfxs(
	size_t bf_num,
	const size_t* bf_pcls,
	const double* bfs)
{
	if (pcl_mem_raw && ori_pcl_num)
	{
		size_t p_id;
		for (size_t bf_id = 0; bf_id < bf_num; ++bf_id)
		{
			p_id = bf_pcls[bf_id];
			Force &bf = pcl_bf[p_id];
			bf.fx += bfs[p_id] * pcl_m[p_id];
		}
	}
	else
	{
		init_bfxs(bf_num);
		for (size_t bf_id = 0; bf_id < bfx_num; ++bf_id)
		{
			BodyForceAtPcl& bf = bfxs[bf_id];
			bf.pcl_id = bf_pcls[bf_id];
			bf.bf = bfs[bf_id];
		}
	}
}

void Model_T2D_ME_mt::init_bfys(
	size_t bf_num,
	const size_t* bf_pcls,
	const double* bfs)
{
	if (pcl_mem_raw && ori_pcl_num)
	{
		size_t p_id;
		for (size_t bf_id = 0; bf_id < bf_num; ++bf_id)
		{
			p_id = bf_pcls[bf_id];
			Force& bf = pcl_bf[p_id];
			bf.fy += bfs[p_id] * pcl_m[p_id];
		}
	}
	else
	{
		init_bfys(bf_num);
		for (size_t bf_id = 0; bf_id < bfy_num; ++bf_id)
		{
			BodyForceAtPcl &bf = bfys[bf_id];
			bf.pcl_id = bf_pcls[bf_id];
			bf.bf = bfs[bf_id];
		}
	}
}

void Model_T2D_ME_mt::init_txs(
	size_t t_num,
	const size_t* t_pcls,
	const double* ts)
{
	if (pcl_mem_raw && ori_pcl_num)
	{
		size_t p_id;
		for (size_t t_id = 0; t_id < t_num; ++t_id)
		{
			p_id = t_pcls[t_id];
			Force &t = pcl_t[p_id];
			t.fx += ts[t_id];
		}
	}
	else
	{
		init_txs(t_num);
		for (size_t t_id = 0; t_id < tx_num; ++t_id)
		{
			TractionBCAtPcl &t = txs[t_id];
			t.pcl_id = t_pcls[t_id];
			t.t = ts[t_id];
		}
	}
}

void Model_T2D_ME_mt::init_tys(
	size_t t_num,
	const size_t* t_pcls,
	const double* ts)
{
	if (pcl_mem_raw && ori_pcl_num)
	{
		for (size_t t_id = 0; t_id < t_num; ++t_id)
		{
			size_t p_id = t_pcls[t_id];
			Force &t = pcl_t[p_id];
			t.fy += ts[t_id];
		}
	}
	else
	{
		init_tys(t_num);
		for (size_t t_id = 0; t_id < ty_num; ++t_id)
		{
			TractionBCAtPcl &t = tys[t_id];
			t.pcl_id = t_pcls[t_id];
			t.t = ts[t_id];
		}
	}
}

void Model_T2D_ME_mt::init_fixed_vx_bc(
	size_t bc_num,
	const size_t* bcs)
{
	for (size_t bc_id = 0; bc_id < bc_num; ++bc_id)
	{
		if (bcs[bc_id] < node_num)
			node_has_vbc[bcs[bc_id]].has_vx_bc = true;
	}
}

void Model_T2D_ME_mt::init_fixed_vy_bc(
	size_t bc_num,
	const size_t* bcs)
{
	for (size_t bc_id = 0; bc_id < bc_num; ++bc_id)
	{
		if (bcs[bc_id] < node_num)
			node_has_vbc[bcs[bc_id]].has_vy_bc = true;
	}
}

void Model_T2D_ME_mt::clear_contact_mem()
{
	if (contact_mem)
	{
		delete[] contact_mem;
		contact_mem = nullptr;
	}
}

void Model_T2D_ME_mt::alloc_contact_mem(size_t pcl_num)
{
	clear_contact_mem();

	contact_mem = new char[(sizeof(size_t)
		+ sizeof(Position) + sizeof(Force)) * pcl_num];

	char* cur_mem = contact_mem;
	contact_substep_id = (size_t*)cur_mem;
	cur_mem += sizeof(size_t) * pcl_num;
	prev_contact_pos = (Position*)cur_mem;
	cur_mem += sizeof(Position) * pcl_num;
	prev_contact_tan_force = (Force*)cur_mem;
}
