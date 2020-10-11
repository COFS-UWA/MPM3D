#include "SimulationsOMP_pcp.h"

#include "TriangleUtils.h"
#include "SearchingGrid2D.hpp"
#include "Model_T2D_ME_mt.h"

Model_T2D_ME_mt::Model_T2D_ME_mt() :
	ori_pcl_num(0), pcl_num(0),
	pcl_mem_raw(nullptr), pcl_mat_model(nullptr),
	mesh_mem_raw(nullptr),
	vx_bc_num(0), vx_bcs(nullptr),
	vy_bc_num(0), vy_bcs(nullptr),
	bfx_num(0), bfxs(nullptr),
	bfy_num(0), bfys(nullptr),
	tx_num(0), txs(nullptr),
	ty_num(0), tys(nullptr),
	grid_x_num(0), grid_y_num(0),
	grids(nullptr), grid_elem_id_array(nullptr) {}

Model_T2D_ME_mt::~Model_T2D_ME_mt()
{
	clear_mesh();
	clear_search_grid();
	clear_pcls();
	clear_bfxs();
	clear_bfys();
	clear_txs();
	clear_tys();
	clear_vx_bcs();
	clear_vy_bcs();
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
	size_t e_num
	)
{
	node_num = uint32_t(n_num);
	elem_num = uint32_t(e_num);

	size_t mem_len = (sizeof(ElemNodeIndex) + sizeof(float)
		+ sizeof(ElemShapeFuncAB) + sizeof(ElemShapeFuncC)
		+ sizeof(uint32_t) * 2 + sizeof(float)
		+ sizeof(ElemStrainInc) + sizeof(ElemStress)
		+ sizeof(float) + sizeof(float)
		+ sizeof(ElemNodeVM) * 3 + sizeof(ElemNodeForce) * 3
		+ sizeof(uint32_t) * 7) * e_num
		+ (sizeof(NodePos) + sizeof(float) * 6) * n_num;
	mesh_mem_raw = new char[mem_len];

	char* cur_mem = mesh_mem_raw;
	elem_node_id = (ElemNodeIndex*)cur_mem;
	cur_mem += sizeof(ElemNodeIndex) * elem_num;
	elem_area = (float*)cur_mem;
	cur_mem += sizeof(float) * elem_num;
	elem_sf_ab = (ElemShapeFuncAB*)cur_mem;
	cur_mem += sizeof(ElemShapeFuncAB) * elem_num;
	elem_sf_c = (ElemShapeFuncC*)cur_mem;
	cur_mem += sizeof(ElemShapeFuncC) * elem_num;

	// elem pcl list 0
	pcl_sorted_var_array[0].elem_has_pcl_num = (uint32_t*)cur_mem;
	cur_mem += sizeof(uint32_t) * elem_num;
	// elem pcl list 1
	pcl_sorted_var_array[1].elem_has_pcl_num = (uint32_t*)cur_mem;
	cur_mem += sizeof(uint32_t) * elem_num;

	elem_density = (float*)cur_mem;
	cur_mem += sizeof(float) * elem_num;
	elem_de = (ElemStrainInc*)cur_mem;
	cur_mem += sizeof(ElemStrainInc) * elem_num;
	elem_stress = (ElemStress*)cur_mem;
	cur_mem += sizeof(ElemStress) * elem_num;

	elem_am = (float*)cur_mem;
	cur_mem += sizeof(float) * elem_num;
	elem_am_de_vol = (float*)cur_mem;
	cur_mem += sizeof(float) * elem_num;

	elem_node_vm = (ElemNodeVM*)cur_mem;
	cur_mem += sizeof(ElemNodeVM) * elem_num * 3;
	elem_node_force = (ElemNodeForce*)cur_mem;
	cur_mem += sizeof(ElemNodeForce) * elem_num * 3;

	elem_id_array = (uint32_t*)cur_mem;
	cur_mem += sizeof(uint32_t) * elem_num * 3;
	node_elem_id_array = (uint32_t*)cur_mem;
	cur_mem += sizeof(uint32_t) * elem_num * 3;
	node_elem_list = (uint32_t*)cur_mem;
	cur_mem += sizeof(uint32_t) * node_num;

	node_pos = (NodePos*)cur_mem;
	cur_mem += sizeof(NodePos) * node_num;
	node_ax = (float*)cur_mem;
	cur_mem += sizeof(float) * node_num;
	node_ay = (float*)cur_mem;
	cur_mem += sizeof(float) * node_num;
	node_vx = (float*)cur_mem;
	cur_mem += sizeof(float) * node_num;
	node_vy = (float*)cur_mem;
	cur_mem += sizeof(float) * node_num;
	node_am = (float*)cur_mem;
	cur_mem += sizeof(float) * node_num;
	node_de_vol = (float*)cur_mem;
	cur_mem += sizeof(float) * node_num;
}

void Model_T2D_ME_mt::init_mesh(const TriangleMesh &mesh)
{
	clear_mesh();

	if (mesh.get_elem_num() == 0 ||
		mesh.get_node_num() == 0)
		return;
	
	alloc_mesh(mesh.get_node_num(), mesh.get_elem_num());

	uint32_t elem_num3 = elem_num * 3;
	uint32_t *elem_id_array_tmp = new uint32_t[size_t(elem_num3) * 2 + size_t(node_num)];
	uint32_t *node_elem_id_array_tmp = elem_id_array_tmp + elem_num3;
	uint32_t *node_bin_tmp = node_elem_id_array_tmp + elem_num3;
	
	memset(node_bin_tmp, 0, sizeof(uint32_t) * node_num);

	PointInTriangle pit;
	const TriangleMesh::Element *elems = mesh.get_elems();
	const TriangleMesh::Node *nodes = mesh.get_nodes();
	uint32_t e_id, n_id, e_id3 = 0;
	for (e_id = 0; e_id < elem_num; ++e_id)
	{
		const TriangleMesh::Element &e = elems[e_id];
		ElemNodeIndex& eni = elem_node_id[e_id];
		// geometry
		eni.n1 = uint32_t(e.n1);
		eni.n2 = uint32_t(e.n2);
		eni.n3 = uint32_t(e.n3);
		elem_area[e_id] = float(e.area);
		// shape functions
		const TriangleMesh::Node& n1 = nodes[e.n1];
		const TriangleMesh::Node& n2 = nodes[e.n2];
		const TriangleMesh::Node& n3 = nodes[e.n3];
		pit.init_triangle(n1, n2, n3, e.area);
		ElemShapeFuncAB &esfab = elem_sf_ab[e_id];
		esfab.dN1_dx = float(pit.dN1_dx());
		esfab.dN1_dy = float(pit.dN1_dy());
		esfab.dN2_dx = float(pit.dN2_dx());
		esfab.dN2_dy = float(pit.dN2_dy());
		esfab.dN3_dx = float(pit.dN3_dx());
		esfab.dN3_dy = float(pit.dN3_dy());
		ElemShapeFuncC& esfc = elem_sf_c[e_id];
		esfc.c1 = float(pit.get_coef1());
		esfc.c2 = float(pit.get_coef2());
		esfc.c3 = float(pit.get_coef3());
		// node-element relation
		elem_id_array_tmp[e_id3] = e_id;
		elem_id_array_tmp[e_id3 + 1] = e_id;
		elem_id_array_tmp[e_id3 + 2] = e_id;
		node_elem_id_array_tmp[e_id3] = e_id3;
		node_elem_id_array_tmp[e_id3 + 1] = e_id3 + 1;
		node_elem_id_array_tmp[e_id3 + 2] = e_id3 + 2;
		e_id3 += 3;
		++node_bin_tmp[e.n1];
		++node_bin_tmp[e.n2];
		++node_bin_tmp[e.n3];
	}

	node_elem_list[0] = node_bin_tmp[0];
	for (n_id = 1; n_id < node_num; ++n_id)
	{
		node_bin_tmp[n_id] += node_bin_tmp[n_id - 1];
		node_elem_list[n_id] = node_bin_tmp[n_id];
	}

	e_id3 = 0;
	for (e_id = elem_num; e_id--;)
	{
		const TriangleMesh::Element& e = elems[e_id];
		--node_bin_tmp[e.n1];
		elem_id_array[node_bin_tmp[e.n1]] = elem_id_array_tmp[e_id3];
		node_elem_id_array[node_bin_tmp[e.n1]] = node_elem_id_array_tmp[e_id3];
		--node_bin_tmp[e.n2];
		elem_id_array[node_bin_tmp[e.n2]] = elem_id_array_tmp[e_id3 + 1];
		node_elem_id_array[node_bin_tmp[e.n2]] = node_elem_id_array_tmp[e_id3 + 1];
		--node_bin_tmp[e.n3];
		elem_id_array[node_bin_tmp[e.n3]] = elem_id_array_tmp[e_id3 + 2];
		node_elem_id_array[node_bin_tmp[e.n3]] = node_elem_id_array_tmp[e_id3 + 2];
		e_id3 += 3;
	}

	delete[] elem_id_array_tmp;

	for (n_id = 0; n_id < node_num; ++n_id)
	{
		const TriangleMesh::Node& n = nodes[n_id];
		NodePos& np = node_pos[n_id];
		np.x = float(n.x);
		np.y = float(n.y);
	}
}

void Model_T2D_ME_mt::clear_search_grid()
{
	if (grids)
	{
		delete[] grids;
		grids = nullptr;
	}
	if (grid_elem_id_array)
	{
		delete[] grid_elem_id_array;
		grid_elem_id_array = nullptr;
	}
	grid_x_num = 0;
	grid_y_num = 0;
}

int Model_T2D_ME_mt::init_search_grid(
	const TriangleMesh& mesh,
	double _hx,
	double _hy
	)
{
	clear_search_grid();
	
	SearchingGrid2D<TriangleMesh> search_grid;
	search_grid.init(mesh, _hx, _hy);

	grid_xl = float(search_grid.get_x_min());
	grid_yl = float(search_grid.get_y_min());
	grid_hx = float(search_grid.get_hx());
	grid_hy = float(search_grid.get_hy());
	grid_x_num = uint32_t(search_grid.get_x_num());
	grid_y_num = uint32_t(search_grid.get_y_num());

	size_t num = size_t(grid_x_num) * size_t(grid_y_num);
	grids = new MeshBgGrid[num];
	typedef SearchingGrid2D<TriangleMesh>::Grid Grid;
	Grid *sg = search_grid.get_grids();

	uint32_t cur_id = 0;
	for (size_t g_id = 0; g_id < num; ++g_id)
	{
		MeshBgGrid& bg_grid = grids[g_id];
		Grid& grid = sg[g_id];
		bg_grid.elem_start_id = cur_id;
		for (auto pe = grid.pelems; pe; pe = pe->next)
			++cur_id;
		bg_grid.elem_end_id = cur_id;
	}

	grid_elem_id_array = new uint32_t[cur_id];
	cur_id = 0;
	for (size_t g_id = 0; g_id < num; ++g_id)
	{
		MeshBgGrid& bg_grid = grids[g_id];
		Grid& grid = sg[g_id];
		for (auto pe = grid.pelems; pe; pe = pe->next)
		{
			grid_elem_id_array[cur_id] = pe->e->id;
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
	if (pcl_mat_model)
	{
		delete[] pcl_mat_model;
		pcl_mat_model = nullptr;
	}
	ori_pcl_num = 0;
	pcl_num = 0;
}

void Model_T2D_ME_mt::alloc_pcls(size_t num)
{
	clear_pcls();

	if (num == 0) return;

	size_t mem_len;
	char* cur_mem;

	ori_pcl_num = uint32_t(num);
	pcl_num = ori_pcl_num;
	mem_len = (sizeof(float) + sizeof(PclBodyForce)
			 + sizeof(PclTraction) + sizeof(PclPos)
			+ (sizeof(uint32_t) + sizeof(float)
			 + sizeof(PclDisp) + sizeof(PclV)
			 + sizeof(PclShapeFunc) + sizeof(PclStress)) * 2
			) * num;
	pcl_mem_raw = new char[mem_len];

	cur_mem = pcl_mem_raw;
	pcl_m = (float *)cur_mem;
	cur_mem += sizeof(float) * num;
	pcl_bf = (PclBodyForce *)(cur_mem);
	cur_mem += sizeof(PclBodyForce) * num;
	pcl_t = (PclTraction *)cur_mem;
	cur_mem += sizeof(PclTraction) * num;
	pcl_pos = (PclPos *)cur_mem;
	cur_mem += sizeof(PclPos) * num;

	PclSortedVarArray &psva0 = pcl_sorted_var_array[0];
	psva0.pcl_index = (uint32_t *)cur_mem;
	cur_mem += sizeof(uint32_t) * num;
	psva0.pcl_density = (float *)cur_mem;
	cur_mem += sizeof(float) * num;
	psva0.pcl_disp = (PclDisp *)cur_mem;
	cur_mem += sizeof(PclDisp) * num;
	psva0.pcl_v = (PclV *)cur_mem;
	cur_mem += sizeof(PclV) * num;
	psva0.pcl_N = (PclShapeFunc *)cur_mem;
	cur_mem += sizeof(PclShapeFunc) * num;
	psva0.pcl_stress = (PclStress*)cur_mem;
	cur_mem += sizeof(PclStress) * num;

	PclSortedVarArray& psva1 = pcl_sorted_var_array[1];
	psva1.pcl_index = (uint32_t *)cur_mem;
	cur_mem += sizeof(uint32_t) * num;
	psva1.pcl_density = (float *)cur_mem;
	cur_mem += sizeof(float) * num;
	psva1.pcl_disp = (PclDisp*)cur_mem;
	cur_mem += sizeof(PclDisp) * num;
	psva1.pcl_v = (PclV*)cur_mem;
	cur_mem += sizeof(PclV) * num;
	psva1.pcl_N = (PclShapeFunc*)cur_mem;
	cur_mem += sizeof(PclShapeFunc) * num;
	psva1.pcl_stress = (PclStress*)cur_mem;
	cur_mem += sizeof(PclStress) * num;

	pcl_mat_model = new MatModel::MaterialModel*[num];
}

void Model_T2D_ME_mt::alloc_pcls(
	size_t num,
	size_t ori_num
	)
{
	clear_pcls();

	if (num == 0 || ori_num == 0)
		return;

	size_t mem_len;
	char* cur_mem;

	ori_pcl_num = uint32_t(ori_num);
	pcl_num = uint32_t(num);
	mem_len = (sizeof(float) + sizeof(PclBodyForce)
			 + sizeof(PclTraction) + sizeof(PclPos)) * ori_num
			 + (sizeof(uint32_t) + sizeof(float)
			  + sizeof(PclDisp) + sizeof(PclV)
			  + sizeof(PclShapeFunc) + sizeof(PclStress)) * 2 * num;
	pcl_mem_raw = new char[mem_len];

	cur_mem = pcl_mem_raw;
	pcl_m = (float*)cur_mem;
	cur_mem += sizeof(float) * ori_num;
	pcl_bf = (PclBodyForce *)(cur_mem);
	cur_mem += sizeof(PclBodyForce) * ori_num;
	pcl_t = (PclTraction *)cur_mem;
	cur_mem += sizeof(PclTraction) * ori_num;
	pcl_pos = (PclPos *)cur_mem;
	cur_mem += sizeof(PclPos) * ori_num;

	PclSortedVarArray& psva0 = pcl_sorted_var_array[0];
	psva0.pcl_index = (uint32_t*)cur_mem;
	cur_mem += sizeof(uint32_t) * num;
	psva0.pcl_density = (float*)cur_mem;
	cur_mem += sizeof(float) * num;
	psva0.pcl_disp = (PclDisp*)cur_mem;
	cur_mem += sizeof(PclDisp) * num;
	psva0.pcl_v = (PclV*)cur_mem;
	cur_mem += sizeof(PclV) * num;
	psva0.pcl_N = (PclShapeFunc*)cur_mem;
	cur_mem += sizeof(PclShapeFunc) * num;
	psva0.pcl_stress = (PclStress*)cur_mem;
	cur_mem += sizeof(PclStress) * num;

	PclSortedVarArray& psva1 = pcl_sorted_var_array[1];
	psva1.pcl_index = (uint32_t*)cur_mem;
	cur_mem += sizeof(uint32_t) * num;
	psva1.pcl_density = (float*)cur_mem;
	cur_mem += sizeof(float) * num;
	psva1.pcl_disp = (PclDisp*)cur_mem;
	cur_mem += sizeof(PclDisp) * num;
	psva1.pcl_v = (PclV*)cur_mem;
	cur_mem += sizeof(PclV) * num;
	psva1.pcl_N = (PclShapeFunc*)cur_mem;
	cur_mem += sizeof(PclShapeFunc) * num;
	psva1.pcl_stress = (PclStress*)cur_mem;
	cur_mem += sizeof(PclStress) * num;

	pcl_mat_model = new MatModel::MaterialModel * [ori_num];
}

int Model_T2D_ME_mt::init_pcls(size_t num, double m, double density)
{
	alloc_pcls(num);

	PclSortedVarArray& psva0 = pcl_sorted_var_array[0];
	PclSortedVarArray& psva1 = pcl_sorted_var_array[1];

	size_t p_id;
	for (p_id = 0; p_id < num; ++p_id)
	{
		pcl_m[p_id] = float(m);
		PclBodyForce &p_bf = pcl_bf[p_id];
		p_bf.bfx = 0.0f;
		p_bf.bfy = 0.0f;
		PclTraction& p_t = pcl_t[p_id];
		p_t.tx = 0.0f;
		p_t.ty = 0.0f;
		psva0.pcl_index[p_id] = uint32_t(p_id);
		psva0.pcl_density[p_id] = float(density);
		PclV& p_v = psva0.pcl_v[p_id];
		p_v.vx = 0.0f;
		p_v.vy = 0.0f;
		PclStress& p_s = psva0.pcl_stress[p_id];
		p_s.s11 = 0.0f;
		p_s.s22 = 0.0f;
		p_s.s12 = 0.0f;
		pcl_mat_model[p_id] = nullptr;
	}
	
	if (bfxs && bfx_num)
	{
		for (size_t bf_id = 0; bf_id < bfx_num; ++bf_id)
		{
			BodyForceAtPcl& bf = bfxs[bf_id];
			p_id = bf.pcl_id;
			pcl_bf[p_id].bfx += pcl_m[p_id] * float(bf.bf);
		}
		clear_bfxs();
	}

	if (bfys && bfy_num)
	{
		for (size_t bf_id = 0; bf_id < bfy_num; ++bf_id)
		{
			BodyForceAtPcl& bf = bfys[bf_id];
			p_id = bf.pcl_id;
			pcl_bf[p_id].bfy += pcl_m[p_id] * float(bf.bf);
		}
		clear_bfys();
	}

	if (txs && tx_num)
	{
		for (size_t t_id = 0; t_id < tx_num; ++t_id)
		{
			TractionBCAtPcl& t = txs[t_id];
			p_id = t.pcl_id;
			pcl_t[p_id].tx += t.t;
		}
		clear_txs();
	}

	if (tys && ty_num)
	{
		for (size_t t_id = 0; t_id < ty_num; ++t_id)
		{
			TractionBCAtPcl& t = tys[t_id];
			p_id = t.pcl_id;
			pcl_t[p_id].ty += t.t;
		}
		clear_tys();
	}

	return 0;
}

int Model_T2D_ME_mt::init_pcls(
	ParticleGenerator2D<TriangleMesh>& pg,
	double density
	)
{
	int res = init_pcls(pg.get_num(), density, density);
	if (res)
		return res;

	typedef ParticleGenerator2D<TriangleMesh>::Particle PgPcl;
	PgPcl* pg_pcl = pg.first();
	for (uint32_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		PclPos &p_p = pcl_pos[p_id];
		p_p.x = float(pg_pcl->x);
		p_p.y = float(pg_pcl->y);
		pcl_m[p_id] *= float(pg_pcl->area);
		pg_pcl = pg.next(pg_pcl);
	}

	return 0;
}

void Model_T2D_ME_mt::init_bfxs(
	size_t bf_num,
	const size_t* bf_pcls,
	const double* bfs
	)
{
	if (pcl_mem_raw && ori_pcl_num)
	{
		size_t p_id;
		for (size_t bf_id = 0; bf_id < bf_num; ++bf_id)
		{
			p_id = bf_pcls[bf_id];
			PclBodyForce &bf = pcl_bf[p_id];
			bf.bfx += float(bfs[p_id]) * pcl_m[p_id];
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
	const double* bfs
	)
{
	if (pcl_mem_raw && ori_pcl_num)
	{
		size_t p_id;
		for (size_t bf_id = 0; bf_id < bf_num; ++bf_id)
		{
			p_id = bf_pcls[bf_id];
			PclBodyForce& bf = pcl_bf[p_id];
			bf.bfy += float(bfs[p_id]) * pcl_m[p_id];
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
	const double* ts
	)
{
	if (pcl_mem_raw && ori_pcl_num)
	{
		size_t p_id;
		for (size_t t_id = 0; t_id < t_num; ++t_id)
		{
			p_id = t_pcls[t_id];
			PclTraction &t = pcl_t[p_id];
			t.tx += float(ts[t_id]);
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
	const double* ts
	)
{
	if (pcl_mem_raw && ori_pcl_num)
	{
		size_t p_id;
		for (size_t t_id = 0; t_id < t_num; ++t_id)
		{
			p_id = t_pcls[t_id];
			PclTraction &t = pcl_t[p_id];
			t.ty += ts[t_id];
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

void Model_T2D_ME_mt::clear_vx_bcs()
{
	if (vx_bcs)
	{
		delete[] vx_bcs;
		vx_bcs = nullptr;
	}
	vx_bc_num = 0;
}

void Model_T2D_ME_mt::alloc_vx_bcs(size_t num)
{
	clear_vx_bcs();
	vx_bc_num = num;
	vx_bcs = new uint32_t[num];
}

void Model_T2D_ME_mt::clear_vy_bcs()
{
	if (vy_bcs)
	{
		delete[] vy_bcs;
		vy_bcs = nullptr;
	}
	vy_bc_num = 0;
}

void Model_T2D_ME_mt::alloc_vy_bcs(size_t num)
{
	clear_vy_bcs();
	vy_bc_num = num;
	vy_bcs = new uint32_t[num];
}

void Model_T2D_ME_mt::init_fixed_vx_bc(
	size_t bc_num,
	const size_t* bcs
	)
{
	alloc_vx_bcs(bc_num);
	for (size_t bc_id = 0; bc_id < bc_num; ++bc_id)
		vx_bcs[bc_id] = uint32_t(bcs[bc_id]);
}

void Model_T2D_ME_mt::init_fixed_vy_bf(
	size_t bc_num,
	const size_t* bcs
	)
{
	alloc_vy_bcs(bc_num);
	for (size_t bc_id = 0; bc_id < bc_num; ++bc_id)
		vy_bcs[bc_id] = uint32_t(bcs[bc_id]);
}
