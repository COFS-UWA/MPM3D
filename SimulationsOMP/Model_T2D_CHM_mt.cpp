#include "SimulationsOMP_pcp.h"

#include <iostream>

#include "TriangleUtils.h"
#include "SearchingGrid2D.hpp"
#include "Model_T2D_CHM_mt.h"

Model_T2D_CHM_mt::Model_T2D_CHM_mt() :
	ori_pcl_num(0), pcl_num(0),
	pcl_mem_raw(nullptr), pcl_mat_model(nullptr),
	node_num(0), elem_num(0),
	mesh_mem_raw(nullptr),
	bfx_num(0), bfxs(nullptr),
	bfy_num(0), bfys(nullptr),
	tx_num(0), txs(nullptr),
	ty_num(0), tys(nullptr),
	grid_x_num(0), grid_y_num(0),
	grid_elem_list(nullptr),
	grid_elem_list_id_array(nullptr),
	rigid_circle_is_valid(false) {}

Model_T2D_CHM_mt::~Model_T2D_CHM_mt()
{
	clear_mesh();
	clear_search_grid();
	clear_pcls();
	clear_bfxs();
	clear_bfys();
	clear_txs();
	clear_tys();
}

Rect Model_T2D_CHM_mt::get_mesh_bbox()
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

void Model_T2D_CHM_mt::clear_mesh()
{
	if (mesh_mem_raw)
	{
		delete[] mesh_mem_raw;
		mesh_mem_raw = nullptr;
	}
}

void Model_T2D_CHM_mt::alloc_mesh(
	size_t n_num,
	size_t e_num
	)
{
	//node_num = n_num;
	//elem_num = e_num;

	//size_t mem_len = (sizeof(ElemNodeIndex) + sizeof(double)
	//	+ sizeof(DShapeFuncAB) + sizeof(DShapeFuncC)
	//	+ sizeof(double) * 4
	//	+ sizeof(ElemStrainInc) + sizeof(ElemStress)
	//	+ sizeof(ElemNodeVM) * 3 + sizeof(ElemNodeForce) * 3
	//	+ sizeof(size_t) * 3 + sizeof(size_t) * 3) * e_num
	//	+ (sizeof(size_t) + sizeof(Position)
	//	+ sizeof(NodeA) + sizeof(NodeV) + sizeof(NodeHasVBC)
	//	+ sizeof(double) * 2) * n_num;
	//mesh_mem_raw = new char[mem_len];

	//char* cur_mem = mesh_mem_raw;
	//elem_node_id = (ElemNodeIndex*)cur_mem;
	//cur_mem += sizeof(ElemNodeIndex) * elem_num;
	//elem_area = (double*)cur_mem;
	//cur_mem += sizeof(double) * elem_num;
	//elem_sf_ab = (DShapeFuncAB*)cur_mem;
	//cur_mem += sizeof(DShapeFuncAB) * elem_num;
	//elem_sf_c = (DShapeFuncC*)cur_mem;
	//cur_mem += sizeof(DShapeFuncC) * elem_num;

	//elem_density = (double*)cur_mem;
	//cur_mem += sizeof(double) * elem_num;
	//elem_pcl_m = (double*)cur_mem;
	//cur_mem += sizeof(double) * elem_num;
	//elem_pcl_vol = (double*)cur_mem;
	//cur_mem += sizeof(double) * elem_num;
	//elem_de = (ElemStrainInc *)cur_mem;
	//cur_mem += sizeof(ElemStrainInc) * elem_num;
	//elem_stress = (ElemStress *)cur_mem;
	//cur_mem += sizeof(ElemStress) * elem_num;
	//elem_m_de_vol = (double*)cur_mem;
	//cur_mem += sizeof(double) * elem_num;

	//elem_node_vm = (ElemNodeVM*)cur_mem;
	//cur_mem += sizeof(ElemNodeVM) * elem_num * 3;
	//elem_node_force = (ElemNodeForce*)cur_mem;
	//cur_mem += sizeof(ElemNodeForce) * elem_num * 3;

	//elem_id_array = (size_t*)cur_mem;
	//cur_mem += sizeof(size_t) * elem_num * 3;
	//node_elem_id_array = (size_t*)cur_mem;
	//cur_mem += sizeof(size_t) * elem_num * 3;
	//node_elem_list = (size_t*)cur_mem;
	//cur_mem += sizeof(size_t) * node_num;

	//node_pos = (Position*)cur_mem;
	//cur_mem += sizeof(Position) * node_num;
	//node_a = (NodeA *)cur_mem;
	//cur_mem += sizeof(NodeA) * node_num;
	//node_v = (NodeV *)cur_mem;
	//cur_mem += sizeof(NodeV) * node_num;
	//node_has_vbc = (NodeHasVBC *)cur_mem;
	//cur_mem += sizeof(NodeHasVBC) * node_num;
	//node_am = (double*)cur_mem;
	//cur_mem += sizeof(double) * node_num;
	//node_de_vol = (double*)cur_mem;
	//cur_mem += sizeof(double) * node_num;
}

void Model_T2D_CHM_mt::init_mesh(const TriangleMesh &mesh)
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
		Position& np = node_pos[n_id];
		np.x = n.x;
		np.y = n.y;
		NodeHasVBC& n_vbc_s = node_has_vbc_s[n_id];
		n_vbc_s.has_vx_bc = false;
		n_vbc_s.has_vy_bc = false;
		NodeHasVBC& n_vbc_f = node_has_vbc_f[n_id];
		n_vbc_f.has_vx_bc = false;
		n_vbc_f.has_vy_bc = false;
	}

	const size_t elem_num3 = elem_num * 3;
	//size_t *elem_id_array_tmp = new size_t[elem_num3 * 2 + node_num];
	//size_t *node_elem_id_array_tmp = elem_id_array_tmp + elem_num3;
	//size_t *node_bin_tmp = node_elem_id_array_tmp + elem_num3;
	//memset(node_bin_tmp, 0, sizeof(size_t) * node_num);

	// init elem connectivity, area and shape functions
	PointInTriangle pit;
	const TriangleMesh::Element *elems = mesh.get_elems();
	size_t e_id3 = 0;
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		const TriangleMesh::Element &e = elems[e_id];
		ElemNodeIndex& eni = elem_node_id[e_id];
		// geometry
		eni.n1 = e.n1;
		eni.n2 = e.n2;
		eni.n3 = e.n3;
		Position &n1_pos = node_pos[e.n1];
		Position &n2_pos = node_pos[e.n2];
		Position &n3_pos = node_pos[e.n3];
		elem_area[e_id] = cal_triangle_area(n1_pos, n2_pos, n3_pos);
		// shape functions
		pit.init_triangle(n1_pos, n2_pos, n3_pos, elem_area[e_id]);
		DShapeFuncAB &esfab = elem_N_ab[e_id];
		esfab.dN1_dx = pit.dN1_dx();
		esfab.dN1_dy = pit.dN1_dy();
		esfab.dN2_dx = pit.dN2_dx();
		esfab.dN2_dy = pit.dN2_dy();
		esfab.dN3_dx = pit.dN3_dx();
		esfab.dN3_dy = pit.dN3_dy();
		DShapeFuncC& esfc = elem_N_c[e_id];
		esfc.c1 = pit.get_coef1();
		esfc.c2 = pit.get_coef2();
		esfc.c3 = pit.get_coef3();
		//// node-element relation
		//elem_id_array_tmp[e_id3] = e_id;
		//elem_id_array_tmp[e_id3 + 1] = e_id;
		//elem_id_array_tmp[e_id3 + 2] = e_id;
		//node_elem_id_array_tmp[e_id3] = e_id3;
		//node_elem_id_array_tmp[e_id3 + 1] = e_id3 + 1;
		//node_elem_id_array_tmp[e_id3 + 2] = e_id3 + 2;
		//e_id3 += 3;
		//++node_bin_tmp[e.n1];
		//++node_bin_tmp[e.n2];
		//++node_bin_tmp[e.n3];
	}

	//node_elem_list[0] = node_bin_tmp[0];
	//for (size_t n_id = 1; n_id < node_num; ++n_id)
	//{
	//	node_bin_tmp[n_id] += node_bin_tmp[n_id - 1];
	//	node_elem_list[n_id] = node_bin_tmp[n_id];
	//}
	//
	//for (size_t e_id = elem_num, e_id3 = elem_num3-1; e_id--; e_id3 -= 3)
	//{
	//	const TriangleMesh::Element& e = elems[e_id];
	//	--node_bin_tmp[e.n3];
	//	elem_id_array[node_bin_tmp[e.n3]] = elem_id_array_tmp[e_id3];
	//	node_elem_id_array[node_bin_tmp[e.n3]] = node_elem_id_array_tmp[e_id3];
	//	--node_bin_tmp[e.n2];
	//	elem_id_array[node_bin_tmp[e.n2]] = elem_id_array_tmp[e_id3 - 1];
	//	node_elem_id_array[node_bin_tmp[e.n2]] = node_elem_id_array_tmp[e_id3 - 1];
	//	--node_bin_tmp[e.n1];
	//	elem_id_array[node_bin_tmp[e.n1]] = elem_id_array_tmp[e_id3 - 2];
	//	node_elem_id_array[node_bin_tmp[e.n1]] = node_elem_id_array_tmp[e_id3 - 2];
	//}

	//delete[] elem_id_array_tmp;
}

void Model_T2D_CHM_mt::clear_search_grid()
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

int Model_T2D_CHM_mt::init_search_grid(
	TriangleMesh& mesh,
	double _hx,
	double _hy
	)
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

	size_t num = grid_x_num * grid_y_num;
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


void Model_T2D_CHM_mt::clear_pcls()
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

void Model_T2D_CHM_mt::alloc_pcls(size_t num)
{
	clear_pcls();

	if (num == 0) return;

	size_t mem_len;
	char* cur_mem;

	ori_pcl_num = num;
	pcl_num = ori_pcl_num;
	mem_len = 0; /*(sizeof(double) + sizeof(Force)
			 + sizeof(Force) + sizeof(PclPos) + sizeof(double)
			+ (sizeof(size_t) + sizeof(double)
			 + sizeof(PclDisp) + sizeof(PclV)
			 + sizeof(PclShapeFunc) + sizeof(PclStress)) * 2
			) * num;*/
	pcl_mem_raw = new char[mem_len];

	//cur_mem = pcl_mem_raw;
	//pcl_m = (double *)cur_mem;
	//cur_mem += sizeof(double) * num;
	//pcl_bf = (PclBodyForce *)(cur_mem);
	//cur_mem += sizeof(PclBodyForce) * num;
	//pcl_t = (Force *)cur_mem;
	//cur_mem += sizeof(Force) * num;
	//pcl_pos = (PclPos *)cur_mem;
	//cur_mem += sizeof(PclPos) * num;
	//pcl_vol = (double *)cur_mem;
	//cur_mem += sizeof(double) * num;

	//PclSortedVarArray &psva0 = pcl_sorted_var_array[0];
	//psva0.pcl_index = (size_t *)cur_mem;
	//cur_mem += sizeof(size_t) * num;
	//psva0.pcl_density = (double *)cur_mem;
	//cur_mem += sizeof(double) * num;
	//psva0.pcl_disp = (PclDisp *)cur_mem;
	//cur_mem += sizeof(PclDisp) * num;
	//psva0.pcl_v = (PclV *)cur_mem;
	//cur_mem += sizeof(PclV) * num;
	//psva0.pcl_N = (PclShapeFunc *)cur_mem;
	//cur_mem += sizeof(PclShapeFunc) * num;
	//psva0.pcl_stress = (PclStress*)cur_mem;
	//cur_mem += sizeof(PclStress) * num;

	//PclSortedVarArray& psva1 = pcl_sorted_var_array[1];
	//psva1.pcl_index = (size_t *)cur_mem;
	//cur_mem += sizeof(size_t) * num;
	//psva1.pcl_density = (double *)cur_mem;
	//cur_mem += sizeof(double) * num;
	//psva1.pcl_disp = (PclDisp*)cur_mem;
	//cur_mem += sizeof(PclDisp) * num;
	//psva1.pcl_v = (PclV*)cur_mem;
	//cur_mem += sizeof(PclV) * num;
	//psva1.pcl_N = (PclShapeFunc*)cur_mem;
	//cur_mem += sizeof(PclShapeFunc) * num;
	//psva1.pcl_stress = (PclStress*)cur_mem;
	//cur_mem += sizeof(PclStress) * num;

	//pcl_mat_model = new MatModel::MaterialModel*[num];
}

void Model_T2D_CHM_mt::alloc_pcls(
	size_t num,
	size_t ori_num
	)
{
	clear_pcls();

	if (num == 0 || ori_num == 0)
		return;

	size_t mem_len;
	char* cur_mem;

	ori_pcl_num = ori_num;
	pcl_num = num;
	mem_len = 0; /*(sizeof(double) + sizeof(PclBodyForce)
			 + sizeof(Force) + sizeof(PclPos)) * ori_num
			 + (sizeof(size_t) + sizeof(double)
			  + sizeof(PclDisp) + sizeof(PclV)
			  + sizeof(PclShapeFunc) + sizeof(PclStress)) * 2 * num;*/
	pcl_mem_raw = new char[mem_len];

	//cur_mem = pcl_mem_raw;
	//pcl_m = (double*)cur_mem;
	//cur_mem += sizeof(double) * ori_num;
	//pcl_bf = (PclBodyForce *)(cur_mem);
	//cur_mem += sizeof(PclBodyForce) * ori_num;
	//pcl_t = (Force *)cur_mem;
	//cur_mem += sizeof(Force) * ori_num;
	//pcl_pos = (PclPos *)cur_mem;
	//cur_mem += sizeof(PclPos) * ori_num;

	//PclSortedVarArray& psva0 = pcl_sorted_var_array[0];
	//psva0.pcl_index = (size_t*)cur_mem;
	//cur_mem += sizeof(size_t) * num;
	//psva0.pcl_density = (double*)cur_mem;
	//cur_mem += sizeof(double) * num;
	//psva0.pcl_disp = (PclDisp*)cur_mem;
	//cur_mem += sizeof(PclDisp) * num;
	//psva0.pcl_v = (PclV*)cur_mem;
	//cur_mem += sizeof(PclV) * num;
	//psva0.pcl_N = (PclShapeFunc*)cur_mem;
	//cur_mem += sizeof(PclShapeFunc) * num;
	//psva0.pcl_stress = (PclStress*)cur_mem;
	//cur_mem += sizeof(PclStress) * num;

	//PclSortedVarArray& psva1 = pcl_sorted_var_array[1];
	//psva1.pcl_index = (size_t*)cur_mem;
	//cur_mem += sizeof(size_t) * num;
	//psva1.pcl_density = (double*)cur_mem;
	//cur_mem += sizeof(double) * num;
	//psva1.pcl_disp = (PclDisp*)cur_mem;
	//cur_mem += sizeof(PclDisp) * num;
	//psva1.pcl_v = (PclV*)cur_mem;
	//cur_mem += sizeof(PclV) * num;
	//psva1.pcl_N = (PclShapeFunc*)cur_mem;
	//cur_mem += sizeof(PclShapeFunc) * num;
	//psva1.pcl_stress = (PclStress*)cur_mem;
	//cur_mem += sizeof(PclStress) * num;

	//pcl_mat_model = new MatModel::MaterialModel * [ori_num];
}

int Model_T2D_CHM_mt::init_pcls(size_t num, double m, double density)
{
	alloc_pcls(num);

	SortedPclVarArrays &spva0 = sorted_pcl_var_arrays[0];
	size_t p_id;
	//for (p_id = 0; p_id < num; ++p_id)
	//{
	//	pcl_m[p_id] = m;
	//	PclBodyForce &p_bf = pcl_bf[p_id];
	//	p_bf.bfx = 0.0;
	//	p_bf.bfy = 0.0;
	//	Force& p_t = pcl_t[p_id];
	//	p_t.tx = 0.0;
	//	p_t.ty = 0.0;
	//	psva0.pcl_index[p_id] = p_id;
	//	psva0.pcl_density[p_id] = density;
	//	PclV& p_v = psva0.pcl_v[p_id];
	//	p_v.vx = 0.0;
	//	p_v.vy = 0.0;
	//	PclStress& p_s = psva0.pcl_stress[p_id];
	//	p_s.s11 = 0.0;
	//	p_s.s22 = 0.0;
	//	p_s.s12 = 0.0;
	//	pcl_mat_model[p_id] = nullptr;
	//}
	//
	//if (bfxs && bfx_num)
	//{
	//	for (size_t bf_id = 0; bf_id < bfx_num; ++bf_id)
	//	{
	//		BodyForceAtPcl& bf = bfxs[bf_id];
	//		p_id = bf.pcl_id;
	//		pcl_bf[p_id].bfx += pcl_m[p_id] * bf.bf;
	//	}
	//	clear_bfxs();
	//}

	//if (bfys && bfy_num)
	//{
	//	for (size_t bf_id = 0; bf_id < bfy_num; ++bf_id)
	//	{
	//		BodyForceAtPcl& bf = bfys[bf_id];
	//		p_id = bf.pcl_id;
	//		pcl_bf[p_id].bfy += pcl_m[p_id] * bf.bf;
	//	}
	//	clear_bfys();
	//}

	//if (txs && tx_num)
	//{
	//	for (size_t t_id = 0; t_id < tx_num; ++t_id)
	//	{
	//		TractionBCAtPcl& t = txs[t_id];
	//		p_id = t.pcl_id;
	//		pcl_t[p_id].tx += t.t;
	//	}
	//	clear_txs();
	//}

	//if (tys && ty_num)
	//{
	//	for (size_t t_id = 0; t_id < ty_num; ++t_id)
	//	{
	//		TractionBCAtPcl& t = tys[t_id];
	//		p_id = t.pcl_id;
	//		pcl_t[p_id].ty += t.t;
	//	}
	//	clear_tys();
	//}

	return 0;
}

int Model_T2D_CHM_mt::init_pcls(
	ParticleGenerator2D<TriangleMesh>& pg,
	double density
	)
{
	int res = init_pcls(pg.get_num(), density, density);
	if (res)
		return res;

	typedef ParticleGenerator2D<TriangleMesh>::Particle PgPcl;
	PgPcl* pg_pcl = pg.first();
	//for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	//{
	//	PclPos &p_p = pcl_pos[p_id];
	//	p_p.x = pg_pcl->x;
	//	p_p.y = pg_pcl->y;
	//	pcl_m[p_id] *= pg_pcl->area;
	//	pg_pcl = pg.next(pg_pcl);
	//}

	return 0;
}

void Model_T2D_CHM_mt::init_bfxs(
	size_t bf_num,
	const size_t* bf_pcls,
	const double* bfs
	)
{
	//if (pcl_mem_raw && ori_pcl_num)
	//{
	//	size_t p_id;
	//	for (size_t bf_id = 0; bf_id < bf_num; ++bf_id)
	//	{
	//		p_id = bf_pcls[bf_id];
	//		PclBodyForce &bf = pcl_bf[p_id];
	//		bf.bfx += bfs[p_id] * pcl_m[p_id];
	//	}
	//}
	//else
	//{
	//	init_bfxs(bf_num);
	//	for (size_t bf_id = 0; bf_id < bfx_num; ++bf_id)
	//	{
	//		BodyForceAtPcl& bf = bfxs[bf_id];
	//		bf.pcl_id = bf_pcls[bf_id];
	//		bf.bf = bfs[bf_id];
	//	}
	//}
}

void Model_T2D_CHM_mt::init_bfys(
	size_t bf_num,
	const size_t* bf_pcls,
	const double* bfs
	)
{
	//if (pcl_mem_raw && ori_pcl_num)
	//{
	//	size_t p_id;
	//	for (size_t bf_id = 0; bf_id < bf_num; ++bf_id)
	//	{
	//		p_id = bf_pcls[bf_id];
	//		PclBodyForce& bf = pcl_bf[p_id];
	//		bf.bfy += bfs[p_id] * pcl_m[p_id];
	//	}
	//}
	//else
	//{
	//	init_bfys(bf_num);
	//	for (size_t bf_id = 0; bf_id < bfy_num; ++bf_id)
	//	{
	//		BodyForceAtPcl &bf = bfys[bf_id];
	//		bf.pcl_id = bf_pcls[bf_id];
	//		bf.bf = bfs[bf_id];
	//	}
	//}
}

void Model_T2D_CHM_mt::init_txs(
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

void Model_T2D_CHM_mt::init_tys(
	size_t t_num,
	const size_t* t_pcls,
	const double* ts
	)
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

void Model_T2D_CHM_mt::init_fixed_vx_s_bc(
	size_t bc_num,
	const size_t* bcs
	)
{
	for (size_t bc_id = 0; bc_id < bc_num; ++bc_id)
	{
		if (bcs[bc_id] < node_num)
			node_has_vbc_s[bcs[bc_id]].has_vx_bc = true;
	}
}

void Model_T2D_CHM_mt::init_fixed_vy_s_bc(
	size_t bc_num,
	const size_t* bcs
	)
{
	for (size_t bc_id = 0; bc_id < bc_num; ++bc_id)
	{
		if (bcs[bc_id] < node_num)
			node_has_vbc_s[bcs[bc_id]].has_vy_bc = true;
	}
}

void Model_T2D_CHM_mt::init_fixed_vx_f_bc(
	size_t bc_num,
	const size_t* bcs
	)
{
	for (size_t bc_id = 0; bc_id < bc_num; ++bc_id)
	{
		if (bcs[bc_id] < node_num)
			node_has_vbc_f[bcs[bc_id]].has_vx_bc = true;
	}
}

void Model_T2D_CHM_mt::init_fixed_vy_f_bc(
	size_t bc_num,
	const size_t* bcs
	)
{
	for (size_t bc_id = 0; bc_id < bc_num; ++bc_id)
	{
		if (bcs[bc_id] < node_num)
			node_has_vbc_f[bcs[bc_id]].has_vy_bc = true;
	}
}
