#include "SimulationsOMP_pcp.h"

#include <iostream>

#include "TetrahedronUtils.h"
#include "SearchingGrid3D.hpp"
#include "Model_T3D_ME_mt.h"

Model_T3D_ME_mt::Model_T3D_ME_mt() :
	ori_pcl_num(0), pcl_num(0),
	pcl_mem_raw(nullptr), pcl_mat_model(nullptr),
	node_num(0), elem_num(0),
	mesh_mem_raw(nullptr),
	bfx_num(0), bfxs(nullptr),
	bfy_num(0), bfys(nullptr),
	bfz_num(0), bfzs(nullptr),
	tx_num(0), txs(nullptr),
	ty_num(0), tys(nullptr),
	tz_num(0), tzs(nullptr),
	grid_x_num(0), grid_y_num(0), grid_z_num(0),
	grid_elem_list(nullptr),
	grid_elem_list_id_array(nullptr) {}

Model_T3D_ME_mt::~Model_T3D_ME_mt()
{
	clear_mesh();
	clear_search_grid();
	clear_pcls();
	clear_bfxs();
	clear_bfys();
	clear_bfzs();
	clear_txs();
	clear_tys();
	clear_tzs();
}

Cube Model_T3D_ME_mt::get_mesh_bbox()
{
	if (!node_num)
		return Cube(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	
	Position& np0 = node_pos[0];
	Cube res(np0.x, np0.x, np0.y, np0.y, np0.z, np0.z);
	for (size_t n_id = 1; n_id < node_num; ++n_id)
	{
		Position& np = node_pos[n_id];
		if (res.xl > np.x)
			res.xl = np.x;
		if (res.xu < np.x)
			res.xu = np.x;
		if (res.yl > np.y)
			res.yl = np.y;
		if (res.yu < np.y)
			res.yu = np.y;
		if (res.zl > np.z)
			res.zl = np.z;
		if (res.zu < np.z)
			res.zu = np.z;
	}
	return res;
}

void Model_T3D_ME_mt::clear_mesh()
{
	if (mesh_mem_raw)
	{
		delete[] mesh_mem_raw;
		mesh_mem_raw = nullptr;
	}
}

void Model_T3D_ME_mt::alloc_mesh(
	size_t n_num,
	size_t e_num
	)
{
	node_num = n_num;
	elem_num = e_num;

	size_t mem_len = (sizeof(ElemNodeIndex) + sizeof(double)
		+ sizeof(DShapeFuncABC) + sizeof(DShapeFuncD)
		+ sizeof(double) * 4
		+ sizeof(StrainInc) + sizeof(Stress)
		+ sizeof(ElemNodeVM) * 4 + sizeof(ElemNodeForce) * 4
		+ sizeof(size_t) * 4 + sizeof(size_t) * 4) * e_num
		+ (sizeof(size_t) + sizeof(Position)
		+ sizeof(Acceleration) + sizeof(Velocity) + sizeof(NodeHasVBC)
		+ sizeof(double) * 2) * n_num;
	mesh_mem_raw = new char[mem_len];

	char* cur_mem = mesh_mem_raw;
	elem_node_id = (ElemNodeIndex*)cur_mem;
	cur_mem += sizeof(ElemNodeIndex) * elem_num;
	elem_vol = (double*)cur_mem;
	cur_mem += sizeof(double) * elem_num;
	elem_dN_abc = (DShapeFuncABC*)cur_mem;
	cur_mem += sizeof(DShapeFuncABC) * elem_num;
	elem_dN_d = (DShapeFuncD*)cur_mem;
	cur_mem += sizeof(DShapeFuncD) * elem_num;

	elem_density = (double*)cur_mem;
	cur_mem += sizeof(double) * elem_num;
	elem_pcl_m = (double*)cur_mem;
	cur_mem += sizeof(double) * elem_num;
	elem_pcl_vol = (double*)cur_mem;
	cur_mem += sizeof(double) * elem_num;
	elem_de = (StrainInc *)cur_mem;
	cur_mem += sizeof(StrainInc) * elem_num;
	elem_stress = (Stress *)cur_mem;
	cur_mem += sizeof(Stress) * elem_num;
	elem_m_de_vol = (double*)cur_mem;
	cur_mem += sizeof(double) * elem_num;

	size_t elem_num4 = elem_num * 4;
	elem_node_vm = (ElemNodeVM*)cur_mem;
	cur_mem += sizeof(ElemNodeVM) * elem_num4;
	elem_node_force = (ElemNodeForce*)cur_mem;
	cur_mem += sizeof(ElemNodeForce) * elem_num4;

	elem_id_array = (size_t*)cur_mem;
	cur_mem += sizeof(size_t) * elem_num4;
	node_elem_id_array = (size_t*)cur_mem;
	cur_mem += sizeof(size_t) * elem_num4;
	node_elem_list = (size_t*)cur_mem;
	cur_mem += sizeof(size_t) * node_num;

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

void Model_T3D_ME_mt::init_mesh(const TetrahedronMesh &mesh)
{
	clear_mesh();

	if (mesh.get_elem_num() == 0 ||
		mesh.get_node_num() == 0)
		return;
	
	alloc_mesh(mesh.get_node_num(), mesh.get_elem_num());

	// init node coordinates
	const TetrahedronMesh::Node* nodes = mesh.get_nodes();
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		const TetrahedronMesh::Node& n = nodes[n_id];
		Position& np = node_pos[n_id];
		np.x = n.x;
		np.y = n.y;
		np.z = n.z;
		NodeHasVBC& n_vbc = node_has_vbc[n_id];
		n_vbc.has_vx_bc = false;
		n_vbc.has_vy_bc = false;
		n_vbc.has_vz_bc = false;
	}

	size_t elem_num4 = elem_num * 4;
	size_t *elem_id_array_tmp = new size_t[elem_num4 * 2 + node_num];
	size_t *node_elem_id_array_tmp = elem_id_array_tmp + elem_num4;
	size_t *node_bin_tmp = node_elem_id_array_tmp + elem_num4;
	
	memset(node_bin_tmp, 0, sizeof(size_t) * node_num);

	// init elem connectivity, area and shape functions
	PointInTetrahedron pit;
	const TetrahedronMesh::Element *elems = mesh.get_elems();
	size_t e_id4 = 0;
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		const TetrahedronMesh::Element &e = elems[e_id];
		ElemNodeIndex& eni = elem_node_id[e_id];
		// geometry
		eni.n1 = e.n1;
		eni.n2 = e.n2;
		eni.n3 = e.n3;
		eni.n4 = e.n4;
		Position &n1_pos = node_pos[e.n1];
		Position &n2_pos = node_pos[e.n2];
		Position &n3_pos = node_pos[e.n3];
		Position &n4_pos = node_pos[e.n4];
		elem_vol[e_id] = cal_tetrahedron_vol(n1_pos, n2_pos, n3_pos, n4_pos);
		// shape functions
		pit.init_tetrahedron(n1_pos, n2_pos, n3_pos, n4_pos, elem_vol[e_id]);
		DShapeFuncABC &edNabc = elem_dN_abc[e_id];
		edNabc.dN1_dx = pit.dN1_dx();
		edNabc.dN1_dy = pit.dN1_dy();
		edNabc.dN1_dz = pit.dN1_dz();
		edNabc.dN2_dx = pit.dN2_dx();
		edNabc.dN2_dy = pit.dN2_dy();
		edNabc.dN2_dz = pit.dN2_dz();
		edNabc.dN3_dx = pit.dN3_dx();
		edNabc.dN3_dy = pit.dN3_dy();
		edNabc.dN3_dz = pit.dN3_dz();
		edNabc.dN4_dx = pit.dN4_dx();
		edNabc.dN4_dy = pit.dN4_dy();
		edNabc.dN4_dz = pit.dN4_dz();
		DShapeFuncD& edNd = elem_dN_d[e_id];
		edNd.d1 = pit.get_coef1();
		edNd.d2 = pit.get_coef2();
		edNd.d3 = pit.get_coef3();
		edNd.d4 = pit.get_coef4();
		// node-element relation
		elem_id_array_tmp[e_id4] = e_id;
		elem_id_array_tmp[e_id4 + 1] = e_id;
		elem_id_array_tmp[e_id4 + 2] = e_id;
		elem_id_array_tmp[e_id4 + 3] = e_id;
		node_elem_id_array_tmp[e_id4] = e_id4;
		node_elem_id_array_tmp[e_id4 + 1] = e_id4 + 1;
		node_elem_id_array_tmp[e_id4 + 2] = e_id4 + 2;
		node_elem_id_array_tmp[e_id4 + 3] = e_id4 + 3;
		e_id4 += 4;
		++node_bin_tmp[e.n1];
		++node_bin_tmp[e.n2];
		++node_bin_tmp[e.n3];
		++node_bin_tmp[e.n4];
	}

	node_elem_list[0] = node_bin_tmp[0];
	for (size_t n_id = 1; n_id < node_num; ++n_id)
	{
		node_bin_tmp[n_id] += node_bin_tmp[n_id - 1];
		node_elem_list[n_id] = node_bin_tmp[n_id];
	}
	
	for (size_t e_id = elem_num, e_id4 = elem_num4 - 1; e_id--; e_id4 -= 4)
	{
		const TetrahedronMesh::Element& e = elems[e_id];
		--node_bin_tmp[e.n4];
		elem_id_array[node_bin_tmp[e.n4]] = elem_id_array_tmp[e_id4];
		node_elem_id_array[node_bin_tmp[e.n4]] = node_elem_id_array_tmp[e_id4];
		--node_bin_tmp[e.n3];
		elem_id_array[node_bin_tmp[e.n3]] = elem_id_array_tmp[e_id4 - 1];
		node_elem_id_array[node_bin_tmp[e.n3]] = node_elem_id_array_tmp[e_id4 - 1];
		--node_bin_tmp[e.n2];
		elem_id_array[node_bin_tmp[e.n2]] = elem_id_array_tmp[e_id4 - 2];
		node_elem_id_array[node_bin_tmp[e.n2]] = node_elem_id_array_tmp[e_id4 - 2];
		--node_bin_tmp[e.n1];
		elem_id_array[node_bin_tmp[e.n1]] = elem_id_array_tmp[e_id4 - 3];
		node_elem_id_array[node_bin_tmp[e.n1]] = node_elem_id_array_tmp[e_id4 - 3];
	}

	delete[] elem_id_array_tmp;
}

void Model_T3D_ME_mt::clear_search_grid()
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
	grid_z_num = 0;
}

int Model_T3D_ME_mt::init_search_grid(
	TetrahedronMesh &mesh,
	double _hx,
	double _hy,
	double _hz
	)
{
	clear_search_grid();
	
	SearchingGrid3D<TetrahedronMesh> search_grid;
	search_grid.init(mesh, _hx, _hy, _hz);

	grid_xl = search_grid.get_xl();
	grid_yl = search_grid.get_yl();
	grid_zl = search_grid.get_zl();
	grid_xu = search_grid.get_xu();
	grid_yu = search_grid.get_yu();
	grid_zu = search_grid.get_zu();
	grid_hx = search_grid.get_hx();
	grid_hy = search_grid.get_hy();
	grid_hz = search_grid.get_hz();
	grid_x_num = search_grid.get_x_num();
	grid_y_num = search_grid.get_y_num();
	grid_z_num = search_grid.get_z_num();
	grid_xy_num = grid_x_num * grid_y_num;

	size_t num = grid_x_num * grid_y_num * grid_z_num;
	typedef SearchingGrid3D<TetrahedronMesh>::Grid Grid;
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


void Model_T3D_ME_mt::clear_pcls()
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

void Model_T3D_ME_mt::alloc_pcls(size_t num)
{
	clear_pcls();

	if (num == 0) return;

	size_t mem_len;
	char* cur_mem;

	ori_pcl_num = num;
	pcl_num = ori_pcl_num;
	mem_len = (sizeof(double) + sizeof(BodyForce)
			+ sizeof(Traction) + sizeof(Position) + sizeof(double)
			+ (sizeof(size_t) + sizeof(double)
			 + sizeof(Displacement) + sizeof(Velocity)
			 + sizeof(Stress) + sizeof(Strain) * 3) * 2
			) * num;
	pcl_mem_raw = new char[mem_len];

	cur_mem = pcl_mem_raw;
	pcl_m = (double *)cur_mem;
	cur_mem += sizeof(double) * num;
	pcl_bf = (BodyForce *)(cur_mem);
	cur_mem += sizeof(BodyForce) * num;
	pcl_t = (Traction *)cur_mem;
	cur_mem += sizeof(Traction) * num;
	pcl_pos = (Position *)cur_mem;
	cur_mem += sizeof(Position) * num;
	pcl_vol = (double *)cur_mem;
	cur_mem += sizeof(double) * num;

	SortedPclVarArray &psva0 = pcl_sorted_var_array[0];
	psva0.pcl_index = (size_t *)cur_mem;
	cur_mem += sizeof(size_t) * num;
	psva0.pcl_density = (double *)cur_mem;
	cur_mem += sizeof(double) * num;
	psva0.pcl_disp = (PclDisp *)cur_mem;
	cur_mem += sizeof(PclDisp) * num;
	psva0.pcl_v = (PclV *)cur_mem;
	cur_mem += sizeof(PclV) * num;
	psva0.pcl_stress = (PclStress*)cur_mem;
	cur_mem += sizeof(PclStress) * num;

	PclSortedVarArray& psva1 = pcl_sorted_var_array[1];
	psva1.pcl_index = (size_t *)cur_mem;
	cur_mem += sizeof(size_t) * num;
	psva1.pcl_density = (double *)cur_mem;
	cur_mem += sizeof(double) * num;
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

void Model_T3D_ME_mt::alloc_pcls(
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
	mem_len = (sizeof(double) + sizeof(PclBodyForce)
			 + sizeof(PclTraction) + sizeof(PclPos)) * ori_num
			 + (sizeof(size_t) + sizeof(double)
			  + sizeof(PclDisp) + sizeof(PclV)
			  + sizeof(PclShapeFunc) + sizeof(PclStress)) * 2 * num;
	pcl_mem_raw = new char[mem_len];

	cur_mem = pcl_mem_raw;
	pcl_m = (double*)cur_mem;
	cur_mem += sizeof(double) * ori_num;
	pcl_bf = (PclBodyForce *)(cur_mem);
	cur_mem += sizeof(PclBodyForce) * ori_num;
	pcl_t = (PclTraction *)cur_mem;
	cur_mem += sizeof(PclTraction) * ori_num;
	pcl_pos = (PclPos *)cur_mem;
	cur_mem += sizeof(PclPos) * ori_num;

	PclSortedVarArray& psva0 = pcl_sorted_var_array[0];
	psva0.pcl_index = (size_t*)cur_mem;
	cur_mem += sizeof(size_t) * num;
	psva0.pcl_density = (double*)cur_mem;
	cur_mem += sizeof(double) * num;
	psva0.pcl_disp = (PclDisp*)cur_mem;
	cur_mem += sizeof(PclDisp) * num;
	psva0.pcl_v = (PclV*)cur_mem;
	cur_mem += sizeof(PclV) * num;
	psva0.pcl_N = (PclShapeFunc*)cur_mem;
	cur_mem += sizeof(PclShapeFunc) * num;
	psva0.pcl_stress = (PclStress*)cur_mem;
	cur_mem += sizeof(PclStress) * num;

	PclSortedVarArray& psva1 = pcl_sorted_var_array[1];
	psva1.pcl_index = (size_t*)cur_mem;
	cur_mem += sizeof(size_t) * num;
	psva1.pcl_density = (double*)cur_mem;
	cur_mem += sizeof(double) * num;
	psva1.pcl_disp = (PclDisp*)cur_mem;
	cur_mem += sizeof(PclDisp) * num;
	psva1.pcl_v = (PclV*)cur_mem;
	cur_mem += sizeof(PclV) * num;
	psva1.pcl_N = (PclShapeFunc*)cur_mem;
	cur_mem += sizeof(PclShapeFunc) * num;
	psva1.pcl_stress = (PclStress*)cur_mem;
	cur_mem += sizeof(PclStress) * num;

	pcl_mat_model = new MatModel::MaterialModel * [ori_num];

	// contact
}

int Model_T3D_ME_mt::init_pcls(size_t num, double m, double density)
{
	alloc_pcls(num);

	PclSortedVarArray& psva0 = pcl_sorted_var_array[0];
	size_t p_id;
	for (p_id = 0; p_id < num; ++p_id)
	{
		pcl_m[p_id] = m;
		PclBodyForce &p_bf = pcl_bf[p_id];
		p_bf.bfx = 0.0;
		p_bf.bfy = 0.0;
		p_bf.bfz = 0.0;
		PclTraction& p_t = pcl_t[p_id];
		p_t.tx = 0.0;
		p_t.ty = 0.0;
		p_t.tz = 0.0;
		psva0.pcl_index[p_id] = p_id;
		psva0.pcl_density[p_id] = density;
		PclV& p_v = psva0.pcl_v[p_id];
		p_v.vx = 0.0;
		p_v.vy = 0.0;
		p_v.vz = 0.0;
		PclStress& p_s = psva0.pcl_stress[p_id];
		p_s.s11 = 0.0;
		p_s.s22 = 0.0;
		p_s.s33 = 0.0;
		p_s.s12 = 0.0;
		p_s.s23 = 0.0;
		p_s.s31 = 0.0;
		pcl_mat_model[p_id] = nullptr;
	}
	
	if (bfxs && bfx_num)
	{
		for (size_t bf_id = 0; bf_id < bfx_num; ++bf_id)
		{
			BodyForceAtPcl& bf = bfxs[bf_id];
			p_id = bf.pcl_id;
			pcl_bf[p_id].bfx += pcl_m[p_id] * bf.bf;
		}
		clear_bfxs();
	}

	if (bfys && bfy_num)
	{
		for (size_t bf_id = 0; bf_id < bfy_num; ++bf_id)
		{
			BodyForceAtPcl& bf = bfys[bf_id];
			p_id = bf.pcl_id;
			pcl_bf[p_id].bfy += pcl_m[p_id] * bf.bf;
		}
		clear_bfys();
	}

	if (bfzs && bfz_num)
	{
		for (size_t bf_id = 0; bf_id < bfz_num; ++bf_id)
		{
			BodyForceAtPcl& bf = bfzs[bf_id];
			p_id = bf.pcl_id;
			pcl_bf[p_id].bfz += pcl_m[p_id] * bf.bf;
		}
		clear_bfzs();
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

	if (tzs && tz_num)
	{
		for (size_t t_id = 0; t_id < tz_num; ++t_id)
		{
			TractionBCAtPcl& t = tzs[t_id];
			p_id = t.pcl_id;
			pcl_t[p_id].tz += t.t;
		}
		clear_tzs();
	}

	return 0;
}

int Model_T3D_ME_mt::init_pcls(
	ParticleGenerator3D<TetrahedronMesh>& pg,
	double density
	)
{
	int res = init_pcls(pg.get_num(), density, density);
	if (res)
		return res;

	typedef ParticleGenerator3D<TetrahedronMesh>::Particle PgPcl;
	PgPcl* pg_pcl = pg.first();
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		PclPos &p_p = pcl_pos[p_id];
		p_p.x = pg_pcl->x;
		p_p.y = pg_pcl->y;
		p_p.z = pg_pcl->z;
		pcl_m[p_id] *= pg_pcl->vol;
		pg_pcl = pg.next(pg_pcl);
	}

	return 0;
}

void Model_T3D_ME_mt::init_bfxs(
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
			bf.bfx += bfs[p_id] * pcl_m[p_id];
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

void Model_T3D_ME_mt::init_bfys(
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
			bf.bfy += bfs[p_id] * pcl_m[p_id];
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

void Model_T3D_ME_mt::init_bfzs(
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
			bf.bfz += bfs[p_id] * pcl_m[p_id];
		}
	}
	else
	{
		init_bfzs(bf_num);
		for (size_t bf_id = 0; bf_id < bfz_num; ++bf_id)
		{
			BodyForceAtPcl &bf = bfzs[bf_id];
			bf.pcl_id = bf_pcls[bf_id];
			bf.bf = bfs[bf_id];
		}
	}
}

void Model_T3D_ME_mt::init_txs(
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
			t.tx += ts[t_id];
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

void Model_T3D_ME_mt::init_tys(
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

void Model_T3D_ME_mt::init_tzs(
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
			PclTraction& t = pcl_t[p_id];
			t.tz += ts[t_id];
		}
	}
	else
	{
		init_tzs(t_num);
		for (size_t t_id = 0; t_id < tz_num; ++t_id)
		{
			TractionBCAtPcl& t = tzs[t_id];
			t.pcl_id = t_pcls[t_id];
			t.t = ts[t_id];
		}
	}
}

void Model_T3D_ME_mt::init_fixed_vx_bc(
	size_t bc_num,
	const size_t* bcs
	)
{
	for (size_t bc_id = 0; bc_id < bc_num; ++bc_id)
	{
		if (bcs[bc_id] < node_num)
			node_has_vbc[bcs[bc_id]].has_vx_bc = true;
	}
}

void Model_T3D_ME_mt::init_fixed_vy_bc(
	size_t bc_num,
	const size_t* bcs
	)
{
	for (size_t bc_id = 0; bc_id < bc_num; ++bc_id)
	{
		if (bcs[bc_id] < node_num)
			node_has_vbc[bcs[bc_id]].has_vy_bc = true;
	}
}

void Model_T3D_ME_mt::init_fixed_vz_bc(
	size_t bc_num,
	const size_t* bcs
)
{
	for (size_t bc_id = 0; bc_id < bc_num; ++bc_id)
	{
		if (bcs[bc_id] < node_num)
			node_has_vbc[bcs[bc_id]].has_vz_bc = true;
	}
}
