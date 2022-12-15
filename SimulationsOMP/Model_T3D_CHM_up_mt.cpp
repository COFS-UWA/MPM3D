#include "SimulationsOMP_pcp.h"

#include "SimulationsOMPUtils.h"
#include "TetrahedronUtils.h"
#include "SearchingGrid3D.hpp"
#include "Model_T3D_CHM_up_mt.h"

//static std::fstream res_file_md_t3d_chm_up_mt;

Model_T3D_CHM_up_mt::Model_T3D_CHM_up_mt() :
	ori_pcl_num(0), pcl_num(0),
	pcl_mem_raw(nullptr),
	node_num(0), elem_num(0),
	mesh_mem_raw(nullptr),
	bfx_s_num(0), bfx_ss(nullptr),
	bfy_s_num(0), bfy_ss(nullptr),
	bfz_s_num(0), bfz_ss(nullptr),
	bfx_f_num(0), bfx_fs(nullptr),
	bfy_f_num(0), bfy_fs(nullptr),
	bfz_f_num(0), bfz_fs(nullptr),
	tx_num(0), txs(nullptr),
	ty_num(0), tys(nullptr),
	tz_num(0), tzs(nullptr),
	grid_x_num(0), grid_y_num(0), grid_z_num(0),
	grid_elem_list(nullptr),
	grid_elem_list_id_array(nullptr),
	rigid_t3d_mesh_is_valid(false),
	contact_mem(nullptr),
	pcm(&smooth_contact),
	m_cav(0.0), f_cav_end(0.05), u_cav(0.0),
	u_cav_x(0.0), u_cav_y(0.0), u_cav_z(0.0)
{
	//res_file_md_t3d_me_mt.open("t3d_mt_model.txt", std::ios::binary | std::ios::out);
}

Model_T3D_CHM_up_mt::~Model_T3D_CHM_up_mt()
{
	clear_mesh();
	clear_search_grid();
	clear_pcls();
	clear_bfx_ss();
	clear_bfy_ss();
	clear_bfz_ss();
	clear_bfx_fs();
	clear_bfy_fs();
	clear_bfz_fs();
	clear_txs();
	clear_tys();
	clear_tzs();
	clear_contact_mem();
}

void Model_T3D_CHM_up_mt::set_cavitation(
	double _m_cav,
	double _u_cav,
	double _f_cav_end,
	double _u_cav_x,
	double _u_cav_y,
	double _u_cav_z) noexcept
{
	m_cav = _m_cav;
	u_cav = _u_cav;
	f_cav_end = _f_cav_end;
	u_cav_x = _u_cav_x;
	u_cav_y = _u_cav_y;
	u_cav_z = _u_cav_z;
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		const ElemNodeIndex& en_ids = elem_node_id[e_id];
		const Position& n1_pos = node_pos[en_ids.n1];
		const Position& n2_pos = node_pos[en_ids.n2];
		const Position& n3_pos = node_pos[en_ids.n3];
		const Position& n4_pos = node_pos[en_ids.n4];
		const double elem_cen_x = (n1_pos.x + n2_pos.x + n3_pos.x + n4_pos.x) * 0.25;
		const double elem_cen_y = (n1_pos.y + n2_pos.y + n3_pos.y + n4_pos.y) * 0.25;
		const double elem_cen_z = (n1_pos.z + n2_pos.z + n3_pos.z + n4_pos.z) * 0.25;
		const double e_u_cav = u_cav + elem_cen_x * u_cav_x + elem_cen_y * u_cav_y + elem_cen_z * u_cav_z;
		elem_u_cav[e_id] = e_u_cav < 0.0 ? e_u_cav : 0.0;
	}
}

Cube Model_T3D_CHM_up_mt::get_mesh_bbox()
{
	if (!node_num)
		return Cube(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
	
	const Position& np0 = node_pos[0];
	Cube res(np0.x, np0.x, np0.y, np0.y, np0.z, np0.z);
	for (size_t n_id = 1; n_id < node_num; ++n_id)
	{
		const Position& np = node_pos[n_id];
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

void Model_T3D_CHM_up_mt::clear_mesh()
{
	if (mesh_mem_raw)
	{
		delete[] mesh_mem_raw;
		mesh_mem_raw = nullptr;
	}
	elem_num = 0;
	node_num = 0;
}

void Model_T3D_CHM_up_mt::alloc_mesh(size_t n_num, size_t e_num)
{
	clear_mesh();

	elem_num = e_num;
	node_num = n_num;

	size_t mem_len = (sizeof(ElemNodeIndex)
		+ sizeof(AdjElemIndex) + sizeof(DShapeFuncABC)
		+ sizeof(DShapeFuncD) + sizeof(StrainInc)
		+ sizeof(ElemNodeVM) * 4 + sizeof(Force) * 4
		+ sizeof(size_t) + sizeof(double) * 17
		+ sizeof(uint16_t) * 4) * e_num
		+ (sizeof(Position) + sizeof(Acceleration) + sizeof(Velocity)
		+ sizeof(NodeHasVBC) + sizeof(NodeVBCVec)
		+ sizeof(double) * 4 + sizeof(uint16_t)) * n_num;
	mesh_mem_raw = new char[mem_len];

	char* cur_mem = mesh_mem_raw;
	elem_node_id = (ElemNodeIndex*)cur_mem; // elem_num
	cur_mem += sizeof(ElemNodeIndex) * elem_num;
	elem_adj_elems = (AdjElemIndex*)cur_mem; // elem_num
	cur_mem += sizeof(AdjElemIndex) * elem_num;
	elem_vol = (double*)cur_mem; // elem_num
	cur_mem += sizeof(double) * elem_num;
	elem_dN_abc = (DShapeFuncABC *)cur_mem; // elem_num
	cur_mem += sizeof(DShapeFuncABC) * elem_num;
	elem_dN_d = (DShapeFuncD *)cur_mem; // elem_num
	cur_mem += sizeof(DShapeFuncD) * elem_num;
	elem_u_cav = (double*)cur_mem; // elem_num
	cur_mem += sizeof(double) * elem_num;
	elem_has_pcls = (size_t *)cur_mem; // elem_num, elem has pcls if == substep_id
	cur_mem += sizeof(size_t) * elem_num;
	elem_pcl_vol = (double *)cur_mem; // elem_num
	cur_mem += sizeof(double) * elem_num;
	elem_pcl_m = (double *)cur_mem; // elem_num
	cur_mem += sizeof(double) * elem_num;
	elem_pcl_n = (double *)cur_mem; // elem_num
	cur_mem += sizeof(double) * elem_num;
	elem_de = (StrainInc *)cur_mem; // elem_num
	cur_mem += sizeof(StrainInc) * elem_num;
	elem_pcl_pm = (double *)cur_mem;
	cur_mem += sizeof(double) * elem_num;
	elem_p = (double *)cur_mem; // elem_num
	cur_mem += sizeof(double) * elem_num;
	elem_density_f = (double *)cur_mem; // elem_num
	cur_mem += sizeof(double) * elem_num;
	elem_m_de_vol_s = (double *)cur_mem; // elem_num, strain enhancement
	cur_mem += sizeof(double) * elem_num;

	elem_node_vm_s = (ElemNodeVM *)cur_mem; // elem_num * 4
	cur_mem += sizeof(ElemNodeVM) * elem_num * 4;
	elem_node_force = (Force *)cur_mem; // elem_num * 4
	cur_mem += sizeof(Force) * elem_num * 4;
	elem_node_p = (double *)cur_mem; // elem_num * 4
	cur_mem += sizeof(double) * elem_num * 4;
	elem_node_p_force = (double *)cur_mem;
	cur_mem += sizeof(double) * elem_num * 4;
	elem_node_at_surface = (uint16_t*)cur_mem; // elem_num * 4
	cur_mem += sizeof(uint16_t) * elem_num * 4;

	node_pos = (Position *)cur_mem; // node_num
	cur_mem += sizeof(Position) * node_num;
	node_has_vbc = (NodeHasVBC*)cur_mem; // node_num
	cur_mem += sizeof(NodeHasVBC) * node_num;
	node_vbc_vec_s = (NodeVBCVec*)cur_mem; // node_num
	cur_mem += sizeof(NodeVBCVec) * node_num;
	node_am = (double *)cur_mem; // node_num
	cur_mem += sizeof(double) * node_num;
	node_a_s = (Acceleration *)cur_mem; // node_num
	cur_mem += sizeof(Acceleration) * node_num;
	node_v_s = (Velocity *)cur_mem; // node_num
	cur_mem += sizeof(Velocity) * node_num;
	node_p = (double*)cur_mem; // node_num
	cur_mem += sizeof(double) * node_num;
	node_dp = (double *)cur_mem; // node_num
	cur_mem += sizeof(double) * node_num;
	node_at_surface = (uint16_t*)cur_mem; // node_num
	cur_mem += sizeof(uint16_t) * node_num;
	node_de_vol_s = (double *)cur_mem; // node_num
}

void Model_T3D_CHM_up_mt::init_mesh(const TetrahedronMesh &mesh)
{
	if (mesh.get_elem_num() == 0 || mesh.get_node_num() == 0)
		return;
	
	alloc_mesh(mesh.get_node_num(), mesh.get_elem_num());

	// init node coordinates
	const TetrahedronMesh::Node* nodes = mesh.get_nodes();
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		const TetrahedronMesh::Node &n = nodes[n_id];
		Position& np = node_pos[n_id];
		np.x = n.x;
		np.y = n.y;
		np.z = n.z;
		NodeHasVBC& n_vbc = node_has_vbc[n_id];
		n_vbc.has_vx_bc = false;
		n_vbc.has_vy_bc = false;
		n_vbc.has_vz_bc = false;
		n_vbc.is_drained = false;
		NodeVBCVec& n_vbcv = node_vbc_vec_s[n_id];
		n_vbcv.x = 0.0;
		n_vbcv.y = 0.0;
		n_vbcv.z = 0.0;
	}

	// init elem connectivity, area and shape functions
	PointInTetrahedron pit;
	const TetrahedronMesh::Element *elems = mesh.get_elems();
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		const TetrahedronMesh::Element &e = elems[e_id];
		// geometry
		ElemNodeIndex& eni = elem_node_id[e_id];
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
	}
	
	TehMeshAdjElementFinder finder;
	finder.find(elem_node_id, elem_num, node_num, elem_adj_elems);
}

void Model_T3D_CHM_up_mt::clear_search_grid()
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

int Model_T3D_CHM_up_mt::init_search_grid(TetrahedronMesh &mesh)
{
	clear_search_grid();
	
	TetrahedronMesh::SearchGrid &search_grid = mesh.get_bg_grid();

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
	typedef TetrahedronMesh::SearchGrid::Grid Grid;
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
	size_t ge_start_id;
	for (size_t g_id = 0; g_id < num; ++g_id)
	{
		Grid& grid = sg[g_id];
		ge_start_id = cur_id;
		for (auto pe = grid.pelems; pe; pe = pe->next)
		{
			grid_elem_list_id_array[cur_id] = pe->e->id;
			++cur_id;
		}
		// sort grid elem list to speed up memory access
		SimulationsOMP::swap_sort_acc(
			grid_elem_list_id_array + ge_start_id,
			cur_id - ge_start_id);
	}

	return 0;
}


void Model_T3D_CHM_up_mt::clear_pcls()
{
	if (pcl_mem_raw)
	{
		delete[] pcl_mem_raw;
		pcl_mem_raw = nullptr;
	}
	ori_pcl_num = 0;
	pcl_num = 0;
}

void Model_T3D_CHM_up_mt::alloc_pcls(size_t num)
{
	clear_pcls();

	if (num == 0)
		return;

	size_t mem_len;
	char* cur_mem;

	ori_pcl_num = num;
	pcl_num = num;
	mem_len = (sizeof(double) * 6 + sizeof(Force) * 3
			+ sizeof(Position) + sizeof(MatModel::MaterialModel*)
			+ (sizeof(size_t) + sizeof(double) * 3
			 + sizeof(Velocity) + sizeof(Displacement)
			 + sizeof(Stress) + sizeof(Strain) * 3
			 + sizeof(ShapeFunc)) * 2
			) * num;
	
	pcl_mem_raw = new char[mem_len];
	cur_mem = pcl_mem_raw;
	pcl_m_s = (double *)cur_mem;
	cur_mem += sizeof(double) * num;
	pcl_density_s = (double*)cur_mem;
	cur_mem += sizeof(double) * num;
	pcl_vol_s = (double*)cur_mem;
	cur_mem += sizeof(double) * num;
	pcl_bf_s = (Force *)cur_mem;
	cur_mem += sizeof(Force) * num;
	pcl_bf_f = (Force *)cur_mem;
	cur_mem += sizeof(Force) * num;
	pcl_t = (Force *)cur_mem;
	cur_mem += sizeof(Force) * num;
	pcl_pos = (Position *)cur_mem;
	cur_mem += sizeof(Position) * num;
	pcl_mat_model = (MatModel::MaterialModel**)cur_mem;
	cur_mem += sizeof(MatModel::MaterialModel*) * num;

	pcl_vol = (double*)cur_mem;
	cur_mem += sizeof(double) * num;
	
	pcl_u_cav = (double*)cur_mem;
	cur_mem += sizeof(double) * num;
	pcl_is_cavitated = (double*)cur_mem;
	cur_mem += sizeof(double) * num;

	SortedPclVarArrays &spva0 = sorted_pcl_var_arrays[0];
	spva0.pcl_index = (size_t *)cur_mem;
	cur_mem += sizeof(size_t) * num;
	spva0.pcl_n = (double*)cur_mem;
	cur_mem += sizeof(double) * num;
	spva0.pcl_density_f = (double *)cur_mem;
	cur_mem += sizeof(double) * num;
	spva0.pcl_v_s = (Velocity*)cur_mem;
	cur_mem += sizeof(Velocity) * num;
	spva0.pcl_u = (Displacement *)cur_mem;
	cur_mem += sizeof(Displacement) * num;
	spva0.pcl_stress = (Stress*)cur_mem;
	cur_mem += sizeof(Stress) * num;
	spva0.pcl_p = (double*)cur_mem;
	cur_mem += sizeof(double) * num;
	spva0.pcl_strain = (Strain*)cur_mem;
	cur_mem += sizeof(Strain) * num;
	spva0.pcl_estrain = (Strain *)cur_mem;
	cur_mem += sizeof(Strain) * num;
	spva0.pcl_pstrain = (Strain*)cur_mem;
	cur_mem += sizeof(Strain) * num;
	spva0.pcl_N = (ShapeFunc*)cur_mem;
	cur_mem += sizeof(ShapeFunc) * num;

	SortedPclVarArrays &spva1 = sorted_pcl_var_arrays[1];
	spva1.pcl_index = (size_t *)cur_mem;
	cur_mem += sizeof(size_t) * num;
	spva1.pcl_n = (double*)cur_mem;
	cur_mem += sizeof(double) * num;
	spva1.pcl_density_f = (double *)cur_mem;
	cur_mem += sizeof(double) * num;
	spva1.pcl_v_s = (Velocity*)cur_mem;
	cur_mem += sizeof(Velocity) * num;
	spva1.pcl_u = (Displacement*)cur_mem;
	cur_mem += sizeof(Displacement) * num;
	spva1.pcl_stress = (Stress*)cur_mem;
	cur_mem += sizeof(Stress) * num;
	spva1.pcl_p = (double*)cur_mem;
	cur_mem += sizeof(double) * num;
	spva1.pcl_strain = (Strain*)cur_mem;
	cur_mem += sizeof(Strain) * num;
	spva1.pcl_estrain = (Strain*)cur_mem;
	cur_mem += sizeof(Strain) * num;
	spva1.pcl_pstrain = (Strain*)cur_mem;
	cur_mem += sizeof(Strain) * num;
	spva1.pcl_N = (ShapeFunc*)cur_mem;
	cur_mem += sizeof(ShapeFunc) * num;

	alloc_contact_mem(num);
}

void Model_T3D_CHM_up_mt::alloc_pcls(
	size_t num,
	size_t ori_num)
{
	clear_pcls();

	if (num == 0 || ori_num == 0)
		return;

	alloc_pcls(ori_num);
	pcl_num = num;
}

int Model_T3D_CHM_up_mt::init_pcls(
	size_t num,
	double n,
	double m_s,
	double density_s,
	double density_f,
	double _Kf,
	double _k,
	double _vis)
{
	if (num == 0)
		return -1;
	alloc_pcls(num);

	Kf = _Kf;
	k = _k;
	dyn_viscosity = _vis;

	SortedPclVarArrays &spva0 = sorted_pcl_var_arrays[0];
	const double vol_s = m_s / density_s;
	const double vol = vol_s / (1.0 - n);
	size_t p_id;
	for (p_id = 0; p_id < num; ++p_id)
	{
		pcl_m_s[p_id] = m_s;
		pcl_density_s[p_id] = density_s;
		pcl_vol_s[p_id] = vol_s;
		pcl_vol[p_id] = vol;
		Force& bf_s = pcl_bf_s[p_id];
		bf_s.fx = 0.0;
		bf_s.fy = 0.0;
		bf_s.fz = 0.0;
		Force& bf_f = pcl_bf_f[p_id];
		bf_f.fx = 0.0;
		bf_f.fy = 0.0;
		bf_f.fz = 0.0;
		Force& p_t = pcl_t[p_id];
		p_t.fx = 0.0;
		p_t.fy = 0.0;
		p_t.fz = 0.0;
		pcl_mat_model[p_id] = nullptr;
		spva0.pcl_index[p_id] = p_id;
		spva0.pcl_n[p_id] = n;
		spva0.pcl_density_f[p_id] = density_f;
		Velocity &p_v = spva0.pcl_v_s[p_id];
		p_v.vx = 0.0;
		p_v.vy = 0.0;
		p_v.vz = 0.0;
		Displacement& p_d = spva0.pcl_u[p_id];
		p_d.ux = 0.0;
		p_d.uy = 0.0;
		p_d.uz = 0.0;
		Stress& p_s = spva0.pcl_stress[p_id];
		p_s.s11 = 0.0;
		p_s.s22 = 0.0;
		p_s.s33 = 0.0;
		p_s.s12 = 0.0;
		p_s.s23 = 0.0;
		p_s.s31 = 0.0;
		spva0.pcl_p[p_id] = 0.0;
		Strain &p_e = spva0.pcl_strain[p_id];
		p_e.e11 = 0.0;
		p_e.e22 = 0.0;
		p_e.e33 = 0.0;
		p_e.e12 = 0.0;
		p_e.e23 = 0.0;
		p_e.e31 = 0.0;
		Strain& p_ee = spva0.pcl_estrain[p_id];
		p_ee.e11 = 0.0;
		p_ee.e22 = 0.0;
		p_ee.e33 = 0.0;
		p_ee.e12 = 0.0;
		p_ee.e23 = 0.0;
		p_ee.e31 = 0.0;
		Strain& p_pe = spva0.pcl_pstrain[p_id];
		p_pe.e11 = 0.0;
		p_pe.e22 = 0.0;
		p_pe.e33 = 0.0;
		p_pe.e12 = 0.0;
		p_pe.e23 = 0.0;
		p_pe.e31 = 0.0;
		// cavitation
		pcl_u_cav[p_id] = 0.0;
		pcl_is_cavitated[p_id] = 1.0;
	}
	
	if (bfx_ss && bfx_s_num)
	{
		for (size_t bf_id = 0; bf_id < bfx_s_num; ++bf_id)
		{
			BodyForceAtPcl& bf = bfx_ss[bf_id];
			p_id = bf.pcl_id;
			pcl_bf_s[p_id].fx += bf.bf;
		}
		clear_bfx_ss();
	}

	if (bfy_ss && bfy_s_num)
	{
		for (size_t bf_id = 0; bf_id < bfy_s_num; ++bf_id)
		{
			BodyForceAtPcl& bf = bfy_ss[bf_id];
			p_id = bf.pcl_id;
			pcl_bf_s[p_id].fy += bf.bf;
		}
		clear_bfy_ss();
	}

	if (bfz_ss && bfz_s_num)
	{
		for (size_t bf_id = 0; bf_id < bfz_s_num; ++bf_id)
		{
			BodyForceAtPcl& bf = bfz_ss[bf_id];
			p_id = bf.pcl_id;
			pcl_bf_s[p_id].fz += bf.bf;
		}
		clear_bfz_ss();
	}

	if (bfx_fs && bfx_f_num)
	{
		for (size_t bf_id = 0; bf_id < bfx_f_num; ++bf_id)
		{
			BodyForceAtPcl& bf = bfx_fs[bf_id];
			p_id = bf.pcl_id;
			pcl_bf_f[p_id].fx += bf.bf;
		}
		clear_bfx_fs();
	}

	if (bfy_fs && bfy_f_num)
	{
		for (size_t bf_id = 0; bf_id < bfy_f_num; ++bf_id)
		{
			BodyForceAtPcl& bf = bfy_fs[bf_id];
			p_id = bf.pcl_id;
			pcl_bf_f[p_id].fy += bf.bf;
		}
		clear_bfy_fs();
	}

	if (bfz_fs && bfz_f_num)
	{
		for (size_t bf_id = 0; bf_id < bfz_f_num; ++bf_id)
		{
			BodyForceAtPcl& bf = bfz_fs[bf_id];
			p_id = bf.pcl_id;
			pcl_bf_f[p_id].fz += bf.bf;
		}
		clear_bfz_fs();
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

	if (tzs && tz_num)
	{
		for (size_t t_id = 0; t_id < tz_num; ++t_id)
		{
			TractionBCAtPcl& t = tzs[t_id];
			p_id = t.pcl_id;
			pcl_t[p_id].fz += t.t;
		}
		clear_tzs();
	}

	return 0;
}

int Model_T3D_CHM_up_mt::init_pcls(
	ParticleGenerator3D<TetrahedronMesh>& pg,
	double n, double density_s, double density_f,
	double _Kf, double _k, double _vis)
{
	int res = init_pcls(pg.get_num(), n, density_s, density_s, density_f, _Kf, _k, _vis);
	if (res)
		return res;

	SortedPclVarArrays& spva = sorted_pcl_var_arrays[0];
	typedef ParticleGenerator3D<TetrahedronMesh>::Particle PgPcl;
	PgPcl* pg_pcl = pg.first();
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		Position& p_p = pcl_pos[p_id];
		p_p.x = pg_pcl->x;
		p_p.y = pg_pcl->y;
		p_p.z = pg_pcl->z;
		pcl_vol[p_id] = pg_pcl->vol;
		pcl_vol_s[p_id] = pg_pcl->vol * (1.0 - spva.pcl_n[p_id]);
		pcl_m_s[p_id] *= pcl_vol_s[p_id];
		pcl_bf_s[p_id].fx *= pcl_vol_s[p_id];
		pcl_bf_s[p_id].fy *= pcl_vol_s[p_id];
		pcl_bf_s[p_id].fz *= pcl_vol_s[p_id];
		pg_pcl = pg.next(pg_pcl);
	}

	return 0;
}

void Model_T3D_CHM_up_mt::init_bfx_ss(
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
			Force &bf = pcl_bf_s[p_id];
			bf.fx += bfs[p_id] * pcl_m_s[p_id];
		}
	}
	else
	{
		init_bfx_ss(bf_num);
		for (size_t bf_id = 0; bf_id < bfx_s_num; ++bf_id)
		{
			BodyForceAtPcl& bf = bfx_ss[bf_id];
			bf.pcl_id = bf_pcls[bf_id];
			bf.bf = bfs[bf_id];
		}
	}
}

void Model_T3D_CHM_up_mt::init_bfy_ss(
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
			Force& bf = pcl_bf_s[p_id];
			bf.fy += bfs[p_id] * pcl_m_s[p_id];
		}
	}
	else
	{
		init_bfy_ss(bf_num);
		for (size_t bf_id = 0; bf_id < bfy_s_num; ++bf_id)
		{
			BodyForceAtPcl &bf = bfy_ss[bf_id];
			bf.pcl_id = bf_pcls[bf_id];
			bf.bf = bfs[bf_id];
		}
	}
}

void Model_T3D_CHM_up_mt::init_bfz_ss(
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
			Force& bf = pcl_bf_s[p_id];
			bf.fz += bfs[p_id] * pcl_m_s[p_id];
		}
	}
	else
	{
		init_bfz_ss(bf_num);
		for (size_t bf_id = 0; bf_id < bfz_s_num; ++bf_id)
		{
			BodyForceAtPcl &bf = bfz_ss[bf_id];
			bf.pcl_id = bf_pcls[bf_id];
			bf.bf = bfs[bf_id];
		}
	}
}

void Model_T3D_CHM_up_mt::init_bfx_fs(
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
			Force& bf = pcl_bf_f[p_id];
			bf.fx += bfs[p_id];
		}
	}
	else
	{
		init_bfx_fs(bf_num);
		for (size_t bf_id = 0; bf_id < bfx_f_num; ++bf_id)
		{
			BodyForceAtPcl& bf = bfx_fs[bf_id];
			bf.pcl_id = bf_pcls[bf_id];
			bf.bf = bfs[bf_id];
		}
	}
}

void Model_T3D_CHM_up_mt::init_bfy_fs(
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
			Force& bf = pcl_bf_f[p_id];
			bf.fy += bfs[p_id];
		}
	}
	else
	{
		init_bfy_fs(bf_num);
		for (size_t bf_id = 0; bf_id < bfy_f_num; ++bf_id)
		{
			BodyForceAtPcl& bf = bfy_fs[bf_id];
			bf.pcl_id = bf_pcls[bf_id];
			bf.bf = bfs[bf_id];
		}
	}
}

void Model_T3D_CHM_up_mt::init_bfz_fs(
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
			Force& bf = pcl_bf_f[p_id];
			bf.fz += bfs[p_id];
		}
	}
	else
	{
		init_bfz_fs(bf_num);
		for (size_t bf_id = 0; bf_id < bfz_f_num; ++bf_id)
		{
			BodyForceAtPcl& bf = bfz_fs[bf_id];
			bf.pcl_id = bf_pcls[bf_id];
			bf.bf = bfs[bf_id];
		}
	}
}

void Model_T3D_CHM_up_mt::init_txs(
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

void Model_T3D_CHM_up_mt::init_tys(
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

void Model_T3D_CHM_up_mt::init_tzs(
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
			t.fz += ts[t_id];
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

void Model_T3D_CHM_up_mt::init_drained_bc(
	size_t bc_num,
	const size_t* bcs)
{
	for (size_t bc_id = 0; bc_id < bc_num; ++bc_id)
	{
		if (bcs[bc_id] < node_num)
			node_has_vbc[bcs[bc_id]].is_drained = true;
	}
}

void Model_T3D_CHM_up_mt::init_fixed_vx_bc(
	size_t bc_num,
	const size_t* bcs)
{
	for (size_t bc_id = 0; bc_id < bc_num; ++bc_id)
	{
		if (bcs[bc_id] < node_num)
			node_has_vbc[bcs[bc_id]].has_vx_bc = true;
	}
}

void Model_T3D_CHM_up_mt::init_fixed_vy_bc(
	size_t bc_num,
	const size_t* bcs)
{
	for (size_t bc_id = 0; bc_id < bc_num; ++bc_id)
	{
		if (bcs[bc_id] < node_num)
			node_has_vbc[bcs[bc_id]].has_vy_bc = true;
	}
}

void Model_T3D_CHM_up_mt::init_fixed_vz_bc(
	size_t bc_num,
	const size_t* bcs)
{
	for (size_t bc_id = 0; bc_id < bc_num; ++bc_id)
	{
		if (bcs[bc_id] < node_num)
			node_has_vbc[bcs[bc_id]].has_vz_bc = true;
	}
}

void Model_T3D_CHM_up_mt::set_vbc_vec(
	size_t n_id,
	double vecx,
	double vecy,
	double vecz)
{
	const double len = sqrt(vecx * vecx + vecy * vecy + vecz * vecz);
	if (len != 0.0)
	{
		NodeVBCVec& nbcv = node_vbc_vec_s[n_id];
		nbcv.x = vecx / len;
		nbcv.y = vecy / len;
		nbcv.z = vecz / len;
	}
}

void Model_T3D_CHM_up_mt::clear_contact_mem()
{
	if (contact_mem)
	{
		delete[] contact_mem;
		contact_mem = nullptr;
	}
}

void Model_T3D_CHM_up_mt::alloc_contact_mem(size_t num)
{
	clear_contact_mem();

	contact_mem = new char[(sizeof(size_t)
		+ sizeof(ContactModel3D::Position)
		+ sizeof(ContactModel3D::Force)
		+ sizeof(double)) * num];
	
	char* cur_mem = contact_mem;
	contact_substep_ids = (size_t *)cur_mem;
	cur_mem += sizeof(size_t) * num;
	prev_contact_poses = (ContactModel3D::Position *)cur_mem;
	cur_mem += sizeof(ContactModel3D::Position) * num;
	prev_contact_tan_forces = (ContactModel3D::Force *)cur_mem;
	cur_mem += sizeof(ContactModel3D::Force) * num;
	prev_contact_dists = (double *)cur_mem;

	// smooth contact
	smooth_contact.contact_substep_ids = contact_substep_ids;
	smooth_contact.prev_contact_dists = prev_contact_dists;
	// smooth contact quad
	smh_contact_quad.contact_substep_ids = contact_substep_ids;
	smh_contact_quad.prev_contact_dists = prev_contact_dists;
	// rough contact
	rough_contact.contact_substep_ids = contact_substep_ids;
	rough_contact.prev_contact_poses = prev_contact_poses;
	rough_contact.prev_contact_tan_forces = prev_contact_tan_forces;
	rough_contact.prev_contact_dists = prev_contact_dists;
	// frictional contact
	fric_contact.contact_substep_ids = contact_substep_ids;
	fric_contact.prev_contact_poses = prev_contact_poses;
	fric_contact.prev_contact_tan_forces = prev_contact_tan_forces;
	fric_contact.prev_contact_dists = prev_contact_dists;
	// sticky contact
	sticky_contact.contact_substep_ids = contact_substep_ids;
	sticky_contact.prev_contact_poses = prev_contact_poses;
	sticky_contact.prev_contact_tan_forces = prev_contact_tan_forces;
	sticky_contact.prev_contact_dists = prev_contact_dists;
}
