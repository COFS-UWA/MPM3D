#include "SimulationsOMP_pcp.h"

#include <iomanip>
#include <fstream>

#include "SimulationsOMPUtils.h"
#include "TetrahedronUtils.h"
#include "SearchingGrid3D.hpp"
#include "Model_T3D_CHM_mt.h"

Model_T3D_CHM_mt::Model_T3D_CHM_mt() :
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
	rigid_cylinder_is_valid(false),
	rigid_t3d_mesh_is_valid(false),
	contact_mem(nullptr),
	pcm_s(&smooth_contact_s),
	pcm_f(&smooth_contact_f)
{

}

Model_T3D_CHM_mt::~Model_T3D_CHM_mt()
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

Cube Model_T3D_CHM_mt::get_mesh_bbox()
{
	if (!node_num)
		return Cube(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

	Position& np0 = node_pos[0];
	Cube res(np0.x, np0.x, np0.y, np0.y, np0.z, np0.z);
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
		if (res.zl > np.z)
			res.zl = np.z;
		if (res.zu < np.z)
			res.zu = np.z;
	}
	return res;
}

void Model_T3D_CHM_mt::clear_mesh()
{
	if (mesh_mem_raw)
	{
		delete[] mesh_mem_raw;
		mesh_mem_raw = nullptr;
	}
	elem_num = 0;
	node_num = 0;
}

void Model_T3D_CHM_mt::alloc_mesh(size_t n_num, size_t e_num)
{
	node_num = n_num;
	elem_num = e_num;

	size_t mem_len = (sizeof(ElemNodeIndex)
		+ sizeof(DShapeFuncABC) + sizeof(DShapeFuncD)
		+ sizeof(double) * 9 + sizeof(StrainInc)
		+ sizeof(ElemNodeVM) * 8 + sizeof(Force) * 9) * e_num
		+ (sizeof(Position) + sizeof(Acceleration) * 2 + sizeof(Velocity) * 2
		+ sizeof(NodeHasVBC) * 2 + sizeof(NodeVBCVec) * 2
		+ sizeof(double) * 4) * n_num;
	mesh_mem_raw = new char[mem_len];

	char* cur_mem = mesh_mem_raw;
	elem_node_id = (ElemNodeIndex*)cur_mem;
	cur_mem += sizeof(ElemNodeIndex) * elem_num;
	elem_N_abc = (DShapeFuncABC*)cur_mem;
	cur_mem += sizeof(DShapeFuncABC) * elem_num;
	elem_N_d = (DShapeFuncD*)cur_mem;
	cur_mem += sizeof(DShapeFuncD) * elem_num;
	elem_vol = (double*)cur_mem;
	cur_mem += sizeof(double) * elem_num;
	elem_density_f = (double*)cur_mem;
	cur_mem += sizeof(double) * elem_num;
	elem_pcl_n = (double*)cur_mem;
	cur_mem += sizeof(double) * elem_num;
	elem_pcl_m_s = (double*)cur_mem;
	cur_mem += sizeof(double) * elem_num;
	elem_pcl_m_f = (double*)cur_mem;
	cur_mem += sizeof(double) * elem_num;
	elem_de = (StrainInc *)cur_mem;
	cur_mem += sizeof(StrainInc) * elem_num;
	elem_p = (double*)cur_mem;
	cur_mem += sizeof(double) * elem_num;
	elem_n2_miu_div_k_vol = (double *)cur_mem;
	cur_mem += sizeof(double) * elem_num;
	elem_seep_force = (Force*)cur_mem;
	cur_mem += sizeof(Force) * elem_num;
	elem_m_de_vol_s = (double*)cur_mem;
	cur_mem += sizeof(double) * elem_num;
	elem_m_de_vol_f = (double*)cur_mem;
	cur_mem += sizeof(double) * elem_num;
	elem_node_vm_s = (ElemNodeVM*)cur_mem;
	cur_mem += sizeof(ElemNodeVM) * elem_num * 4;
	elem_node_vm_f = (ElemNodeVM*)cur_mem;
	cur_mem += sizeof(ElemNodeVM) * elem_num * 4;
	elem_node_force_s = (Force*)cur_mem;
	cur_mem += sizeof(Force) * elem_num * 4;
	elem_node_force_f = (Force*)cur_mem;
	cur_mem += sizeof(Force) * elem_num * 4;

	node_pos = (Position*)cur_mem;
	cur_mem += sizeof(Position) * node_num;
	node_a_s = (Acceleration *)cur_mem;
	cur_mem += sizeof(Acceleration) * node_num;
	node_a_f = (Acceleration*)cur_mem;
	cur_mem += sizeof(Acceleration) * node_num;
	node_v_s = (Velocity *)cur_mem;
	cur_mem += sizeof(Velocity) * node_num;
	node_v_f = (Velocity*)cur_mem;
	cur_mem += sizeof(Velocity) * node_num;
	node_has_vbc_s = (NodeHasVBC *)cur_mem;
	cur_mem += sizeof(NodeHasVBC) * node_num;
	node_has_vbc_f = (NodeHasVBC*)cur_mem;
	cur_mem += sizeof(NodeHasVBC) * node_num;
	node_vbc_vec_s = (NodeVBCVec*)cur_mem;
	cur_mem += sizeof(NodeVBCVec) * node_num;
	node_vbc_vec_f = (NodeVBCVec*)cur_mem;
	cur_mem += sizeof(NodeVBCVec) * node_num;
	node_am_s = (double*)cur_mem;
	cur_mem += sizeof(double) * node_num;
	node_am_f = (double*)cur_mem;
	cur_mem += sizeof(double) * node_num;
	node_de_vol_s = (double*)cur_mem;
	cur_mem += sizeof(double) * node_num;
	node_de_vol_f = (double*)cur_mem;
	cur_mem += sizeof(double) * node_num;
}

void Model_T3D_CHM_mt::init_mesh(const TetrahedronMesh &mesh)
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
		NodeHasVBC& n_vbc_s = node_has_vbc_s[n_id];
		n_vbc_s.has_vx_bc = false;
		n_vbc_s.has_vy_bc = false;
		n_vbc_s.has_vz_bc = false;
		NodeHasVBC& n_vbc_f = node_has_vbc_f[n_id];
		n_vbc_f.has_vx_bc = false;
		n_vbc_f.has_vy_bc = false;
		n_vbc_f.has_vz_bc = false;
		NodeVBCVec& n_vbcv_s = node_vbc_vec_s[n_id];
		n_vbcv_s.x = 0.0;
		n_vbcv_s.y = 0.0;
		n_vbcv_s.z = 0.0;
		NodeVBCVec& n_vbcv_f = node_vbc_vec_f[n_id];
		n_vbcv_f.x = 0.0;
		n_vbcv_f.y = 0.0;
		n_vbcv_f.z = 0.0;
	}

	// init elem connectivity, area and shape functions
	PointInTetrahedron pit;
	const TetrahedronMesh::Element *elems = mesh.get_elems();
	size_t e_id3 = 0;
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
		Position& n4_pos = node_pos[e.n4];
		elem_vol[e_id] = cal_tetrahedron_vol(n1_pos, n2_pos, n3_pos, n4_pos);
		// shape functions
		pit.init_tetrahedron(n1_pos, n2_pos, n3_pos, n4_pos, elem_vol[e_id]);
		DShapeFuncABC& edNabc = elem_N_abc[e_id];
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
		DShapeFuncD& edNd = elem_N_d[e_id];
		edNd.d1 = pit.get_coef1();
		edNd.d2 = pit.get_coef2();
		edNd.d3 = pit.get_coef3();
		edNd.d4 = pit.get_coef4();
	}
}

void Model_T3D_CHM_mt::clear_search_grid()
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

int Model_T3D_CHM_mt::init_search_grid(TetrahedronMesh &mesh)
{
	clear_search_grid();
	
	TetrahedronMesh::SearchGrid& search_grid = mesh.get_bg_grid();

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

	const size_t num = grid_x_num * grid_y_num * grid_z_num;
	typedef TetrahedronMesh::SearchGrid::Grid Grid;
	Grid* sg = search_grid.get_grids();
	grid_elem_list = new size_t[num + 1];
	grid_elem_list[0] = 0;
	size_t cur_id = 0;
	for (size_t g_id = 0; g_id < num; ++g_id)
	{
		Grid& grid = sg[g_id];
		for (auto pe = grid.pelems; pe; pe = pe->next)
			++cur_id;
		grid_elem_list[g_id + 1] = cur_id;
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

void Model_T3D_CHM_mt::clear_pcls()
{
	if (pcl_mem_raw)
	{
		delete[] pcl_mem_raw;
		pcl_mem_raw = nullptr;
	}
	ori_pcl_num = 0;
	pcl_num = 0;
}

void Model_T3D_CHM_mt::alloc_pcls(size_t num)
{
	clear_pcls();

	if (num == 0)
		return;

	ori_pcl_num = num;
	pcl_num = ori_pcl_num;
	size_t mem_len = (sizeof(double)  * 4
		+ sizeof(Force) * 3	+ sizeof(Position)
		+ sizeof(MatModel::MaterialModel *)
		+ (sizeof(size_t) + sizeof(double) * 3
		+ sizeof(Velocity) * 2 + sizeof(Displacement) * 2
		+ sizeof(Stress) + sizeof(Strain) * 3
		+ sizeof(ShapeFunc)) * 2) * num;
	pcl_mem_raw = new char[mem_len];

	char* cur_mem = pcl_mem_raw;
	pcl_m_s = (double *)cur_mem;
	cur_mem += sizeof(double) * num;
	pcl_density_s = (double *)cur_mem;
	cur_mem += sizeof(double) * num;
	pcl_vol_s = (double *)cur_mem;
	cur_mem += sizeof(double) * num;
	pcl_bf_s = (Force *)cur_mem;
	cur_mem += sizeof(Force) * num;
	pcl_bf_f = (Force *)cur_mem;
	cur_mem += sizeof(Force) * num;
	pcl_t = (Force *)cur_mem;
	cur_mem += sizeof(Force) * num;
	pcl_pos = (Position *)cur_mem;
	cur_mem += sizeof(Position) * num;
	pcl_vol = (double *)cur_mem;
	cur_mem += sizeof(double) * num;
	pcl_mat_model = (MatModel::MaterialModel **)cur_mem;
	cur_mem += sizeof(MatModel::MaterialModel *) * num;

	SortedPclVarArrays &spva0 = sorted_pcl_var_arrays[0];
	spva0.pcl_index = (size_t*)cur_mem;
	cur_mem += sizeof(size_t) * num;
	spva0.pcl_n = (double *)cur_mem;
	cur_mem += sizeof(double) * num;
	spva0.pcl_density_f = (double *)cur_mem;
	cur_mem += sizeof(double) * num;
	spva0.pcl_v_s = (Velocity *)cur_mem;
	cur_mem += sizeof(Velocity) * num;
	spva0.pcl_v_f = (Velocity *)cur_mem;
	cur_mem += sizeof(Velocity) * num;
	spva0.pcl_u_s = (Displacement *)cur_mem;
	cur_mem += sizeof(Displacement) * num;
	spva0.pcl_u_f = (Displacement *)cur_mem;
	cur_mem += sizeof(Displacement) * num;
	spva0.pcl_stress = (Stress *)cur_mem;
	cur_mem += sizeof(Stress) * num;
	spva0.pcl_p = (double *)cur_mem;
	cur_mem += sizeof(double) * num;
	spva0.pcl_strain = (Strain *)cur_mem;
	cur_mem += sizeof(Strain) * num;
	spva0.pcl_estrain = (Strain *)cur_mem;
	cur_mem += sizeof(Strain) * num;
	spva0.pcl_pstrain = (Strain *)cur_mem;
	cur_mem += sizeof(Strain) * num;
	spva0.pcl_N = (ShapeFunc *)cur_mem;
	cur_mem += sizeof(ShapeFunc) * num;

	SortedPclVarArrays& spva1 = sorted_pcl_var_arrays[1];
	spva1.pcl_index = (size_t*)cur_mem;
	cur_mem += sizeof(size_t) * num;
	spva1.pcl_n = (double*)cur_mem;
	cur_mem += sizeof(double) * num;
	spva1.pcl_density_f = (double*)cur_mem;
	cur_mem += sizeof(double) * num;
	spva1.pcl_v_s = (Velocity*)cur_mem;
	cur_mem += sizeof(Velocity) * num;
	spva1.pcl_v_f = (Velocity*)cur_mem;
	cur_mem += sizeof(Velocity) * num;
	spva1.pcl_u_s = (Displacement*)cur_mem;
	cur_mem += sizeof(Displacement) * num;
	spva1.pcl_u_f = (Displacement*)cur_mem;
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

	alloc_contact_mem(num);
}

void Model_T3D_CHM_mt::alloc_pcls(size_t num, size_t ori_num)
{
	clear_pcls();
	if (num == 0 || ori_num == 0)
		return;
	alloc_pcls(ori_num);
	pcl_num = num;
}

int Model_T3D_CHM_mt::init_pcls(
	size_t num,
	double n,
	double m_s,
	double density_s,
	double density_f,
	double _Kf,
	double _k,
	double _miu)
{
	if (num == 0)
		return -1;
	alloc_pcls(num);
	SortedPclVarArrays &spva = sorted_pcl_var_arrays[0];
	const double vol_s = m_s / density_s;
	const double vol = vol_s / (1.0 - n);
	size_t p_id;
	for (p_id = 0; p_id < num; ++p_id)
	{
		pcl_m_s[p_id] = m_s;
		pcl_density_s[p_id] = density_s;
		pcl_vol_s[p_id] = vol_s;
		pcl_vol[p_id] = vol;
		Force &bf_s = pcl_bf_s[p_id];
		bf_s.fx = 0.0;
		bf_s.fy = 0.0;
		bf_s.fz = 0.0;
		Force &bf_f = pcl_bf_f[p_id];
		bf_f.fx = 0.0;
		bf_f.fy = 0.0;
		bf_f.fz = 0.0;
		Force &t = pcl_t[p_id];
		t.fx = 0.0;
		t.fy = 0.0;
		t.fz = 0.0;
		spva.pcl_index[p_id] = p_id;
		spva.pcl_n[p_id] = n;
		spva.pcl_density_f[p_id] = density_f;
		Displacement& u_s = spva.pcl_u_s[p_id];
		u_s.ux = 0.0;
		u_s.uy = 0.0;
		u_s.uz = 0.0;
		Displacement &u_f = spva.pcl_u_f[p_id];
		u_f.ux = 0.0;
		u_f.uy = 0.0;
		u_f.uz = 0.0;
		Velocity& v_s = spva.pcl_v_s[p_id];
		v_s.vx = 0.0;
		v_s.vy = 0.0;
		v_s.vz = 0.0;
		Velocity& v_f = spva.pcl_v_f[p_id];
		v_f.vx = 0.0;
		v_f.vy = 0.0;
		v_f.vz = 0.0;
		Stress& p_s = spva.pcl_stress[p_id];
		p_s.s11 = 0.0;
		p_s.s22 = 0.0;
		p_s.s33 = 0.0;
		p_s.s12 = 0.0;
		p_s.s23 = 0.0;
		p_s.s31 = 0.0;
		spva.pcl_p[p_id] = 0.0;
		Strain& p_e = spva.pcl_strain[p_id];
		p_e.e11 = 0.0;
		p_e.e22 = 0.0;
		p_e.e33 = 0.0;
		p_e.e12 = 0.0;
		p_e.e23 = 0.0;
		p_e.e31 = 0.0;
		Strain& p_ee = spva.pcl_estrain[p_id];
		p_ee.e11 = 0.0;
		p_ee.e22 = 0.0;
		p_ee.e33 = 0.0;
		p_ee.e12 = 0.0;
		p_ee.e23 = 0.0;
		p_ee.e31 = 0.0;
		Strain& p_pe = spva.pcl_pstrain[p_id];
		p_pe.e11 = 0.0;
		p_pe.e22 = 0.0;
		p_pe.e33 = 0.0;
		p_pe.e12 = 0.0;
		p_pe.e23 = 0.0;
		p_pe.e31 = 0.0;
		pcl_mat_model[p_id] = nullptr;
	}

	Kf = _Kf; k = _k; miu = _miu;
	
	if (bfx_ss && bfx_s_num)
	{
		for (size_t bf_id = 0; bf_id < bfx_s_num; ++bf_id)
		{
			const BodyForceAtPcl& bf = bfx_ss[bf_id];
			p_id = bf.pcl_id;
			pcl_bf_s[p_id].fx += pcl_m_s[p_id] * bf.bf;
		}
		clear_bfx_ss();
	}

	if (bfy_ss && bfy_s_num)
	{
		for (size_t bf_id = 0; bf_id < bfy_s_num; ++bf_id)
		{
			const BodyForceAtPcl& bf = bfy_ss[bf_id];
			p_id = bf.pcl_id;
			pcl_bf_s[p_id].fy += pcl_m_s[p_id] * bf.bf;
		}
		clear_bfy_ss();
	}

	if (bfz_ss && bfz_s_num)
	{
		for (size_t bf_id = 0; bf_id < bfz_s_num; ++bf_id)
		{
			const BodyForceAtPcl& bf = bfz_ss[bf_id];
			p_id = bf.pcl_id;
			pcl_bf_s[p_id].fz += pcl_m_s[p_id] * bf.bf;
		}
		clear_bfz_ss();
	}
	
	double pcl_m_f;
	if (bfx_fs && bfx_f_num)
	{
		for (size_t bf_id = 0; bf_id < bfx_f_num; ++bf_id)
		{
			const BodyForceAtPcl& bf = bfx_fs[bf_id];
			p_id = bf.pcl_id;
			pcl_m_f = pcl_vol_s[p_id] / (1.0 - spva.pcl_n[p_id])
				* spva.pcl_n[p_id] * spva.pcl_density_f[p_id];
			pcl_bf_f[p_id].fx += pcl_m_f * bf.bf;
		}
		clear_bfx_fs();
	}

	if (bfy_fs && bfy_f_num)
	{
		for (size_t bf_id = 0; bf_id < bfy_f_num; ++bf_id)
		{
			const BodyForceAtPcl& bf = bfy_fs[bf_id];
			p_id = bf.pcl_id;
			pcl_m_f = pcl_vol_s[p_id] / (1.0 - spva.pcl_n[p_id])
				* spva.pcl_n[p_id] * spva.pcl_density_f[p_id];
			pcl_bf_f[p_id].fy += pcl_m_f * bf.bf;
		}
		clear_bfy_fs();
	}

	if (bfz_fs && bfz_f_num)
	{
		for (size_t bf_id = 0; bf_id < bfz_f_num; ++bf_id)
		{
			const BodyForceAtPcl& bf = bfz_fs[bf_id];
			p_id = bf.pcl_id;
			pcl_m_f = pcl_vol_s[p_id] / (1.0 - spva.pcl_n[p_id])
				* spva.pcl_n[p_id] * spva.pcl_density_f[p_id];
			pcl_bf_f[p_id].fz += pcl_m_f * bf.bf;
		}
		clear_bfz_fs();
	}
	
	if (txs && tx_num)
	{
		for (size_t t_id = 0; t_id < tx_num; ++t_id)
		{
			const TractionBCAtPcl& t = txs[t_id];
			pcl_t[t.pcl_id].fx += t.t;
		}
		clear_txs();
	}

	if (tys && ty_num)
	{
		for (size_t t_id = 0; t_id < ty_num; ++t_id)
		{
			const TractionBCAtPcl& t = tys[t_id];
			pcl_t[t.pcl_id].fy += t.t;
		}
		clear_tys();
	}

	if (tzs && tz_num)
	{
		for (size_t t_id = 0; t_id < tz_num; ++t_id)
		{
			const TractionBCAtPcl& t = tzs[t_id];
			pcl_t[t.pcl_id].fz += t.t;
		}
		clear_tzs();
	}
	
	return 0;
}

int Model_T3D_CHM_mt::init_pcls(
	ParticleGenerator3D<TetrahedronMesh>& pg,
	double n, double density_s, double density_f,
	double _Kf, double _k, double _miu)
{
	int res = init_pcls(pg.get_num(), n, density_s, density_s, density_f, _Kf, _k, _miu);
	if (res)
		return res;

	SortedPclVarArrays& spva = sorted_pcl_var_arrays[0];
	auto *pg_pcl = pg.first();
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		Position &p_p = pcl_pos[p_id];
		p_p.x = pg_pcl->x;
		p_p.y = pg_pcl->y;
		p_p.z = pg_pcl->z;
		pcl_vol[p_id] = pg_pcl->vol;
		pcl_vol_s[p_id] = pg_pcl->vol * (1.0 - spva.pcl_n[p_id]);
		pcl_m_s[p_id] *= pcl_vol_s[p_id];
		Force &bf_s = pcl_bf_s[p_id];
		bf_s.fx *= pcl_vol_s[p_id];
		bf_s.fy *= pcl_vol_s[p_id];
		bf_s.fz *= pcl_vol_s[p_id];
		const double pcl_vol_f = pcl_vol[p_id] * spva.pcl_n[p_id];
		Force& bf_f = pcl_bf_f[p_id];
		bf_f.fx *= pcl_vol_f;
		bf_f.fy *= pcl_vol_f;
		bf_f.fz *= pcl_vol_f;
		pg_pcl = pg.next(pg_pcl);
	}
	return 0;
}

void Model_T3D_CHM_mt::init_bfx_ss(
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
			pcl_bf_s[p_id].fx += pcl_m_s[p_id] * bfs[p_id];
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

void Model_T3D_CHM_mt::init_bfy_ss(
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
			pcl_bf_s[p_id].fy += pcl_m_s[p_id] * bfs[p_id];
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

void Model_T3D_CHM_mt::init_bfz_ss(
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
			pcl_bf_s[p_id].fz += pcl_m_s[p_id] * bfs[p_id];
		}
	}
	else
	{
		init_bfz_ss(bf_num);
		for (size_t bf_id = 0; bf_id < bfz_s_num; ++bf_id)
		{
			BodyForceAtPcl& bf = bfz_ss[bf_id];
			bf.pcl_id = bf_pcls[bf_id];
			bf.bf = bfs[bf_id];
		}
	}
}

void Model_T3D_CHM_mt::init_bfx_fs(
	size_t bf_num,
	const size_t* bf_pcls,
	const double* bfs)
{
	if (pcl_mem_raw && ori_pcl_num)
	{
		SortedPclVarArrays& spva = sorted_pcl_var_arrays[0];
		for (size_t bf_id = 0; bf_id < bf_num; ++bf_id)
		{
			const size_t p_id = bf_pcls[bf_id];
			const double pcl_m_f = pcl_vol_s[p_id] / (1.0 - spva.pcl_n[p_id])
				* spva.pcl_n[p_id] * spva.pcl_density_f[p_id];
			pcl_bf_f[p_id].fx += pcl_m_f * bfs[p_id];
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

void Model_T3D_CHM_mt::init_bfy_fs(
	size_t bf_num,
	const size_t* bf_pcls,
	const double* bfs)
{
	if (pcl_mem_raw && ori_pcl_num)
	{
		SortedPclVarArrays& spva = sorted_pcl_var_arrays[0];
		for (size_t bf_id = 0; bf_id < bf_num; ++bf_id)
		{
			const size_t p_id = bf_pcls[bf_id];
			const double pcl_m_f = pcl_vol_s[p_id] / (1.0 - spva.pcl_n[p_id])
				* spva.pcl_n[p_id] * spva.pcl_density_f[p_id];
			pcl_bf_f[p_id].fy += pcl_m_f * bfs[p_id];
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

void Model_T3D_CHM_mt::init_bfz_fs(
	size_t bf_num,
	const size_t* bf_pcls,
	const double* bfs)
{
	if (pcl_mem_raw && ori_pcl_num)
	{
		SortedPclVarArrays& spva = sorted_pcl_var_arrays[0];
		for (size_t bf_id = 0; bf_id < bf_num; ++bf_id)
		{
			const size_t p_id = bf_pcls[bf_id];
			const double pcl_m_f = pcl_vol_s[p_id] / (1.0 - spva.pcl_n[p_id])
				* spva.pcl_n[p_id] * spva.pcl_density_f[p_id];
			pcl_bf_f[p_id].fz += pcl_m_f * bfs[p_id];
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

void Model_T3D_CHM_mt::init_txs(
	size_t t_num,
	const size_t* t_pcls,
	const double* ts)
{
	if (pcl_mem_raw && ori_pcl_num)
	{
		for (size_t t_id = 0; t_id < t_num; ++t_id)
			pcl_t[t_pcls[t_id]].fx += ts[t_id];
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

void Model_T3D_CHM_mt::init_tys(
	size_t t_num,
	const size_t* t_pcls,
	const double* ts)
{
	if (pcl_mem_raw && ori_pcl_num)
	{
		for (size_t t_id = 0; t_id < t_num; ++t_id)
			pcl_t[t_pcls[t_id]].fy += ts[t_id];
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

void Model_T3D_CHM_mt::init_tzs(
	size_t t_num,
	const size_t* t_pcls,
	const double* ts)
{
	if (pcl_mem_raw && ori_pcl_num)
	{
		for (size_t t_id = 0; t_id < t_num; ++t_id)
			pcl_t[t_pcls[t_id]].fz += ts[t_id];
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

void Model_T3D_CHM_mt::init_fixed_vx_s_bc(
	size_t bc_num, const size_t* bcs)
{
	for (size_t bc_id = 0; bc_id < bc_num; ++bc_id)
	{
		if (bcs[bc_id] < node_num)
			node_has_vbc_s[bcs[bc_id]].has_vx_bc = true;
	}
}

void Model_T3D_CHM_mt::init_fixed_vy_s_bc(
	size_t bc_num, const size_t* bcs)
{
	for (size_t bc_id = 0; bc_id < bc_num; ++bc_id)
	{
		if (bcs[bc_id] < node_num)
			node_has_vbc_s[bcs[bc_id]].has_vy_bc = true;
	}
}

void Model_T3D_CHM_mt::init_fixed_vz_s_bc(
	size_t bc_num, const size_t* bcs)
{
	for (size_t bc_id = 0; bc_id < bc_num; ++bc_id)
	{
		if (bcs[bc_id] < node_num)
			node_has_vbc_s[bcs[bc_id]].has_vz_bc = true;
	}
}

void Model_T3D_CHM_mt::init_fixed_vx_f_bc(
	size_t bc_num, const size_t* bcs)
{
	for (size_t bc_id = 0; bc_id < bc_num; ++bc_id)
	{
		if (bcs[bc_id] < node_num)
			node_has_vbc_f[bcs[bc_id]].has_vx_bc = true;
	}
}

void Model_T3D_CHM_mt::init_fixed_vy_f_bc(
	size_t bc_num, const size_t* bcs)
{
	for (size_t bc_id = 0; bc_id < bc_num; ++bc_id)
	{
		if (bcs[bc_id] < node_num)
			node_has_vbc_f[bcs[bc_id]].has_vy_bc = true;
	}
}

void Model_T3D_CHM_mt::init_fixed_vz_f_bc(
	size_t bc_num, const size_t* bcs)
{
	for (size_t bc_id = 0; bc_id < bc_num; ++bc_id)
	{
		if (bcs[bc_id] < node_num)
			node_has_vbc_f[bcs[bc_id]].has_vz_bc = true;
	}
}

void Model_T3D_CHM_mt::set_vbc_vec_s(
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

void Model_T3D_CHM_mt::set_vbc_vec_f(
	size_t n_id,
	double vecx,
	double vecy,
	double vecz)
{
	const double len = sqrt(vecx * vecx + vecy * vecy + vecz * vecz);
	if (len != 0.0)
	{
		NodeVBCVec& nbcv = node_vbc_vec_f[n_id];
		nbcv.x = vecx / len;
		nbcv.y = vecy / len;
		nbcv.z = vecz / len;
	}
}

void Model_T3D_CHM_mt::clear_contact_mem()
{
	if (contact_mem)
	{
		delete[] contact_mem;
		contact_mem = nullptr;
	}
}

void Model_T3D_CHM_mt::alloc_contact_mem(size_t num)
{
	clear_contact_mem();

	if (num == 0)
		return;

	contact_mem = new char[(sizeof(size_t)
		+ sizeof(ContactModel3D::Position)
		+ sizeof(ContactModel3D::Force)
		+ sizeof(double)) * num * 2];
	char* cur_mem = contact_mem;
	// solid
	contact_substep_ids_s = (size_t*)cur_mem;
	cur_mem += sizeof(size_t) * num;
	prev_contact_poses_s = (ContactModel3D::Position*)cur_mem;
	cur_mem += sizeof(ContactModel3D::Position) * num;
	prev_contact_tan_forces_s = (ContactModel3D::Force *)cur_mem;
	cur_mem += sizeof(ContactModel3D::Force) * num;
	prev_contact_dists_s = (double *)cur_mem;
	cur_mem += sizeof(double) * num;
	// fluid
	contact_substep_ids_f = (size_t*)cur_mem;
	cur_mem += sizeof(size_t) * num;
	prev_contact_poses_f = (ContactModel3D::Position*)cur_mem;
	cur_mem += sizeof(ContactModel3D::Position) * num;
	prev_contact_tan_forces_f = (ContactModel3D::Force*)cur_mem;
	cur_mem += sizeof(ContactModel3D::Force) * num;
	prev_contact_dists_f = (double *)cur_mem;

	smooth_contact_s.contact_substep_ids = contact_substep_ids_s;
	smooth_contact_s.prev_contact_dists = prev_contact_dists_s;
	rough_contact_s.contact_substep_ids = contact_substep_ids_s;
	rough_contact_s.prev_contact_poses = prev_contact_poses_s;
	rough_contact_s.prev_contact_tan_forces = prev_contact_tan_forces_s;
	rough_contact_s.prev_contact_dists = prev_contact_dists_s;
	fric_contact_s.contact_substep_ids = contact_substep_ids_s;
	fric_contact_s.prev_contact_poses = prev_contact_poses_s;
	fric_contact_s.prev_contact_tan_forces = prev_contact_tan_forces_s;
	fric_contact_s.prev_contact_dists = prev_contact_dists_s;
	sticky_contact_s.contact_substep_ids = contact_substep_ids_s;
	sticky_contact_s.prev_contact_poses = prev_contact_poses_s;
	sticky_contact_s.prev_contact_tan_forces = prev_contact_tan_forces_s;
	sticky_contact_s.prev_contact_dists = prev_contact_dists_s;

	smooth_contact_f.contact_substep_ids = contact_substep_ids_f;
	smooth_contact_f.prev_contact_dists = prev_contact_dists_f;
	rough_contact_f.contact_substep_ids = contact_substep_ids_f;
	rough_contact_f.prev_contact_poses = prev_contact_poses_f;
	rough_contact_f.prev_contact_tan_forces = prev_contact_tan_forces_f;
	rough_contact_f.prev_contact_dists = prev_contact_dists_f;
}
