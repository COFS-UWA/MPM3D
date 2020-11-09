#include "SimulationsOMP_pcp.h"

#include "Model_R2D_ME_mt.h"

Model_R2D_ME_mt::Model_R2D_ME_mt() :
	ori_pcl_num(0), pcl_num(0),
	pcl_mem_raw(nullptr), pcl_mat_model(nullptr),
	elem_x_num(0), elem_y_num(0), elem_num(0),
	node_x_num(0), node_y_num(0), node_num(0),
	actual_elem_x_num(0), actual_elem_num(0),
	mesh_mem_raw(nullptr),
	bfx_num(0), bfxs(nullptr),
	bfy_num(0), bfys(nullptr),
	tx_num(0), txs(nullptr),
	ty_num(0), tys(nullptr),
	rigid_rect_is_valid(false) {}

Model_R2D_ME_mt::~Model_R2D_ME_mt()
{
	clear_mesh();
	clear_pcls();
	clear_bfxs();
	clear_bfys();
	clear_txs();
	clear_tys();
}

void Model_R2D_ME_mt::clear_mesh()
{
	if (mesh_mem_raw)
	{
		delete[] mesh_mem_raw;
		mesh_mem_raw = nullptr;
	}
	elem_x_num = 0;
	elem_y_num = 0;
	elem_num = 0;
	node_x_num = 0;
	node_y_num = 0;
	node_num = 0;
	actual_elem_x_num = 0;
	actual_elem_num = 0;
}

void Model_R2D_ME_mt::alloc_mesh(
	size_t n_num,
	size_t e_num
	)
{
	clear_mesh();

	node_num = n_num;
	elem_num = e_num;

	size_t mem_len = (sizeof(double) * 2 + sizeof(size_t)
		+ sizeof(ElemNodeVM) * 4 + sizeof(ElemNodeForce) * 4) * e_num
		+ (sizeof(NodeA) + sizeof(NodeV) + sizeof(NodeHasVBC)) * n_num;
	mesh_mem_raw = new char[mem_len];

	char *cur_mem = mesh_mem_raw;
	elem_density = (double *)(cur_mem); // elem_num
	cur_mem += sizeof(double) * elem_num;
	elem_pcl_m = (double *)cur_mem; // elem_num
	cur_mem += sizeof(double) * elem_num;
	elem_substp_id = (size_t *)cur_mem; // elem_num
	cur_mem += sizeof(size_t) * elem_num;

	elem_node_vm = (ElemNodeVM*)cur_mem; // elem_num * 4
	cur_mem += sizeof(ElemNodeVM) * elem_num * 4;
	elem_node_force = (ElemNodeForce *)cur_mem; // elem_num * 4
	cur_mem += sizeof(ElemNodeForce) * elem_num * 4;

	node_a = (NodeA *)cur_mem; // node_num
	cur_mem += sizeof(NodeA) * node_num;
	node_v = (NodeV *)cur_mem; // node_num
	cur_mem += sizeof(NodeV) * node_num;
	node_has_vbc = (NodeHasVBC *)cur_mem; // node_num
}

#define N_low(xi)  ((1.0 - (xi)) * 0.5)
#define N_high(xi) ((1.0 + (xi)) * 0.5)
#define dN_dxi_low -0.5
#define dN_dxi_high 0.5

int Model_R2D_ME_mt::init_mesh(
	double _xl, double _yl,
	double _xu,	double _yu,
	size_t _x_num, size_t _y_num
	)
{
	elem_x_num = _x_num;
	elem_y_num = _y_num;
	elem_num = elem_x_num * elem_y_num;
	node_x_num = elem_x_num + 1;
	node_y_num = elem_y_num + 1;
	node_num = node_x_num * node_y_num;
	actual_elem_x_num = elem_num + 2;
	actual_elem_num = actual_elem_x_num * (elem_y_num + 2);
	mh_xl = _xl;
	mh_yl = _yl;
	mh_xu = _xu;
	mh_yu = _yu;
	elem_hx = (mh_xu - mh_xl) / double(node_x_num);
	elem_hy = (mh_yu - mh_yl) / double(node_y_num);
	mh_xl -= elem_hx;
	mh_xu += elem_hx;
	mh_yl -= elem_hy;
	mh_yu += elem_hy;
	elem_area = elem_hx * elem_hy;
	inv_elem_hx = 1.0 / elem_hx;
	inv_elem_hy = 1.0 / elem_hy;
	dxi_dx = 2.0 * inv_elem_hx;
	deta_dy = 2.0 * inv_elem_hy;
	elem_dN.dN1_dx = dN_dxi_low * N_low(0.5) * dxi_dx;
	elem_dN.dN1_dy = N_low(0.5) * dN_dxi_low * deta_dy;
	elem_dN.dN2_dx = dN_dxi_high * N_low(0.5) * dxi_dx;
	elem_dN.dN2_dy = N_high(0.5) * dN_dxi_low * deta_dy;
	elem_dN.dN3_dx = dN_dxi_low * N_high(0.5) * dxi_dx;
	elem_dN.dN3_dy = N_low(0.5) * dN_dxi_high * deta_dy;
	elem_dN.dN4_dx = dN_dxi_high * N_high(0.5) * dxi_dx;
	elem_dN.dN4_dy = N_high(0.5) * dN_dxi_high * deta_dy;

	alloc_mesh(node_num, actual_elem_num);

	memset(elem_density, 0, sizeof(double) * elem_num);
	memset(elem_pcl_m, 0, sizeof(double) * elem_num);
	memset(elem_substp_id, 0xff, sizeof(double) * elem_num);
	memset(elem_node_vm, 0, sizeof(ElemNodeVM) * elem_num * 4);
	memset(elem_node_force, 0, sizeof(ElemNodeForce) * elem_num * 4);
	memset(node_has_vbc, 0, sizeof(NodeHasVBC) * node_num);
	return 0;
}

#undef N_low
#undef N_high
#undef dN_dxi_low
#undef dN_dxi_high

void Model_R2D_ME_mt::clear_pcls()
{
	if (pcl_mem_raw)
	{
		delete[] pcl_mem_raw;
		pcl_mem_raw = nullptr;
	}
	ori_pcl_num = 0;
	pcl_num = 0;
}

void Model_R2D_ME_mt::alloc_pcls(size_t num)
{
	clear_pcls();

	ori_pcl_num = num;
	pcl_num = num;

	size_t mem_len = (sizeof(double)
		+ sizeof(PclBodyForce) + sizeof(PclTraction)
		+ sizeof(PclPos) + sizeof(double)
		+ sizeof(ShapeFunc) + sizeof(DShapeFunc)
		+ sizeof(MatModel::MaterialModel *)
		+ (sizeof(size_t) + sizeof(double)
		+ sizeof(PclDisp) + sizeof(PclV)
		+ sizeof(Stress) + sizeof(Strain) * 3) * 2
		+ sizeof(size_t) + sizeof(ContPos)) * ori_pcl_num;
	pcl_mem_raw = new char[mem_len];

	char* cur_mem = pcl_mem_raw;
	pcl_m = (double *)cur_mem; // ori_pcl_num
	cur_mem += sizeof(double) * ori_pcl_num;
	pcl_bf = (PclBodyForce *)cur_mem; // ori_pcl_num
	cur_mem += sizeof(PclBodyForce) * ori_pcl_num;
	pcl_t = (PclTraction *)cur_mem; // ori_pcl_num
	cur_mem += sizeof(PclTraction) * ori_pcl_num;
	pcl_pos = (PclPos *)cur_mem; // ori_pcl_num
	cur_mem += sizeof(PclPos) * ori_pcl_num;
	pcl_vol = (double *)cur_mem; // ori_pcl_num
	cur_mem += sizeof(double) * ori_pcl_num;
	pcl_N = (ShapeFunc *)cur_mem; // ori_pcl_num
	cur_mem += sizeof(ShapeFunc) * ori_pcl_num;
	pcl_dN = (DShapeFunc *)cur_mem; // ori_pcl_num
	cur_mem += sizeof(DShapeFunc) * ori_pcl_num;
	pcl_mat_model = (MatModel::MaterialModel **)cur_mem; // ori_pcl_num
	cur_mem += sizeof(MatModel::MaterialModel *) * ori_pcl_num;

	PclSortedVarArray& psva0 = pcl_sorted_var_array[0];
	psva0.pcl_index = (size_t *)cur_mem; // ori_pcl_num
	cur_mem += sizeof(size_t) * ori_pcl_num;
	psva0.pcl_density = (double *)cur_mem; // ori_pcl_num
	cur_mem += sizeof(double) * ori_pcl_num;
	psva0.pcl_disp = (PclDisp *)cur_mem; // ori_pcl_num
	cur_mem += sizeof(PclDisp) * ori_pcl_num;
	psva0.pcl_v = (PclV *)cur_mem; // ori_pcl_num
	cur_mem += sizeof(PclV) * ori_pcl_num;
	psva0.pcl_stress = (Stress *)cur_mem; // ori_pcl_num
	cur_mem += sizeof(Stress) * ori_pcl_num;
	psva0.pcl_strain = (Strain *)cur_mem; // ori_pcl_num
	cur_mem += sizeof(Strain) * ori_pcl_num;
	psva0.pcl_estrain = (Strain *)cur_mem; // ori_pcl_num
	cur_mem += sizeof(Strain) * ori_pcl_num;
	psva0.pcl_pstrain = (Strain *)cur_mem; // ori_pcl_num
	cur_mem += sizeof(Strain) * ori_pcl_num;

	PclSortedVarArray& psva1 = pcl_sorted_var_array[1];
	psva1.pcl_index = (size_t*)cur_mem; // ori_pcl_num
	cur_mem += sizeof(size_t) * ori_pcl_num;
	psva1.pcl_density = (double*)cur_mem; // ori_pcl_num
	cur_mem += sizeof(double) * ori_pcl_num;
	psva1.pcl_disp = (PclDisp*)cur_mem; // ori_pcl_num
	cur_mem += sizeof(PclDisp) * ori_pcl_num;
	psva1.pcl_v = (PclV*)cur_mem; // ori_pcl_num
	cur_mem += sizeof(PclV) * ori_pcl_num;
	psva1.pcl_stress = (Stress*)cur_mem; // ori_pcl_num
	cur_mem += sizeof(Stress) * ori_pcl_num;
	psva1.pcl_strain = (Strain*)cur_mem; // ori_pcl_num
	cur_mem += sizeof(Strain) * ori_pcl_num;
	psva1.pcl_estrain = (Strain*)cur_mem; // ori_pcl_num
	cur_mem += sizeof(Strain) * ori_pcl_num;
	psva1.pcl_pstrain = (Strain*)cur_mem; // ori_pcl_num
}

void Model_R2D_ME_mt::alloc_pcls(
	size_t num,
	size_t ori_num
	)
{
	alloc_pcls(num);
	ori_pcl_num = ori_num;
}

int Model_R2D_ME_mt::init_pcls(
	size_t num,
	double m,
	double density
	)
{
	alloc_pcls(num);

	PclSortedVarArray &psva = pcl_sorted_var_array[0];
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		pcl_m[p_id] = m;
		PclBodyForce &p_bf = pcl_bf[p_id];
		p_bf.bfx = 0.0;
		p_bf.bfy = 0.0;
		PclTraction& p_t = pcl_t[p_id];
		p_t.tx = 0.0;
		p_t.ty = 0.0;
		psva.pcl_index[p_id] = p_id;
		psva.pcl_density[p_id] = density;
		PclV &p_v = psva.pcl_v[p_id];
		p_v.vx = 0.0;
		p_v.vy = 0.0;
		Stress& p_s = psva.pcl_stress[p_id];
		p_s.s11 = 0.0;
		p_s.s22 = 0.0;
		p_s.s12 = 0.0;
		Strain& p_e = psva.pcl_strain[p_id];
		p_e.e11 = 0.0;
		p_e.e22 = 0.0;
		p_e.e12 = 0.0;
		Strain &p_ee = psva.pcl_estrain[p_id];
		p_ee.e11 = 0.0;
		p_ee.e22 = 0.0;
		p_ee.e12 = 0.0;
		Strain& p_pe = psva.pcl_pstrain[p_id];
		p_pe.e11 = 0.0;
		p_pe.e22 = 0.0;
		p_pe.e12 = 0.0;
	}

	return 0;
}

void Model_R2D_ME_mt::init_pcls(
	double _xl, double _yl,
	double _xu, double _yu,
	size_t _x_num, size_t _y_num,
	double density
	)
{
	double p_hx = (_xu - _xl) / double(_x_num);
	double p_hy = (_yu - _yl) / double(_y_num);
	double p_vol = p_hx * p_hy;

	init_pcls(_x_num * _y_num, density, density);
	
	size_t p_id = 0;
	double p_x, p_y = _yl + p_hy * 0.5;
	for (size_t y_id = 0; y_id < _y_num; ++y_id)
	{
		p_x = _xl + p_hx * 0.5;
		for (size_t x_id = 0; x_id < _x_num; ++x_id)
		{
			pcl_m[p_id] *= p_vol;
			pcl_vol[p_id] = p_vol;
			PclPos& p_p = pcl_pos[p_id];
			p_p.x = p_x;
			p_p.y = p_y;
			++p_id;
			p_x += p_hx;
		}
		p_y += p_hy;
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
}

void Model_R2D_ME_mt::init_bfxs(
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

void Model_R2D_ME_mt::init_bfys(
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
			BodyForceAtPcl& bf = bfys[bf_id];
			bf.pcl_id = bf_pcls[bf_id];
			bf.bf = bfs[bf_id];
		}
	}
}

void Model_R2D_ME_mt::init_txs(
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
			PclTraction& t = pcl_t[p_id];
			t.tx += ts[t_id];
		}
	}
	else
	{
		init_txs(t_num);
		for (size_t t_id = 0; t_id < tx_num; ++t_id)
		{
			TractionBCAtPcl& t = txs[t_id];
			t.pcl_id = t_pcls[t_id];
			t.t = ts[t_id];
		}
	}
}

void Model_R2D_ME_mt::init_tys(
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
			TractionBCAtPcl& t = tys[t_id];
			t.pcl_id = t_pcls[t_id];
			t.t = ts[t_id];
		}
	}
}

void Model_R2D_ME_mt::init_fixed_vx_bc(
	size_t vx_bc_num,
	const size_t* vx_bcs
	)
{
	for (size_t bc_id = 0; bc_id < vx_bc_num; ++bc_id)
	{
		if (vx_bcs[bc_id] < node_num)
			node_has_vbc[vx_bcs[bc_id]].has_vx_bc = true;
	}
}

void Model_R2D_ME_mt::init_fixed_vy_bc(
	size_t vy_bc_num,
	const size_t* vy_bcs
	)
{
	for (size_t bc_id = 0; bc_id < vy_bc_num; ++bc_id)
	{
		if (vy_bcs[bc_id] < node_num)
			node_has_vbc[vy_bcs[bc_id]].has_vy_bc = true;
	}
}
