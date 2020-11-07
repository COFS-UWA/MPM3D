#include "SimulationsOMP_pcp.h"

#include "Model_R2D_ME_mt.h"

Model_R2D_ME_mt::Model_R2D_ME_mt() :
	ori_pcl_num(0), pcl_num(0),
	pcl_mem_raw(nullptr), pcl_mat_model(nullptr),
	elem_x_num(0), elem_y_num(0), elem_num(0),
	node_x_num(0), node_y_num(0), node_num(0),
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
}

void Model_R2D_ME_mt::alloc_mesh(
	size_t n_num,
	size_t e_num
	)
{
	clear_mesh();

	node_num = n_num;
	elem_num = e_num;

	size_t mem_len = (sizeof(ElemShapeFunc)
		+ sizeof(double) * 4
		+ sizeof(ElemStrainInc) + sizeof(ElemStress)
		+ sizeof(ElemNodeVM) * 3 + sizeof(ElemNodeForce) * 3) * e_num
		+ (sizeof(NodeA) + sizeof(NodeV)
		+ sizeof(NodeHasVBC) + sizeof(double) * 2) * n_num;
	mesh_mem_raw = new char[mem_len];

	char *cur_mem = mesh_mem_raw;
	elem_sf = (ElemShapeFunc*)cur_mem;
	cur_mem += sizeof(ElemShapeFunc) * elem_num;

	elem_density = (double*)cur_mem;
	cur_mem += sizeof(double) * elem_num;
	elem_pcl_m = (double*)cur_mem;
	cur_mem += sizeof(double) * elem_num;
	elem_pcl_vol = (double*)cur_mem;
	cur_mem += sizeof(double) * elem_num;
	elem_de = (ElemStrainInc*)cur_mem;
	cur_mem += sizeof(ElemStrainInc) * elem_num;
	elem_stress = (ElemStress*)cur_mem;
	cur_mem += sizeof(ElemStress) * elem_num;
	elem_m_de_vol = (double*)cur_mem;
	cur_mem += sizeof(double) * elem_num;

	elem_node_vm = (ElemNodeVM*)cur_mem;
	cur_mem += sizeof(ElemNodeVM) * elem_num * 3;
	elem_node_force = (ElemNodeForce*)cur_mem;
	cur_mem += sizeof(ElemNodeForce) * elem_num * 3;

	node_a = (NodeA*)cur_mem;
	cur_mem += sizeof(NodeA) * node_num;
	node_v = (NodeV*)cur_mem;
	cur_mem += sizeof(NodeV) * node_num;
	node_has_vbc = (NodeHasVBC*)cur_mem;
	cur_mem += sizeof(NodeHasVBC) * node_num;
	node_am = (double*)cur_mem;
	cur_mem += sizeof(double) * node_num;
	node_de_vol = (double*)cur_mem;
	cur_mem += sizeof(double) * node_num;
}

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
	mh_xl = _xl;
	mh_yl = _yl;
	mh_xu = _xu;
	mh_yu = _yu;
	mh_hx = (_xu - _xl) / double(_x_num);
	mh_hy = (_yu - _yl) / double(_y_num);
	mh_area = mh_hx * mh_hy;
	inv_mh_hx = 1.0 / mh_hx;
	inv_mh_hy = 1.0 / mh_hy;
	dxi_dx = 2.0 * inv_mh_hx;
	deta_dy = 2.0 * inv_mh_hy;

	alloc_mesh(node_num, elem_num);
	
	node_num;

	return 0;
}

void Model_R2D_ME_mt::clear_pcls()
{

}

void Model_R2D_ME_mt::alloc_pcls(size_t num)
{

}

void Model_R2D_ME_mt::alloc_pcls(
	size_t num,
	size_t ori_num
	)
{

}

int Model_R2D_ME_mt::init_pcls(
	size_t num,
	double m,
	double density
	)
{

}

void Model_R2D_ME_mt::init_pcls(
	double _xl, double _yl,
	double _xn, double _yn,
	size_t _x_num, size_t _y_num,
	double density
	)
{

}

void Model_R2D_ME_mt::init_bfxs(
	size_t bf_num,
	const size_t* bf_pcls,
	const double* bfs
	)
{

}

void Model_R2D_ME_mt::init_bfys(
	size_t bf_num,
	const size_t* bf_pcls,
	const double* bfs
	)
{

}

void Model_R2D_ME_mt::init_txs(
	size_t t_num,
	const size_t* t_pcls,
	const double* ts
	)
{

}

void Model_R2D_ME_mt::init_tys(
	size_t t_num,
	const size_t* t_pcls,
	const double* ts
	)
{

}

void Model_R2D_ME_mt::init_fixed_vx_bc(
	size_t vx_bc_num,
	const size_t* vx_bcs
	)
{

}

void Model_R2D_ME_mt::init_fixed_vy_bc(
	size_t vy_bc_num,
	const size_t* vy_bcs
	)
{

}
