#include "Simulations_pcp.h"

#include "Geometry2D.h"
#include "Model_T2D_CHM_d.h"

Model_T2D_CHM_d::Model_T2D_CHM_d() :
	Model("Model_T2D_CHM_d"),
	spcl_num(0), spcls(nullptr),
	fpcl_num(0), fpcls(nullptr),
	bfx_num(0), bfxs(nullptr),
	bfy_num(0), bfys(nullptr),
	tx_num(0), txs(nullptr),
	ty_num(0), tys(nullptr),
	asx_num(0), asxs(nullptr),
	asy_num(0), asys(nullptr),
	afx_num(0), afxs(nullptr),
	afy_num(0), afys(nullptr),
	vsx_num(0), vsxs(nullptr),
	vsy_num(0), vsys(nullptr),
	vfx_num(0), vfxs(nullptr),
	vfy_num(0), vfys(nullptr),
	rigid_circle_is_init(false) {}

Model_T2D_CHM_d::~Model_T2D_CHM_d()
{
	clear_solid_pcls();
	clear_fluid_pcls();

	clear_bfxs();
	clear_bfys();
	clear_txs();
	clear_tys();
	clear_asxs();
	clear_asys();
	clear_afxs();
	clear_afys();
	clear_vsxs();
	clear_vsys();
	clear_vfxs();
	clear_vfys();
}

void Model_T2D_CHM_d::init_mesh_shape_funcs()
{
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		Element& e = elems[e_id];
		e.pt_in_tri.init_triangle<Node>(
			nodes[e.n1],
			nodes[e.n2],
			nodes[e.n3]
			);
	}
}

int Model_T2D_CHM_d::init_mesh(
	double *node_coords, size_t n_num,
	size_t *elem_n_ids,  size_t e_num
	)
{
	int res = ParentMesh::init_mesh(node_coords, n_num, elem_n_ids, e_num);
	if (res < 0)
		return res;
	init_mesh_shape_funcs();
	return 0;
}

int Model_T2D_CHM_d::load_mesh_from_hdf5(const char* file_name)
{
	int res = ParentMesh::load_mesh_from_hdf5(file_name);
	if (res < 0)
		return res;
	init_mesh_shape_funcs();
	return 0;
}

int Model_T2D_CHM_d::init_search_grid(double _hx, double _hy)
{
	return search_bg_grid.init(*this, _hx, _hy);
}

int Model_T2D_CHM_d::init_solid_pcls(
	size_t num,
	double m,
	double density,
	double n,
	double _k
	)
{
	clear_solid_pcls();
	if (num == 0)
		return -1;

	alloc_solid_pcls(num);
	for (size_t pcl_id = 0; pcl_id < num; ++pcl_id)
	{
		SolidParticle &pcl = spcls[pcl_id];
		pcl.id = pcl_id;
		pcl.ux = 0.0;
		pcl.uy = 0.0;
		pcl.vx = 0.0;
		pcl.vy = 0.0;
		pcl.n = n;
		pcl.m = m;
		pcl.density = density;
		pcl.e11 = 0.0;
		pcl.e22 = 0.0;
		pcl.e12 = 0.0;
		pcl.s11 = 0.0;
		pcl.s12 = 0.0;
		pcl.s22 = 0.0;
	}

	k = _k;
	return 0;
}

int Model_T2D_CHM_d::init_fluid_pcls(
	size_t num,
	double m,
	double density,
	double _Kf,
	double _miu
	)
{
	clear_fluid_pcls();
	if (num == 0)
		return -1;

	alloc_fluid_pcls(num);
	for (size_t pcl_id = 0; pcl_id < num; ++pcl_id)
	{
		FluidParticle& pcl = fpcls[pcl_id];
		pcl.id = pcl_id;
		pcl.ux = 0.0;
		pcl.uy = 0.0;
		pcl.vx = 0.0;
		pcl.vy = 0.0;
		pcl.m = m;
		pcl.density = density;
		pcl.p;
	}

	Kf = _Kf;
	miu = _miu;
	return 0;
}

int Model_T2D_CHM_d::init_solid_pcls(
	ParticleGenerator2D<Model_T2D_CHM_d>& pg,
	double density,
	double n,
	double _k
	)
{
	int res = init_solid_pcls(pg.get_num(), density, density, n, _k);
	if (res)
		return res;

	typedef ParticleGenerator2D<Model_T2D_CHM_d>::Particle PgPcl;
	PgPcl* pg_pcl = pg.first();
	for (size_t pcl_id = 0; pcl_id < spcl_num; ++pcl_id)
	{
		SolidParticle &pcl = spcls[pcl_id];
		pcl.x = pg_pcl->x;
		pcl.y = pg_pcl->y;
		pcl.m *= pg_pcl->area * (1.0 - pcl.n);
		pg_pcl = pg.next(pg_pcl);
	}

	return 0;
}

int Model_T2D_CHM_d::init_fluid_pcls(
	ParticleGenerator2D<Model_T2D_CHM_d>& pg,
	double density,
	double _Kf,
	double _miu
	)
{
	int res = init_fluid_pcls(pg.get_num(), density, density, _Kf, _miu);
	if (res)
		return res;

	typedef ParticleGenerator2D<Model_T2D_CHM_d>::Particle PgPcl;
	PgPcl* pg_pcl = pg.first();
	for (size_t pcl_id = 0; pcl_id < fpcl_num; ++pcl_id)
	{
		FluidParticle& pcl = fpcls[pcl_id];
		pcl.x = pg_pcl->x;
		pcl.y = pg_pcl->y;
		pcl.m *= pg_pcl->area;
		pg_pcl = pg.next(pg_pcl);
	}

	return 0;
}

void Model_T2D_CHM_d::alloc_solid_pcls(size_t num)
{
	spcls = new SolidParticle[num];
	spcl_num = num;
}

void Model_T2D_CHM_d::clear_solid_pcls()
{
	if (spcls)
	{
		delete[] spcls;
		spcls = nullptr;
	}
	spcl_num = 0;
}

void Model_T2D_CHM_d::alloc_fluid_pcls(size_t num)
{
	fpcls = new FluidParticle[num];
	fpcl_num = num;
}

void Model_T2D_CHM_d::clear_fluid_pcls()
{
	if (fpcls)
	{
		delete[] fpcls;
		fpcls = nullptr;
	}
	fpcl_num = 0;
}

//int Model_T2D_CHM_d::apply_rigid_circle(double dt)
//{
	//rigid_circle.reset_rf();

	//double dist, norm_x, norm_y;
	//double fs_cont, fsx_cont, fsy_cont;
	//double ff_cont, ffx_cont, ffy_cont;
	//double nfsx_cont, nfsy_cont, ndasx, ndasy;
	//double nffx_cont, nffy_cont, ndafx, ndafy;
	//for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	//{
	//	Particle& pcl = pcls[p_id];
	//	if (pcl.pe && rigid_circle.detect_collision_with_point(
	//		pcl.x, pcl.y, pcl.vol, dist, norm_x, norm_y))
	//	{
	//		fs_cont = Ks_cont * dist;
	//		fsx_cont = fs_cont * norm_x;
	//		fsy_cont = fs_cont * norm_y;
	//		ff_cont = Kf_cont * dist;
	//		ffx_cont = ff_cont * norm_x;
	//		ffy_cont = ff_cont * norm_y;
	//		// reaction force by the rigid object
	//		rigid_circle.add_rf(pcl.x, pcl.y,
	//			-(fsx_cont + ffx_cont), -(fsy_cont + ffy_cont));
	//		// adjust velocity at nodes
	//		Element& e = *pcl.pe;
	//		// node 1
	//		Node& n1 = nodes[e.n1];
	//		nfsx_cont = pcl.N1 * fsx_cont;
	//		nfsy_cont = pcl.N1 * fsy_cont;
	//		ndasx = nfsx_cont / n1.m_s;
	//		ndasy = nfsy_cont / n1.m_s;
	//		n1.ax_s += ndasx;
	//		n1.ay_s += ndasy;
	//		n1.vx_s += ndasx * dt;
	//		n1.vy_s += ndasy * dt;
	//		nffx_cont = pcl.N1 * ffx_cont;
	//		nffy_cont = pcl.N1 * ffy_cont;
	//		ndafx = nffx_cont / n1.m_f;
	//		ndafy = nffy_cont / n1.m_f;
	//		n1.ax_f += ndafx;
	//		n1.ay_f += ndafy;
	//		n1.vx_f += ndafx * dt;
	//		n1.vy_f += ndafy * dt;
	//		// node 2
	//		Node& n2 = nodes[e.n2];
	//		nfsx_cont = pcl.N2 * fsx_cont;
	//		nfsy_cont = pcl.N2 * fsy_cont;
	//		ndasx = nfsx_cont / n2.m_s;
	//		ndasy = nfsy_cont / n2.m_s;
	//		n2.ax_s += ndasx;
	//		n2.ay_s += ndasy;
	//		n2.vx_s += ndasx * dt;
	//		n2.vy_s += ndasy * dt;
	//		nffx_cont = pcl.N2 * ffx_cont;
	//		nffy_cont = pcl.N2 * ffy_cont;
	//		ndafx = nffx_cont / n2.m_f;
	//		ndafy = nffy_cont / n2.m_f;
	//		n2.ax_f += ndafx;
	//		n2.ay_f += ndafy;
	//		n2.vx_f += ndafx * dt;
	//		n2.vy_f += ndafy * dt;
	//		// node 3
	//		Node& n3 = nodes[e.n3];
	//		nfsx_cont = pcl.N3 * fsx_cont;
	//		nfsy_cont = pcl.N3 * fsy_cont;
	//		ndasx = nfsx_cont / n3.m_s;
	//		ndasy = nfsy_cont / n3.m_s;
	//		n3.ax_s += ndasx;
	//		n3.ay_s += ndasy;
	//		n3.vx_s += ndasx * dt;
	//		n3.vy_s += ndasy * dt;
	//		nffx_cont = pcl.N3 * ffx_cont;
	//		nffy_cont = pcl.N3 * ffy_cont;
	//		ndafx = nffx_cont / n3.m_f;
	//		ndafy = nffy_cont / n3.m_f;
	//		n3.ax_f += ndafx;
	//		n3.ay_f += ndafy;
	//		n3.vx_f += ndafx * dt;
	//		n3.vy_f += ndafy * dt;
	//	}
	//}

	//rigid_circle.update_motion(dt);
//	return 0;
//}
