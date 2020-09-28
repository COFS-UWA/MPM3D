#include "Simulations_pcp.h"

#include "Geometry2D.h"
#include "Model_T2D_CHM_d.h"

Model_T2D_CHM_d::Model_T2D_CHM_d() :
	Model("Model_T2D_CHM_d"),
	spcl_num(0), spcls(nullptr),
	fpcl_num(0), fpcls(nullptr),
	bfsx_num(0), bfsxs(nullptr),
	bfsy_num(0), bfsys(nullptr),
	bffx_num(0), bffxs(nullptr),
	bffy_num(0), bffys(nullptr),
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
	sgrid_xl(0.0), sgrid_yl(0.0),
	sgrid_hx(0.0), sgrid_hy(0.0),
	sg_x_num(0), sg_y_num(0),
	sgrids(nullptr),
	rigid_circle_is_init(false) {}

Model_T2D_CHM_d::~Model_T2D_CHM_d()
{
	clear_solid_pcls();
	clear_fluid_pcls();

	clear_bfsxs();
	clear_bfsys();
	clear_bffxs();
	clear_bffys();
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

	clear_spcl_grids();
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

int Model_T2D_CHM_d::init_search_grid(
	double _hx,
	double _hy
	)
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
		pcl.p = 0.0;
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

void Model_T2D_CHM_d::init_spcl_grids(double _hx, double _hy)
{

}

void Model_T2D_CHM_d::clear_spcl_grids()
{
	if (sgrids)
	{
		delete[] sgrids;
		sgrids = nullptr;
	}
	sgrid_xl = 0.0;
	sgrid_yl = 0.0;
	sgrid_hx = 0.0;
	sgrid_hy = 0.0;
	sg_x_num = 0;
	sg_y_num = 0;
}

void Model_T2D_CHM_d::reset_spcl_grids()
{
	spcl_pt_buffer.reset_and_optimize();
	
	size_t x_id, y_id;
	GridSPclList *cur_g = sgrids;
	for (y_id = 0; y_id < sg_y_num; ++y_id)
		for (x_id = 0; x_id < sg_x_num; ++x_id)
		{
			cur_g->spcls = nullptr;
			++cur_g;
		}

	double pcl_hlen;
	for (size_t pcl_id = 0; pcl_id < spcl_num; ++pcl_id)
	{
		SolidParticle& spcl = spcls[pcl_id];
		pcl_hlen = 0.5 * sqrt(spcl.vol);
		Rect &pcl_bbox = spcl.bbox;
		pcl_bbox.xl = spcl.x - pcl_hlen;
		pcl_bbox.xu = spcl.x + pcl_hlen;
		pcl_bbox.yl = spcl.y - pcl_hlen;
		pcl_bbox.yu = spcl.y + pcl_hlen;
		//IdRect;

		x_id = (spcl.x - sgrid_xl) / sgrid_hx;
		y_id = (spcl.y - sgrid_yl) / sgrid_hy;
		GridSPclList& g = sgrid_by_id(x_id, y_id);
		add_pcl(g, spcl);
	}
}

bool Model_T2D_CHM_d::is_in_solid(double x, double y)
{
	size_t x_id = (x - sgrid_xl) / sgrid_hx;
	size_t y_id = (y - sgrid_yl) / sgrid_hy;
	GridSPclList& g = sgrid_by_id(x_id, y_id);
	for (SPclPointer *sp = g.spcls; sp; sp = sp->next)
	{
		if (sp->pcl->bbox.is_in_box(x, y))
			return true;
	}
	return false;
}
