#include "Simulations_pcp.h"

#include "Geometry2D.h"
#include "Model_T2D_CHM_s.h"

Model_T2D_CHM_s::Model_T2D_CHM_s() :
	Model("Model_T2D_CHM_s"),
	pcls(nullptr), pcl_num(0),
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

Model_T2D_CHM_s::~Model_T2D_CHM_s()
{
	clear_pcls();

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

void Model_T2D_CHM_s::init_mesh_shape_funcs()
{
	double area2;
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		Element& e = elems[e_id];
		Node& n1 = nodes[e.n1];
		Node& n2 = nodes[e.n2];
		Node& n3 = nodes[e.n3];
		area2 = n1.x * n2.y - n2.x * n1.y
			  + n2.x * n3.y - n3.x * n2.y
			  + n3.x * n1.y - n1.x * n3.y;
		e.a1 = (n2.y - n3.y) / area2;
		e.b1 = (n3.x - n2.x) / area2;
		e.coef1 = (n2.x * n3.y - n3.x * n2.y) / area2;
		e.a2 = (n3.y - n1.y) / area2;
		e.b2 = (n1.x - n3.x) / area2;
		e.coef2 = (n3.x * n1.y - n1.x * n3.y) / area2;
		e.a3 = (n1.y - n2.y) / area2;
		e.b3 = (n2.x - n1.x) / area2;
		e.coef3 = (n1.x * n2.y - n2.x * n1.y) / area2;
		// shape function derivatives
		e.dN1_dx = e.a1;
		e.dN1_dy = e.b1;
		e.dN2_dx = e.a2;
		e.dN2_dy = e.b2;
		e.dN3_dx = e.a3;
		e.dN3_dy = e.b3;
	}
}

int Model_T2D_CHM_s::init_mesh(
	double *node_coords, size_t n_num,
	size_t *elem_n_ids,  size_t e_num
	)
{
	int res = BgMesh::init_mesh(node_coords, n_num, elem_n_ids, e_num);
	if (res < 0)
		return res;
	init_mesh_shape_funcs();
	return 0;
}

int Model_T2D_CHM_s::load_mesh_from_hdf5(const char* file_name)
{
	int res = BgMesh::load_mesh_from_hdf5(file_name);
	if (res < 0)
		return res;
	init_mesh_shape_funcs();
	return 0;
}

int Model_T2D_CHM_s::init_search_grid(double _hx, double _hy)
{
	return search_bg_grid.init(*this, _hx, _hy);
}

int Model_T2D_CHM_s::init_pcls(size_t num,
	double n, double m_s, double density_s, double density_f,
	double _Kf, double _k, double _miu)
{
	clear_pcls();
	if (num == 0)
		return -1;

	alloc_pcls(num);
	for (size_t pcl_id = 0; pcl_id < num; ++pcl_id)
	{
		Particle &pcl = pcls[pcl_id];
		pcl.id = pcl_id;
		pcl.ux_s = 0.0;
		pcl.uy_s = 0.0;
		pcl.ux_f = 0.0;
		pcl.uy_f = 0.0;
		pcl.vx_s = 0.0;
		pcl.vy_s = 0.0;
		pcl.vx_f = 0.0;
		pcl.vy_f = 0.0;
		pcl.n = n;
		pcl.m_s = m_s;
		pcl.density_s = density_s;
		pcl.density_f = density_f;
		pcl.e11 = 0.0;
		pcl.e22 = 0.0;
		pcl.e12 = 0.0;
		pcl.s11 = 0.0;
		pcl.s12 = 0.0;
		pcl.s22 = 0.0;
		pcl.p = 0.0;
	}
	
	Kf = _Kf;
	k = _k;
	miu = _miu;
	return 0;
}

int Model_T2D_CHM_s::init_pcls(
	ParticleGenerator2D<Model_T2D_CHM_s>& pg,
	double n, double density_s, double density_f,
	double _Kf, double _k, double _miu)
{
	int res = init_pcls(pg.get_num(), n, density_s,
						density_s, density_f, _Kf, _k, _miu);
	if (res)
		return res;

	typedef ParticleGenerator2D<Model_T2D_CHM_s>::Particle PgPcl;
	PgPcl* pg_pcl = pg.first();
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		Particle &pcl = pcls[pcl_id];
		pcl.x = pg_pcl->x;
		pcl.y = pg_pcl->y;
		pcl.m_s *= pg_pcl->area * (1.0 - pcl.n);
		pg_pcl = pg.next(pg_pcl);
	}

	return 0;
}

void Model_T2D_CHM_s::alloc_pcls(size_t num)
{
	pcls = new Particle[num];
	pcl_num = num;
}

void Model_T2D_CHM_s::clear_pcls()
{
	if (pcls)
	{
		delete[] pcls;
		pcls = nullptr;
	}
	pcl_num = 0;
}

int Model_T2D_CHM_s::apply_rigid_circle(double dt)
{
	rigid_circle.reset_rf();

	double dist, norm_x, norm_y;
	double fs_cont, fsx_cont, fsy_cont;
	double ff_cont, ffx_cont, ffy_cont;
	double nfsx_cont, nfsy_cont, ndasx, ndasy;
	double nffx_cont, nffy_cont, ndafx, ndafy;
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		Particle& pcl = pcls[p_id];
		if (pcl.pe && rigid_circle.detect_collision_with_point(
			pcl.x, pcl.y, pcl.vol, dist, norm_x, norm_y))
		{
			fs_cont = Ks_cont * dist;
			fsx_cont = fs_cont * norm_x;
			fsy_cont = fs_cont * norm_y;
			ff_cont = Kf_cont * dist;
			ffx_cont = ff_cont * norm_x;
			ffy_cont = ff_cont * norm_y;
			// reaction force by the rigid object
			rigid_circle.add_rf(pcl.x, pcl.y,
				-(fsx_cont + ffx_cont), -(fsy_cont + ffy_cont));
			// adjust velocity at nodes
			Element& e = *pcl.pe;
			// node 1
			Node& n1 = nodes[e.n1];
			nfsx_cont = pcl.N1 * fsx_cont;
			nfsy_cont = pcl.N1 * fsy_cont;
			ndasx = nfsx_cont / n1.m_s;
			ndasy = nfsy_cont / n1.m_s;
			n1.ax_s += ndasx;
			n1.ay_s += ndasy;
			n1.vx_s += ndasx * dt;
			n1.vy_s += ndasy * dt;
			nffx_cont = pcl.N1 * ffx_cont;
			nffy_cont = pcl.N1 * ffy_cont;
			ndafx = nffx_cont / n1.m_f;
			ndafy = nffy_cont / n1.m_f;
			n1.ax_f += ndafx;
			n1.ay_f += ndafy;
			n1.vx_f += ndafx * dt;
			n1.vy_f += ndafy * dt;
			// node 2
			Node& n2 = nodes[e.n2];
			nfsx_cont = pcl.N2 * fsx_cont;
			nfsy_cont = pcl.N2 * fsy_cont;
			ndasx = nfsx_cont / n2.m_s;
			ndasy = nfsy_cont / n2.m_s;
			n2.ax_s += ndasx;
			n2.ay_s += ndasy;
			n2.vx_s += ndasx * dt;
			n2.vy_s += ndasy * dt;
			nffx_cont = pcl.N2 * ffx_cont;
			nffy_cont = pcl.N2 * ffy_cont;
			ndafx = nffx_cont / n2.m_f;
			ndafy = nffy_cont / n2.m_f;
			n2.ax_f += ndafx;
			n2.ay_f += ndafy;
			n2.vx_f += ndafx * dt;
			n2.vy_f += ndafy * dt;
			// node 3
			Node& n3 = nodes[e.n3];
			nfsx_cont = pcl.N3 * fsx_cont;
			nfsy_cont = pcl.N3 * fsy_cont;
			ndasx = nfsx_cont / n3.m_s;
			ndasy = nfsy_cont / n3.m_s;
			n3.ax_s += ndasx;
			n3.ay_s += ndasy;
			n3.vx_s += ndasx * dt;
			n3.vy_s += ndasy * dt;
			nffx_cont = pcl.N3 * ffx_cont;
			nffy_cont = pcl.N3 * ffy_cont;
			ndafx = nffx_cont / n3.m_f;
			ndafy = nffy_cont / n3.m_f;
			n3.ax_f += ndafx;
			n3.ay_f += ndafy;
			n3.vx_f += ndafx * dt;
			n3.vy_f += ndafy * dt;
		}
	}

	rigid_circle.update_motion(dt);

	return 0;
}

//int Model_T2D_CHM_s::apply_rigid_circle(double dt)
//{
//	rigid_circle.reset_rf();
//
//	double dist, norm_x, norm_y;
//	double fs_cont, fsx_cont, fsy_cont;
//	double ff_cont, ffx_cont, ffy_cont;
//	double nfsx_cont, nfsy_cont, ndasx, ndasy;
//	double nffx_cont, nffy_cont, ndafx, ndafy;
//	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
//	{
//		Particle &pcl = pcls[p_id];
//		if (!pcl.pe)
//			continue;
//
//		// solid particle
//		if (rigid_circle.detect_collision_with_point(
//				pcl.x, pcl.y, pcl.vol, dist, norm_x, norm_y))
//		{
//			fs_cont = Ks_cont * dist;
//			fsx_cont = fs_cont * norm_x;
//			fsy_cont = fs_cont * norm_y;
//			// reaction force by the rigid object
//			rigid_circle.add_rf(pcl.x, pcl.y, -fsx_cont, -fsy_cont);
//			// adjust velocity at nodes
//			Element &e = *pcl.pe;
//			// node 1
//			Node &n1 = nodes[e.n1];
//			nfsx_cont = pcl.N1 * fsx_cont;
//			nfsy_cont = pcl.N1 * fsy_cont;
//			ndasx = nfsx_cont / n1.m_s;
//			ndasy = nfsy_cont / n1.m_s;
//			n1.ax_s += ndasx;
//			n1.ay_s += ndasy;
//			n1.vx_s += ndasx * dt;
//			n1.vy_s += ndasy * dt;
//			// node 2
//			Node &n2 = nodes[e.n2];
//			nfsx_cont = pcl.N2 * fsx_cont;
//			nfsy_cont = pcl.N2 * fsy_cont;
//			ndasx = nfsx_cont / n2.m_s;
//			ndasy = nfsy_cont / n2.m_s;
//			n2.ax_s += ndasx;
//			n2.ay_s += ndasy;
//			n2.vx_s += ndasx * dt;
//			n2.vy_s += ndasy * dt;
//			// node 3
//			Node &n3 = nodes[e.n3];
//			nfsx_cont = pcl.N3 * fsx_cont;
//			nfsy_cont = pcl.N3 * fsy_cont;
//			ndasx = nfsx_cont / n3.m_s;
//			ndasy = nfsy_cont / n3.m_s;
//			n3.ax_s += ndasx;
//			n3.ay_s += ndasy;
//			n3.vx_s += ndasx * dt;
//			n3.vy_s += ndasy * dt;
//		}
//
//		// fluid particle
//		if (rigid_circle.detect_collision_with_point(
//				pcl.x_f, pcl.y_f, pcl.vol, dist, norm_x, norm_y))
//		{
//			ff_cont = Kf_cont * dist;
//			ffx_cont = ff_cont * norm_x;
//			ffy_cont = ff_cont * norm_y;
//			// reaction force by the rigid object
//			rigid_circle.add_rf(pcl.x, pcl.y, -ffx_cont, -ffy_cont);
//			// adjust velocity at nodes
//			Element& e = *pcl.pe;
//			// node 1
//			Node& n1 = nodes[e.n1];
//			nffx_cont = pcl.N1 * ffx_cont;
//			nffy_cont = pcl.N1 * ffy_cont;
//			ndafx = nffx_cont / n1.m_f;
//			ndafy = nffy_cont / n1.m_f;
//			n1.ax_f += ndafx;
//			n1.ay_f += ndafy;
//			n1.vx_f += ndafx * dt;
//			n1.vy_f += ndafy * dt;
//			// node 2
//			Node& n2 = nodes[e.n2];
//			nffx_cont = pcl.N2 * ffx_cont;
//			nffy_cont = pcl.N2 * ffy_cont;
//			ndafx = nffx_cont / n2.m_f;
//			ndafy = nffy_cont / n2.m_f;
//			n2.ax_f += ndafx;
//			n2.ay_f += ndafy;
//			n2.vx_f += ndafx * dt;
//			n2.vy_f += ndafy * dt;
//			// node 3
//			Node& n3 = nodes[e.n3];
//			nffx_cont = pcl.N3 * ffx_cont;
//			nffy_cont = pcl.N3 * ffy_cont;
//			ndafx = nffx_cont / n3.m_f;
//			ndafy = nffy_cont / n3.m_f;
//			n3.ax_f += ndafx;
//			n3.ay_f += ndafy;
//			n3.vx_f += ndafx * dt;
//			n3.vy_f += ndafy * dt;
//		}
//	}
//
//	rigid_circle.update_motion(dt);
//
//	return 0;
//}

#include <fstream>

void Model_T2D_CHM_s::sum_vol_for_each_elements()
{
	// init elements
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		Element &e = elems[e_id];
		e.pcls = nullptr;
		e.pcl_vol = 0.0;
	}

	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		Particle &pcl = pcls[p_id];
		if (!(pcl.pe = find_in_which_element(pcl)))
			continue;
		pcl.pe->add_pcl(pcl);

		pcl.vol_s = pcl.m_s / pcl.density_s;
		pcl.vol = pcl.vol_s / (1.0 - pcl.n);
		pcl.pe->pcl_vol += pcl.vol;
	}

	std::fstream out_file;
	out_file.open("area.txt", std::ios::binary | std::ios::out);
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		Element &e = elems[e_id];
		if (e.pcls)
		{
			if (e.pcl_vol > e.area)
			{
				out_file << "id: " << e.id << ", elem_a: " << e.area
					<< ", pcl_a:" << e.pcl_vol 
					<< " * +" << (e.pcl_vol - e.area) / e.area * 100.0 << "% *\n";
			}
			else if (e.pcl_vol < e.area)
			{
				out_file << "id: " << e.id << ", elem_a: " << e.area
					<< ", pcl_a:" << e.pcl_vol
					<< " * -" << (e.area - e.pcl_vol) / e.area * 100.0 <<"% *\n";
			}
			else
			{
				out_file << "id: " << e.id << ", elem_a: " << e.area
					<< ", pcl_a:" << e.pcl_vol << "\n";
			}
		}
	}
	out_file.close();
}
