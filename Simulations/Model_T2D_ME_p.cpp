#include "Simulations_pcp.h"

#include <set>
#include <exception>

#include "Geometry.h"
#include "Model_T2D_ME_p.h"

Model_T2D_ME_p::Model_T2D_ME_p() :
	Model("Model_T2D_ME_p"),
	node2elems(nullptr),
	pcls(nullptr), pcl_num(0),
	bfx_num(0), bfxs(nullptr),
	bfy_num(0), bfys(nullptr),
	tx_num(0), txs(nullptr),
	ty_num(0), tys(nullptr),
	ax_num(0), axs(nullptr),
	ay_num(0), ays(nullptr),
	vx_num(0), vxs(nullptr),
	vy_num(0), vys(nullptr),
	rigid_circle_is_init(false) {}

Model_T2D_ME_p::~Model_T2D_ME_p()
{
	clear_n2e_info();

	clear_pcls();

	clear_bfxs();
	clear_bfys();
	clear_txs();
	clear_tys();
	clear_axs();
	clear_ays();
	clear_vxs();
	clear_vys();
}

void Model_T2D_ME_p::init_mesh_shape_funcs()
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

void Model_T2D_ME_p::alloc_n2e_info(size_t num)
{
	clear_n2e_info();
	node2elems = new NodeToElem[num];
}

void Model_T2D_ME_p::clear_n2e_info()
{
	if (node2elems)
	{
		delete[] node2elems;
		node2elems = nullptr;
	}
}

void Model_T2D_ME_p::init_node_to_elem_info()
{
	struct N2E
	{
		Element* pe;
		unsigned char n_id;
		N2E* next;
	};

	struct N2EStack
	{
		N2E *top;
		size_t num;
		inline void init()
		{
			top = nullptr;
			num = 0;
		}
		inline void add(N2E &n2e)
		{
			n2e.next = top;
			top = &n2e;
			++num;
		}
	};

	N2EStack *n2e_stacks = new N2EStack[node_num];
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		N2EStack& n2es = n2e_stacks[n_id];
		n2es.init();
	}

	N2E* n2es = new N2E[elem_num * 3];
	size_t cur_n2e_id = 0;
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		Element& e = elems[e_id];
		// node 1
		N2E &n2e1 = n2es[cur_n2e_id];
		n2e1.pe = &e;
		n2e1.n_id = 0;
		N2EStack &ns1 = n2e_stacks[e.n1];
		ns1.add(n2e1);
		++cur_n2e_id;
		// node 2
		N2E& n2e2 = n2es[cur_n2e_id];
		n2e2.pe = &e;
		n2e2.n_id = 1;
		N2EStack &ns2 = n2e_stacks[e.n2];
		ns2.add(n2e2);
		++cur_n2e_id;
		// node 3
		N2E& n2e3 = n2es[cur_n2e_id];
		n2e3.pe = &e;
		n2e3.n_id = 2;
		N2EStack &ns3 = n2e_stacks[e.n3];
		ns3.add(n2e3);
		++cur_n2e_id;
	}

	alloc_n2e_info(elem_num * 3);
	NodeToElem *cur_node2elems = node2elems;
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node& n = nodes[n_id];
		N2EStack &ns = n2e_stacks[n_id];
		n.n2e_num = ns.num;
		n.n2es = cur_node2elems;
		size_t i = 0;
		for (N2E *n2e = ns.top; n2e; n2e = n2e->next)
		{
			n.n2es[i].e_id = n2e->pe->id;
			n.n2es[i].n_id = n2e->n_id;
			++i;
		}
		cur_node2elems += n.n2e_num;
	}

	delete[] n2e_stacks;
	delete[] n2es;
}

namespace
{
	void swap_abc(AccelerationBC& abc1, AccelerationBC& abc2)
	{
		size_t n_id_tmp;
		n_id_tmp = abc1.node_id;
		abc1.node_id = abc2.node_id;
		abc2.node_id = abc1.node_id;

		double a_tmp;
		a_tmp = abc1.a;
		abc1.a = abc2.a;
		abc2.a = a_tmp;
	}
}

void Model_T2D_ME_p::init_bcs()
{
	// traction bc and body force
	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		Particle& pcl = pcls[p_id];
		pcl.has_fx_ext = false;
		pcl.has_fy_ext = false;
		pcl.fx_ext = 0.0;
		pcl.fy_ext = 0.0;
	}

	for (size_t t_id = 0; t_id < tx_num; ++t_id)
	{
		TractionBCAtPcl& tbc = txs[t_id];
		Particle &pcl = pcls[tbc.pcl_id];
		pcl.has_fx_ext = true;
		pcl.fx_ext += tbc.t;
	}

	for (size_t t_id = 0; t_id < ty_num; ++t_id)
	{
		TractionBCAtPcl& tbc = tys[t_id];
		Particle& pcl = pcls[tbc.pcl_id];
		pcl.has_fy_ext = true;
		pcl.fy_ext += tbc.t;
	}

	for (size_t bf_id = 0; bf_id < bfx_num; ++bf_id)
	{
		BodyForceAtPcl &bf = bfxs[bf_id];
		Particle& pcl = pcls[bf.pcl_id];
		pcl.has_fx_ext = true;
		pcl.fx_ext += pcl.m * bf.bf;
	}

	for (size_t bf_id = 0; bf_id < bfy_num; ++bf_id)
	{
		BodyForceAtPcl &bf = bfys[bf_id];
		Particle& pcl = pcls[bf.pcl_id];
		pcl.has_fy_ext = true;
		pcl.fy_ext += pcl.m * bf.bf;
	}

	// acceleration and velocity bcs
	std::pair<std::set<size_t>::iterator, bool> res;
	char error_str[100];

	std::set<size_t> ax_n_id_set;
	for (size_t a_id = 0; a_id < ax_num; ++a_id)
	{
		AccelerationBC& abc = axs[a_id];
		if (abc.node_id >= node_num)
		{
			snprintf(error_str, 100,
				"ax bc %zu has node id %zu out of range.\n",
				a_id, abc.node_id);
			throw std::exception(error_str);
		}
		res = ax_n_id_set.insert(abc.node_id);
		if (!res.second)
		{
			snprintf(error_str, 100,
				"Already axist ax bc %zu at node %zu.\n",
				a_id, abc.node_id);
			throw std::exception(error_str);
		}
	}

	std::set<size_t> ay_n_id_set;
	for (size_t a_id = 0; a_id < ay_num; ++a_id)
	{
		AccelerationBC& abc = ays[a_id];
		if (abc.node_id >= node_num)
		{
			snprintf(error_str, 100,
				"ay bc %zu has node id %zu out of range.\n",
				a_id, abc.node_id);
			throw std::exception(error_str);
		}
		res = ay_n_id_set.insert(abc.node_id);
		if (!res.second)
		{
			snprintf(error_str, 100,
				"Already axist ay bc %zu at node %zu.\n",
				a_id, abc.node_id);
			throw std::exception(error_str);
		}
	}

	std::set<size_t> vx_n_id_set;
	for (size_t v_id = 0; v_id < vx_num; ++v_id)
	{
		VelocityBC& vbc = vxs[v_id];
		if (vbc.node_id >= node_num)
		{
			snprintf(error_str, 100,
				"vx bc %zu has node id %zu out of range.\n",
				v_id, vbc.node_id);
			throw std::exception(error_str);
		}
		res = vx_n_id_set.insert(vbc.node_id);
		if (!res.second)
		{
			snprintf(error_str, 100,
				"Already axist vx bc %zu at node %zu.\n",
				v_id, vbc.node_id);
			throw std::exception(error_str);
		}
	}

	std::set<size_t> vy_n_id_set;
	for (size_t v_id = 0; v_id < vy_num; ++v_id)
	{
		VelocityBC& vbc = vys[v_id];
		if (vbc.node_id >= node_num)
		{
			snprintf(error_str, 100,
				"vy bc %zu has node id %zu out of range.\n",
				v_id, vbc.node_id);
			throw std::exception(error_str);
		}
		res = vy_n_id_set.insert(vbc.node_id);
		if (!res.second)
		{
			snprintf(error_str, 100,
				"Already axist vy bc %zu at node %zu.\n",
				v_id, vbc.node_id);
			throw std::exception(error_str);
		}
	}

	// ax and vx, ay and vy bcs should be mutually exclusive
	for (size_t a_id = 0; a_id < ax_num;)
	{
		AccelerationBC& abc = axs[a_id];
		if (vx_n_id_set.find(abc.node_id) != vx_n_id_set.end())
		{
			--ax_num;
			swap_abc(abc, axs[ax_num]);
		}
		else
		{
			++a_id;
		}
	}

	for (size_t a_id = 0; a_id < ay_num;)
	{
		AccelerationBC& abc = ays[a_id];
		if (vy_n_id_set.find(abc.node_id) != vy_n_id_set.end())
		{
			--ay_num;
			swap_abc(abc, ays[ay_num]);
		}
		else
		{
			++a_id;
		}
	}
}

int Model_T2D_ME_p::init_mesh(
	double *node_coords, size_t n_num,
	size_t *elem_n_ids,  size_t e_num
	)
{
	int res = BgMesh::init_mesh(node_coords, n_num, elem_n_ids, e_num);
	if (res < 0)
		return res;
	init_mesh_shape_funcs();
	init_node_to_elem_info();
	return 0;
}

int Model_T2D_ME_p::load_mesh_from_hdf5(const char* file_name)
{
	int res = BgMesh::load_mesh_from_hdf5(file_name);
	if (res < 0)
		return res;
	init_mesh_shape_funcs();
	init_node_to_elem_info();
	return 0;
}

int Model_T2D_ME_p::init_search_grid(double _hx, double _hy)
{
	return search_bg_grid.init(*this, _hx, _hy);
}

int Model_T2D_ME_p::init_pcls(size_t num, double m, double density)
{
	clear_pcls();

	if (num == 0)
		return -1;
	alloc_pcls(num);
	for (size_t pcl_id = 0; pcl_id < num; ++pcl_id)
	{
		Particle &pcl = pcls[pcl_id];
		pcl.id = pcl_id;
		pcl.ux = 0.0;
		pcl.uy = 0.0;
		pcl.vx = 0.0;
		pcl.vy = 0.0;
		pcl.m = m;
		pcl.density = density;
		pcl.e11 = 0.0;
		pcl.e22 = 0.0;
		pcl.e12 = 0.0;
		pcl.s11 = 0.0;
		pcl.s12 = 0.0;
		pcl.s22 = 0.0;
	}
	return 0;
}

int Model_T2D_ME_p::init_pcls(
	ParticleGenerator2D<Model_T2D_ME_p>& pg,
	double density
	)
{
	int res = init_pcls(pg.get_num(), density, density);
	if (res)
		return res;

	typedef ParticleGenerator2D<Model_T2D_ME_p>::Particle PgPcl;
	PgPcl* pg_pcl = pg.first();
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		Particle& pcl = pcls[pcl_id];
		pcl.x = pg_pcl->x;
		pcl.y = pg_pcl->y;
		pcl.m *= pg_pcl->area;
		pg_pcl = pg.next(pg_pcl);
	}

	return 0;
}


void Model_T2D_ME_p::alloc_pcls(size_t num)
{
	pcls = new Particle[num];
	pcl_num = num;
}

void Model_T2D_ME_p::clear_pcls()
{
	if (pcls)
	{
		delete[] pcls;
		pcls = nullptr;
	}
	pcl_num = 0;
}

#include <fstream>

void Model_T2D_ME_p::sum_vol_for_all_elements()
{
	struct ElemTmpData
	{
		Particle *pcls;
		double pcl_vol;
		inline void add_pcl(Particle &pcl)
		{
			pcl.next_pcl_in_elem = pcls;
			pcls = &pcl;
		}
	};

	// init
	ElemTmpData* etds = new ElemTmpData[elem_num];
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		ElemTmpData &etd = etds[e_id];
		etd.pcls = nullptr;
		etd.pcl_vol = 0.0;
	}

	for (size_t p_id = 0; p_id < pcl_num; ++p_id)
	{
		Particle& pcl = pcls[p_id];
		if (!(pcl.pe = find_in_which_element(pcl)))
			continue;
		ElemTmpData& etd = etds[pcl.pe->id];
		etd.add_pcl(pcl);
		pcl.vol = pcl.m / pcl.density;
		etd.pcl_vol += pcl.vol;
	}

	// write to file
	std::fstream out_file;
	out_file.open("area.txt", std::ios::binary | std::ios::out);
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		Element& e = elems[e_id];
		ElemTmpData& etd = etds[e_id];
		if (etd.pcls)
		{
			if (etd.pcl_vol > e.area)
			{
				out_file << "id: " << e.id << ", elem_a: " << e.area
					<< ", pcl_a:" << etd.pcl_vol
					<< " * +" << (etd.pcl_vol - e.area) / e.area * 100.0 << "% *\n";
			}
			else if (e.mi_pcl_vol < e.area)
			{
				out_file << "id: " << e.id << ", elem_a: " << e.area
					<< ", pcl_a:" << etd.pcl_vol
					<< " * -" << (e.area - etd.pcl_vol) / e.area * 100.0 << "% *\n";
			}
			else
			{
				out_file << "id: " << e.id << ", elem_a: " << e.area
					<< ", pcl_a:" << etd.pcl_vol << "\n";
			}
		}
	}

	// complete
	out_file.close();
	delete[] etds;
}
