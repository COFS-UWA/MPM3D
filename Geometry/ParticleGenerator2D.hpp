#ifndef __Particle_Generator_2D_hpp__
#define __Particle_Generator_2D_hpp__

#include "ItemArray.hpp"
#include "ItemBuffer.hpp"
#include "LinkList.hpp"
#include "Geometry.h"
#include "PDSBgGrid2D.hpp"
#include "RandNumGenerator.h"
#include "RandomPointQueue.h"

template <typename TriangleMesh>
class ParticleGenerator2D
{
public:
	typedef typename TriangleMesh::Edge Edge;
	typedef typename TriangleMesh::Node Node;
	typedef typename TriangleMesh::Element Element;

	struct Particle
	{
		friend ParticleGenerator2D;
	public:
		double x, y, area;
		typename TriangleMesh::Element* pe;
	protected:
		LinkListPointer<Particle> pointer;
		Particle *next2;
		Particle* next3; // for pds bg grid
	};

protected:
	typedef MemoryUtils::ItemBuffer<Particle> ParticleBuffer;
	typedef LinkList<Particle, offsetof(Particle, pointer)> ParticleList;
	ParticleBuffer pcl_buf;
	ParticleList pcl_list;

public:
	inline Particle* first() noexcept { return pcl_list.first(); }
	inline Particle* next(Particle* cur) noexcept { return pcl_list.next(cur); }
	inline bool is_end(Particle* cur) noexcept { return pcl_list.is_end(cur); }
	inline bool is_not_end(Particle* cur) noexcept { return pcl_list.is_not_end(cur); }
	inline bool is_empty() noexcept { return pcl_list.is_empty(); }
	inline size_t get_num() noexcept { return pcl_list.get_num(); }

public:
	ParticleGenerator2D()
	{
		pg_param.line_div_num = 1;
		pg_param.area_div_num = 1;
	}
	~ParticleGenerator2D() { clear(); }

	inline void clear()
	{
		pcl_buf.clear();
		pcl_list.reset();
	}

	// single particle operation
	inline Particle& add_pcl(Particle& pcl)
	{
		Particle* res = pcl_buf.alloc();
		res->x = pcl.x;
		res->y = pcl.y;
		res->area = pcl.area;
		pcl_list.append(res);
		return *res;
	}

	inline void del_pcl(Particle& pcl)
	{
		pcl_list.del(pcl);
		pcl_buf.del(&pcl);
	}

// ============ generate particle in grid layout ============ 
	void clear_pcls_in_rect(Rect &range);
	void generate_pcls_in_grid_layout(Rect &range, double pcl_dx, double pcl_dy);
	void replace_with_pcls_in_grid_layout(Rect &range, double pcl_dx, double pcl_dy);

// ============ generate particle at gauss points ============ 
protected:
	// parameters
	TriangleMesh* mesh;
	struct
	{
		Point2D p1, p2, p3;
		double area;
		size_t line_div_num; // evenly distributed generator
		size_t area_div_num; // randomly distributed generator
	} pg_param;
	Particle pcl_res;

	// particle generators
	void FirstOrderGaussPoint2DGenerator();
	void SecondOrderGaussPoint2DGenerator();

	typedef void (ParticleGenerator2D::* GeneratorFunc)();
	void generate_pcls(GeneratorFunc cur_generator_func);

public:
	inline void generate_pcls_at_1st_gauss(TriangleMesh& _mesh)
	{
		mesh = &_mesh;
		generate_pcls(&ParticleGenerator2D::FirstOrderGaussPoint2DGenerator);
	}
	inline void generate_pcls_at_2nd_gauss(TriangleMesh& _mesh)
	{
		mesh = &_mesh;
		generate_pcls(&ParticleGenerator2D::SecondOrderGaussPoint2DGenerator);
	}

// ============ generate randomly distributed particle ============ 
// using poisson disk sampling method
protected:
	double g_pt_dist; // global point distance
	struct PDSNodeInfo
	{
		size_t id;
		double l_pt_dist; // local point distance
	};
	size_t pds_node_info_num;
	PDSNodeInfo *pds_node_infos;
	double min_pt_dist;
	void init_node_info();
	double get_l_pt_size_at_pcl(Particle& pcl);
	Particle gen_rand_point_around(Particle& pcl, double dist);

public:
	void generate_randomly_distributed_pcls(TriangleMesh& _mesh, double g_dist);

// ============ adjust particle size to fit elements ===========
protected: // helper functions and data structures
	struct ElemInfo
	{
		size_t id;
		double pcl_area;
		size_t pcl_num;
		Particle* pcls;
		inline void init()
		{
			pcl_area = 0.0;
			pcls = nullptr;
			pcl_num = 0;
		}
		inline void add_pcl(Particle* pcl)
		{
			pcl_area += pcl->area;
			pcl->next2 = pcls;
			pcls = pcl;
			++pcl_num;
		}
	};

public:
	int adjust_pcl_size_to_fit_elems(TriangleMesh &mesh);
};

template <typename TriangleMesh>
void ParticleGenerator2D<TriangleMesh>::clear_pcls_in_rect(Rect &range)
{
	for (Particle* p_iter = first();
		 is_not_end(p_iter); p_iter = next(p_iter))
	{
		Particle& pcl = *p_iter;
		if (pcl.x >= range.xl && pcl.x <= range.xu &&
			pcl.y >= range.yl && pcl.y <= range.yu)
			del_pcl(pcl);
	}
}

template <typename TriangleMesh>
void ParticleGenerator2D<TriangleMesh>::generate_pcls_in_grid_layout(Rect& range, double pcl_dx, double pcl_dy)
{
	double x_len, y_len;
	double start_x, start_y;
	size_t x_num, y_num;
	// x direction
	x_len = range.xu - range.xl;
	x_num = size_t(ceil(x_len / pcl_dx));
	pcl_dx = x_len / double(x_num);
	start_x = range.xl + pcl_dx * 0.5f;
	// y direction
	y_len = range.yu - range.yl;
	y_num = size_t(ceil(y_len / pcl_dy));
	pcl_dy = y_len / double(y_num);
	start_y = range.yl + pcl_dy * 0.5f;
	// generate particle
	Particle pcl;
	pcl.area = pcl_dx * pcl_dy;
	pcl_buf.set_page_size(x_num * y_num);
		for (size_t y_id = 0; y_id < y_num; ++y_id)
			for (size_t x_id = 0; x_id < x_num; ++x_id)
			{
				pcl.x = start_x + pcl_dx * double(x_id);
				pcl.y = start_y + pcl_dy * double(y_id);
				add_pcl(pcl);
			}
}

template <typename TriangleMesh>
void ParticleGenerator2D<TriangleMesh>::replace_with_pcls_in_grid_layout(
	Rect& range, double pcl_dx, double pcl_dy)
{
	clear_pcls_in_rect(range);
	generate_pcls_in_grid_layout(range, pcl_dx, pcl_dy);
}

/* =========== Generate particles from tetrahedron mesh =========== */
template <typename TriangleMesh>
void ParticleGenerator2D<TriangleMesh>::FirstOrderGaussPoint2DGenerator()
{
	Point2D& p1 = pg_param.p1;
	Point2D& p2 = pg_param.p2;
	Point2D& p3 = pg_param.p3;
	pcl_res.x = (p1.x + p2.x + p3.x) * 0.33333333333;
	pcl_res.y = (p1.y + p2.y + p3.y) * 0.33333333333;
	pcl_res.area = pg_param.area;
	add_pcl(pcl_res);
}

template <typename TriangleMesh>
void ParticleGenerator2D<TriangleMesh>::SecondOrderGaussPoint2DGenerator()
{
#define alpha (2.0/3.0)
#define belta (1.0/6.0)
	pcl_res.area = pg_param.area * 0.3333333333;
	Point2D& p1 = pg_param.p1;
	Point2D& p2 = pg_param.p2;
	Point2D& p3 = pg_param.p3;
	// pcl 1
	pcl_res.x = alpha * p1.x + belta * p2.x + belta * p3.x;
	pcl_res.y = alpha * p1.y + belta * p2.y + belta * p3.y;
	add_pcl(pcl_res);
	// pcl 2
	pcl_res.x = belta * p1.x + alpha * p2.x + belta * p3.x;
	pcl_res.y = belta * p1.y + alpha * p2.y + belta * p3.y;
	add_pcl(pcl_res);
	// pcl 3
	pcl_res.x = belta * p1.x + belta * p2.x + alpha * p3.x;
	pcl_res.y = belta * p1.y + belta * p2.y + alpha * p3.y;
	add_pcl(pcl_res);
#undef alpha
#undef belta
}

/* =========== Generate particles in grid =========== */
template <typename TriangleMesh>
void ParticleGenerator2D<TriangleMesh>::generate_pcls(GeneratorFunc cur_generator_func)
{
	size_t elem_num = mesh->get_elem_num();
	typename TriangleMesh::Element* elems = mesh->get_elems();
	typename TriangleMesh::Node* nodes = mesh->get_nodes();
	pcl_buf.set_page_size(elem_num);

	for (size_t elem_id = 0; elem_id < elem_num; ++elem_id)
	{
		typename TriangleMesh::Element& elem = elems[elem_id];
		typename TriangleMesh::Node& n1 = nodes[elem.n1];
		typename TriangleMesh::Node& n2 = nodes[elem.n2];
		typename TriangleMesh::Node& n3 = nodes[elem.n3];
		Point2D& p1 = pg_param.p1;
		p1.x = n1.x;
		p1.y = n1.y;
		Point2D& p2 = pg_param.p2;
		p2.x = n2.x;
		p2.y = n2.y;
		Point2D& p3 = pg_param.p3;
		p3.x = n3.x;
		p3.y = n3.y;
		pg_param.area = elem.area;
		(this->*cur_generator_func)();
	}
}

template <typename TriangleMesh>
int ParticleGenerator2D<TriangleMesh>::
	adjust_pcl_size_to_fit_elems(TriangleMesh& mesh)
{
	typedef typename TriangleMesh::Element Element;

	// init element info
	size_t elem_num = mesh.get_elem_num();
	if (!elem_num)
		return -1;
	ElemInfo *elem_infos = new ElemInfo[elem_num];
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		ElemInfo& ei = elem_infos[e_id];
		ei.id = e_id;
		ei.init();
	}

	Element* pe;
	Particle* pcl_tmp;
	Particle *pcl = first();
	size_t pid = 0;
	while (is_not_end(pcl))
	{
		++pid;
		if (pid == 871)
			pe = mesh.find_in_which_element_bf(*pcl);
		pe = mesh.find_in_which_element(*pcl);
		if (!pe)
		{
			pcl_tmp = pcl;
			pcl = next(pcl);
			del_pcl(*pcl_tmp);
			continue;
		}
		elem_infos[pe->id].add_pcl(pcl);
		pcl = next(pcl);
	}

	Element* elems = mesh.get_elems();
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		Element& e = elems[e_id];
		ElemInfo& ei = elem_infos[e_id];

		if (ei.pcl_num == 0)
			continue;

		double area_ratio = e.area / ei.pcl_area;
		for (Particle* pcl = ei.pcls; pcl; pcl = pcl->next2)
			pcl->area *= area_ratio;
	}

	delete[] elem_infos;
	
	return 0;
}


template <typename TriangleMesh>
void ParticleGenerator2D<TriangleMesh>::init_node_info()
{
	pds_node_info_num = mesh->get_node_num();
	pds_node_infos = new PDSNodeInfo[pds_node_info_num];
	for (size_t n_id = 0; n_id < pds_node_info_num; ++n_id)
	{
		PDSNodeInfo& ni = pds_node_infos[n_id];
		ni.id = n_id;
		ni.l_pt_size = std::numeric_limits<double>::max();
	}

	size_t edge_num = mesh->get_edge_num();
	Edge* edges = mesh->get_edges();
	Node* nodes = mesh->get_nodes();
	double edge_len;
	min_pt_dist = std::numeric_limits<double>::max();
	for (size_t e_id = 0; e_id < edge_num; ++e_id)
	{
		Edge &e = edges[e_id];
		Node& n1 = nodes[e.n1];
		Node& n2 = nodes[e.n2];
		edge_len = cal_distance_2D(n1, n2);
		PDSNodeInfo& ni1 = pds_node_infos[e.n1];
		if (ni1.l_pt_size > edge_len)
			ni1.l_pt_size = edge_len;
		PDSNodeInfo& ni2 = pds_node_infos[e.n2];
		if (ni2.l_pt_size > edge_len)
			ni2.l_pt_size = edge_len;
		if (min_pt_dist > edge_len)
			min_pt_dist = edge_len;
	}

	for (size_t n_id = 0; n_id < pds_node_info_num; ++n_id)
	{
		PDSNodeInfo& ni = pds_node_infos[n_id];
		ni.l_pt_size *= g_pt_dist;
	}
	min_pt_size *= g_pt_dist;
}

template <typename TriangleMesh>
double ParticleGenerator2D<TriangleMesh>::get_l_pt_size_at_pcl(Particle &pcl)
{
	Element& e = *pcl.pe;
	double N1 = mesh->cal_N1(e, pcl.x, pcl.y);
	double N2 = mesh->cal_N2(e, pcl.x, pcl.y);
	double N3 = mesh->cal_N3(e, pcl.x, pcl.y);
	PDSNodeInfo& ni1 = nodes[e.n1];
	PDSNodeInfo& ni2 = nodes[e.n2];
	PDSNodeInfo& ni3 = nodes[e.n3];
	return N1 * ni1.l_pt_size + N2 * ni2.l_pt_size + N3 * ni3.l_pt_size;
}

template <typename TriangleMesh>
ParticleGenerator2D<TriangleMesh>::Particle
	ParticleGenerator2D<TriangleMesh>::gen_rand_point_around(Particle& pcl, double dist)
{
	double radius = dist * (1.0 + RandNumGenerator::get_double());
	double angle = 2.0 * 3.14159265359 * RandNumGenerator::get_double();
	Particle res;
	res.x = pcl.x + radius * cos(angle);
	res.y = pcl.y + radius * sin(angle);
	return res;
}

#define NEW_POINTS_COUNT 30
template <typename TriangleMesh>
void ParticleGenerator2D<TriangleMesh>
	::generate_randomly_distributed_pcls(
	TriangleMesh& _mesh,
	double g_dist
	)
{
	mesh = &_mesh;
	g_pt_dist = g_dist;
	init_node_info();
	
	double cell_size = min_pt_size / sqrt(2.0);
	double bg_grid_xl = mesh->get_xl();
	double bg_grid_xu = mesh->get_xu();
	double bg_grid_yl = mesh->get_yl();
	double bg_grid_yu = mesh->get_yu();
	size_t bg_grid_x_num = size_t(ceil((bg_grid_xu - bg_grid_xl) / cell_size));
	size_t bg_grid_y_num = size_t(ceil((bg_grid_yu - bg_grid_yl) / cell_size));
	PDSBgGrid2D<Particle, offsetof(Particle, next3)> bg_grid;
	bg_grid.init(
		bg_grid_xl,
		bg_grid_xu,
		bg_grid_yl,
		bg_grid_yu,
		bg_grid_x_num,
		bg_grid_y_num
		);

	RandomPointQueue random_point_queue;

	// generate the first point
	Particle pcl_tmp;
	double pcl_min_dist;
	Element* pe;
	do
	{
		pcl_tmp.x = RandNumGenerator::get_double(bg_grid_xl, bg_grid_xu);
		pcl_tmp.y = RandNumGenerator::get_double(bg_grid_yl, bg_grid_yu);
		pe = mesh->find_in_which_element(pcl_tmp);
	} while (!pe);
	pcl_tmp.pe = pe;
	pcl_tmp.area = get_l_pt_size_at_pcl(pcl_tmp);
	random_point_queue.add_point(pcl_tmp);
	add_pcl(pcl_tmp);
	bg_grid.add_point(pcl_tmp);
	
	Particle cur_pcl;
	while (random_point_queue.get_point(cur_pcl))
	{
		for (size_t i = 0; i < NEW_POINTS_COUNT; ++i)
		{
			pcl_tmp = gen_rand_point_around(cur_pcl, pcl_min_dist.area);
			if (!(pe = mesh->find_in_which_element(pcl_tmp))
				continue;

			pcl_tmp.pe = pe;
			pcl_tmp.area = get_l_pt_size_at_pcl(pcl_tmp);
			if (!bg_grid.has_point_nearby(pcl_tmp, pcl_tmp.area))
			{
				random_point_queue.add_point(pcl_tmp);
				add_pcl(pcl_tmp);
				bg_grid.add_point(pcl_tmp);
			}
		}
	}

	double pcl_len;
	for (Particle *piter = first(); is_not_end(piter); piter = next(piter))
	{
		pcl_len = piter->area;
		piter->area = pcl_len;
	}
	adjust_pcl_size_to_fit_elems(); // adjust pcl size

	delete[] pds_nodes_infos;
}
#undef NEW_POINTS_COUNT

#endif