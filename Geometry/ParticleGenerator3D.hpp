#ifndef __Particle_Generator_3D_h__
#define __Particle_Generator_3D_h__

#include "ItemArray.hpp"
#include "ItemBuffer.hpp"
#include "LinkList.hpp"
#include "Geometry3D.h"

template <typename TetrahedronMesh>
class ParticleGenerator3D
{
public:
	struct Particle
	{
	public:
		double x, y, z, vol;
	protected:
		friend ParticleGenerator3D;
		LinkListPointer<Particle> pointer;
	};

protected:
	typedef MemoryUtils::ItemBuffer<Particle> ParticleBuffer;
	typedef LinkList<Particle, offsetof(Particle, pointer)> ParticleList;
	ParticleBuffer pcl_buf;
	ParticleList pcl_list;
	
public:
	inline Particle *first() noexcept { return pcl_list.first(); }
	inline Particle *next(Particle *cur) noexcept { return pcl_list.next(cur); }
	inline bool is_end(Particle *cur) noexcept { return pcl_list.is_end(cur); }
	inline bool is_not_end(Particle *cur) noexcept { return pcl_list.is_not_end(cur); }
	inline bool is_empty() noexcept { return pcl_list.is_empty(); }
	inline size_t get_num() noexcept { return pcl_list.get_num(); }

public:
	ParticleGenerator3D()
	{
		pg_param.line_div_num = 1;
		pg_param.vol_div_num = 1;
	}
	~ParticleGenerator3D() { clear(); }

	inline void clear()
	{
		pcl_buf.clear();
		pcl_list.reset();
	}

	// single particle operation
	inline Particle &add_pcl(Particle &pcl)
	{
		Particle *res = pcl_buf.alloc();
		res->x = pcl.x;
		res->y = pcl.y;
		res->z = pcl.z;
		res->vol = pcl.vol;
		pcl_list.append(res);
		return *res;
	}

	inline void del_pcl(Particle &pcl)
	{
		pcl_list.del(pcl);
		pcl_buf.del(&pcl);
	}

	// generate particle in grid layout
	void clear_pcls_in_cube(Cube &range);
	void generate_pcls_grid(Cube &range, double pcl_dx, double pcl_dy, double pcl_dz);
	void replace_with_pcls_grid(Cube &range, double pcl_dx, double pcl_dy, double pcl_dz);
	
protected: // generate particle from tetrahedron mesh
	// parameters
	TetrahedronMesh *mesh;
	struct
	{
		Point3D p1, p2, p3, p4;
		double vol;
		size_t line_div_num; // evenly distributed generator
		size_t vol_div_num; // randomly distributed generator
	} pg_param;
	Particle pcl_res;

	// particle generators
	void FirstOrderGaussPoint3DGenerator();
	void SecondOrderGaussPoint3DGenerator();
	void RandomlyDistributedPoint3DGenerator();

	typedef void (ParticleGenerator3D::*GeneratorFunc)();
	void generate_pcls(GeneratorFunc cur_generator_func);

public:
	inline void generate_pcls_first_order_gauss(TetrahedronMesh &_mesh)
	{
		mesh = &_mesh;
		generate_pcls(&ParticleGenerator3D::FirstOrderGaussPoint3DGenerator);
	}
	inline void generate_pcls_second_order_gauss(TetrahedronMesh &_mesh)
	{
		mesh = &_mesh;
		generate_pcls(&ParticleGenerator3D::SecondOrderGaussPoint3DGenerator);
	}
	inline void generate_pcls_random(TetrahedronMesh &_mesh, size_t _vol_div_num)
	{
		mesh = &_mesh;
		pg_param.vol_div_num = _vol_div_num;
		generate_pcls(&ParticleGenerator3D::RandomlyDistributedPoint3DGenerator);
	}
};

template <typename TetrahedronMesh>
void ParticleGenerator3D<TetrahedronMesh>::clear_pcls_in_cube(Cube &range)
{
	for (Particle *p_iter = first(); is_not_end(p_iter); p_iter = next(p_iter))
	{
		Particle &pcl = *p_iter;
		if (pcl.x >= range.xl && pcl.x <= range.xu &&
			pcl.y >= range.yl && pcl.y <= range.yu &&
			pcl.z >= range.zl && pcl.z <= range.zu)
			del_pcl(pcl);
	}
}

template <typename TetrahedronMesh>
void ParticleGenerator3D<TetrahedronMesh>::generate_pcls_grid(
	Cube &range,
	double pcl_dx,
	double pcl_dy,
	double pcl_dz
	)
{
	double x_len, y_len, z_len;
	double start_x, start_y, start_z;
	size_t x_num, y_num, z_num;
	// x direction
	x_len = range.xu - range.xl;
	x_num = size_t(ceil(x_len / pcl_dx));
	pcl_dx = x_len / double(x_num);
	start_x = range.xl + pcl_dx * 0.5;
	// y direction
	y_len = range.yu - range.yl;
	y_num = size_t(ceil(y_len / pcl_dy));
	pcl_dy = y_len / double(y_num);
	start_y = range.yl + pcl_dy * 0.5;
	// z direction
	z_len = range.zu - range.zl;
	z_num = size_t(ceil(z_len / pcl_dz));
	pcl_dz = z_len / double(z_num);
	start_z = range.zl + pcl_dz * 0.5;
	// generate particle
	Particle pcl;
	pcl.vol = pcl_dx * pcl_dy * pcl_dz;
	pcl_buf.set_page_size(x_num * y_num * z_num);
	for (size_t z_id = 0; z_id < z_num; ++z_id)
		for (size_t y_id = 0; y_id < y_num; ++y_id)
			for (size_t x_id = 0; x_id < x_num; ++x_id)
			{
				pcl.x = start_x + pcl_dx * double(x_id);
				pcl.y = start_y + pcl_dy * double(y_id);
				pcl.z = start_z + pcl_dz * double(z_id);
				add_pcl(pcl);
			}
}

template <typename TetrahedronMesh>
void ParticleGenerator3D<TetrahedronMesh>::replace_with_pcls_grid(
	Cube &range,
	double pcl_dx,
	double pcl_dy,
	double pcl_dz
	)
{
	clear_pcls_in_cube(range);
	generate_pcls_grid(range, pcl_dx, pcl_dy, pcl_dz);
}

/* =========== Generate particles from tetrahedron mesh =========== */
template <typename TetrahedronMesh>
void ParticleGenerator3D<TetrahedronMesh>::FirstOrderGaussPoint3DGenerator()
{
	Point3D &p1 = pg_param.p1;
	Point3D &p2 = pg_param.p2;
	Point3D &p3 = pg_param.p3;
	Point3D &p4 = pg_param.p4;
	pcl_res.x = (p1.x + p2.x + p3.x + p4.x) * 0.25;
	pcl_res.y = (p1.y + p2.y + p3.y + p4.y) * 0.25;
	pcl_res.z = (p1.z + p2.z + p3.z + p4.z) * 0.25;
	pcl_res.vol = pg_param.vol;
	add_pcl(pcl_res);
}

template <typename TetrahedronMesh>
void ParticleGenerator3D<TetrahedronMesh>::SecondOrderGaussPoint3DGenerator()
{
#define alpha 0.58541020
#define belta 0.13819660
	pcl_res.vol = pg_param.vol * 0.25;
	Point3D &p1 = pg_param.p1;
	Point3D &p2 = pg_param.p2;
	Point3D &p3 = pg_param.p3;
	Point3D &p4 = pg_param.p4;
	// pcl 1
	pcl_res.x = alpha * p1.x + belta * p2.x + belta * p3.x + belta * p4.x;
	pcl_res.y = alpha * p1.y + belta * p2.y + belta * p3.y + belta * p4.y;
	pcl_res.z = alpha * p1.z + belta * p2.z + belta * p3.z + belta * p4.z;
	add_pcl(pcl_res);
	// pcl 2
	pcl_res.x = belta * p1.x + alpha * p2.x + belta * p3.x + belta * p4.x;
	pcl_res.y = belta * p1.y + alpha * p2.y + belta * p3.y + belta * p4.y;
	pcl_res.z = belta * p1.z + alpha * p2.z + belta * p3.z + belta * p4.z;
	add_pcl(pcl_res);
	// pcl 3
	pcl_res.x = belta * p1.x + belta * p2.x + alpha * p3.x + belta * p4.x;
	pcl_res.y = belta * p1.y + belta * p2.y + alpha * p3.y + belta * p4.y;
	pcl_res.z = belta * p1.z + belta * p2.z + alpha * p3.z + belta * p4.z;
	add_pcl(pcl_res);
	// pcl 4
	pcl_res.x = belta * p1.x + belta * p2.x + belta * p3.x + alpha * p4.x;
	pcl_res.y = belta * p1.y + belta * p2.y + belta * p3.y + alpha * p4.y;
	pcl_res.z = belta * p1.z + belta * p2.z + belta * p3.z + alpha * p4.z;
	add_pcl(pcl_res);
#undef alpha
#undef belta
}

template <typename TetrahedronMesh>
void ParticleGenerator3D<TetrahedronMesh>::RandomlyDistributedPoint3DGenerator()
{
	double xi, eta, zeta, N4;
	Point3D &p1 = pg_param.p1;
	Point3D &p2 = pg_param.p2;
	Point3D &p3 = pg_param.p3;
	Point3D &p4 = pg_param.p4;

	pcl_res.vol = pg_param.vol / double(pg_param.vol_div_num);
	for (size_t p_id = 0; p_id < pg_param.vol_div_num; ++p_id)
	{
		do
		{
			xi = (double)rand() / (double)(RAND_MAX + 1);
			eta = (double)rand() / (double)(RAND_MAX + 1);
			zeta = (double)rand() / (double)(RAND_MAX + 1);
		} while (xi + eta + zeta > 1.0);
		N4 = 1.0 - xi - eta - zeta;
		pcl_res.x = xi * p1.x + eta * p2.x + zeta * p3.x + N4 * p4.x;
		pcl_res.y = xi * p1.y + eta * p2.y + zeta * p3.y + N4 * p4.y;
		pcl_res.z = xi * p1.z + eta * p2.z + zeta * p3.z + N4 * p4.z;
		add_pcl(pcl_res);
	}
}

/* =========== Generate particles in grid =========== */
template <typename TetrahedronMesh>
void ParticleGenerator3D<TetrahedronMesh>::generate_pcls(GeneratorFunc cur_generator_func)
{
	size_t elem_num = mesh->get_elem_num();
	typename TetrahedronMesh::Element *elems = mesh->get_elems();
	typename TetrahedronMesh::Node *nodes = mesh->get_nodes();
	pcl_buf.set_page_size(elem_num);

	for (size_t elem_id = 0; elem_id < elem_num; ++elem_id)
	{
		typename TetrahedronMesh::Element &elem = elems[elem_id];
		typename TetrahedronMesh::Node &n1 = nodes[elem.n1];
		typename TetrahedronMesh::Node &n2 = nodes[elem.n2];
		typename TetrahedronMesh::Node &n3 = nodes[elem.n3];
		typename TetrahedronMesh::Node &n4 = nodes[elem.n4];
		Point3D &p1 = pg_param.p1;
		Point3D &p2 = pg_param.p2;
		Point3D &p3 = pg_param.p3;
		Point3D &p4 = pg_param.p4;
		p1.x = n1.x;
		p1.y = n1.y;
		p1.z = n1.z;
		p2.x = n2.x;
		p2.y = n2.y;
		p2.z = n2.z;
		p3.x = n3.x;
		p3.y = n3.y;
		p3.z = n3.z;
		p4.x = n4.x;
		p4.y = n4.y;
		p4.z = n4.z;
		pg_param.vol = elem.vol;
		(this->*cur_generator_func)();
	}
}

#endif