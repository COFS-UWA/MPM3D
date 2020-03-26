#ifndef __TRIANGLE_MESH_TO_PARTICLES_HPP__
#define __TRIANGLE_MESH_TO_PARTICLES_HPP__

#include "ItemArray.hpp"
#include "ItemBuffer.hpp"

#include "geometry.h"
#include "TriangleMesh.h"

class TriangleMeshToParticles
{
public:
	struct Particle
	{
		friend TriangleMeshToParticles;
		double x, y, vol;
	protected:
		Particle *prev;
		Particle *next;
	};
	typedef MemoryUtilities::ItemBuffer<Particle> ParticleBuffer;
	enum class GeneratorType : unsigned int
	{
		FirstOrderGaussPoint = 0,
		SecondOrderGaussPoint = 1,
		EvenlyDistributedPoint = 2,
		RandomlyDistributedPoint = 3
	};
protected: // particle generator functions
	typedef void(TriangleMeshToParticles::*GeneratorFunc)(Point &p1, Point &p2, Point &p3, double vol);

public:
	inline Particle *first(void) noexcept { return head.next; }
	inline Particle *next(Particle *cur) noexcept { return cur->next; }
	inline bool not_end_yet(Particle *cur) noexcept { return cur != &head; }
	inline bool is_end(Particle *cur) noexcept { return cur == &head; }

public:
	TriangleMeshToParticles(TriangleMesh &_mesh,
		GeneratorType _type = GeneratorType::FirstOrderGaussPoint) : 
		mesh(_mesh), pcl_num(0), evenly_div_num(1)
	{
		head.next = &head;
		head.prev = &head;

		if (unsigned int(_type) < generator_num)
			type = _type;
		else
			type = GeneratorType::FirstOrderGaussPoint; // default value
		cur_generator_func = generator_funcs[unsigned int(_type)];
		cur_generator_pcl_num = generator_pcl_num[unsigned int(_type)];
	}

	~TriangleMeshToParticles() { clear(); }

	inline void clear(void) noexcept
	{
		head.next = &head;
		head.prev = &head;
		pcl_num = 0;
		particle_buffer.clear();
	}

	inline int set_generator(GeneratorType _type) noexcept
	{
		if (unsigned int(_type) < generator_num)
		{
			type = _type;
			cur_generator_func = generator_funcs[unsigned int(_type)];
			cur_generator_pcl_num = generator_pcl_num[unsigned int(_type)];
			return 0;
		}
		return -1;
	}

	// main function
	// max_pcl_size == 0.0 means no restriction on particle size
	void generate_pcls(double max_pcl_area = 0.0);
	inline size_t get_pcl_num_per_elem(void) const noexcept { return generator_pcl_num[unsigned int(type)]; }
	inline size_t get_pcl_num(void) const noexcept { return pcl_num; }

protected:
	TriangleMesh &mesh;
	GeneratorType type;
	GeneratorFunc cur_generator_func;
	size_t cur_generator_pcl_num;

	Particle head;
	size_t pcl_num;
	ParticleBuffer particle_buffer;
	inline Particle *add_pcl(void)
	{
		Particle *res = particle_buffer.alloc();
		res->next = &head;
		res->prev = head.prev;
		head.prev->next = res;
		head.prev = res;
		++pcl_num;
		return res;
	}

public:
	inline Particle *add_pcl(Particle &pcl)
	{
		Particle *res = add_pcl();
		res->x = pcl.x;
		res->y = pcl.y;
		res->vol = pcl.vol;
	}

	inline void del_pcl(Particle &pcl)
	{
		pcl.prev->next = pcl.next;
		pcl.next->prev = pcl.prev;
		particle_buffer.del(&pcl);
		--pcl_num;
	}

protected:
	static const GeneratorFunc generator_funcs[];
	static const size_t generator_pcl_num[];
	static const size_t generator_num;
	// particle generator functions
	void FirstOrderGaussPointGenerator(Point &p1, Point &p2, Point &p3, double vol);
	void SecondOrderGaussPointGenerator(Point &p1, Point &p2, Point &p3, double vol);

	size_t evenly_div_num;
	void EvenlyDistributedPointGenerator(Point &p1, Point &p2, Point &p3, double vol);
public:
	inline void set_even_div_num(size_t num) noexcept { evenly_div_num = num; }

protected:
	void RandomlyDistributedPointGenerator(Point &p1, Point &p2, Point &p3, double vol);

public:
	void clear_points_in_rect(double xl, double xu, double yl, double yu);
	void generate_grid_points(double xl, double xu, double yl, double yu, double pcl_w, double pcl_h);
	void replace_with_grid_points(double xl, double xu, double yl, double yu, double pcl_w, double pcl_h);
};

#endif