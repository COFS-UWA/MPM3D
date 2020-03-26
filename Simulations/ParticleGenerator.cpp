#include "SimulationCore_pcp.h"

#include "TriangleMeshToParticles.h"

const typename TriangleMeshToParticles::GeneratorFunc
	TriangleMeshToParticles::generator_funcs[] = {
	&FirstOrderGaussPointGenerator,
	&SecondOrderGaussPointGenerator,
	&EvenlyDistributedPointGenerator,
	&RandomlyDistributedPointGenerator
};

const size_t TriangleMeshToParticles::generator_pcl_num[] = {
	1,
	3,
	1, // ad hoc, need to be fixed
	1  // ad hoc, need to be fixed
};

const size_t TriangleMeshToParticles::generator_num =
	sizeof(TriangleMeshToParticles::generator_funcs) /
	sizeof(TriangleMeshToParticles::generator_funcs[0]);

void TriangleMeshToParticles::FirstOrderGaussPointGenerator(Point &p1, Point &p2, Point &p3, double vol)
{
	Particle &pcl = *add_pcl();
	pcl.x = (p1.x + p2.x + p3.x) / 3.0;
	pcl.y = (p1.y + p2.y + p3.y) / 3.0;
	pcl.vol = vol;
}

void TriangleMeshToParticles::SecondOrderGaussPointGenerator(Point &p1, Point &p2, Point &p3, double vol)
{
	vol /= 3.0;
	Particle &pcl1 = *add_pcl();
	pcl1.x = 2.0 / 3.0 * p1.x + 1.0 / 6.0 * p2.x + 1.0 / 6.0 * p3.x;
	pcl1.y = 2.0 / 3.0 * p1.y + 1.0 / 6.0 * p2.y + 1.0 / 6.0 * p3.y;
	pcl1.vol = vol;
	Particle &pcl2 = *add_pcl();
	pcl2.x = 1.0 / 6.0 * p1.x + 2.0 / 3.0 * p2.x + 1.0 / 6.0 * p3.x;
	pcl2.y = 1.0 / 6.0 * p1.y + 2.0 / 3.0 * p2.y + 1.0 / 6.0 * p3.y;
	pcl2.vol = vol;
	Particle &pcl3 = *add_pcl();
	pcl3.x = 1.0 / 6.0 * p1.x + 1.0 / 6.0 * p2.x + 2.0 / 3.0 * p3.x;
	pcl3.y = 1.0 / 6.0 * p1.y + 1.0 / 6.0 * p2.y + 2.0 / 3.0 * p3.y;
	pcl3.vol = vol;
}

void TriangleMeshToParticles::EvenlyDistributedPointGenerator(Point &p1, Point &p2, Point &p3, double vol)
{
	vol /= double(evenly_div_num * evenly_div_num);
	double dx21 = (p2.x - p1.x) / double(evenly_div_num);
	double dy21 = (p2.y - p1.y) / double(evenly_div_num);
	double dx32 = (p3.x - p2.x) / double(evenly_div_num);
	double dy32 = (p3.y - p2.y) / double(evenly_div_num);
	double row_inc, col_inc;
	for (size_t row_id = 0; row_id < evenly_div_num; ++row_id)
	{
		for (size_t col_id = 0; col_id < row_id; ++col_id)
		{
			Particle &pcl1 = *add_pcl();
			row_inc = double(row_id) + 2.0 / 3.0;
			col_inc = double(col_id) + 1.0 / 3.0;
			pcl1.x = p1.x + row_inc * dx21 + col_inc * dx32;
			pcl1.y = p1.y + row_inc * dy21 + col_inc * dy32;
			pcl1.vol = vol;
			Particle &pcl2 = *add_pcl();
			row_inc = double(row_id) + 1.0 / 3.0;
			col_inc = double(col_id) + 2.0 / 3.0;
			pcl2.x = p1.x + row_inc * dx21 + col_inc * dx32;
			pcl2.y = p1.y + row_inc * dy21 + col_inc * dy32;
			pcl2.vol = vol;
		}
		Particle &pcln = *add_pcl(); // last pcl of each row
		row_inc = double(row_id) + 2.0 / 3.0;
		col_inc = double(row_id) + 1.0 / 3.0;
		pcln.x = p1.x + row_inc * dx21 + col_inc * dx32;
		pcln.y = p1.y + row_inc * dy21 + col_inc * dy32;
		pcln.vol = vol;
	}
}

namespace
{
	inline void get_rand_num(double &xi, double &eta)
	{
		// xi  - [0,        0.5]
		// eta - [0.5 - xi, 0.5]
		xi  = (double)rand() / (double)(RAND_MAX + 1) * 0.5;
		eta = (double)rand() / (double)(RAND_MAX + 1) * xi + 0.5 - xi;
	}
}

void TriangleMeshToParticles::RandomlyDistributedPointGenerator(Point &p1, Point &p2, Point &p3, double vol)
{
	size_t evenly_div_num2 = evenly_div_num * evenly_div_num;
	vol /= double(evenly_div_num2);
	double dx21 = (p2.x - p1.x) / double(evenly_div_num);
	double dy21 = (p2.y - p1.y) / double(evenly_div_num);
	double dx32 = (p3.x - p2.x) / double(evenly_div_num);
	double dy32 = (p3.y - p2.y) / double(evenly_div_num);
	double row_inc, col_inc;
	double xi, eta;
	for (size_t row_id = 0; row_id < evenly_div_num; ++row_id)
	{
		for (size_t col_id = 0; col_id < row_id; ++col_id)
		{
			Particle &pcl1 = *add_pcl();
			get_rand_num(xi, eta);
			row_inc = double(row_id) + 1.0 - xi;
			col_inc = double(col_id) + 1.0 - xi - eta;
			pcl1.x = p1.x + row_inc * dx21 + col_inc * dx32;
			pcl1.y = p1.y + row_inc * dy21 + col_inc * dy32;
			pcl1.vol = vol;
			Particle &pcl2 = *add_pcl();
			get_rand_num(xi, eta);
			row_inc = double(row_id) + xi;
			col_inc = double(col_id) + xi + eta;
			pcl2.x = p1.x + row_inc * dx21 + col_inc * dx32;
			pcl2.y = p1.y + row_inc * dy21 + col_inc * dy32;
			pcl2.vol = vol;
		}
		Particle &pcln = *add_pcl(); // last pcl of each row
		get_rand_num(xi, eta);
		row_inc = double(row_id) + 1.0 - xi;
		col_inc = double(row_id) + 1.0 - xi - eta;
		pcln.x = p1.x + row_inc * dx21 + col_inc * dx32;
		pcln.y = p1.y + row_inc * dy21 + col_inc * dy32;
		pcln.vol = vol;
	}
}

void TriangleMeshToParticles::generate_pcls(double max_pcl_area)
{
	size_t elem_num = mesh.get_elem_num();
	TriangleMesh::Element *elems = mesh.get_elems();
	TriangleMesh::Node *nodes = mesh.get_nodes();
	Point elem_n1, elem_n2, elem_n3;
	double elem_area;

	if (max_pcl_area == 0.0)
	{
		particle_buffer.set_page_size(elem_num);
		for (size_t elem_id = 0; elem_id < elem_num; ++elem_id)
		{
			TriangleMesh::Element &elem = elems[elem_id];
			TriangleMesh::Node &n1 = nodes[elem.n1];
			elem_n1.x = n1.x;
			elem_n1.y = n1.y;
			TriangleMesh::Node &n2 = nodes[elem.n2];
			elem_n2.x = n2.x;
			elem_n2.y = n2.y;
			TriangleMesh::Node &n3 = nodes[elem.n3];
			elem_n3.x = n3.x;
			elem_n3.y = n3.y;
			elem_area = abs(((n2.x - n1.x) * (n3.y - n1.y) - (n3.x - n1.x) * (n2.y - n1.y)) / 2.0);
			(this->*cur_generator_func)(elem_n1, elem_n2, elem_n3, elem_area);
		}
	}
	else
	{
		particle_buffer.set_page_size(mesh.get_area() / max_pcl_area / 2.0);
		MemoryUtilities::ItemArray<Point, 2, 22> pt_mem;
		Point *low_pts, *upp_pts;
		double dx21, dy21, dx32, dy32;
		for (size_t elem_id = 0; elem_id < elem_num; ++elem_id)
		{
			TriangleMesh::Element elem = elems[elem_id];
			TriangleMesh::Node &n1 = nodes[elem.n1];
			elem_n1.x = n1.x;
			elem_n1.y = n1.y;
			TriangleMesh::Node &n2 = nodes[elem.n2];
			elem_n2.x = n2.x;
			elem_n2.y = n2.y;
			TriangleMesh::Node &n3 = nodes[elem.n3];
			elem_n3.x = n3.x;
			elem_n3.y = n3.y;
			elem_area = abs(((n2.x - n1.x) * (n3.y - n1.y) - (n3.x - n1.x) * (n2.y - n1.y)) / 2.0);
			size_t div_num = size_t(ceil(sqrt(elem_area / (max_pcl_area * get_pcl_num_per_elem()))));
			if (div_num == 0) div_num = 1;
			elem_area /= double(div_num * div_num);
			dx21 = (elem_n2.x - elem_n1.x) / double(div_num);
			dy21 = (elem_n2.y - elem_n1.y) / double(div_num);
			dx32 = (elem_n3.x - elem_n2.x) / double(div_num);
			dy32 = (elem_n3.y - elem_n2.y) / double(div_num);
			pt_mem.reserve(div_num + div_num + 2);
			upp_pts = pt_mem.get_mem();
			low_pts = upp_pts + div_num + 1;
			low_pts[0].x = elem_n1.x;
			low_pts[0].y = elem_n1.y;
			for (size_t line_id = 0; line_id < div_num; ++line_id)
			{
				Point *tmp = upp_pts;
				upp_pts = low_pts;
				low_pts = tmp;
				low_pts[0].x = upp_pts[0].x + dx21;
				low_pts[0].y = upp_pts[0].y + dy21;
				size_t col_num = line_id + 1;
				for (size_t col_id = 0; col_id < col_num; ++col_id)
				{
					low_pts[col_id + 1].x = low_pts[col_id].x + dx32;
					low_pts[col_id + 1].y = low_pts[col_id].y + dy32;
					(this->*cur_generator_func)(upp_pts[col_id], low_pts[col_id], low_pts[col_id + 1], elem_area);
				}
				for (size_t col_id = 0; col_id < line_id; ++col_id)
				{
					(this->*cur_generator_func)(upp_pts[col_id], low_pts[col_id], upp_pts[col_id + 1], elem_area);
				}
			}
		}
	}
}

void TriangleMeshToParticles::clear_points_in_rect(
	double xl, double xu, double yl, double yu)
{
	for (Particle *ppcl = first(); not_end_yet(ppcl); ppcl = next(ppcl))
	{
		Particle &pcl = *ppcl;
		if (pcl.x >= xl && pcl.x <= xu &&
			pcl.y >= yl && pcl.y <= yu)
			del_pcl(pcl);
	}
}

void TriangleMeshToParticles::generate_grid_points(
	double xl, double xu,
	double yl, double yu,
	double pcl_w, double pcl_h
	)
{
	double width, height, pcl_start_x, pcl_start_y;
	size_t pcl_x_num, pcl_y_num;
	// x direction
	width = xu - xl;
	pcl_x_num = size_t(ceil(width / pcl_w));
	pcl_w = width / double(pcl_x_num);
	pcl_start_x = xl + pcl_w * 0.5;
	// y direction
	height = yu - yl;
	pcl_y_num = size_t(ceil(height / pcl_h));
	pcl_h = height / double(pcl_y_num);
	pcl_start_y = yl + pcl_h * 0.5;
	for (size_t row_id = 0; row_id < pcl_y_num; ++row_id)
		for (size_t col_id = 0; col_id < pcl_x_num; ++col_id)
		{
			Particle &pcl = *add_pcl();
			pcl.x = pcl_start_x + pcl_w * double(col_id);
			pcl.y = pcl_start_y + pcl_h * double(row_id);
			pcl.vol = pcl_w * pcl_h;
		}
}

void TriangleMeshToParticles::replace_with_grid_points(
	double xl, double xu,
	double yl, double yu,
	double pcl_w, double pcl_h
	)
{
	clear_points_in_rect(xl, xu, yl, yu);
	generate_grid_points(xl, xu, yl, yu, pcl_w, pcl_h);
}
