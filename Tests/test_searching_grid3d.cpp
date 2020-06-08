#include "Tests_pcp.h"

#include "Model_T3D_ME_s.h"

#include "test_geometry.h"

namespace
{

class Model_test : public Model_T3D_ME_s
{
public:
	int alloc_grids(Cube box, double _hx, double _hy, double _hz)
	{
		int res = search_bg_grid.alloc_grids(box, _hx, _hy, _hz);
		search_bg_grid.set_mesh_info(*this);
		return res;
	}

	bool detect_AABB_tetrahedron_collision(Cube &aabb, Element &e)
	{
		search_bg_grid.cal_seperating_axes(e);
		bool res = search_bg_grid.detect_AABB_tetrahedron_collision(aabb, e);
		if (res)
			std::cout << "contact \n";
		else
			std::cout << "not contact\n";
		return res;
	}

};

}

// test aabb tetrahedron collision
void test_searching_grid3d1()
{
	Model_test model;

	double node_coords[] = {
		// te1
		0.5, 0.5, 0.5,
		3.0, 0.5, 0.5,
		0.5, 3.0, 0.5,
		0.5, 0.5, 3.0,
		// te2
		0.75, 0.75, -0.5,
		5.0,  0.75, -0.5,
		0.75,  5.0, -0.5,
		0.75, 0.75, 5.0,
		// te3
		1.25, 0.25, 0.25,
		1.25, 0.75, 0.25,
		-0.25, 0.25, 0.25,
		1.25, 0.25, 0.75,
		// te4
		2.0, 0.25, 0.25,
		3.0, 0.25, 0.25,
		2.0, 0.75, 0.25,
		2.0, 0.25, 1.0,
		// te5
		2.0, 0.5, 0.25,
		3.0, 0.25, 0.3,
		4.0, 1.0, 0.5,
		3.0, 0.58, 0.75
	};
	size_t elem_n_ids[] = {
		0, 1, 2, 3,
		4, 5, 6, 7,
		8, 9, 10, 11,
		12, 13, 14, 15,
		16, 17, 18, 19
	};
	model.init_mesh(node_coords, sizeof(node_coords) / (sizeof(double) * 3),
		elem_n_ids, sizeof(elem_n_ids) / (sizeof(size_t) * 4));

	Cube grid_box = { 0.0, 1.0, 0.0, 1.0, 0.0, 1.0 };
	model.alloc_grids(grid_box, 1.0, 1.0, 1.0);

	Model_test::Element *elems = model.get_elems();
	// in contact
	model.detect_AABB_tetrahedron_collision(grid_box, elems[0]);
	// in contact 
	model.detect_AABB_tetrahedron_collision(grid_box, elems[1]);
	// in contact
	model.detect_AABB_tetrahedron_collision(grid_box, elems[2]);
	// not contact
	model.detect_AABB_tetrahedron_collision(grid_box, elems[3]);
	// not contact
	model.detect_AABB_tetrahedron_collision(grid_box, elems[4]);
}

void test_searching_grid3d2()
{
	Model_T3D_ME_s model;

	// load mesh
	model.load_mesh_from_hdf5("..\\..\\Asset\\brick_mesh.h5");

	//model.init_search_grid(2.0, 2.0, 2.0);
	model.init_search_grid(1.5, 1.5, 1.5);

	typedef SearchingGrid3D<Model_T3D_ME_s> SG;
	SG &sg = model.search_bg_grid;
	size_t grid_num = sg.get_grid_num();
	SG::Grid *grids = sg.get_grids();
	std::cout << "grid num: " << grid_num << "\n";
	for (size_t g_id = 0; g_id < grid_num; ++g_id)
	{
		SG::Grid &grd = grids[g_id];
		std::cout << "grid (" << grd.x_id << ", " << grd.y_id << ", " << grd.z_id << "): ";
		for (SG::ElemPointer *ep = grd.pelems; ep; ep = ep->next)
		{
			Model_T3D_ME_s::Element &e = *ep->e;
			std::cout << e.id << ", ";
		}
		std::cout << "\n";
	}
}

void test_searching_grid3d3()
{
	Model_T3D_ME_s model;

	// load mesh
	model.load_mesh_from_hdf5("..\\..\\Asset\\column_mesh.h5");

	// searching grid
	model.init_search_grid(5.0, 5.0, 5.0);

	// generate particles
	ParticleGenerator3D<Model_T3D_ME_s> pcl_gen;
	pcl_gen.generate_pcls_first_order_gauss(model);
	model.init_pcls(pcl_gen, 10.0);

	size_t pcl_num = model.get_pcl_num();
	Model_T3D_ME_s::Particle *pcls = model.get_pcls();
	size_t e_id1 = 0, e_id2 = 0;
	for (size_t pcl_id = 0; pcl_id < pcl_num; ++pcl_id)
	{
		Model_T3D_ME_s::Particle &pcl = pcls[pcl_id];
		e_id1 = model.find_in_which_element(pcl)->id;
		e_id2 = model.find_in_which_element_bf(pcl)->id;
		if (e_id1 == e_id2)
			std::cout << "pcl " << pcl.id << " in elem: " << e_id1 << ", " << e_id2 << "\n";
		else
			std::cout << "pcl " << pcl.id << " in elem: " << e_id1 << ", " << e_id2 << " error!\n";
	}
}

