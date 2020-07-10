#include "Tests_pcp.h"

#include "Model_T2D_ME_s.h"

#include "test_geometry.h"

namespace
{

class Model_test : public Model_T2D_ME_s
{
public:
	int alloc_grids(Rect &box, double _hx, double _hy)
	{
		int res = search_bg_grid.alloc_grids(box, _hx, _hy);
		search_bg_grid.set_mesh_info(*this);
		return res;
	}

	bool detect_AABB_triangle_collision(Rect &aabb, Element &e)
	{
		bool res = search_bg_grid.detect_AABB_triangle_collision(aabb, e);
		if (res)
			std::cout << "contact \n";
		else
			std::cout << "not contact\n";
		return res;
	}

};

}

// test aabb tetrahedron collision
void test_searching_grid2d1()
{
	Model_test model;

	double node_coords[] = {
		// te1
		0.4,  0.5,
		0.6, -1.0,
		1.5, -0.1,
		// te2
		0.1, 1.1,
		1.1, 0.1,
		1.1, 1.1,
		// te3
		1.0, 2.0,
		2.0, 1.0,
		2.0, 2.0,
		// te4
		0.1, 1.1,
		1.1, 0.1,
		1.2, 0.3
	};
	size_t elem_n_ids[] = {
		0, 1, 2,
		3, 4, 5,
		6, 7, 8,
		9, 10, 11
	};
	model.init_mesh(node_coords, sizeof(node_coords) / (sizeof(double) * 2),
					elem_n_ids, sizeof(elem_n_ids) / (sizeof(size_t) * 3));

	Rect grid_box = { 0.0, 1.0, 0.0, 1.0 };
	model.alloc_grids(grid_box, 1.0, 1.0);

	Model_test::Element *elems = model.get_elems();
	// in contact
	model.detect_AABB_triangle_collision(grid_box, elems[0]);
	// in contact 
	model.detect_AABB_triangle_collision(grid_box, elems[1]);
	// not contacct
	model.detect_AABB_triangle_collision(grid_box, elems[2]);
	// in contacct
	model.detect_AABB_triangle_collision(grid_box, elems[3]);
}
