#ifndef __Triangle_Mesh_h__
#define __Triangle_Mesh_h__

#include "TriangleMeshTemplate.hpp"
#include "SearchingGrid2D.hpp"

namespace TriangleMesh_Internal
{
	struct Node
	{
		size_t id;
		double x, y;
	};

	struct Element
	{
		size_t id;
		size_t n1, n2, n3;
		double area;
	};

	struct Edge { size_t n1, n2; };
};

struct TriangleMesh : 
	public TriangleMeshTemplate<TriangleMesh_Internal::Node,
			TriangleMesh_Internal::Element, TriangleMesh_Internal::Edge>
{
public:
	typedef TriangleMesh_Internal::Node Node;
	typedef TriangleMesh_Internal::Element Element;
	typedef TriangleMesh_Internal::Edge Edge;

	TriangleMesh() {}
	~TriangleMesh() { clear(); }

	inline SearchingGrid2D<TriangleMesh>& get_bg_grid() { return search_bg_grid; }
	inline double get_bg_grid_xl() { return search_bg_grid.get_x_min(); }
	inline double get_bg_grid_xu() { return search_bg_grid.get_x_max(); }
	inline double get_bg_grid_yl() { return search_bg_grid.get_y_min(); }
	inline double get_bg_grid_yu() { return search_bg_grid.get_y_max(); }
	inline double get_bg_grid_hx() { return search_bg_grid.get_hx(); }
	inline double get_bg_grid_hy() { return search_bg_grid.get_hy(); }

	void clear()
	{
		TriangleMeshTemplate<TriangleMesh_Internal::Node,
			TriangleMesh_Internal::Element,
			TriangleMesh_Internal::Edge>::clear();
		search_bg_grid.clear();
	}

	int init_search_grid(double _hx, double _hy);

	// search using background grid
	template <typename Point2D>
	inline Element* find_in_which_element(Point2D& pt)
	{
		return search_bg_grid.find_in_which_element<Point2D>(pt);
	}

protected:
	SearchingGrid2D<TriangleMesh> search_bg_grid;
};

#endif