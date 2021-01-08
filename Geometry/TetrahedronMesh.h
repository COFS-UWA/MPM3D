#ifndef __Tetrahedron_Mesh_H__
#define __Tetrahedron_Mesh_H__

#include "TetrahedronMeshTemplate.hpp"
#include "SearchingGrid3D.hpp"

namespace TetrahedronMesh_Internal
{
	struct Node
	{
		size_t id;
		double x, y, z;
	};

	struct Element
	{
		size_t id;
		size_t n1, n2, n3, n4;
		double vol;
	};

	struct Edge { size_t n1, n2; };
};

struct TetrahedronMesh : 
	public TetrahedronMeshTemplate<TetrahedronMesh_Internal::Node,
			TetrahedronMesh_Internal::Element, TetrahedronMesh_Internal::Edge>
{
public:
	typedef TetrahedronMesh_Internal::Node Node;
	typedef TetrahedronMesh_Internal::Element Element;
	typedef TetrahedronMesh_Internal::Edge Edge;
	typedef SearchingGrid3D<TetrahedronMesh> SearchGrid;

protected:
	SearchGrid search_bg_grid;

public:
	inline SearchGrid &get_bg_grid() { return search_bg_grid; }
	inline double get_bg_grid_xl() { return search_bg_grid.get_xl(); }
	inline double get_bg_grid_xu() { return search_bg_grid.get_xu(); }
	inline double get_bg_grid_yl() { return search_bg_grid.get_yl(); }
	inline double get_bg_grid_yu() { return search_bg_grid.get_yu(); }
	inline double get_bg_grid_zl() { return search_bg_grid.get_zl(); }
	inline double get_bg_grid_zu() { return search_bg_grid.get_zu(); }
	inline double get_bg_grid_hx() { return search_bg_grid.get_hx(); }
	inline double get_bg_grid_hy() { return search_bg_grid.get_hy(); }
	inline double get_bg_grid_hz() { return search_bg_grid.get_hz(); }

	int init_search_grid(double _hx, double _hy, double _hz);

	void rotate_mesh(double dx_ang, double dy_ang, double dz_ang) noexcept;
	void translate_mesh(double dx, double dy, double dz) noexcept;

	// search using background grid
	template <typename Point3D>
	inline Element* find_in_which_element(Point3D& pt)
	{ return search_bg_grid.find_in_which_element<Point3D>(pt); }
};

#endif