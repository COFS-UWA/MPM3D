#ifndef __Searching_Grid_3D_hpp__
#define __Searching_Grid_3D_hpp__

#include "ItemBuffer.hpp"
#include "Geometry.h"

// Accelerate spatial searching of tetrahedron mesh
// Assumptions:
//	1. Tetrahedron mesh has is_in_tetrahedron(MeshElement &e, double x, double y, double z)
//	2. Tetrahedron mesh has get_elems() and get_elem_num()
//  3. Tetrahedron mesh has get_bounding_box()
template <typename TetrahedronMesh>
class SearchingGrid3D
{
public:
	typedef typename TetrahedronMesh::Element MeshElement;
	typedef typename TetrahedronMesh::Node MeshNode;
	
	struct ElemPointer
	{
		MeshElement *e;
		ElemPointer *next;
	};

	struct Grid
	{
		size_t x_id, y_id, z_id;
		ElemPointer *pelems;
	};

protected:
	size_t x_num, y_num, z_num;
	size_t xy_num, num;
	double hx, hy, hz;
	double x_min, x_max;
	double y_min, y_max;
	double z_min, z_max;
	Grid *grids;

	TetrahedronMesh *mesh;
	size_t node_num;
	MeshNode *nodes;
	size_t elem_num;
	MeshElement *elems;

	MemoryUtils::ItemBuffer<ElemPointer> pe_buffer;

	// seperating axises: 4 face normal + 3 * 6 edge cross products
	Vector3D seperating_axes[22];

public:
	SearchingGrid3D() :
		x_num(0), y_num(0), z_num(0),
		xy_num(0), num(0),
		grids(nullptr), mesh(nullptr) {}
	~SearchingGrid3D() { clear(); }

	inline double get_hx() { return hx; }
	inline double get_hy() { return hy; }
	inline double get_hz() { return hz; }
	inline size_t get_grid_num() { return num; }
	inline Grid *get_grids() { return grids; }

	int init(TetrahedronMesh &_mesh, double hx, double hy, double hz)
	{
		Cube &mesh_bbox = _mesh.get_bounding_box();
		if (alloc_grids(mesh_bbox, hx, hy, hz) < 0)
			return -1;

		set_mesh_info(_mesh);
		
		for (size_t e_id = 0; e_id < elem_num; ++e_id)
			add_elem_to_grids(elems[e_id]);

		return 0;
	}

	void clear()
	{
		if (grids)
		{
			delete[] grids;
			grids = nullptr;
		}
		x_num = 0;
		y_num = 0;
		z_num = 0;
		xy_num = 0;
		num = 0;

		mesh = nullptr;
		node_num = 0;
		nodes = nullptr;
		elem_num = 0;
		elems = nullptr;

		pe_buffer.clear();
	}

	template <typename Point3D>
	inline MeshElement *find_in_which_element(Point3D &point)
	{
		if (point.x < x_min || point.x > x_max ||
			point.y < y_min || point.y > y_max ||
			point.z < z_min || point.z > z_max)
			return nullptr;

		size_t x_id = size_t((point.x - x_min) / hx);
		size_t y_id = size_t((point.y - y_min) / hy);
		size_t z_id = size_t((point.z - z_min) / hz);
		Grid &g = get_grid(x_id, y_id, z_id);
		MeshElement *elem;
		for (ElemPointer *pelem = g.pelems; pelem; pelem = pelem->next)
		{
			elem = pelem->e;
			if (mesh->is_in_tetrahedron(*elem, point))
				return elem;
		}
		return nullptr;
	}

protected: // helper functions
public:
	int alloc_grids(Cube box, double _hx, double _hy, double _hz)
	{
		clear();
		double x_len = box.xu - box.xl;
		double y_len = box.yu - box.yl;
		double z_len = box.zu - box.zl;
		if (x_len <= 0.0 || y_len <= 0.0 || z_len <= 0.0)
			return -1;

		hx = _hx;
		hy = _hy;
		hz = _hz;
		x_num = size_t(ceil(x_len / hx));
		y_num = size_t(ceil(y_len / hy));
		z_num = size_t(ceil(z_len / hz));
		double x_pad = (double(x_num) * hx - x_len) * 0.5;
		double y_pad = (double(y_num) * hy - y_len) * 0.5;
		double z_pad = (double(z_num) * hz - z_len) * 0.5;
		x_min = box.xl - x_pad;
		x_max = box.xu + x_pad;
		y_min = box.yl - y_pad;
		y_max = box.yu + y_pad;
		z_min = box.zl - z_pad;
		z_max = box.zu + z_pad;
		xy_num = x_num * y_num;
		num = xy_num * z_num;
		grids = new Grid[num];
		Grid *pg = grids;
		for (size_t z_id = 0; z_id < z_num; ++z_id)
			for (size_t y_id = 0; y_id < y_num; ++y_id)
				for (size_t x_id = 0; x_id < x_num; ++x_id)
				{
					pg->x_id = x_id;
					pg->y_id = y_id;
					pg->z_id = z_id;
					pg->pelems = nullptr;
					++pg;
				}

		return 0;
	}

	void set_mesh_info(TetrahedronMesh &_mesh)
	{
		mesh = &_mesh;
		node_num = mesh->get_node_num();
		nodes = mesh->get_nodes();
		elem_num = mesh->get_elem_num();
		elems = mesh->get_elems();
	}

	void add_elem_to_grids(MeshElement &elem)
	{
		cal_seperating_axes(elem);
		
		Cube elem_bbox;
		get_elem_bbox(elem, elem_bbox);

		size_t xl_id = size_t(floor((elem_bbox.xl - x_min) / hx));
		size_t xu_id = size_t(ceil( (elem_bbox.xu - x_min) / hx));
		size_t yl_id = size_t(floor((elem_bbox.yl - y_min) / hy));
		size_t yu_id = size_t(ceil( (elem_bbox.yu - y_min) / hy));
		size_t zl_id = size_t(floor((elem_bbox.zl - z_min) / hz));
		size_t zu_id = size_t(ceil( (elem_bbox.zu - z_min) / hz));
		Cube grid_box;
		for (size_t z_id = zl_id; z_id < zu_id; ++z_id)
		{
			grid_box.zl = z_min + double(z_id) * hz;
			grid_box.zu = grid_box.zl + hz;
			for (size_t y_id = yl_id; y_id < yu_id; ++y_id)
			{
				grid_box.yl = y_min + double(y_id) * hy;
				grid_box.yu = grid_box.yl + hy;
				for (size_t x_id = xl_id; x_id < xu_id; ++x_id)
				{
					grid_box.xl = x_min + double(x_id) * hx;
					grid_box.xu = grid_box.xl + hx;
					if (detect_AABB_tetrahedron_collision(grid_box, elem))
					{
						Grid &cur_grid = get_grid(x_id, y_id, z_id);
						add_elem_to_grid(cur_grid, elem);
					}
				}
			}
		}
	}

	inline Grid &get_grid(size_t x_id, size_t y_id, size_t z_id)
	{
		return grids[xy_num * z_id + x_num * y_id + x_id];
	}

	inline void add_elem_to_grid(Grid &g, MeshElement &e)
	{
		ElemPointer *pe = pe_buffer.alloc();
		pe->e = &e;
		pe->next = g.pelems;
		g.pelems = pe;
	}

	inline void clear_elems_in_grid(Grid &g)
	{
		ElemPointer *pe = g.pelems;
		ElemPointer *pe_tmp;
		while (pe)
		{
			pe_tmp = pe;
			pe = pe->next;
			pe_buffer.del(pe_tmp);
		}
	}

	inline void get_elem_bbox(MeshElement& elem, Cube &elem_bbox)
	{
		MeshNode &n1 = nodes[elem.n1];
		MeshNode &n2 = nodes[elem.n2];
		MeshNode &n3 = nodes[elem.n3];
		MeshNode &n4 = nodes[elem.n4];

		double xl, xu, yl, yu, zl, zu;
		xl = n1.x;
		xu = xl;
		if (xl > n2.x)
			xl = n2.x;
		if (xu < n2.x)
			xu = n2.x;
		if (xl > n3.x)
			xl = n3.x;
		if (xu < n3.x)
			xu = n3.x;
		if (xl > n4.x)
			xl = n4.x;
		if (xu < n4.x)
			xu = n4.x;

		yl = n1.y;
		yu = yl;
		if (yl > n2.y)
			yl = n2.y;
		if (yu < n2.y)
			yu = n2.y;
		if (yl > n3.y)
			yl = n3.y;
		if (yu < n3.y)
			yu = n3.y;
		if (yl > n4.y)
			yl = n4.y;
		if (yu < n4.y)
			yu = n4.y;

		zl = n1.z;
		zu = zl;
		if (zl > n2.z)
			zl = n2.z;
		if (zu < n2.z)
			zu = n2.z;
		if (zl > n3.z)
			zl = n3.z;
		if (zu < n3.z)
			zu = n3.z;
		if (zl > n4.z)
			zl = n4.z;
		if (zu < n4.z)
			zu = n4.z;

		elem_bbox.xl = xl;
		elem_bbox.xu = xu;
		elem_bbox.yl = yl;
		elem_bbox.yu = yu;
		elem_bbox.zl = zl;
		elem_bbox.zu = zu;
	}

	// test if aligned-axis bounding box intersect tetrahedron
	bool detect_AABB_tetrahedron_collision(Cube &aabb, MeshElement &e)
	{
		MeshNode &n1 = nodes[e.n1];
		MeshNode &n2 = nodes[e.n2];
		MeshNode &n3 = nodes[e.n3];
		MeshNode &n4 = nodes[e.n4];

		// whether tetrahedron nodes locate in box
		// efficient when tetrahedron is much smaller than grid
		if (aabb.is_in_box(n1) || aabb.is_in_box(n2) || 
			aabb.is_in_box(n3) || aabb.is_in_box(n4))
			return true;

		// whether box corners locate in tetrahedron
		// efficient when grid is much smaller than tetrahedron
		if (mesh->is_in_tetrahedron(e, aabb.xl, aabb.yl, aabb.zl) ||
			mesh->is_in_tetrahedron(e, aabb.xl, aabb.yl, aabb.zu) ||
			mesh->is_in_tetrahedron(e, aabb.xl, aabb.yu, aabb.zl) ||
			mesh->is_in_tetrahedron(e, aabb.xl, aabb.yu, aabb.zu) ||
			mesh->is_in_tetrahedron(e, aabb.xu, aabb.yl, aabb.zl) ||
			mesh->is_in_tetrahedron(e, aabb.xu, aabb.yl, aabb.zu) ||
			mesh->is_in_tetrahedron(e, aabb.xu, aabb.yu, aabb.zl) ||
			mesh->is_in_tetrahedron(e, aabb.xu, aabb.yu, aabb.zu))
			return true;

		// applied separating axis theory
		// for case when grid and tetrahedron is of comparable size
		// move origin to centre of the Cube
		double box_xc = (aabb.xl + aabb.xu) * 0.5;
		double box_yc = (aabb.yl + aabb.yu) * 0.5;
		double box_zc = (aabb.zl + aabb.zu) * 0.5;
		Point3D n1_m, n2_m, n3_m, n4_m; // moved tetrahedron nodes
		n1_m.x = n1.x - box_xc;
		n1_m.y = n1.y - box_yc;
		n1_m.z = n1.z - box_zc;
		n2_m.x = n2.x - box_xc;
		n2_m.y = n2.y - box_yc;
		n2_m.z = n2.z - box_zc;
		n3_m.x = n3.x - box_xc;
		n3_m.y = n3.y - box_yc;
		n3_m.z = n3.z - box_zc;
		n4_m.x = n4.x - box_xc;
		n4_m.y = n4.y - box_yc;
		n4_m.z = n4.z - box_zc;
		// if there is one seperating axis, there is no collision
		if (is_seperating_axis(seperating_axes[0], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[1], n1_m, n2_m, n3_m, n4_m) || 
			is_seperating_axis(seperating_axes[2], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[3], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[4], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[5], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[6], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[7], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[8], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[9], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[10], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[11], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[12], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[13], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[14], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[15], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[16], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[17], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[18], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[19], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[20], n1_m, n2_m, n3_m, n4_m) ||
			is_seperating_axis(seperating_axes[21], n1_m, n2_m, n3_m, n4_m))
			return false;
		return true;
	}

	// return ture if this is the seperating axis
	// assume that the origin is the centre of the box
	inline bool is_seperating_axis(Vector3D &sa,
		Point3D &p1, Point3D &p2, Point3D &p3, Point3D &p4)
	{
#define Norm_tol 1.0e-5
		if (sa.norm() < Norm_tol)
			return false;
		double box_range = 0.5 * (hx * abs(sa.x) + hy * abs(sa.y) + hz * abs(sa.z));
		double p1_proj = p1.x * sa.x + p1.y * sa.y + p1.z * sa.z;
		double p2_proj = p2.x * sa.x + p2.y * sa.y + p2.z * sa.z;
		double p3_proj = p3.x * sa.x + p3.y * sa.y + p3.z * sa.z;
		double p4_proj = p4.x * sa.x + p4.y * sa.y + p4.z * sa.z;
		if ((p1_proj >  box_range && p2_proj >  box_range && p3_proj >  box_range && p4_proj >  box_range) ||
			(p1_proj < -box_range && p2_proj < -box_range && p3_proj < -box_range && p4_proj < -box_range))
			return true;
		return false;
#undef Norm_tol
	}

	void cal_seperating_axes(MeshElement &e)
	{
		MeshNode &n1 = nodes[e.n1];
		MeshNode &n2 = nodes[e.n2];
		MeshNode &n3 = nodes[e.n3];
		MeshNode &n4 = nodes[e.n4];
		double e12_x = n1.x - n2.x;
		double e12_y = n1.y - n2.y;
		double e12_z = n1.z - n2.z;
		double e13_x = n1.x - n3.x;
		double e13_y = n1.y - n3.y;
		double e13_z = n1.z - n3.z;
		double e14_x = n1.x - n4.x;
		double e14_y = n1.y - n4.y;
		double e14_z = n1.z - n4.z;
		double e23_x = n2.x - n3.x;
		double e23_y = n2.y - n3.y;
		double e23_z = n2.z - n3.z;
		double e24_x = n2.x - n4.x;
		double e24_y = n2.y - n4.y;
		double e24_z = n2.z - n4.z;
		double e34_x = n3.x - n4.x;
		double e34_y = n3.y - n4.y;
		double e34_z = n3.z - n4.z;
		// need normalization?
		// 4 face normal
		seperating_axes[0].cross(e13_x, e13_y, e13_z, e12_x, e12_y, e12_z);
		seperating_axes[1].cross(e12_x, e12_y, e12_z, e14_x, e14_y, e14_z);
		seperating_axes[2].cross(e14_x, e14_y, e14_z, e13_x, e13_y, e13_z);
		seperating_axes[3].cross(e23_x, e23_y, e23_z, e24_x, e24_y, e24_z);
		// 3 * 6 edge cross product
		seperating_axes[4].cross(1.0, 0.0, 0.0, e12_x, e12_y, e12_z);
		seperating_axes[5].cross(0.0, 1.0, 0.0, e12_x, e12_y, e12_z);
		seperating_axes[6].cross(0.0, 0.0, 1.0, e12_x, e12_y, e12_z);
		seperating_axes[7].cross(1.0, 0.0, 0.0, e13_x, e13_y, e13_z);
		seperating_axes[8].cross(0.0, 1.0, 0.0, e13_x, e13_y, e13_z);
		seperating_axes[9].cross(0.0, 0.0, 1.0, e13_x, e13_y, e13_z);
		seperating_axes[10].cross(1.0, 0.0, 0.0, e14_x, e14_y, e14_z);
		seperating_axes[11].cross(0.0, 1.0, 0.0, e14_x, e14_y, e14_z);
		seperating_axes[12].cross(0.0, 0.0, 1.0, e14_x, e14_y, e14_z);
		seperating_axes[13].cross(1.0, 0.0, 0.0, e23_x, e23_y, e23_z);
		seperating_axes[14].cross(0.0, 1.0, 0.0, e23_x, e23_y, e23_z);
		seperating_axes[15].cross(0.0, 0.0, 1.0, e23_x, e23_y, e23_z);
		seperating_axes[16].cross(1.0, 0.0, 0.0, e24_x, e24_y, e24_z);
		seperating_axes[17].cross(0.0, 1.0, 0.0, e24_x, e24_y, e24_z);
		seperating_axes[18].cross(0.0, 0.0, 1.0, e24_x, e24_y, e24_z);
		seperating_axes[19].cross(1.0, 0.0, 0.0, e34_x, e34_y, e34_z);
		seperating_axes[20].cross(0.0, 1.0, 0.0, e34_x, e34_y, e34_z);
		seperating_axes[21].cross(0.0, 0.0, 1.0, e34_x, e34_y, e34_z);
	}
};

#endif