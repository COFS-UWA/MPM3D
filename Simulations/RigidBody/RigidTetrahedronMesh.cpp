#include "Simulations_pcp.h"

#include <unordered_map>
#include "Geometry.h"
#include "TetrahedronUtils.h"

#include "RigidTetrahedronMesh.h"

namespace
{
	typedef RigidTetrahedronMesh_Internal::Face Face;

	struct FaceMapItem
	{
		size_t n1, n2, n3, elem_num;
		FaceMapItem(size_t _n1, size_t _n2, size_t _n3) :
			n1(_n1), n2(_n2), n3(_n3), elem_num(0) {}
	};

	class FaceMap
	{
	public:
		FaceMap() {}
		~FaceMap() {}

		void add_face(size_t n1, size_t n2, size_t n3);
		Face* output_boundary_face(size_t& bface_num);

	protected:
		typedef std::unordered_map<std::string, FaceMapItem> Map;
		Map map;
	};

	void FaceMap::add_face(size_t n1, size_t n2, size_t n3)
	{
		size_t ns[3];
		ns[0] = n1;
		ns[1] = n2;
		ns[2] = n3;
		sort_array_3_acc(ns);
		if (ns[0] == 0 && ns[1] == 3 && ns[2] == 6)
		{
			int efef = 0;
		}
		char nid_str[50];
		snprintf(nid_str, sizeof(nid_str), "%zu_%zu_%zu", ns[0], ns[1], ns[2]);
		auto res = map.emplace(std::string(nid_str), FaceMapItem(n1, n2, n3));
		if (!res.second)
			++(res.first->second.elem_num);
	}

	Face* FaceMap::output_boundary_face(size_t& bface_num)
	{
		size_t fnum = 0;
		for (Map::iterator iter = map.begin(); iter != map.end(); ++iter)
		{
			FaceMapItem& fi = iter->second;
			if (fi.elem_num == 0)
				++fnum;
		}

		bface_num = fnum;
		if (fnum == 0)
			return nullptr;
		Face* res = new Face[fnum];
		Face* cur_face = res;
		size_t fid = 0;
		for (Map::iterator iter = map.begin(); iter != map.end(); ++iter)
		{
			FaceMapItem& fi = iter->second;
			if (fi.elem_num == 0)
			{
				cur_face->id = fid;
				cur_face->n1 = fi.n1;
				cur_face->n2 = fi.n2;
				cur_face->n3 = fi.n3;
				++cur_face;
				++fid;
			}
		}
		return res;
	}

}

RigidTetrahedronMesh::RigidTetrahedronMesh() :
	density(1.0),
	ax(0.0), ay(0.0), az(0.0),
	vx(0.0), vy(0.0), vz(0.0),
	pax(&ax), pay(&ay), paz(&az),
	pvx(&vx), pvy(&vy), pvz(&vz),
	fx_ext(0.0), fy_ext(0.0), fz_ext(0.0),
	grids(nullptr),
	bface_num(0), bfaces(nullptr)
{
	init_cal_var();
}

RigidTetrahedronMesh::~RigidTetrahedronMesh()
{
	clear_bfaces();
	clear_bg_grids();
}

int RigidTetrahedronMesh::init_mesh(
	const char* file_name,
	double dx,
	double dy,
	double dz
	)
{
	int res = load_mesh_from_hdf5(file_name);
	if (res)
		return res;
	move_mesh(dx, dy, dz);

	extract_bfaces();

	// init bg_mesh

	return res;
}

void RigidTetrahedronMesh::move_mesh(
	double dx,
	double dy,
	double dz
	)
{
	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node &n = nodes[n_id];
		n.x += dx;
		n.y += dy;
		n.z += dz;
	}

	bounding_box.xl += dx;
	bounding_box.xu += dx;
	bounding_box.yl += dy;
	bounding_box.yu += dy;
	bounding_box.zl += dz;
	bounding_box.zu += dz;

	centre.x += dx;
	centre.y += dy;
	centre.z += dz;
}

void RigidTetrahedronMesh::extract_bfaces()
{
	clear_bfaces();

	FaceMap face_map;
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		Element& e = elems[e_id];
		face_map.add_face(e.n1, e.n2, e.n3);
		face_map.add_face(e.n1, e.n4, e.n2);
		face_map.add_face(e.n2, e.n4, e.n3);
		face_map.add_face(e.n1, e.n3, e.n4);
	}

	bfaces = face_map.output_boundary_face(bface_num);
}

int RigidTetrahedronMesh::init_bg_grids(
	double _g_h,
	double expand_size
	)
{
	clear_bg_grids();

	g_h = _g_h;

	double xlen, ylen, zlen;
	g_bbox.xl = bounding_box.xl - expand_size;
	g_bbox.xu = bounding_box.xu + expand_size;
	xlen = g_bbox.xu - g_bbox.xl;
	g_bbox.yl = bounding_box.yl - expand_size;
	g_bbox.yu = bounding_box.yu + expand_size;
	ylen = g_bbox.yu - g_bbox.yl;
	g_bbox.zl = bounding_box.zl - expand_size;
	g_bbox.zu = bounding_box.zu + expand_size;
	zlen = g_bbox.zu - g_bbox.zl;

	double pad_len;
	// x
	g_x_num = size_t(ceil(xlen / g_h));
	pad_len = (double(g_x_num) * g_h - xlen) * 0.5;
	g_bbox.xl -= pad_len;
	g_bbox.xu += pad_len;
	// y
	g_y_num = size_t(ceil(ylen / g_h));
	pad_len = (double(g_y_num) * g_h - ylen) * 0.5;
	g_bbox.yl -= pad_len;
	g_bbox.yu += pad_len;
	// z
	g_z_num = size_t(ceil(zlen / g_h));
	pad_len = (double(g_z_num) * g_h - zlen) * 0.5;
	g_bbox.zl -= pad_len;
	g_bbox.zu += pad_len;
	//
	g_xy_num = g_x_num * g_y_num;
	g_num = g_xy_num * g_z_num;
	grids = new Grid[g_num];
	Grid* cur_pg = grids;
	for (size_t z_id = 0; z_id < g_z_num; ++z_id)
		for (size_t y_id = 0; y_id < g_y_num; ++y_id)
			for (size_t x_id = 0; x_id < g_x_num; ++x_id)
			{
				cur_pg->x_id = x_id;
				cur_pg->y_id = y_id;
				cur_pg->z_id = z_id;
				cur_pg->pos_type = PosType::Outside;
				cur_pg->bfaces = nullptr;
				++cur_pg;
			}

	Cube e_bbox;
	size_t xl_id, xu_id, yl_id, yu_id, zl_id, zu_id;
	Cube grid_box;
	// apply all elements to grids
	for (size_t e_id = 0; e_id < elem_num; ++e_id)
	{
		Element& e = elems[e_id];
		Node& n1 = nodes[e.n1];
		Node& n2 = nodes[e.n2];
		Node& n3 = nodes[e.n3];
		Node& n4 = nodes[e.n4];
		e_bbox = get_tetrahedron_bounding_box(n1, n2, n3, n4);
		xl_id = size_t(floor((e_bbox.xl - g_bbox.xl) / g_h));
		xu_id = size_t(ceil((e_bbox.xu - g_bbox.xl) / g_h));
		yl_id = size_t(floor((e_bbox.yl - g_bbox.yl) / g_h));
		yu_id = size_t(ceil((e_bbox.yu - g_bbox.yl) / g_h));
		zl_id = size_t(floor((e_bbox.zl - g_bbox.zl) / g_h));
		zu_id = size_t(ceil((e_bbox.zu - g_bbox.zl) / g_h));
		for (size_t z_id = zl_id; z_id < zu_id; ++z_id)
			for (size_t y_id = yl_id; y_id < yu_id; ++y_id)
				for (size_t x_id = xl_id; x_id < xu_id; ++x_id)
				{
					Grid& g = grid_by_id(x_id, y_id, z_id);
					grid_box = grid_box_by_id(x_id, y_id, z_id);
					if (test_tetrahedron_aabb_intersection(
						n1, n2, n3, n4, grid_box))
						g.pos_type = PosType::Inside;
				}
	}

	// apply all faces to grids
	for (size_t f_id = 0; f_id < bface_num; ++f_id)
	{
		Face& f = bfaces[f_id];
		Node& n1 = nodes[f.n1];
		Node& n2 = nodes[f.n2];
		Node& n3 = nodes[f.n3];
		e_bbox = get_3Dtriangle_bounding_box(n1, n2, n3);
		xl_id = size_t(floor((e_bbox.xl - g_bbox.xl) / g_h));
		xu_id = size_t(ceil((e_bbox.xu - g_bbox.xl) / g_h));
		yl_id = size_t(floor((e_bbox.yl - g_bbox.yl) / g_h));
		yu_id = size_t(ceil((e_bbox.yu - g_bbox.yl) / g_h));
		zl_id = size_t(floor((e_bbox.zl - g_bbox.zl) / g_h));
		zu_id = size_t(ceil((e_bbox.zu - g_bbox.zl) / g_h));
		for (size_t z_id = zl_id; z_id < zu_id; ++z_id)
			for (size_t y_id = yl_id; y_id < yu_id; ++y_id)
				for (size_t x_id = xl_id; x_id < xu_id; ++x_id)
				{
					Grid& g = grid_by_id(x_id, y_id, z_id);
					grid_box = grid_box_by_id(x_id, y_id, z_id);
					if (test_3Dtriangle_aabb_intersection(
						n1, n2, n3, grid_box))
					{
						g.pos_type = PosType::AtBoundary;
						add_bface_to_grid(g, f);
					}
				}
	}

	return 0;
}
