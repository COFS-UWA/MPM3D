#include "Simulations_pcp.h"

#include <unordered_map>
#include "GeometryUtils.h"
#include "Geometry3D.h"

#include "RigidTetrahedronMesh.h"

namespace
{
	typedef RigidTetrahedronMesh_Internal::Face Face;

	class FaceMap
	{
	public:
		FaceMap() {}
		~FaceMap() {}

		void add_face(size_t n1, size_t n2, size_t n3);
		Face* output_boundary_face(size_t& bface_num);

	protected:
		struct FaceMapItem
		{
			size_t n1, n2, n3, elem_num;
			FaceMapItem(size_t _n1, size_t _n2, size_t _n3) :
				n1(_n1), n2(_n2), n3(_n3), elem_num(0) {}
		};

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

	x = centre.x;
	y = centre.y;
	z = centre.z;

	for (size_t n_id = 0; n_id < node_num; ++n_id)
	{
		Node& n = nodes[n_id];
		n.x -= x;
		n.y -= y;
		n.z -= z;
	}

	bounding_box.xl -= x;
	bounding_box.xu -= x;
	bounding_box.yl -= y;
	bounding_box.yu -= y;
	bounding_box.zl -= z;
	bounding_box.zu -= z;

	centre.x = 0.0;
	centre.y = 0.0;
	centre.z = 0.0;

	x += dx;
	y += dy;
	z += dz;

	extract_bfaces();

	// init bg_mesh

	return res;
}

void RigidTetrahedronMesh::set_init_state(
	double _density,
	double _fx_contact,
	double _fy_contact,
	double _fz_contact,
	double _ax, double _ay, double _az,
	double _vx, double _vy, double _vz,
	double _x, double _y, double _z
	)
{
	density = _density;
	fx_con = _fx_contact;
	fy_con = _fy_contact;
	fz_con = _fz_contact;
	ax = _ax;
	ay = _ay;
	az = _az;
	vx = _vx;
	vy = _vy;
	vz = _vz;
	x = _x;
	y = _y;
	z = _z;
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
	for (size_t bf_id = 0; bf_id < bface_num; ++bf_id)
	{
		Face &f = bfaces[bf_id];
		Node& n1 = nodes[f.n1];
		Node& n2 = nodes[f.n2];
		Node& n3 = nodes[f.n3];
		f.pt_tri_dist.init_triangle(n1, n2, n3);
	}

}

void RigidTetrahedronMesh::clear_bg_grids()
{
	if (grids)
	{
		delete[] grids;
		grids = nullptr;
	}
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
				Grid& g = *cur_pg;
				g.x_id = x_id;
				g.y_id = y_id;
				g.z_id = z_id;
				g.pos_type = PosType::Outside;
				g.close_to_boundary = true;
				g.bfaces = nullptr;
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
		init_teh_aabb_collision(e);
		for (size_t z_id = zl_id; z_id < zu_id; ++z_id)
			for (size_t y_id = yl_id; y_id < yu_id; ++y_id)
				for (size_t x_id = xl_id; x_id < xu_id; ++x_id)
				{
					Grid& g = grid_by_id(x_id, y_id, z_id);
					grid_box = grid_box_by_id(x_id, y_id, z_id);
					if (detect_teh_aabb_collision(grid_box))
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
		init_tri_aabb_collision(f);
		for (size_t z_id = zl_id; z_id < zu_id; ++z_id)
			for (size_t y_id = yl_id; y_id < yu_id; ++y_id)
				for (size_t x_id = xl_id; x_id < xu_id; ++x_id)
				{
					Grid& g = grid_by_id(x_id, y_id, z_id);
					grid_box = grid_box_by_id(x_id, y_id, z_id);
					if (detect_tri_aabb_collision(grid_box))
					{
						g.pos_type = PosType::AtBoundary;
						add_bface_to_grid(g, f);
					}
				}
	}

	dist_max = g_bbox.xu - g_bbox.xl;
	id_dist_max = g_x_num;
	if (id_dist_max < g_y_num)
	{
		dist_max = g_bbox.yu - g_bbox.yl;
		id_dist_max = g_x_num;
	}
	if (id_dist_max < g_z_num)
	{
		dist_max = g_bbox.zu - g_bbox.zl;
		id_dist_max = g_z_num;
	}
	const long long id_range = id_dist_max + id_dist_max + 1;
	height_max = 0;
	id_stride_max = 1;
	while (id_stride_max < id_range)
	{
		++height_max;
		id_stride_max = id_stride_max << 1;
	}
	stride_max = double(id_stride_max) * g_h;
	
	return 0;
}

bool RigidTetrahedronMesh::detect_teh_aabb_collision(Cube& box)
{
	// whether tetrahedron nodes locate in box
	// efficient when tetrahedron is much smaller than grid
	if (box.is_in_box(teh_aabb_collision.get_n1()) ||
		box.is_in_box(teh_aabb_collision.get_n2()) ||
		box.is_in_box(teh_aabb_collision.get_n3()) ||
		box.is_in_box(teh_aabb_collision.get_n4()))
		return true;

	// whether box corners locate in tetrahedron
	// efficient when grid is much smaller than tetrahedron
	if (pt_in_teh.is_in_tetrahedron(box.xl, box.yl, box.zl) ||
		pt_in_teh.is_in_tetrahedron(box.xl, box.yl, box.zu) ||
		pt_in_teh.is_in_tetrahedron(box.xl, box.yu, box.zl) ||
		pt_in_teh.is_in_tetrahedron(box.xl, box.yu, box.zu) ||
		pt_in_teh.is_in_tetrahedron(box.xu, box.yl, box.zl) ||
		pt_in_teh.is_in_tetrahedron(box.xu, box.yl, box.zu) ||
		pt_in_teh.is_in_tetrahedron(box.xu, box.yu, box.zl) ||
		pt_in_teh.is_in_tetrahedron(box.xu, box.yu, box.zu))
		return true;

	return teh_aabb_collision.detect_collision_with_cube(box);
}

bool RigidTetrahedronMesh::detect_tri_aabb_collision(Cube& box)
{
	// whether tetrahedron nodes locate in box
	// efficient when tetrahedron is much smaller than grid
	if (box.is_in_box(tri_aabb_collision.get_n1()) ||
		box.is_in_box(tri_aabb_collision.get_n2()) ||
		box.is_in_box(tri_aabb_collision.get_n3()))
		return true;
	
	return tri_aabb_collision.detect_collision_with_cube(box);
}

void RigidTetrahedronMesh::set_dist_max(double _dist_max)
{
	dist_max = _dist_max > g_h ? _dist_max : g_h;
	id_dist_max = long long(ceil(_dist_max / g_h));
	const long long id_range = id_dist_max + id_dist_max + 1;
	height_max = 0;
	id_stride_max = 1;
	while (id_stride_max < id_range)
	{
		++height_max;
		id_stride_max = id_stride_max << 1;
	}
	stride_max = double(id_stride_max) * g_h;

	init_close_enough_to_boundary();
}

void RigidTetrahedronMesh::init_close_enough_to_boundary()
{
	for (size_t g_id = 0; g_id < g_num; ++g_id)
	{
		Grid& g = grids[g_id];
		g.close_to_boundary = true;
	}

	struct GridNode { bool close_to_boundary; };

	size_t gn_x_num = g_x_num + 1;
	size_t gn_y_num = g_y_num + 1;
	size_t gn_z_num = g_z_num + 1;
	size_t gn_xy_num = gn_x_num * gn_y_num;
	GridNode *gns = new GridNode[gn_x_num * gn_y_num * gn_z_num];
	GridNode* cur_pgn = gns;
	Point3D gn_pt;
	double dist, nx, ny, nz;
	gn_pt.z = g_bbox.zl;
	for (size_t z_id = 0; z_id < gn_z_num; ++z_id)
	{
		gn_pt.y = g_bbox.yl;
		for (size_t y_id = 0; y_id < gn_y_num; ++y_id)
		{
			gn_pt.x = g_bbox.xl;
			for (size_t x_id = 0; x_id < gn_x_num; ++x_id)
			{
				cur_pgn->close_to_boundary
					= cal_dist_and_dir_to_pt_internal(
						gn_pt,
						dist,
						nx,
						ny,
						nz
						);
				++cur_pgn;
				gn_pt.x += g_h;
			}
			gn_pt.y += g_h;
		}
		gn_pt.z += g_h;
	}

	Grid* cur_pg = grids;
	GridNode *gn_z_start = gns, *gn_y_start;
	for (size_t z_id = 0; z_id < g_z_num; ++z_id)
	{
		gn_y_start = gn_z_start;
		for (size_t y_id = 0; y_id < g_y_num; ++y_id)
		{
			cur_pgn = gn_y_start;
			for (size_t x_id = 0; x_id < g_x_num; ++x_id)
			{
				GridNode& n1 = *cur_pgn;
				GridNode& n2 = *(cur_pgn + 1);
				GridNode& n3 = *(cur_pgn + gn_x_num);
				GridNode& n4 = *(&n3 + 1);
				GridNode& n5 = *(cur_pgn + gn_xy_num);
				GridNode& n6 = *(&n5 + 1);
				GridNode& n7 = *(&n5 + gn_x_num);
				GridNode& n8 = *(&n7 + 1);
				if (n1.close_to_boundary == false &&
					n2.close_to_boundary == false &&
					n3.close_to_boundary == false &&
					n4.close_to_boundary == false &&
					n5.close_to_boundary == false &&
					n6.close_to_boundary == false &&
					n7.close_to_boundary == false &&
					n8.close_to_boundary == false)
					cur_pg->close_to_boundary = false;
				++cur_pg;
				++cur_pgn;
			}
			gn_y_start += gn_x_num;
		}
		gn_z_start += gn_xy_num;
	}

	delete[] gns;
}

bool RigidTetrahedronMesh::cal_dist_and_dir_to_pt_internal(
	Point3D &pt,
	double& dist, 
	double& nx,
	double& ny,
	double& nz
	)
{
	if (!g_bbox.is_in_box(pt))
		return false;

	long long ptx_id = long long((pt.x - g_bbox.xl) / g_h);
	long long pty_id = long long((pt.y - g_bbox.yl) / g_h);
	long long ptz_id = long long((pt.z - g_bbox.zl) / g_h);
	Grid& g = grid_by_id(ptx_id, pty_id, ptz_id);
	if (!g.close_to_boundary)
		return false;
	
	cur_pt = pt;
	cur_dist = dist_max;
	cur_height = height_max;
	cur_id_stride = id_stride_max;
	cur_face = nullptr;

	id_range.xl_id = ptx_id - id_dist_max;
	id_range.yl_id = pty_id - id_dist_max;
	id_range.zl_id = ptz_id - id_dist_max;
	id_range.xu_id = ptx_id + id_dist_max + 1;
	id_range.yu_id = pty_id + id_dist_max + 1;
	id_range.zu_id = ptz_id + id_dist_max + 1;
	SearchClosestFaceParam param;
	param.xl_id = id_range.xl_id;
	param.yl_id = id_range.yl_id;
	param.zl_id = id_range.zl_id;
	if (id_range.xl_id < 0)
		id_range.xl_id = 0;
	if (id_range.xu_id > g_x_num)
		id_range.xu_id = g_x_num;
	if (id_range.yl_id < 0)
		id_range.yl_id = 0;
	if (id_range.yu_id > g_y_num)
		id_range.yu_id = g_y_num;
	if (id_range.zl_id < 0)
		id_range.zl_id = 0;
	if (id_range.zu_id > g_z_num)
		id_range.zu_id = g_z_num;
	Cube &box = param.box;
	box.xl = g_bbox.xl + double(param.xl_id) * g_h;
	box.xu = box.xl + stride_max;
	box.yl = g_bbox.yl + double(param.yl_id) * g_h;
	box.yu = box.yl + stride_max;
	box.zl = g_bbox.zl + double(param.zl_id) * g_h;
	box.zu = box.zl + stride_max;

	search_closest_face(param);

	if (cur_face)
	{
		dist = cur_dist;
		Vector3D fnormal;
		cur_face->pt_tri_dist.cal_normal_to_point(
			cur_pt,
			cur_norm_type,
			fnormal
			);
		if (dist >= 0.0)
			fnormal.reverse();
		nx = fnormal.x;
		ny = fnormal.y;
		nz = fnormal.z;
		return true;
	}
	return false;
}

void RigidTetrahedronMesh::search_closest_face(SearchClosestFaceParam &param)
{
	IdCube my_id_range;
	my_id_range.xl_id = param.xl_id;
	my_id_range.yl_id = param.yl_id;
	my_id_range.zl_id = param.zl_id;
	my_id_range.xu_id = my_id_range.xl_id + cur_id_stride;
	my_id_range.yu_id = my_id_range.yl_id + cur_id_stride;
	my_id_range.zu_id = my_id_range.zl_id + cur_id_stride;
	if (id_range.does_not_overlap(my_id_range))
		return;

	if (cur_height)
	{
		--cur_height;
		cur_id_stride = cur_id_stride >> 1;
		long long xm_id = param.xl_id + cur_id_stride;
		long long ym_id = param.yl_id + cur_id_stride;
		long long zm_id = param.zl_id + cur_id_stride;
		Cube& box = param.box;
		double box_xm = (box.xl + box.xu) * 0.5;
		double box_ym = (box.yl + box.yu) * 0.5;
		double box_zm = (box.zl + box.zu) * 0.5;
		IdDistPair id_pairs[8];
		SearchClosestFaceParam params[8];
		// child 0
		SearchClosestFaceParam& param0 = params[0];
		param0.xl_id = param.xl_id;
		param0.yl_id = param.yl_id;
		param0.zl_id = param.zl_id;
		Cube &box0 = param0.box;
		box0.xl = box.xl;
		box0.xu = box_xm;
		box0.yl = box.yl;
		box0.yu = box_ym;
		box0.zl = box.zl;
		box0.zu = box_zm;
		IdDistPair& pair0 = id_pairs[0];
		pair0.id = 0;
		pair0.dist = cal_cube_point_distance(box0, cur_pt);
		// child 1
		SearchClosestFaceParam& param1 = params[1];
		param1.xl_id = xm_id;
		param1.yl_id = param.yl_id;
		param1.zl_id = param.zl_id;
		Cube &box1 = param1.box;
		box1.xl = box_xm;
		box1.xu = box.xu;
		box1.yl = box.yl;
		box1.yu = box_ym;
		box1.zl = box.zl;
		box1.zu = box_zm;
		IdDistPair& pair1 = id_pairs[1];
		pair1.id = 1;
		pair1.dist = cal_cube_point_distance(box1, cur_pt);
		// child 2
		SearchClosestFaceParam& param2 = params[2];
		param2.xl_id = param.xl_id;
		param2.yl_id = ym_id;
		param2.zl_id = param.zl_id;
		Cube &box2 = param2.box;
		box2.xl = box.xl;
		box2.xu = box_xm;
		box2.yl = box_ym;
		box2.yu = box.yu;
		box2.zl = box.zl;
		box2.zu = box_zm;
		IdDistPair& pair2 = id_pairs[2];
		pair2.id = 2;
		pair2.dist = cal_cube_point_distance(box2, cur_pt);
		// child 3
		SearchClosestFaceParam& param3 = params[3];
		param3.xl_id = xm_id;
		param3.yl_id = ym_id;
		param3.zl_id = param.zl_id;
		Cube& box3 = param3.box;
		box3.xl = box_xm;
		box3.xu = box.xu;
		box3.yl = box_ym;
		box3.yu = box.yu;
		box3.zl = box.zl;
		box3.zu = box_zm;
		IdDistPair& pair3 = id_pairs[3];
		pair3.id = 3;
		pair3.dist = cal_cube_point_distance(box3, cur_pt);
		// child 4
		SearchClosestFaceParam& param4 = params[4];
		param4.xl_id = param.xl_id;
		param4.yl_id = param.yl_id;
		param4.zl_id = zm_id;
		Cube& box4 = param4.box;
		box4.xl = box.xl;
		box4.xu = box_xm;
		box4.yl = box.yl;
		box4.yu = box_ym;
		box4.zl = box_zm;
		box4.zu = box.zu;
		IdDistPair& pair4 = id_pairs[4];
		pair4.id = 4;
		pair4.dist = cal_cube_point_distance(box4, cur_pt);
		// child 5
		SearchClosestFaceParam& param5 = params[5];
		param5.xl_id = xm_id;
		param5.yl_id = param.yl_id;
		param5.zl_id = zm_id;
		Cube& box5 = param5.box;
		box5.xl = box_xm;
		box5.xu = box.xu;
		box5.yl = box.yl;
		box5.yu = box_ym;
		box5.zl = box_zm;
		box5.zu = box.zu;
		IdDistPair& pair5 = id_pairs[5];
		pair5.id = 5;
		pair5.dist = cal_cube_point_distance(box5, cur_pt);
		// child 6
		SearchClosestFaceParam& param6 = params[6];
		param6.xl_id = param.xl_id;
		param6.yl_id = ym_id;
		param6.zl_id = zm_id;
		Cube& box6 = param6.box;
		box6.xl = box.xl;
		box6.xu = box_xm;
		box6.yl = box_ym;
		box6.yu = box.yu;
		box6.zl = box_zm;
		box6.zu = box.zu;
		IdDistPair& pair6 = id_pairs[6];
		pair6.id = 6;
		pair6.dist = cal_cube_point_distance(box6, cur_pt);
		// child 7
		SearchClosestFaceParam& param7 = params[7];
		param7.xl_id = xm_id;
		param7.yl_id = ym_id;
		param7.zl_id = zm_id;
		Cube& box7 = param7.box;
		box7.xl = box_xm;
		box7.xu = box.xu;
		box7.yl = box_ym;
		box7.yu = box.yu;
		box7.zl = box_zm;
		box7.zu = box.zu;
		IdDistPair& pair7 = id_pairs[7];
		pair7.id = 7;
		pair7.dist = cal_cube_point_distance(box7, cur_pt);

		sort_acc_id_dist_pairs_8(id_pairs);

		// 0
		if (pair0.dist > abs(cur_dist))
			goto complete_searching_child;
		search_closest_face(params[pair0.id]);
		// 1
		if (pair1.dist > abs(cur_dist))
			goto complete_searching_child;
		search_closest_face(params[pair1.id]);
		// 2
		if (pair2.dist > abs(cur_dist))
			goto complete_searching_child;
		search_closest_face(params[pair2.id]);
		// 3
		if (pair3.dist > abs(cur_dist))
			goto complete_searching_child;
		search_closest_face(params[pair3.id]);
		// 4
		if (pair4.dist > abs(cur_dist))
			goto complete_searching_child;
		search_closest_face(params[pair4.id]);
		// 5
		if (pair5.dist > abs(cur_dist))
			goto complete_searching_child;
		search_closest_face(params[pair5.id]);
		// 6
		if (pair6.dist > abs(cur_dist))
			goto complete_searching_child;
		search_closest_face(params[pair6.id]);
		// 7
		if (pair7.dist > abs(cur_dist))
			goto complete_searching_child;
		search_closest_face(params[pair7.id]);

	complete_searching_child:
		++cur_height;
		cur_id_stride = cur_id_stride << 1;
		return;
	}

	double dist_tmp;
	unsigned char norm_type_tmp;
	Grid& g = grid_by_id(param.xl_id, param.yl_id, param.zl_id);
	for (FacePointer* fp = g.bfaces; fp; fp = fp->next)
	{
		Face& f = *(fp->pface);
		// cal dist and dist_tye
		norm_type_tmp = f.pt_tri_dist.cal_distance_to_point(cur_pt, dist_tmp);
		if (abs(cur_dist) > abs(dist_tmp))
		{
			cur_dist = dist_tmp;
			cur_face = &f;
			cur_norm_type = norm_type_tmp;
		}
	}
}
