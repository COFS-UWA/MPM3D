#ifndef __Rigid_Tetrahedron_Mesh_h__
#define __Rigid_Tetrahedron_Mesh_h__

#include "ItemBuffer.hpp"
#include "ItemStack.hpp"
#include "TetrahedronUtils.h"
#include "TetrahedronMeshTemplate.hpp"

// Usage:
// 1. init_mesh
// 2. init_bg_grids
namespace RigidTetrahedronMesh_Internal
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

	struct Face
	{
		size_t id;
		size_t n1, n2, n3;
		PointToTriangleDistance<Node> pt_tri_dist;
	};
}

class RigidTetrahedronMesh :
	public TetrahedronMeshTemplate<RigidTetrahedronMesh_Internal::Node,
								   RigidTetrahedronMesh_Internal::Element,
								   RigidTetrahedronMesh_Internal::Edge>
{
public:
	typedef RigidTetrahedronMesh_Internal::Node Node;
	typedef RigidTetrahedronMesh_Internal::Element Element;
	typedef RigidTetrahedronMesh_Internal::Edge Edge;
	typedef RigidTetrahedronMesh_Internal::Face Face;

	enum class PosType : unsigned char
	{
		Inside = 0,
		AtBoundary = 1,
		Outside = 2
	};

	struct FacePointer
	{
		Face* pface;
		FacePointer *next;
	};

	struct Grid
	{
		size_t x_id, y_id, z_id;
		PosType pos_type;
		bool close_to_boundary;
		FacePointer *bfaces;
	};

protected:
	double density;

	double x, y, z;
	double m;
	
	double ax, ay, az;
	double vx, vy, vz;
	//double ax_ang, ay_ang, az_ang;
	//double vx_ang, vy_ang, vz_ang;

	double fx_ext, fy_ext, fz_ext;
	double ax_bc, ay_bc, az_bc;
	double vx_bc, vy_bc, vz_bc;

	double *pax, *pay, *paz;
	double *pvx, *pvy, *pvz;

	// for precision, x = x_ori + ux ...
	double x_ori, ux;
	double y_ori, uy;
	double z_ori, uz;
	// contact force
	double fx_con, fy_con, fz_con;

	inline void init_cal_var() noexcept
	{
		m = get_vol() * density;
		// moi;
		// moi_inv;
	}
	
	// boundary faces
	size_t bface_num;
	Face* bfaces;

	void clear_bfaces()
	{
		if (bfaces)
		{
			delete[] bfaces;
			bfaces = nullptr;
		}
		bface_num = 0;
	}

	void extract_bfaces();

public:
	RigidTetrahedronMesh();
	~RigidTetrahedronMesh();

	inline double get_density() { return density; }
	inline double get_m() { return m; }
	inline double get_x() { return x; }
	inline double get_y() { return y; }
	inline double get_z() { return z; }
	inline double get_ax() { return ax; }
	inline double get_ay() { return ay; }
	inline double get_az() { return az; }
	inline double get_vx() { return vx; }
	inline double get_vy() { return vy; }
	inline double get_vz() { return vz; }
	inline double get_fx_contact() { return fx_con; }
	inline double get_fy_contact() { return fy_con; }
	inline double get_fz_contact() { return fz_con; }
	inline double get_fx_ext() { return fx_ext; }
	inline double get_fy_ext() { return fy_ext; }
	inline double get_fz_ext() { return fz_ext; }
	inline bool has_ax_bc() { return pax == &ax_bc; }
	inline bool has_ay_bc() { return pay == &ay_bc; }
	inline bool has_az_bc() { return paz == &az_bc; }
	inline double get_ax_bc() { return ax_bc; }
	inline double get_ay_bc() { return ay_bc; }
	inline double get_az_bc() { return az_bc; }
	inline bool has_vx_bc() { return pvx == &vx_bc; }
	inline bool has_vy_bc() { return pvy == &vy_bc; }
	inline bool has_vz_bc() { return pvz == &vz_bc; }
	inline double get_vx_bc() { return vx_bc; }
	inline double get_vy_bc() { return vy_bc; }
	inline double get_vz_bc() { return vz_bc; }

	inline void set_ax_bc(double _a) { pax = &ax_bc; ax_bc = _a; }
	inline void set_ay_bc(double _a) { pay = &ay_bc; ay_bc = _a; }
	inline void set_az_bc(double _a) { paz = &az_bc; az_bc = _a; }
	inline void set_a_bc(double _ax, double _ay, double _az)
	{ set_ax_bc(_ax); set_ay_bc(_ay); set_az_bc(_az); }

	inline void set_vx_bc(double _v) { pvx = &vx_bc; vx_bc = _v; }
	inline void set_vy_bc(double _v) { pvy = &vy_bc; vy_bc = _v; }
	inline void set_vz_bc(double _v) { pvz = &vz_bc; vz_bc = _v; }
	inline void set_v_bc(double _vx, double _vy, double _vz)
	{ set_vx_bc(_vx); set_vy_bc(_vy); set_vz_bc(_vz); }

	inline void add_fx_ext(double _f) { fx_ext += _f; }
	inline void add_fy_ext(double _f) { fy_ext += _f; }
	inline void add_fz_ext(double _f) { fz_ext += _f; }

	inline size_t get_bface_num() const { return bface_num; }
	inline Face* get_bfaces() { return bfaces; }

	int init_mesh(const char *file_name, double dx, double dy, double dz);

	inline void init_calculation() noexcept
	{
		x_ori = x;
		y_ori = y;
		z_ori = z;
		ux = 0.0;
		uy = 0.0;
		uz = 0.0;
	}

	inline void reset() noexcept
	{
		fx_con = 0.0;
		fy_con = 0.0;
		fz_con = 0.0;
	}

	inline void update_motion(double dt)
	{
		// update a
		ax = (fx_ext + fx_con) / m;
		ay = (fy_ext + fy_con) / m;
		az = (fz_ext + fz_con) / m;
		// apply abc
		ax = *pax;
		ay = *pay;
		az = *paz;
		// update velocity
		vx += ax * dt;
		vy += ay * dt;
		vz += az * dt;
		// apply vbc
		vx = *pvx;
		vy = *pvy;
		vz = *pvz;
		// update position
		x += vx * dt;
		y += vy * dt;
		z += vz * dt;
	}

	inline void add_con_force(double fx, double fy, double fz,
					double posx, double posy, double posz) noexcept
	{
		fx_con += fx;
		fy_con += fy;
		fz_con += fz;
	}
	
	inline Point3D to_local_coord(Point3D &gp) const noexcept
	{ return Point3D(gp.x - x, gp.y - y, gp.z - z); }
	inline Point3D to_global_coord(Point3D &lp) const noexcept
	{ return Point3D(lp.x + x, lp.y + y, lp.z + z); }

	inline double get_grid_h() const noexcept { return g_h; }
	inline const Cube& get_grid_bbox() const noexcept { return g_bbox; }
	inline size_t get_grid_x_num() const noexcept { return g_x_num; }
	inline size_t get_grid_y_num() const noexcept { return g_y_num; }
	inline size_t get_grid_z_num() const noexcept { return g_z_num; }

	void clear_bg_grids();
	int init_bg_grids(double _g_h, double expand_size);

	void set_dist_max(double _dist_max);
	bool cal_distance_to_boundary(Point3D& pt,
		double& dist, double& nx, double& ny, double& nz);
	
protected: // background grid
	double g_h;
	Cube g_bbox;
	size_t g_x_num, g_y_num, g_z_num;
	size_t g_xy_num, g_num;
	Grid* grids;

	inline Grid& grid_by_id(size_t x_id, size_t y_id, size_t z_id)
	{ return grids[g_xy_num * z_id + g_x_num * y_id + x_id]; }
	inline Cube grid_box_by_id(size_t x_id, size_t y_id, size_t z_id)
	{
		Cube res;
		res.xl = g_bbox.xl + double(x_id) * g_h;
		res.xu = res.xl + g_h;
		res.yl = g_bbox.yl + double(y_id) * g_h;
		res.yu = res.yl + g_h;
		res.zl = g_bbox.zl + double(z_id) * g_h;
		res.zu = res.zl + g_h;
		return res;
	}
	
	MemoryUtils::ItemBuffer<FacePointer> fp_array;
	inline void add_bface_to_grid(Grid& g, Face& f)
	{
		FacePointer* fp = fp_array.alloc();
		fp->pface = &f;
		fp->next = g.bfaces;
		g.bfaces = fp;
	}
	
	PointInTetrahedron<Node> pt_in_teh;
	TetrahedronAABBCollisionSAT<Node> teh_aabb_collision;
	inline void init_teh_aabb_collision(Element& e)
	{
		Node& n1 = nodes[e.n1];
		Node& n2 = nodes[e.n2];
		Node& n3 = nodes[e.n3];
		Node& n4 = nodes[e.n4];
		pt_in_teh.init_tetrahedron(n1, n2, n3, n4);
		teh_aabb_collision.init_tetrahedron(n1, n2, n3, n4);
	}
	bool detect_teh_aabb_collision(Cube &box);

	TriangleAABBCollisionSAT<Node> tri_aabb_collision;
	void init_tri_aabb_collision(Face& f)
	{
		Node& n1 = nodes[f.n1];
		Node& n2 = nodes[f.n2];
		Node& n3 = nodes[f.n3];
		tri_aabb_collision.init_triangle(n1, n2, n3);
	}
	bool detect_tri_aabb_collision(Cube& box);

	// search cloest face with bg grids
	struct IdDistPair
	{
		size_t id;
		double dist;
	};

	inline void swap_id_dist_pair(IdDistPair &idp1, IdDistPair &idp2)
	{
		IdDistPair tmp;
		tmp.id = idp1.id;
		tmp.dist = idp1.dist;
		idp1.id = idp2.id;
		idp1.dist = idp2.dist;
		idp2.id = tmp.id;
		idp2.dist = tmp.dist;
	}

	void sort_acc_id_dist_pairs_8(IdDistPair *id_pairs)
	{
		size_t min_id;
		for (size_t i = 0; i < 7; ++i)
		{
			min_id = i;
			for (size_t j = i + 1; j < 8; ++j)
			{
				if (id_pairs[j].dist < id_pairs[min_id].dist)
					min_id = j;
			}
			if (min_id != i)
				swap_id_dist_pair(id_pairs[i], id_pairs[min_id]);
		}
	}

	double dist_max;
	unsigned char height_max;
	long long id_dist_max, id_stride_max;
	double stride_max;

	// variables for search_closest_face
	struct SearchClosestFaceParam
	{
		long long xl_id, yl_id, zl_id;
		Cube box;
	};

	IdCube id_range;
	Point3D cur_pt;
	double cur_dist;
	unsigned char cur_height;
	long long cur_id_stride;
	Face* cur_face;
	unsigned char cur_norm_type;

	// acceleration searching
	void init_close_enough_to_boundary();
	void search_closest_face(SearchClosestFaceParam &param);
};

#endif