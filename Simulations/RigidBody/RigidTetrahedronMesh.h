#ifndef __Rigid_Tetrahedron_Mesh_h__
#define __Rigid_Tetrahedron_Mesh_h__

#include "ItemArray.hpp"
#include "TetrahedronUtils.h"
#include "TetrahedronMeshTemplate.hpp"

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

	void move_mesh(double dx, double dy, double dz);

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

	inline int distance_from_boundary(Point3D &p, double dist_max,
				double &dist, double &nx, double &ny, double &nz)
	{
		Point3D lp;

		return 0;
	}
	inline Point3D to_local_coord(Point3D &gp) const noexcept
	{
		Point3D lp;
		lp.x = gp.x - centre.x;
		lp.y = gp.y - centre.y;
		lp.z = gp.z - centre.z;
		return lp;
	}
	inline Point3D to_global_coord(Point3D &lp) const noexcept
	{
		Point3D gp;
		gp.x = lp.x + centre.x;
		gp.y = lp.y + centre.y;
		gp.z = lp.z + centre.z;
		return gp;
	}

protected: // background grid
	double g_h;
	Cube g_bbox;
	size_t g_x_num, g_y_num, g_z_num;
	size_t g_xy_num, g_num;
	Grid* grids;

	inline double get_g_h() { return g_h; }
	inline size_t get_grid_num() { return g_num; }
	inline Grid* get_grids() { return grids; }
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

	MemoryUtils::ItemArray<FacePointer> fp_array;
	inline void add_bface_to_grid(Grid& g, Face& f)
	{
		FacePointer* fp = fp_array.alloc();
		fp->pface = &f;
		fp->next = g.bfaces;
		g.bfaces = fp;
	}

	TetrahedronAABBCollisionSAT<Node> teh_aabb_collision;
	inline void init_teh_aabb_collision(Element& e)
	{
		Node& n1 = nodes[e.n1];
		Node& n2 = nodes[e.n2];
		Node& n3 = nodes[e.n3];
		Node& n4 = nodes[e.n4];
		teh_aabb_collision.init(n1, n2, n3, n4);
	}
	bool detect_teh_aabb_collision(Cube &box);

	TriangleAABBCollisionSAT<Node> tri_aabb_collision;
	void init_tri_aabb_collision(Face& f)
	{
		Node& n1 = nodes[f.n1];
		Node& n2 = nodes[f.n2];
		Node& n3 = nodes[f.n3];
		tri_aabb_collision.init(n1, n2, n3);
	}
	bool detect_tri_aabb_collision(Cube& box);

	void clear_bg_grids();
	int init_bg_grids(double _g_h, double expand_size);
};

#endif