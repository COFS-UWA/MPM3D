#ifndef __Rigid_Tetrahedron_Mesh_h__
#define __Rigid_Tetrahedron_Mesh_h__

#include <Eigen/Dense>

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

	inline void clear_bfaces()
	{
		if (bfaces)
		{
			delete[] bfaces;
			bfaces = nullptr;
		}
		bface_num = 0;
	}

	void extract_bfaces();

	// rotation
	double ax_ang, ay_ang, az_ang;
	double vx_ang, vy_ang, vz_ang;

	double ax_ang_bc, ay_ang_bc, az_ang_bc;
	union
	{
		struct { double vx_ang_bc, vy_ang_bc, vz_ang_bc; };
		Vector3D v_ang;
	};

	double* pax_ang, * pay_ang, * paz_ang;
	double* pvx_ang, * pvy_ang, * pvz_ang;

	double mx_ext, my_ext, mz_ext;
	double mx_con, my_con, mz_con;
	Vector3D ix, iy, iz;

	typedef Eigen::Matrix<double, 3, 3> Matrix3x3;
	typedef Eigen::Matrix<double, 3, 1> Vector3;

	Matrix3x3 moi_mat; // moment of inertia

public:
	RigidTetrahedronMesh();
	~RigidTetrahedronMesh();

	inline double get_density() { return density; }
	inline double get_m() const noexcept { return m; }
	inline double get_x() const noexcept { return x; }
	inline double get_y() const noexcept { return y; }
	inline double get_z() const noexcept { return z; }
	inline double get_ax() const noexcept { return ax; }
	inline double get_ay() const noexcept { return ay; }
	inline double get_az() const noexcept { return az; }
	inline double get_vx() const noexcept { return vx; }
	inline double get_vy() const noexcept { return vy; }
	inline double get_vz() const noexcept { return vz; }
	inline double get_fx_ext() const noexcept { return fx_ext; }
	inline double get_fy_ext() const noexcept { return fy_ext; }
	inline double get_fz_ext() const noexcept { return fz_ext; }
	inline double get_fx_contact() const noexcept { return fx_con; }
	inline double get_fy_contact() const noexcept { return fy_con; }
	inline double get_fz_contact() const noexcept { return fz_con; }
	inline bool has_ax_bc() const noexcept { return pax == &ax_bc; }
	inline bool has_ay_bc() const noexcept { return pay == &ay_bc; }
	inline bool has_az_bc() const noexcept { return paz == &az_bc; }
	inline double get_ax_bc() const noexcept { return ax_bc; }
	inline double get_ay_bc() const noexcept { return ay_bc; }
	inline double get_az_bc() const noexcept { return az_bc; }
	inline bool has_vx_bc() const noexcept { return pvx == &vx_bc; }
	inline bool has_vy_bc() const noexcept { return pvy == &vy_bc; }
	inline bool has_vz_bc() const noexcept { return pvz == &vz_bc; }
	inline double get_vx_bc() const noexcept { return vx_bc; }
	inline double get_vy_bc() const noexcept { return vy_bc; }
	inline double get_vz_bc() const noexcept { return vz_bc; }

	inline double get_ax_ang() const noexcept { return ax_ang; }
	inline double get_ay_ang() const noexcept { return ay_ang; }
	inline double get_az_ang() const noexcept { return az_ang; }
	inline double get_vx_ang() const noexcept { return vx_ang; }
	inline double get_vy_ang() const noexcept { return vy_ang; }
	inline double get_vz_ang() const noexcept { return vz_ang; }
	inline double get_mx_ext() const noexcept { return mx_ext; }
	inline double get_my_ext() const noexcept { return my_ext; }
	inline double get_mz_ext() const noexcept { return mz_ext; }
	inline double get_mx_contact() const noexcept { return mx_con; }
	inline double get_my_contact() const noexcept { return my_con; }
	inline double get_mz_contact() const noexcept { return mz_con; }
	inline bool has_ax_ang_bc() const noexcept { return pax_ang == &ax_ang_bc; }
	inline bool has_ay_ang_bc() const noexcept { return pay_ang == &ay_ang_bc; }
	inline bool has_az_ang_bc() const noexcept { return paz_ang == &az_ang_bc; }
	inline double get_ax_ang_bc() const noexcept { return ax_ang_bc; }
	inline double get_ay_ang_bc() const noexcept { return ay_ang_bc; }
	inline double get_az_ang_bc() const noexcept { return az_ang_bc; }
	inline bool has_vx_ang_bc() const noexcept { return pvx_ang == &vx_ang_bc; }
	inline bool has_vy_ang_bc() const noexcept { return pvy_ang == &vy_ang_bc; }
	inline bool has_vz_ang_bc() const noexcept { return pvz_ang == &vz_ang_bc; }
	inline double get_vx_ang_bc() const noexcept { return vx_ang_bc; }
	inline double get_vy_ang_bc() const noexcept { return vy_ang_bc; }
	inline double get_vz_ang_bc() const noexcept { return vz_ang_bc; }

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

	inline void set_ax_ang_bc(double _a_ang) { pax_ang = &ax_ang_bc; ax_ang_bc = _a_ang; }
	inline void set_ay_ang_bc(double _a_ang) { pay_ang = &ay_ang_bc; ay_ang_bc = _a_ang; }
	inline void set_az_ang_bc(double _a_ang) { paz_ang = &az_ang_bc; az_ang_bc = _a_ang; }
	inline void set_a_ang_bc(double _ax_ang, double _ay_ang, double _az_ang)
	{ set_ax_bc(_ax_ang); set_ay_bc(_ay_ang); set_az_bc(_az_ang); }

	inline void set_vx_ang_bc(double _v_ang) { pvx_ang = &vx_ang_bc; vx_bc = _v_ang; }
	inline void set_vy_ang_bc(double _v_ang) { pvy_ang = &vy_ang_bc; vy_bc = _v_ang; }
	inline void set_vz_ang_bc(double _v_ang) { pvz_ang = &vz_ang_bc; vz_bc = _v_ang; }
	inline void set_v_ang_bc(double _vx_ang, double _vy_ang, double _vz_ang)
	{ set_vx_bc(_vx_ang); set_vy_bc(_vy_ang); set_vz_bc(_vz_ang); }

	inline void add_fx_ext(double _f) { fx_ext += _f; }
	inline void add_fy_ext(double _f) { fy_ext += _f; }
	inline void add_fz_ext(double _f) { fz_ext += _f; }
	inline void add_mx_ext(double _m) { mx_ext += _m; }
	inline void add_my_ext(double _m) { my_ext += _m; }
	inline void add_mz_ext(double _m) { mz_ext += _m; }
	inline void add_f_ext(double _fx, double _fy, double _fz,
						  double _x,  double _y,  double _z)
	{
		fx_ext += _fx;
		fy_ext += _fy;
		fz_ext += _fz;
		Point3D gp(_x, _y, _z);
		Point3D lp = to_local_coord(gp);
		mx_ext += lp.y * _fz - lp.z * _fy;
		my_ext += lp.z * _fx - lp.x * _fz;
		mz_ext += lp.x * _fy - lp.y * _fx;
	}

	inline size_t get_bface_num() const { return bface_num; }
	inline Face* get_bfaces() { return bfaces; }

	int init_mesh(const char *file_name, double dx, double dy, double dz);
	
	void set_init_state(double density,
		double fx_contact, double fy_contact, double fz_contact,
		double ax, double ay, double az,
		double vx, double vy, double vz,
		double x,  double y,  double z);
	
	inline Face* alloc_bfaces(size_t num)
	{
		clear_bfaces();
		if (num == 0)
			return nullptr;
		bfaces = new Face[num];
		bface_num = num;
		return bfaces;
	}

	inline void init_calculation() noexcept
	{
		x_ori = x;
		y_ori = y;
		z_ori = z;
		ux = 0.0;
		uy = 0.0;
		uz = 0.0;
	}

	inline void reset_substep() noexcept
	{
		fx_con = 0.0;
		fy_con = 0.0;
		fz_con = 0.0;
		mx_con = 0.0;
		my_con = 0.0;
		mz_con = 0.0;
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
		
		// 3D rotation
		// transform I
		Matrix3x3 T_mat;
		T_mat << ix.x, ix.y, ix.z,
				 iy.x, iy.y, iy.z,
				 iz.x, iz.y, iz.z;
		Matrix3x3 cur_moi = T_mat.transpose() * moi_mat * T_mat;
		Vector3 v_ang_vec(vx_ang, vy_ang, vz_ang);
		Vector3 m_vec(mx_ext + mx_con, my_ext + my_con, mz_ext + mz_con);
		Vector3 a_ang = cur_moi.ldlt().solve(m_vec - v_ang_vec.cross(cur_moi * v_ang_vec));
		// cal angular velocity
		ax_ang = a_ang.x();
		ay_ang = a_ang.y();
		az_ang = a_ang.z();
		ax_ang = *pax_ang;
		ay_ang = *pay_ang;
		az_ang = *paz_ang;
		vx_ang += ax_ang * dt;
		vy_ang += ay_ang * dt;
		vz_ang += az_ang * dt;
		vx_ang = *pvx_ang;
		vy_ang = *pvy_ang;
		vz_ang = *pvz_ang;
		// adjust local axises
		double tan_ang = tan(v_ang.norm() * dt);
		Vector3D tmp;
		tmp.cross(v_ang, ix);
		tmp.scale(tan_ang);
		ix.add(tmp).normalize();
		tmp.cross(v_ang, iy);
		tmp.scale(tan_ang);
		iy.add(tmp).normalize();
		iz.cross(ix, iy);
		iy.cross(iz, ix);
	}

	inline void add_con_force(double _fx, double _fy, double _fz,
							  double _x,  double  _y, double  _z) noexcept
	{
		fx_con += _fx;
		fy_con += _fy;
		fz_con += _fz;
		Point3D gp(_x, _y, _z);
		Point3D lp = to_local_coord(gp);
		mx_ext += lp.y * _fz - lp.z * _fy;
		my_ext += lp.z * _fx - lp.x * _fz;
		mz_ext += lp.x * _fy - lp.y * _fx;
	}
	
	template <typename Point3DType>
	inline Point3D to_local_coord(Point3DType&gp) const noexcept
	{
		double dx = gp.x - x;
		double dy = gp.y - y;
		double dz = gp.z - z;
		return Point3D (ix.x * dx + ix.y * dy + ix.z * dz,
						iy.x * dx + iy.y * dy + iy.z * dz,
						iz.x * dx + iz.y * dy + iz.z * dz);
	}
	inline Point3D to_global_coord(Point3D &lp) const noexcept
	{
		double dx = ix.x * lp.x + iy.x * lp.y + iz.x * lp.z;
		double dy = ix.y * lp.x + iy.y * lp.y + iz.y * lp.z;
		double dz = ix.z * lp.x + iy.z * lp.y + iz.z * lp.z;
		return Point3D(dx + x, dy + y, dz + z);
	}

	inline double get_grid_h() const noexcept { return g_h; }
	inline size_t get_grid_x_num() const noexcept { return g_x_num; }
	inline size_t get_grid_y_num() const noexcept { return g_y_num; }
	inline size_t get_grid_z_num() const noexcept { return g_z_num; }

	void clear_bg_grids();
	int init_bg_grids(double _g_h, double expand_size);
	void set_dist_max(double _dist_max);

	inline Cube get_cur_bbox()
	{
		Cube res;
		res.xl = g_bbox.xl + x;
		res.xu = g_bbox.xu + x;
		res.yl = g_bbox.yl + y;
		res.yu = g_bbox.yu + y;
		res.zl = g_bbox.zl + z;
		res.zu = g_bbox.zu + z;
		return res;
	}
	template <typename Point3DType>
	inline bool cal_dist_and_dir_to_pt(Point3DType& pt,
		double& dist, double& nx, double& ny, double& nz)
	{
		Point3D lpt = to_local_coord(pt);
		return cal_dist_and_dir_to_pt_internal(lpt, dist, nx, ny, nz);
	}

protected: // background grid
	double g_h;
	Cube g_bbox;
	size_t g_x_num, g_y_num, g_z_num;
	size_t g_xy_num, g_num;
	Grid* grids;

	inline const Cube& get_grid_bbox() const noexcept { return g_bbox; }
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
	bool cal_dist_and_dir_to_pt_internal(Point3D& pt,
		double& dist, double& nx, double& ny, double& nz);
	void search_closest_face(SearchClosestFaceParam &param);
};

#endif