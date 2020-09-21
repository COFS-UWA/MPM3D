#ifndef __Rigid_Tetrahedron_Mesh_h__
#define __Rigid_Tetrahedron_Mesh_h__

#include <Eigen/Dense>

#include "ItemBuffer.hpp"
#include "ItemStack.hpp"
#include "GeometryUtils.h"
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
		PointToTriangleDistance pt_tri_dist;
	};

	typedef TetrahedronMeshTemplate<Node, Element, Edge> ParentTehClass;
}

class RigidTetrahedronMesh :
	protected RigidTetrahedronMesh_Internal::ParentTehClass
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

protected:
	typedef RigidTetrahedronMesh_Internal::ParentTehClass ParentTehClass;

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

	union
	{
		struct { double x, y, z; };
		Point3D cen_pos;
	};
	double m;
	
	double ax, ay, az;
	double vx, vy, vz;

	double ax_bc, ay_bc, az_bc;
	double vx_bc, vy_bc, vz_bc;

	double *pax, *pay, *paz;
	double *pvx, *pvy, *pvz;

	// external force
	double fx_ext, fy_ext, fz_ext;
	// contact force
	double fx_con, fy_con, fz_con;
	
	// for precision
	double x_ori, ux;
	double y_ori, uy;
	double z_ori, uz;

	// rotation
	double ax_ang, ay_ang, az_ang;
	double vx_ang, vy_ang, vz_ang;
	union
	{
		struct { double x_ang, y_ang, z_ang; };
		Vector3D cen_ang;
	};

	double ax_ang_bc, ay_ang_bc, az_ang_bc;
	double vx_ang_bc, vy_ang_bc, vz_ang_bc;

	double* pax_ang, * pay_ang, * paz_ang;
	double* pvx_ang, * pvy_ang, * pvz_ang;

	double mx_ext, my_ext, mz_ext;
	double mx_con, my_con, mz_con;
	Vector3D ix, iy, iz;

	typedef Eigen::Matrix<double, 3, 3> Matrix3x3;
	typedef Eigen::Matrix<double, 3, 1> Vector3;

	Matrix3x3 moi_mat; // moment of inertia

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

	void cal_m_and_moi();

public:
	RigidTetrahedronMesh();
	~RigidTetrahedronMesh();

	inline double get_density() { return density; }
	inline double get_m() const noexcept { return m; }
	inline const Matrix3x3 &get_moi() const noexcept { return moi_mat; }

	inline double get_ax() const noexcept { return ax; }
	inline double get_ay() const noexcept { return ay; }
	inline double get_az() const noexcept { return az; }
	inline double get_vx() const noexcept { return vx; }
	inline double get_vy() const noexcept { return vy; }
	inline double get_vz() const noexcept { return vz; }
	inline double get_x() const noexcept { return x; }
	inline double get_y() const noexcept { return y; }
	inline double get_z() const noexcept { return z; }
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
	inline double get_x_ang() const noexcept { return x_ang; }
	inline double get_y_ang() const noexcept { return y_ang; }
	inline double get_z_ang() const noexcept { return z_ang; }
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
	
	template <typename Point3DType>
	inline void get_velocity(Point3DType &pos, Vector3D &v) const noexcept
	{
		double dx = pos.x - x;
		double dy = pos.y - y;
		double dz = pos.z - z;
		v.x = vx + vy_ang * dz - vz_ang * dy;
		v.y = vy + vz_ang * dx - vx_ang * dz;
		v.z = vz + vx_ang * dy - vy_ang * dx;
	}

	inline void set_density(double den)
	{
		density = den;
		cal_m_and_moi();
	}

	inline void set_ax_bc(double _a) { pax = &ax_bc; ax_bc = _a; }
	inline void set_ay_bc(double _a) { pay = &ay_bc; ay_bc = _a; }
	inline void set_az_bc(double _a) { paz = &az_bc; az_bc = _a; }

	inline void set_vx_bc(double _v) { pvx = &vx_bc; vx_bc = _v; }
	inline void set_vy_bc(double _v) { pvy = &vy_bc; vy_bc = _v; }
	inline void set_vz_bc(double _v) { pvz = &vz_bc; vz_bc = _v; }

	inline void set_ax_ang_bc(double _a_ang) { pax_ang = &ax_ang_bc; ax_ang_bc = _a_ang; }
	inline void set_ay_ang_bc(double _a_ang) { pay_ang = &ay_ang_bc; ay_ang_bc = _a_ang; }
	inline void set_az_ang_bc(double _a_ang) { paz_ang = &az_ang_bc; az_ang_bc = _a_ang; }

	inline void set_vx_ang_bc(double _v_ang) { pvx_ang = &vx_ang_bc; vx_ang_bc = _v_ang; }
	inline void set_vy_ang_bc(double _v_ang) { pvy_ang = &vy_ang_bc; vy_ang_bc = _v_ang; }
	inline void set_vz_ang_bc(double _v_ang) { pvz_ang = &vz_ang_bc; vz_ang_bc = _v_ang; }

	inline void set_v_bc(
		double _vx,
		double _vy,
		double _vz,
		double _vx_ang = 0.0,
		double _vy_ang = 0.0,
		double _vz_ang = 0.0
		)
	{
		set_vx_bc(_vx); set_vy_bc(_vy); set_vz_bc(_vz);
		set_vx_ang_bc(_vx_ang); set_vy_ang_bc(_vy_ang); set_vz_ang_bc(_vz_ang);
	}

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
		Point3D lp;
		to_local_coord(gp, lp);
		mx_ext += lp.y * _fz - lp.z * _fy;
		my_ext += lp.z * _fx - lp.x * _fz;
		mz_ext += lp.x * _fy - lp.y * _fx;
	}

	inline const Point3D& get_centre() const noexcept { return cen_pos; }
	inline const Vector3D& get_ix() const noexcept { return ix; }
	inline const Vector3D& get_iy() const noexcept { return iy; }
	inline const Vector3D& get_iz() const noexcept { return iz; }

	inline size_t get_node_num() const noexcept { return ParentTehClass::get_node_num(); }
	inline const Node* get_nodes() const noexcept { return ParentTehClass::get_nodes(); }
	inline Node* get_nodes() noexcept { return ParentTehClass::get_nodes(); }
	inline size_t get_elem_num() const noexcept { return ParentTehClass::get_elem_num(); }
	inline const Element* get_elems() const noexcept { return ParentTehClass::get_elems(); }
	inline Element* get_elems() noexcept { return ParentTehClass::get_elems(); }
	inline size_t get_bface_num() const noexcept { return bface_num; }
	inline const Face* get_bfaces() const noexcept { return bfaces; }
	inline Face* get_bfaces() noexcept { return bfaces; }

	int init(double _density, const char* file_name,
			 double dx, double dy, double dz,
			 double dx_ang, double dy_ang, double dz_ang);

	int init_mesh(const char *file_name, double dx, double dy, double dz,
				  double dx_ang, double dy_ang, double dz_ang);
	
	void init_mesh_properties_after_loading();

public: // for reading object from file
	void set_init_state(double _density,
		double _fx_contact, double _fy_contact, double _fz_contact,
		double _ax, double _ay, double _az,
		double _vx, double _vy, double _vz,
		double _x,  double _y,  double _z,
		double _mx_contact, double _my_contact, double _mz_contact,
		double _ax_ang, double _ay_ang, double _az_ang,
		double _vx_ang, double _vy_ang, double _vz_ang,
		double _x_ang, double _y_ang, double _z_ang);
	
	inline Node* alloc_nodes(size_t num) { return ParentTehClass::alloc_nodes(num); }
	inline Element* alloc_elements(size_t num) { return ParentTehClass::alloc_elements(num); }
	inline Face* alloc_bfaces(size_t num)
	{
		clear_bfaces();
		if (num == 0)
			return nullptr;
		bfaces = new Face[num];
		bface_num = num;
		return bfaces;
	}

public: // calculation functions
	inline void init_calculation() noexcept
	{
		x_ori = x;
		y_ori = y;
		z_ori = z;
		ux = 0.0;
		uy = 0.0;
		uz = 0.0;
	}

	inline void get_cur_bbox(Cube &bbox) const noexcept
	{
		double hx = (g_bbox.xu - g_bbox.xl) * 0.5;
		double hy = (g_bbox.yu - g_bbox.yl) * 0.5;
		double hz = (g_bbox.zu - g_bbox.zl) * 0.5;
		double rx = abs(ix.x * hx) + abs(iy.x * hy) + abs(iz.x * hz);
		double ry = abs(ix.y * hx) + abs(iy.y * hy) + abs(iz.y * hz);
		double rz = abs(ix.z * hx) + abs(iy.z * hy) + abs(iz.z * hz);
		bbox.xl = x - rx;
		bbox.xu = x + rx;
		bbox.yl = y - ry;
		bbox.yu = y + ry;
		bbox.zl = z - rz;
		bbox.zu = z + rz;
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

	inline void add_contact_force(double _fx, double _fy, double _fz,
								  double _x,  double  _y, double  _z) noexcept
	{
		fx_con += _fx;
		fy_con += _fy;
		fz_con += _fz;
		double dx = _x - x;
		double dy = _y - y;
		double dz = _z - z;
		mx_ext += dy * _fz - dz * _fy;
		my_ext += dz * _fx - dx * _fz;
		mz_ext += dx * _fy - dy * _fx;
	}
	
	inline void update_motion(double dt) noexcept
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
		Vector3 a_ang = cur_moi.partialPivLu().solve(m_vec - v_ang_vec.cross(cur_moi * v_ang_vec));
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
		x_ang += vx_ang * dt;
		trim_to_pi(x_ang);
		y_ang += vy_ang * dt;
		trim_to_pi(y_ang);
		z_ang += vz_ang * dt;
		trim_to_pi(z_ang);
		// update ix, iy, iz
		ix.x = 1.0, ix.y = 0.0, ix.z = 0.0;
		iy.x = 0.0, iy.y = 1.0, iy.z = 0.0;
		iz.x = 0.0, iz.y = 0.0, iz.z = 1.0;
		rotate_axses_by_angle(cen_ang, ix, iy, iz);
	}

	template <typename Point3DType>
	inline void to_local_coord(const Point3DType&gp, Point3D &lp) const noexcept
	{ point_from_global_to_local_coordinate<Point3DType, Point3D>(cen_pos, ix, iy, iz, gp, lp); }
	template <typename Point3DType>
	inline void to_global_coord(const Point3DType &lp, Point3D &gp) const noexcept
	{ point_from_local_to_global_coordinate<Point3DType, Point3D>(cen_pos, ix, iy, iz, lp, gp);	}

public: // bg search grid
	inline double get_grid_h() const noexcept { return g_h; }
	inline size_t get_grid_x_num() const noexcept { return g_x_num; }
	inline size_t get_grid_y_num() const noexcept { return g_y_num; }
	inline size_t get_grid_z_num() const noexcept { return g_z_num; }

	void clear_bg_grids();
	int init_bg_grids(double _g_h, double expand_size);
	void set_dist_max(double _dist_max);

	template <typename Point3DType>
	inline bool cal_dist_and_dir_to_pt(const Point3DType& gpt,
		double& dist, double& nx, double& ny, double& nz)
	{
		Point3D lpt;
		to_local_coord<Point3DType>(gpt, lpt);
		Vector3D lnorm, gnorm;
		if (cal_dist_and_dir_to_pt_internal(lpt, dist, lnorm.x, lnorm.y, lnorm.z))
		{
			vector_from_local_to_global_coordinate<Vector3D, Vector3D>(ix, iy, iz, lnorm, gnorm);
			nx = gnorm.x;
			ny = gnorm.y;
			nz = gnorm.z;
			return true;
		}
		return false;
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
	
	PointInTetrahedron pt_in_teh;
	TetrahedronAABBCollisionSAT teh_aabb_collision;
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

	TriangleAABBCollisionSAT tri_aabb_collision;
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