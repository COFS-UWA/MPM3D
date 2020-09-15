#ifndef __Model_T2D_CHM_d_h__
#define __Model_T2D_CHM_d_h__

#include "TriangleUtils.h"
#include "BCs.h"
#include "TriangleMeshTemplate.hpp"
#include "MatModelContainer.h"
#include "SearchingGrid2D.hpp"
#include "ParticleGenerator2D.hpp"
#include "RigidBody/RigidCircle.h"
#include "macro_utils.h"
#include "Model.h"

namespace Model_T2D_CHM_d_Internal
{
struct Node
{
	size_t id;
	double x, y;

	// solid phase
	bool has_mp_s;
	double m_s;
	double ax_s, ay_s;
	double vx_s, vy_s;
	double dux_s, duy_s;
	double fx_ext_s, fy_ext_s;
	double fx_int_s, fy_int_s;
	// strain enhancement
	double pcl_vol_s, de_vol_s;

	// fluid phase
	bool has_mp_f;
	double m_f;
	double ax_f, ay_f;
	double vx_f, vy_f;
	double dux_f, duy_f;
	double fx_ext_f, fy_ext_f;
	double fx_int_f, fy_int_f;
	// strain enhancement
	double pcl_vol_f, de_vol_f;

	// solid - fluid interaction
	double fx_drag, fy_drag;
};

struct Element;

struct SolidParticle
{
	size_t id;
	double x, y;
	double vx, vy;
	double m, density, n;

	double e11, e22, e12;
	double s11, s22, s12;

	double x_ori, y_ori;
	double ux, uy;
	
	double vol_s, vol;
	inline double get_vol() { return m / (density * (1.0 - n)); }

	Element* pe;
	double N1, N2, N3;

	MatModel::MaterialModel* mm;
	inline void set_mat_model(MatModel::MaterialModel& _mm)
	{ _mm.ext_data_pt = this; mm = &_mm; }
};

struct FluidParticle
{
	size_t id;
	double x, y;
	double vx, vy;
	double m, density;
	
	double p;

	double x_ori, y_ori;
	double ux, uy;

	double vol;
	inline double get_vol() { return m / density; }

	Element* pe;
	double N1, N2, N3;
};

struct Element
{
	size_t id;
	size_t n1, n2, n3;

	double area;
	PointInTriangle pt_in_tri;

	// solid
	bool has_mp_s;
	double vol_s, n, s11, s22, s12;
	double dde11, dde22, de12, de_vol_s;

	// seepage fluid
	bool has_mp_fs;
	double vol_fs, p_s;
	// pure fluid
	bool has_mp_fp;
	double vol_fp, p_p;
	double de_vol_f;
};

struct Edge { size_t n1, n2; };

typedef TriangleMeshTemplate<Node, Element, Edge> ParentMesh;

}

class Step_T2D_CHM_d;
int solve_substep_T2D_CHM_d(void* _self);

class Model_T2D_CHM_d;
class ResultFile_hdf5;
namespace Model_T2D_CHM_d_hdf5_utilities
{
	int output_background_mesh_to_hdf5_file(Model_T2D_CHM_d& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_background_mesh_from_hdf5_file(Model_T2D_CHM_d& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_boundary_condition_to_hdf5_file(Model_T2D_CHM_d& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_boundary_condition_from_hdf5_file(Model_T2D_CHM_d& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_pcl_data_to_hdf5_file(Model_T2D_CHM_d& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_pcl_data_from_hdf5_file(Model_T2D_CHM_d& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_material_model_to_hdf5_file(Model_T2D_CHM_d& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_material_model_from_hdf5_file(Model_T2D_CHM_d& md, ResultFile_hdf5& rf, hid_t grp_id);
	int output_rigid_circle_to_hdf5_file(Model_T2D_CHM_d& md, ResultFile_hdf5& rf, hid_t grp_id);
	int load_rigid_circle_from_hdf5_file(Model_T2D_CHM_d& md, ResultFile_hdf5& rf, hid_t grp_id);
}

struct Model_T2D_CHM_d : public Model,
	public Model_T2D_CHM_d_Internal::ParentMesh,
	public MatModel::MatModelContainer
{
	friend class Step_T2D_CHM_d;
	friend int solve_substep_T2D_CHM_d(void *_self);
	friend int Model_T2D_CHM_d_hdf5_utilities::output_background_mesh_to_hdf5_file(Model_T2D_CHM_d& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_CHM_d_hdf5_utilities::load_background_mesh_from_hdf5_file(Model_T2D_CHM_d& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_CHM_d_hdf5_utilities::output_boundary_condition_to_hdf5_file(Model_T2D_CHM_d& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_CHM_d_hdf5_utilities::load_boundary_condition_from_hdf5_file(Model_T2D_CHM_d& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_CHM_d_hdf5_utilities::output_pcl_data_to_hdf5_file(Model_T2D_CHM_d& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_CHM_d_hdf5_utilities::load_pcl_data_from_hdf5_file(Model_T2D_CHM_d& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_CHM_d_hdf5_utilities::output_material_model_to_hdf5_file(Model_T2D_CHM_d& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_CHM_d_hdf5_utilities::load_material_model_from_hdf5_file(Model_T2D_CHM_d& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_CHM_d_hdf5_utilities::output_rigid_circle_to_hdf5_file(Model_T2D_CHM_d& md, ResultFile_hdf5& rf, hid_t grp_id);
	friend int Model_T2D_CHM_d_hdf5_utilities::load_rigid_circle_from_hdf5_file(Model_T2D_CHM_d& md, ResultFile_hdf5& rf, hid_t grp_id);

public:
	typedef Model_T2D_CHM_d_Internal::ParentMesh ParentMesh;
	typedef Model_T2D_CHM_d_Internal::Node Node;
	typedef Model_T2D_CHM_d_Internal::Element Element;
	typedef Model_T2D_CHM_d_Internal::Edge Edge;
	typedef Model_T2D_CHM_d_Internal::SolidParticle SolidParticle;
	typedef Model_T2D_CHM_d_Internal::FluidParticle FluidParticle;

protected:
	double Kf;  // Bulk modulus of water
	double k;   // Permeability
	double miu; // Dynamic viscosity

	size_t spcl_num;
	SolidParticle *spcls;
	size_t fpcl_num;
	FluidParticle *fpcls;

	// boundary conditions
	size_t bfx_num, bfy_num;
	BodyForceAtPcl *bfxs, *bfys;
	size_t tx_num, ty_num;
	TractionBCAtPcl *txs, *tys;
	size_t asx_num, asy_num;
	AccelerationBC *asxs, *asys;
	size_t afx_num, afy_num;
	AccelerationBC *afxs, *afys;
	size_t vsx_num, vsy_num;
	VelocityBC* vsxs, * vsys;
	size_t vfx_num, vfy_num;
	VelocityBC *vfxs, *vfys;
	
	// background grid to accelerate searching
	SearchingGrid2D<Model_T2D_CHM_d> search_bg_grid;

	void init_mesh_shape_funcs();
	
public:
	Model_T2D_CHM_d();
	~Model_T2D_CHM_d();
	
	inline double get_bg_grid_xl() { return search_bg_grid.get_x_min(); }
	inline double get_bg_grid_xu() { return search_bg_grid.get_x_max(); }
	inline double get_bg_grid_yl() { return search_bg_grid.get_y_min(); }
	inline double get_bg_grid_yu() { return search_bg_grid.get_y_max(); }
	inline double get_bg_grid_hx() { return search_bg_grid.get_hx(); }
	inline double get_bg_grid_hy() { return search_bg_grid.get_hy(); }

	inline size_t get_solid_pcl_num() const noexcept { return spcl_num; }
	inline SolidParticle* get_solid_pcls() { return spcls; }
	inline size_t get_fluid_pcl_num() const noexcept { return fpcl_num; }
	inline FluidParticle* get_fluid_pcls() { return fpcls; }

	int init_mesh(double *node_coords, size_t n_num,
				  size_t *elem_n_ids,  size_t e_num);
	int load_mesh_from_hdf5(const char* file_name);

	int init_search_grid(double _hx, double _hy);

	int init_solid_pcls(size_t num, double m, double density,
						double n, double _k);
	int init_fluid_pcls(size_t num, double m, double density,
						double _Kf, double _miu);

	int init_solid_pcls(ParticleGenerator2D<Model_T2D_CHM_d>& pg,
						double density, double n, double _k);
	int init_fluid_pcls(ParticleGenerator2D<Model_T2D_CHM_d>& pg,
						double density, double _Kf, double _miu);
	
	void alloc_solid_pcls(size_t num);
	void clear_solid_pcls();
	void alloc_fluid_pcls(size_t num);
	void clear_fluid_pcls();

	INIT_BC_TEMPLATE(bfx, BodyForceAtPcl)
	INIT_BC_TEMPLATE(bfy, BodyForceAtPcl)
	INIT_BC_TEMPLATE(tx, TractionBCAtPcl)
	INIT_BC_TEMPLATE(ty, TractionBCAtPcl)
	INIT_BC_TEMPLATE(asx, AccelerationBC)
	INIT_BC_TEMPLATE(asy, AccelerationBC)
	INIT_BC_TEMPLATE(afx, AccelerationBC)
	INIT_BC_TEMPLATE(afy, AccelerationBC)
	INIT_BC_TEMPLATE(vsx, VelocityBC)
	INIT_BC_TEMPLATE(vsy, VelocityBC)
	INIT_BC_TEMPLATE(vfx, VelocityBC)
	INIT_BC_TEMPLATE(vfy, VelocityBC)

	inline bool is_in_triangle(Element& e, double x, double y)
	{ return e.pt_in_tri.is_in_triangle(x, y); }

	template <typename Point2D>
	inline bool is_in_triangle(Element &e, Point2D &pt)
	{ return e.pt_in_tri.is_in_triangle<Point2D>(pt); }

	// search using background grid
	template <typename Point2D>
	inline Element* find_in_which_element(Point2D& pt)
	{
		return search_bg_grid.find_in_which_element<Point2D>(pt);
	}

// ===================== rigid circle =====================
protected:
	bool rigid_circle_is_init;
	RigidCircle rigid_circle;
	// contact stiffness
	double Ks_cont, Kfs_cont, Kfp_cont;

public: // interaction with rigid circle
	inline bool rigid_circle_is_valid() { return rigid_circle_is_init; }
	inline RigidCircle& get_rigid_circle() { return rigid_circle; }

	inline void init_rigid_circle(
		double _Ks_cont, double _Kfs_cont, double _Kfp_cont,
		double x, double y, double r, double density = 1.0)
	{
		rigid_circle_is_init = true;
		Ks_cont = _Ks_cont;
		Kfs_cont = _Kfs_cont;
		Kfp_cont = _Kfp_cont;
		rigid_circle.init(x, y, r, density);
	}
	inline void set_rigid_circle_velocity(double vx, double vy, double w)
	{ rigid_circle.set_v_bc(vx, vy, w); }

};

#endif