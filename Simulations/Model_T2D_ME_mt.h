#ifndef __Model_T2D_ME_mt_h__
#define __Model_T2D_ME_mt_h__

#include <stdint.h>

#include "BCs.h"
#include "macro_utils.h"
#include "Model.h"
#include "MatModelContainer.h"

class Step_T2D_ME_mt;
int solve_substep_T2D_ME_mt(void* _self);

struct Model_T2D_ME_mt : public Model
{
	friend class Step_T2D_ME_mt;
	friend int solve_substep_T2D_ME_mt(void* _self);
	
protected:
	struct PclMass { float m; };
	struct PclBodyForce { float bfx, bfy; };
	struct PclTraction { float tx, ty; };
	struct PclPos { float x, y; };

	struct PclIndex { uint32_t id; };
	struct PclDensity { float density; };
	struct PclDisp { float ux, uy; };
	struct PclV { float vx, vy; };
	struct PclShapeFunc { float N1, N2, N3; };
	struct PclStress { float s11, s22, s12; };

	struct ElemPclList { uint32_t start_id, end_id; };
	struct ElemNodeVarOffset { uint32_t n1, n2, n3; };
	struct ElemDensity { float density; };
	struct ElemStrainInc { float de11, de22, de33; };
	struct ElemStress { float s11, s22, s12; };
	struct ElemShapeFuncAB
	{
		union { float a1; float dN1_dx; };
		union { float b1; float dN1_dy; };
		union { float a2; float dN2_dx; };
		union { float b2; float dN2_dy; };
		union { float a3; float dN3_dx; };
		union { float b3; float dN3_dy; };
	};
	struct ElemShapeFuncC { float c1, c2, c3; };

	struct NodeElemVM { float vm, vmx, vmy; };
	struct NodeElemAF { float am, fx, fy; };
	struct NodeElemDeVol { float am_de_vol; };

	struct NodeElemList { uint32_t start_id, end_id; };
	struct NodeMotion { float ax, ay, vx, vy; };

	uint32_t pcl_num;
	uint32_t node_num;
	uint32_t elem_num;

	PclMass *pcl_m;
	PclBodyForce* pcl_bf;
	PclTraction* pcl_t;
	PclPos* pcl_pos;

	PclIndex* pcl_index;
	PclDensity *pcl_density;
	PclDisp* pcl_disp;
	PclV* pcl_v;
	PclShapeFunc* pcl_N;
	PclStress* pcl_stress;

	ElemPclList *elem_pcl_list;
	ElemNodeVarOffset *elem_node_var_offset;
	ElemDensity* elem_density;
	ElemStrainInc *elem_de;
	ElemStress *elem_stress;
	ElemShapeFuncAB *elem_sf_ab;
	ElemShapeFuncC *elem_sf_c;

	NodeElemAF* ne_af;
	NodeElemVM* ne_vm;
	NodeElemDeVol* ne_de_vol;

	uint32_t* node_elem_offset;
	NodeElemList* node_elem_list;
	NodeMotion *node_motion;

	// boundary conditions
	size_t bfx_num, bfy_num;
	BodyForceAtPcl* bfxs, * bfys;
	size_t tx_num, ty_num;
	TractionBCAtPcl* txs, * tys;
	size_t ax_num, ay_num;
	AccelerationBC* axs, * ays;
	size_t vx_num, vy_num;
	VelocityBC* vxs, * vys;

public:
	Model_T2D_ME_mt();
	~Model_T2D_ME_mt();

	inline size_t get_pcl_num() { return pcl_num; }

	int init_mesh(double* node_coords, size_t n_num,
		size_t* elem_n_ids, size_t e_num);
	int load_mesh_from_hdf5(const char* file_name);

	int init_search_grid(double _hx, double _hy);

	int init_pcls(size_t num, double m, double density);

	void alloc_pcls(size_t num);
	void clear_pcls();

	INIT_BC_TEMPLATE(bfx, BodyForceAtPcl)
	INIT_BC_TEMPLATE(bfy, BodyForceAtPcl)
	INIT_BC_TEMPLATE(tx, TractionBCAtPcl)
	INIT_BC_TEMPLATE(ty, TractionBCAtPcl)
	INIT_BC_TEMPLATE(ax, AccelerationBC)
	INIT_BC_TEMPLATE(ay, AccelerationBC)
	INIT_BC_TEMPLATE(vx, VelocityBC)
	INIT_BC_TEMPLATE(vy, VelocityBC)

	//inline double cal_N1(Element& e, double x, double y)
	//{
	//	return e.a1 * x + e.b1 * y + e.coef1;
	//}
	//inline double cal_N2(Element& e, double x, double y)
	//{
	//	return e.a2 * x + e.b2 * y + e.coef2;
	//}
	//inline double cal_N3(Element& e, double x, double y)
	//{
	//	return e.a3 * x + e.b3 * y + e.coef3;
	//}

//protected: // helper functions
//	void init_mesh_shape_funcs();
};

#endif