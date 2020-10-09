#ifndef __Model_T2D_ME_mt_h__
#define __Model_T2D_ME_mt_h__

#include <stdint.h>

#include "BCs.h"
#include "macro_utils.h"
#include "Model.h"
#include "MatModelContainer.h"
#include "ParticleGenerator2D.hpp"
#include "TriangleMesh.h"

class Step_T2D_ME_mt;
int solve_substep_T2D_ME_mt(void* _self);

struct Model_T2D_ME_mt : public Model,
	public MatModel::MatModelContainer
{
	friend class Step_T2D_ME_mt;
	friend int solve_substep_T2D_ME_mt(void* _self);
	
public:
	struct PclBodyForce { float bfx, bfy; };
	struct PclTraction { float tx, ty; };
	struct PclPos { float x, y; };

	struct PclDisp { float ux, uy; };
	struct PclV { float vx, vy; };
	struct PclShapeFunc { float N1, N2, N3; };
	struct PclStress { float s11, s22, s12; };

	struct ElemNodeIndex { uint32_t n1, n2, n3; };
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

	struct ElemStrainInc { float de11, de22, de12; };
	struct ElemStress { float s11, s22, s12; };

	struct ElemNodeVM { float vm, vmx, vmy; };
	struct ElemNodeForce { float fx, fy; };

	struct PclSortedVarArray
	{
		uint32_t* pcl_index; // ori_pcl_num
		float* pcl_density; // ori_pcl_num
		PclDisp* pcl_disp; // ori_pcl_num
		PclV* pcl_v; // ori_pcl_num
		PclShapeFunc* pcl_N; // ori_pcl_num
		PclStress* pcl_stress; // ori_pcl_num
		uint32_t* elem_has_pcl_num; // elem_num
	};
	
protected:
	uint32_t pcl_num;
	uint32_t elem_num;
	uint32_t node_num;
	uint32_t vx_bc_num;
	uint32_t vy_bc_num;
	
	float *pcl_m; // ori_pcl_num
	PclBodyForce* pcl_bf; // ori_pcl_num
	PclTraction* pcl_t; // ori_pcl_num
	PclPos* pcl_pos; // ori_pcl_num
	MatModel::MaterialModel **pcl_mat_model; // ori_pcl_num

	PclSortedVarArray pcl_sorted_var_array[2];

	// element data
	ElemNodeIndex *elem_node_id; // elem_num
	float *elem_area; // elem_num
	ElemShapeFuncAB* elem_sf_ab; // elem_num
	ElemShapeFuncC* elem_sf_c; // elem_num

	// element calculation data
	float *elem_density; // elem_num
	ElemStrainInc *elem_de; // elem_num
	ElemStress *elem_stress; // elem_num
	float *elem_am; // elem_num
	float *elem_am_de_vol; // elem_num

	// element-node data
	ElemNodeVM* elem_node_vm; // elem_num * 3
	ElemNodeForce* elem_node_force; // elem_num * 3
	
	uint32_t *elem_id_array; // elem_num * 3 
	uint32_t *node_elem_id_array; // elem_num * 3
	uint32_t *node_elem_list;  // node_num

	float *node_ax; // node_num
	float *node_ay; // node_num
	float *node_vx; // node_num
	float *node_vy; // node_num
	float *node_am; // node_num
	float *node_de_vol; // node_num

	// nodal velocity bcs
	uint32_t* vx_bcs; // vx_bc_num
	uint32_t* vy_bcs; // vy_bc_num

// ======== Non calculation data ========
	uint32_t ori_pcl_num;
	char *pcl_mem_raw;

	char *mesh_mem_raw;

	// boundary conditions
	size_t bfx_num, bfy_num;
	BodyForceAtPcl* bfxs, * bfys;
	size_t tx_num, ty_num;
	TractionBCAtPcl* txs, * tys;

public:
	Model_T2D_ME_mt();
	~Model_T2D_ME_mt();

	inline size_t get_pcl_num() const noexcept { return size_t(pcl_num); }
	inline size_t get_elem_num() const noexcept { return size_t(elem_num); }
	inline size_t get_node_num() const noexcept { return size_t(node_num); }
	inline MatModel::MaterialModel **get_mat_models() noexcept { return pcl_mat_model; }

	void clear_mesh();
	void init_mesh(const TriangleMesh& mesh);
	
	void clear_search_grid();
	int init_search_grid(const TriangleMesh& mesh, double _hx, double _hy);
	
	void alloc_pcls(size_t num);
	void clear_pcls();
	int init_pcls(size_t num, double m, double density);
	int Model_T2D_ME_mt::init_pcls(
		ParticleGenerator2D<TriangleMesh>& pg,
		double density);

	INIT_BC_TEMPLATE(bfx, BodyForceAtPcl)
	INIT_BC_TEMPLATE(bfy, BodyForceAtPcl)
	INIT_BC_TEMPLATE(tx, TractionBCAtPcl)
	INIT_BC_TEMPLATE(ty, TractionBCAtPcl)
	void init_bfxs(size_t bf_num, const size_t *bf_pcls, const double *bfs);
	void init_bfys(size_t bf_num, const size_t *bf_pcls, const double *bfs);
	void init_txs(size_t t_num, const size_t *t_pcls, const double *ts);
	void init_tys(size_t t_num, const size_t *t_pcls, const double *ts);

	void alloc_vx_bcs(size_t vx_bc_num);
	void clear_vx_bcs();
	void alloc_vy_bcs(size_t vy_bc_num);
	void clear_vy_bcs();
	void init_fixed_vx_bc(size_t vx_bc_num, const size_t *vx_bcs);
	void init_fixed_vy_bf(size_t vy_bc_num, const size_t* vy_bcs);

protected:
	inline bool is_in_element(
		float pcl_x, float pcl_y,
		uint32_t elem_id,
		PclShapeFunc &pcl_sf
		) noexcept
	{
		ElemShapeFuncAB& e_sfab = elem_sf_ab[elem_id];
		ElemShapeFuncC& e_sfc = elem_sf_c[elem_id];
		pcl_sf.N1 = e_sfab.a1 * pcl_x + e_sfab.b1 * pcl_y + e_sfc.c1;
		pcl_sf.N2 = e_sfab.a2 * pcl_x + e_sfab.b2 * pcl_y + e_sfc.c2;
		pcl_sf.N3 = 1.0f - pcl_sf.N1 - pcl_sf.N2;
		return pcl_sf.N1 >= 0.0f && pcl_sf.N1 <= 1.0f
			&& pcl_sf.N2 >= 0.0f && pcl_sf.N2 <= 1.0f
			&& pcl_sf.N3 >= 0.0f && pcl_sf.N3 <= 1.0f;
	}
	
	// background grid for mesh
	struct MeshBgGrid
	{
		uint32_t elem_start_id, elem_end_id;
	};

	float grid_xl, grid_yl;
	float grid_hx, grid_hy;
	uint32_t grid_x_num, grid_y_num;
	MeshBgGrid *grids;
	uint32_t* grid_elem_id_array;

	inline uint32_t find_pcl_in_which_elem(
		float pcl_x, float pcl_y,
		PclShapeFunc& pcl_sf
		) noexcept
	{
		uint32_t x_id, y_id, elem_id;
		x_id = (pcl_x - grid_xl) / grid_hx;
		y_id = (pcl_y - grid_yl) / grid_hy;
		MeshBgGrid &g = grids[grid_x_num * y_id + x_id];
		for (uint32_t e_id = g.elem_start_id; e_id < g.elem_end_id; ++e_id)
		{
			elem_id = grid_elem_id_array[e_id];
			if (is_in_element(pcl_x, pcl_y, elem_id, pcl_sf))
				return elem_id;
		}
		return UINT32_MAX;
	}
};

#endif