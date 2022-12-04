#ifndef __Step_T3D_CHM_up_TBB_Task_h__
#define __Step_T3D_CHM_up_TBB_Task_h__

#define TBB_SUPPRESS_DEPRECATED_MESSAGES 1
#include <tbb/tbb.h>

#include "RigidObject/RigidObjectByT3DMesh.h"

class Step_T3D_CHM_up_TBB;

namespace Step_T3D_CHM_up_TBB_Task
{
	constexpr size_t init_pcl_task_num_per_thread = 100;
	constexpr size_t map_pcl_to_mesh_task_num_per_thread = 100;
	constexpr size_t cal_find_soil_surface_task_num_per_thread = 20;
	constexpr size_t update_node_av_task_num_per_thread = 20;
	constexpr size_t cal_elem_de_task_num_per_thread = 20;
	constexpr size_t cal_node_de_task_num_per_thread = 20;
	constexpr size_t map_mesh_to_pcl_task_num_per_thread = 100;

	constexpr size_t min_pcl_num_per_task = 100;
	constexpr size_t min_node_elem_num_per_task = 20;
	constexpr size_t min_elem_num_per_task = 20;

	constexpr size_t cache_line_size = 64;
	
	struct InitPclRes
	{
		size_t pcl_num;
		double max_pcl_vol;
		inline void join(const InitPclRes& other)
		{
			pcl_num += other.pcl_num;
			if (max_pcl_vol < other.max_pcl_vol)
				max_pcl_vol = other.max_pcl_vol;
		}
	};
	
	class InitPcl
	{
	protected:
		typedef Model_T3D_CHM_up_mt::ShapeFunc ShapeFunc;
		typedef Model_T3D_CHM_up_mt::Displacement Displacement;
		typedef Model_T3D_CHM_up_mt::Position Position;
		
		Step_T3D_CHM_up_TBB &stp;
		const double* pcl_vol_s;
		const double* pcl_n0;
		size_t* in_pcl_in_elems, * in_prev_pcl_ids;
		size_t task_num;
	
	public:
		InitPcl(Step_T3D_CHM_up_TBB &_stp) : stp(_stp) {}
		void init(size_t thread_num) noexcept;
		inline size_t get_task_num() const { return task_num; }
		void work(size_t tsk_id, InitPclRes &res) const;
		inline tbb::task* operator() (tbb::task& parent, size_t tsk_id, InitPclRes &res) const
		{ work(tsk_id, res); return nullptr; }
	};
	
	struct MapPclToBgMeshRes
	{
		Force3D react_force;
		inline void join(const MapPclToBgMeshRes &res) noexcept
		{ react_force.combine(res.react_force); }
	};

	static constexpr double max_Kf_ratio_divider = 1.0e10;
	
	class MapPclToBgMesh
	{
	protected:
		typedef Model_T3D_CHM_up_mt::Force Force;
		typedef Model_T3D_CHM_up_mt::Acceleration Acceleration;
		typedef Model_T3D_CHM_up_mt::Velocity Velocity;
		typedef Model_T3D_CHM_up_mt::Displacement Displacement;
		typedef Model_T3D_CHM_up_mt::Position Position;
		typedef Model_T3D_CHM_up_mt::Stress Stress;
		typedef Model_T3D_CHM_up_mt::ShapeFunc ShapeFunc;
		typedef Model_T3D_CHM_up_mt::DShapeFuncABC DShapeFuncABC;
		typedef Model_T3D_CHM_up_mt::DShapeFuncD DShapeFuncD;
		typedef Model_T3D_CHM_up_mt::ElemNodeVM ElemNodeVM;

		Step_T3D_CHM_up_TBB& stp;

		double m_cav, f_cav_end;
		double u_cav_off, u_div_u_cav_lim;
		double Kf0, k, dyn_viscosity;
		
		// pcl range
		const size_t* pcl_in_elems;
		const size_t* prev_pcl_ids;
		
		// pcl data
		const double* pcl_m_s;
		const double *pcl_vol_s;
		double* pcl_vol;
		const Force* pcl_bf_s, * pcl_bf_f, * pcl_t;

		// bg mesh data
		const DShapeFuncABC *elem_dN_abc;
		const double *elem_vol;
		size_t* elem_has_pcls;
		double* elem_pcl_m, *elem_pcl_pm;
		double* elem_pcl_n, * elem_density_f;
		double *elem_p, * elem_pcl_vol;
		ElemNodeVM* elem_node_vm_s;
		double *elem_node_p;
		Force* elem_node_force;
		double* elem_node_p_force;
		uint16_t *elem_node_at_surface;
		double* elem_u_cav;

		// pcl_vars0
		size_t* pcl_index0;
		double* pcl_n0;
		double* pcl_density_f0;
		Velocity* pcl_v_s0;
		Displacement* pcl_u0;
		Stress* pcl_stress0;
		//double* pcl_p0;
		ShapeFunc* pcl_N0;
		// pcl_vars1
		const size_t* pcl_index1;
		double* pcl_n1;
		const double* pcl_density_f1;
		const Velocity* pcl_v_s1;
		const Displacement* pcl_u1;
		const Stress* pcl_stress1;
		const double* pcl_p1;
		const ShapeFunc* pcl_N1;
		
		size_t substep_index;
		double dtime;
		size_t pcl_num, task_num;

	public:
		MapPclToBgMesh(Step_T3D_CHM_up_TBB& _stp) : stp(_stp) {}
		void init() noexcept;
		void update(size_t tsk_num) noexcept;
		void work(size_t tsk_id, MapPclToBgMeshRes &res) const;
		inline tbb::task* operator() (tbb::task& parent, size_t tsk_id, MapPclToBgMeshRes &res) const
		{ work(tsk_id, res); return nullptr; }
	};
	
	class FindSoilSurface
	{
	protected:
		typedef Model_T3D_CHM_up_mt::AdjElemIndex AdjElemIndex;

		Step_T3D_CHM_up_TBB& stp;

		const AdjElemIndex* elem_adj_elems;
		const size_t* elem_has_pcls;
		uint16_t *elem_node_at_surface;

		// elem ranges
		const size_t* elem_ids;

		size_t substep_index;
		size_t elem_num, task_num;

	public:
		FindSoilSurface(Step_T3D_CHM_up_TBB& _stp) : stp(_stp) {}
		void init() noexcept;
		void update(size_t tsk_num) noexcept;
		void work(size_t tsk_id) const;
		inline tbb::task* operator() (tbb::task& parent, size_t tsk_id) const { work(tsk_id); return nullptr; }
	};

	class UpdateAccelerationAndVelocity
	{
	protected:
		typedef Model_T3D_CHM_up_mt::Force Force;
		typedef Model_T3D_CHM_up_mt::ElemNodeVM ElemNodeVM;
		typedef Model_T3D_CHM_up_mt::Acceleration Acceleration;
		typedef Model_T3D_CHM_up_mt::Velocity Velocity;
		typedef Model_T3D_CHM_up_mt::NodeHasVBC NodeHasVBC;
		typedef Model_T3D_CHM_up_mt::NodeVBCVec NodeVBCVec;

		Step_T3D_CHM_up_TBB& stp;

		const double* elem_pcl_m;
		const Force* elem_node_force;
		const ElemNodeVM *elem_node_vm_s;
		const double* elem_node_p;
		const uint16_t *elem_node_at_surface;

		const NodeHasVBC *node_has_vbc;
		const NodeVBCVec *node_vbc_vec_s;
		double* node_am;
		Acceleration* node_a_s;
		Velocity* node_v_s;
		double* node_p;
		uint16_t *node_at_surface;

		// node ranges
		const size_t* node_ids;
		const size_t* node_elem_offs;
		
		double dtime;
		size_t four_elem_num, task_num;

	public:
		UpdateAccelerationAndVelocity(Step_T3D_CHM_up_TBB &_stp) : stp(_stp) {}
		void init() noexcept;
		void update(size_t tsk_num) noexcept;
		void work(size_t tsk_id) const;
		inline tbb::task* operator() (tbb::task& parent, size_t tsk_id) const { work(tsk_id); return nullptr; }
	};

	class CalElemDeAndMapToNode
	{
	protected:
		typedef Model_T3D_CHM_up_mt::ElemNodeIndex ElemNodeIndex;
		typedef Model_T3D_CHM_up_mt::DShapeFuncABC DShapeFuncABC;
		typedef Model_T3D_CHM_up_mt::Acceleration Acceleration;
		typedef Model_T3D_CHM_up_mt::Velocity Velocity;
		typedef Model_T3D_CHM_up_mt::Force Force;
		typedef Model_T3D_CHM_up_mt::StrainInc StrainInc;
		typedef Model_T3D_CHM_up_mt::AdjElemIndex AdjElemIndex;

		Step_T3D_CHM_up_TBB &stp;
	
		double k, dyn_viscosity;

		const ElemNodeIndex* elem_node_id;
		const DShapeFuncABC* elem_dN_abc;
		const double* elem_pcl_m;
		const double* elem_pcl_vol;
		const double *elem_density_f;
		StrainInc* elem_de;
		double* elem_m_de_vol_s;

		double* elem_node_p_force;
		uint16_t* elem_node_at_surface;

		const Acceleration* node_a_s;
		const Velocity* node_v_s;
		const double* node_p;

		// elem ranges
		const size_t *elem_ids;

		double dtime;
		size_t elem_num, task_num;

	public:
		CalElemDeAndMapToNode(Step_T3D_CHM_up_TBB &_stp) : stp(_stp) {}
		void init() noexcept;
		void update(size_t tsk_num) noexcept;
		void work(size_t tsk_id) const;
		inline tbb::task* operator() (tbb::task& parent, size_t tsk_id) const { work(tsk_id); return nullptr; }
	};

	class CalNodeDe
	{
	protected:
		typedef Model_T3D_CHM_up_mt::NodeHasVBC NodeHasVBC;

		Step_T3D_CHM_up_TBB &stp;

		const double* elem_pcl_pm;
		const double* elem_node_p_force;
		const NodeHasVBC* node_has_vbc;
		const uint16_t* node_at_surface;
		double* node_dp;

		// strain enhancement
		const double* elem_m_de_vol_s;
		const double* node_am;
		double* node_de_vol_s;

		// node ranges
		const size_t* node_ids;
		const size_t* node_elem_offs;
		
		size_t four_elem_num, task_num;

	public:
		CalNodeDe(Step_T3D_CHM_up_TBB &_stp) : stp(_stp) {}
		void init() noexcept;
		void update(size_t tsk_num) noexcept;
		void work(size_t tsk_id) const;
		inline tbb::task* operator() (tbb::task& parent, size_t tsk_id) const { work(tsk_id); return nullptr; }
	};

	struct MapBgMeshToPclRes
	{
		size_t pcl_num;
		inline void join(const MapBgMeshToPclRes &res)
		{ pcl_num += res.pcl_num; }
	};

	class MapBgMeshToPcl
	{
	protected:
		typedef Model_T3D_CHM_up_mt::ElemNodeIndex ElemNodeIndex;
		typedef Model_T3D_CHM_up_mt::Acceleration Acceleration;
		typedef Model_T3D_CHM_up_mt::Velocity Velocity;
		typedef Model_T3D_CHM_up_mt::Displacement Displacement;
		typedef Model_T3D_CHM_up_mt::Position Position;
		typedef Model_T3D_CHM_up_mt::Stress Stress;
		typedef Model_T3D_CHM_up_mt::Strain Strain;
		typedef Model_T3D_CHM_up_mt::StrainInc StrainInc;
		typedef Model_T3D_CHM_up_mt::ShapeFunc ShapeFunc;
		
		Step_T3D_CHM_up_TBB& stp;

		double Kf0;
		double m_cav, f_cav_end;
		double u_cav_off, u_div_u_cav_lim;

		// pcls
		const Position* pcl_pos;
		MatModel::MaterialModel** pcl_mat_model;
		// bg mesh
		const ElemNodeIndex *elem_node_id;
		const Acceleration *node_a_s;
		const Velocity *node_v_s;
		double* node_dp;
		double* node_p;
		double* elem_p;
		double* elem_pcl_n;
		double * elem_density_f;
		double *elem_u_cav;
		StrainInc *elem_de;
		const double *node_de_vol_s;

		// pcl ranges
		const size_t* pcl_in_elems, *prev_pcl_ids;
		size_t* in_pcl_in_elems, *in_prev_pcl_ids;

		// pcl vars
		const size_t* pcl_index0;
		double* pcl_n0;
		double* pcl_density_f0;
		Velocity* pcl_v_s0;
		Displacement* pcl_u0;
		ShapeFunc* pcl_N0;
		Stress* pcl_stress0;
		double* pcl_p0;
		Strain* pcl_strain0;
		Strain* pcl_estrain0;
		Strain* pcl_pstrain0;
		//
		const Strain* pcl_strain1;
		const Strain* pcl_estrain1;
		const Strain* pcl_pstrain1;

		double dtime;
		size_t pcl_num, task_num;
		
	public:
		MapBgMeshToPcl(Step_T3D_CHM_up_TBB &_stp) : stp(_stp) {}
		void init() noexcept;
		void update(size_t tsk_num) noexcept;
		void work(size_t tsk_id, MapBgMeshToPclRes &res) const;
		inline tbb::task* operator() (tbb::task& parent, size_t tsk_id, MapBgMeshToPclRes &res) const
		{ work(tsk_id, res); return nullptr; }
	};
	
	class ContactRigidBody
	{
	protected:
		typedef Model_T3D_CHM_up_mt::Force Force;
		typedef Model_T3D_CHM_up_mt::Displacement Displacement;
		typedef Model_T3D_CHM_up_mt::Position Position;
		typedef Model_T3D_CHM_up_mt::ShapeFunc ShapeFunc;
		typedef Model_T3D_CHM_up_mt::SortedPclVarArrays SortedPclVarArrays;

		Step_T3D_CHM_up_TBB& stp;

		RigidObjectByT3DMesh* prmesh;

		ContactModel3D* pcm;

		const size_t* pcl_in_elems;

		const Position* pcl_pos;
		const double* pcl_vol;
		Force* elem_node_force;
		uint16_t *elem_node_at_surface;

		const size_t* pcl_index;
		const Displacement* pcl_u;
		const ShapeFunc* pcl_N;

		size_t substep_index;

	public:
		ContactRigidBody(Step_T3D_CHM_up_TBB& _stp) : stp(_stp) {}
		void init(double max_pcl_vol, bool is_first_step) noexcept;
		void update() noexcept;
		inline bool has_rigid_mesh() const noexcept { return prmesh != nullptr; }
		void apply_t3d_rigid_object(size_t p_id0, size_t p_id1, Force3D& rc_cf) const noexcept;
	};
}

#endif