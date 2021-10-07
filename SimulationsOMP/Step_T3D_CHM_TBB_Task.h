#ifndef __Step_T3D_CHM_TBB_Task_h__
#define __Step_T3D_CHM_TBB_Task_h__

#define TBB_SUPPRESS_DEPRECATED_MESSAGES 1
#include <tbb/tbb.h>

#include "RigidObject/RigidCylinder.h"

class Step_T3D_CHM_TBB;

namespace Step_T3D_CHM_TBB_Task
{
	constexpr size_t task_num_per_thread = 4;

	constexpr size_t min_pcl_num_per_init_pcl_task = 100;
	constexpr size_t min_pcl_num_per_task = 100;
	constexpr size_t min_node_elem_num_per_task = 20;
	constexpr size_t min_elem_num_per_task = 20;

	constexpr size_t cache_line_size = 64;

	struct PclRange
	{
		size_t p_id0, p_id1;
		char padding[cache_line_size];
	};

	struct NodeElemRange
	{
		size_t ve_id0, ve_id1;
		char padding[cache_line_size];
	};
	
	class InitPcl;
	struct InitPclRes
	{
		size_t pcl_num;
		InitPcl* init_pcl;
		// for initializeing
		InitPclRes(InitPcl &_ip) : pcl_num(0), init_pcl(&_ip) {}
		// tbb::parallel_for
		InitPclRes(InitPclRes &other, tbb::split) :
			pcl_num(0), init_pcl(other.init_pcl) {}
		void operator() (tbb::blocked_range<size_t>& range);
		inline void join(const InitPclRes &other)
		{ pcl_num += other.pcl_num; }
		// parallel_for
		InitPclRes() : pcl_num(0), init_pcl(nullptr) {}
	};

	class InitPcl
	{
	protected:
		typedef Model_T3D_CHM_mt::ShapeFunc ShapeFunc;
		typedef Model_T3D_CHM_mt::Displacement Displacement;
		typedef Model_T3D_CHM_mt::Position Position;
		
		Step_T3D_CHM_TBB& stp;
		size_t *in_pcl_in_elems, *in_prev_pcl_ids;
		size_t task_num;

	public:
		InitPcl(Step_T3D_CHM_TBB &_stp) : stp(_stp) {}
		void init(size_t thread_num) noexcept;
		size_t work(size_t wk_id) const;
		inline tbb::task* operator() (tbb::task& parent, size_t wk_id, InitPclRes &res) const
		{ res.pcl_num = work(wk_id); return nullptr; }
		inline size_t get_task_num() const noexcept { return task_num; }
	};
	
	class MapPclToBgMesh
	{
	protected:
		typedef Model_T3D_CHM_mt::Force Force;
		typedef Model_T3D_CHM_mt::Acceleration Acceleration;
		typedef Model_T3D_CHM_mt::Velocity Velocity;
		typedef Model_T3D_CHM_mt::Displacement Displacement;
		typedef Model_T3D_CHM_mt::Stress Stress;
		typedef Model_T3D_CHM_mt::ShapeFunc ShapeFunc;
		typedef Model_T3D_CHM_mt::DShapeFuncABC DShapeFuncABC;
		typedef Model_T3D_CHM_mt::DShapeFuncD DShapeFuncD;
		typedef Model_T3D_CHM_mt::ElemNodeVM ElemNodeVM;

		Step_T3D_CHM_TBB &stp;

		// pcl data
		const double* pcl_m_s;
		const double* pcl_vol_s;
		const Force* pcl_bf_s;
		const Force* pcl_bf_f;
		const Force* pcl_t;
		double* pcl_vol;
		// bg mesh data
		const DShapeFuncABC *elem_N_abc;
		const double *elem_vol;
		double* elem_pcl_m_s;
		double* elem_pcl_m_f;
		double* elem_pcl_n;
		double* elem_density_f;
		double* elem_p;
		ElemNodeVM* elem_node_vm_s;
		ElemNodeVM* elem_node_vm_f;
		Force* elem_node_force_s;
		Force* elem_node_force_f;
		// pcl range
		const size_t* pcl_in_elems;
		const size_t* prev_pcl_ids;
		PclRange* pcl_ranges;

		// pcl_vars0
		size_t *pcl_index0;
		double* pcl_n0;
		double *pcl_density_f0;
		Velocity* pcl_v_s0;
		Velocity *pcl_v_f0;
		Displacement* pcl_u_s0;
		Displacement* pcl_u_f0;
		Stress *pcl_stress0;
		ShapeFunc *pcl_N0;
		// pcl_vars1
		const size_t *pcl_index1;
		const double* pcl_n1;
		const double *pcl_density_f1;
		const Velocity* pcl_v_s1;
		const Velocity* pcl_v_f1;
		const Displacement *pcl_u_s1;
		const Displacement* pcl_u_f1;
		const Stress *pcl_stress1;
		double* pcl_p1;
		const ShapeFunc *pcl_N1;

		size_t task_num;

	public:
		MapPclToBgMesh(Step_T3D_CHM_TBB &_stp) : stp(_stp) {}
		void init() noexcept;
		void update(size_t tsk_num) noexcept;
		void work(size_t wk_id) const;
		inline void operator() (const tbb::blocked_range<size_t>& range) const { work(range.begin()); }
		inline tbb::task* operator() (tbb::task& parent, size_t wk_id) const { work(wk_id); return nullptr; }
	};

	class ContactRigidCylinder;
	struct ContactForceRes
	{
		Force3D react_force;
		ContactRigidCylinder* contact_rigid_cylinder;
		// initializing
		ContactForceRes(ContactRigidCylinder &_cont) :
			contact_rigid_cylinder(&_cont) { react_force.reset(); }
		// tbb::parallel_reduce
		ContactForceRes(ContactForceRes &other, tbb::split) :
			contact_rigid_cylinder(other.contact_rigid_cylinder) { react_force.reset(); }
		inline void join(ContactForceRes& res)
		{ react_force.combine(res.react_force); }
		void operator() (const tbb::blocked_range<size_t>& range);
		// parallel_reduce
		ContactForceRes() : contact_rigid_cylinder(nullptr) { react_force.reset(); }
	};

	class ContactRigidCylinder
	{
	protected:
		typedef Model_T3D_CHM_mt::Force Force;
		typedef Model_T3D_CHM_mt::Position Position;
		typedef Model_T3D_CHM_mt::Displacement Displacement;
		typedef Model_T3D_CHM_mt::ShapeFunc ShapeFunc;

		Step_T3D_CHM_TBB& stp;

		RigidCylinder *prcy;
		ContactModel3D* pcm_s, * pcm_f;

		const Position *pcl_pos;
		const double *pcl_vol;
		Force *elem_node_force_s;
		Force *elem_node_force_f;

		const size_t* pcl_index;
		const ShapeFunc* pcl_N;
		const Displacement* pcl_u_s;

		const size_t* pcl_in_elems;
		PclRange* pcl_ranges;

	public:
		ContactRigidCylinder(Step_T3D_CHM_TBB &_stp) : stp(_stp) {}
		void init(Model_T3D_CHM_mt& md) noexcept;
		void update() noexcept;
		Force3D work(size_t wk_id) const;
		inline tbb::task *operator() (tbb::task &parent, size_t wk_id, ContactForceRes& rcy_cf) const
		{ rcy_cf.react_force = work(wk_id); return nullptr; }
	};
	
	class UpdateAccelerationAndVelocity
	{
	protected:
		typedef Model_T3D_CHM_mt::Force Force;
		typedef Model_T3D_CHM_mt::ElemNodeVM ElemNodeVM;
		typedef Model_T3D_CHM_mt::Acceleration Acceleration;
		typedef Model_T3D_CHM_mt::Velocity Velocity;
		typedef Model_T3D_CHM_mt::NodeHasVBC NodeHasVBC;

		Step_T3D_CHM_TBB& stp;

		const double* elem_pcl_m_s;
		const double* elem_pcl_m_f;
		const Force* elem_node_force_s;
		const Force* elem_node_force_f;
		const ElemNodeVM *elem_node_vm_s;
		const ElemNodeVM* elem_node_vm_f;
		double *node_am_s;
		double* node_am_f;
		Acceleration* node_a_s;
		Acceleration* node_a_f;
		Velocity *node_v_s;
		Velocity* node_v_f;
		NodeHasVBC *node_has_vbc_s;
		NodeHasVBC* node_has_vbc_f;

		const size_t* node_ids;
		const size_t* node_elem_offs;
		NodeElemRange *node_elem_ranges;
		
		size_t four_elem_num, task_num;

	public:
		UpdateAccelerationAndVelocity(Step_T3D_CHM_TBB& _stp) : stp(_stp) {}
		void init() noexcept;
		void update(size_t tsk_num);
		void work(size_t wk_id) const;
		inline void operator() (const tbb::blocked_range<size_t>& range) const
		{ work(range.begin()); }
		inline tbb::task* operator() (tbb::task& parent, size_t wk_id) const
		{ work(wk_id); return nullptr; }
	};

	class CalElemDeAndMapToNode
	{
	protected:
		typedef Model_T3D_CHM_mt::ElemNodeIndex ElemNodeIndex;
		typedef Model_T3D_CHM_mt::Velocity Velocity;
		typedef Model_T3D_CHM_mt::DShapeFuncABC DShapeFuncABC;
		typedef Model_T3D_CHM_mt::StrainInc StrainInc;

		Step_T3D_CHM_TBB& stp;

		const ElemNodeIndex* elem_node_id;
		const DShapeFuncABC *elem_N_abc;
		const double* elem_pcl_m_s, *elem_pcl_m_f;
		const double* elem_pcl_n;
		const Velocity *node_v_s, *node_v_f;
		StrainInc* elem_de;
		double* elem_m_de_vol_s, *elem_m_de_vol_f;

		const size_t *elem_ids;

		size_t elem_num, task_num;

	public:
		CalElemDeAndMapToNode(Step_T3D_CHM_TBB &_stp) : stp(_stp) {}
		void init() noexcept;
		void update(size_t tsk_num);
		void work(size_t wk_id) const;
		inline void operator() (const tbb::blocked_range<size_t>& range) const
		{ work(range.begin()); }
		inline tbb::task* operator() (tbb::task& parent, size_t wk_id) const
		{ work(wk_id); return nullptr; }
	};

	class CalNodeDe
	{
	protected:
		Step_T3D_CHM_TBB& stp;
		
		const double* elem_m_de_vol_s, *elem_m_de_vol_f;
		const double* node_am_s, *node_am_f;
		double* node_de_vol_s, *node_de_vol_f;
		
		const size_t* node_ids, * node_elem_offs;
		NodeElemRange* node_elem_ranges;

	public:
		CalNodeDe(Step_T3D_CHM_TBB &_stp) : stp(_stp) {}
		void init() noexcept;
		void work(size_t wk_id) const;
		inline void operator() (const tbb::blocked_range<size_t>& range) const
		{ work(range.begin()); }
		inline tbb::task* operator() (tbb::task& parent, size_t wk_id) const
		{ work(wk_id); return nullptr; }
	};

	class MapBgMeshToPcl;
	struct MapBgMeshToPclRes
	{
		size_t pcl_num;
		MapBgMeshToPcl* map_bg_mesh_to_pcl;
		// initializing
		MapBgMeshToPclRes(MapBgMeshToPcl &_map) :
			pcl_num(0), map_bg_mesh_to_pcl(&_map) {}
		// parallel_reduce
		MapBgMeshToPclRes() : pcl_num(0), map_bg_mesh_to_pcl(nullptr) {}
		// tbb::parallel_reduce
		MapBgMeshToPclRes(MapBgMeshToPclRes& other, tbb::split)
			: pcl_num(0), map_bg_mesh_to_pcl(other.map_bg_mesh_to_pcl) {}
		inline void operator() (const tbb::blocked_range<size_t>& range);
		inline void join(const MapBgMeshToPclRes &other)
		{ pcl_num += other.pcl_num; }
	};

	class MapBgMeshToPcl
	{
	protected:
		typedef Model_T3D_CHM_mt::ElemNodeIndex ElemNodeIndex;
		typedef Model_T3D_CHM_mt::Acceleration Acceleration;
		typedef Model_T3D_CHM_mt::Velocity Velocity;
		typedef Model_T3D_CHM_mt::Displacement Displacement;
		typedef Model_T3D_CHM_mt::Position Position;
		typedef Model_T3D_CHM_mt::Stress Stress;
		typedef Model_T3D_CHM_mt::Strain Strain;
		typedef Model_T3D_CHM_mt::StrainInc StrainInc;
		typedef Model_T3D_CHM_mt::ShapeFunc ShapeFunc;
		
		Step_T3D_CHM_TBB& stp;

		const Position* pcl_pos;
		MatModel::MaterialModel** pcl_mat_model;
		
		const ElemNodeIndex *elem_node_id;
		const Acceleration *node_a_s, *node_a_f;
		const Velocity *node_v_s, *node_v_f;
		double* elem_density_f;
		double* elem_pcl_n;
		double* elem_p;
		StrainInc* elem_de;
		const double *node_de_vol_s, *node_de_vol_f;

		const size_t *pcl_index0;
		double* pcl_density_f0;
		double* pcl_n0;
		double* pcl_p0;
		ShapeFunc* pcl_N0;
		Velocity* pcl_v_s0;
		Velocity* pcl_v_f0;
		Displacement* pcl_u_s0;
		Displacement* pcl_u_f0;
		Stress* pcl_stress0;
		Strain* pcl_strain0;
		Strain* pcl_estrain0;
		Strain* pcl_pstrain0;
		const Strain* pcl_strain1;
		const Strain* pcl_estrain1;
		const Strain* pcl_pstrain1;

		const size_t* pcl_in_elems;
		const size_t* prev_pcl_ids;
		PclRange* pcl_ranges;
		size_t* in_pcl_in_elems;
		size_t* in_prev_pcl_ids;
		
	public:
		MapBgMeshToPcl(Step_T3D_CHM_TBB &_stp) : stp(_stp) {}
		void init() noexcept;
		void update(size_t thread_num) noexcept;
		size_t work(size_t wk_id) const;
		inline tbb::task *operator() (tbb::task &parent, size_t wk_id, MapBgMeshToPclRes &res) const
		{ res.pcl_num = work(wk_id); return nullptr; }
	};
}

#endif