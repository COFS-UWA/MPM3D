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
	
	struct InitPclRes
	{
		size_t pcl_num;
		inline void join(const InitPclRes &other)
		{ pcl_num += other.pcl_num; }
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
		void work(size_t tsk_id, InitPclRes &res) const;
		inline tbb::task* operator() (tbb::task& parent, size_t tsk_id, InitPclRes &res) const
		{ work(tsk_id, res); return nullptr; }
		inline size_t get_task_num() const noexcept { return task_num; }
	};
	
	struct MapPclToBgMeshRes
	{
		Force3D react_force;
		inline void join(const MapPclToBgMeshRes &res)
		{ react_force.combine(res.react_force); }
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
		void work(size_t tsk_id, MapPclToBgMeshRes &res) const;
		inline tbb::task* operator() (tbb::task& parent, size_t tsk_id, MapPclToBgMeshRes &res) const
		{ work(tsk_id, res); return nullptr; }
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
		void work(size_t tsk_id) const;
		inline void operator() (const tbb::blocked_range<size_t>& range) const
		{ work(range.begin()); }
		inline tbb::task* operator() (tbb::task& parent, size_t tsk_id) const
		{ work(tsk_id); return nullptr; }
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
		void work(size_t tsk_id) const;
		inline void operator() (const tbb::blocked_range<size_t>& range) const
		{ work(range.begin()); }
		inline tbb::task* operator() (tbb::task& parent, size_t tsk_id) const
		{ work(tsk_id); return nullptr; }
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
		void work(size_t tsk_id) const;
		inline void operator() (const tbb::blocked_range<size_t>& range) const
		{ work(range.begin()); }
		inline tbb::task* operator() (tbb::task& parent, size_t tsk_id) const
		{ work(tsk_id); return nullptr; }
	};

	struct MapBgMeshToPclRes
	{
		size_t pcl_num;
		inline void join(const MapBgMeshToPclRes &other) noexcept
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
		void work(size_t tsk_id, MapBgMeshToPclRes &res) const;
		inline tbb::task *operator() (tbb::task &parent, size_t tsk_id, MapBgMeshToPclRes &res) const
		{ work(tsk_id, res); return nullptr; }
	};

	class ContactRigidBody
	{
	protected:
		typedef Model_T3D_CHM_mt::Force Force;
		typedef Model_T3D_CHM_mt::Displacement Displacement;
		typedef Model_T3D_CHM_mt::Position Position;
		typedef Model_T3D_CHM_mt::ShapeFunc ShapeFunc;
		typedef Model_T3D_CHM_mt::SortedPclVarArrays SortedPclVarArrays;

		Step_T3D_CHM_TBB& stp;

		RigidCylinder* prcy;
		RigidObjectByT3DMesh* prmesh;

		ContactModel3D* pcm_s, *pcm_f;

		const Position* pcl_pos;
		const double* pcl_vol;
		Force* elem_node_force_s, *elem_node_force_f;

		const size_t* pcl_index;
		const Displacement* pcl_u_s, *pcl_u_f;
		const ShapeFunc* pcl_N;

		const PclRange* pcl_ranges;
		const size_t* pcl_in_elems;

	public:
		ContactRigidBody(Step_T3D_CHM_TBB& _stp) : stp(_stp) {}
		void init() noexcept;
		void update() noexcept;
		inline bool has_rigid_cylinder() const noexcept { return prcy != nullptr; }
		inline bool has_rigid_mesh() const noexcept { return prmesh != nullptr; }
		void apply_rigid_cylinder(size_t p_id0, size_t p_id1, Force3D& rc_cf) const noexcept;
		void apply_t3d_rigid_object(size_t p_id0, size_t p_id1, Force3D& rc_cf) const noexcept;
	};

	// tbb::parallel_reduce
	struct InitPclTbb
	{
		InitPclRes res;
		InitPcl& worker;
		inline void reset() { res.pcl_num = 0; }
		InitPclTbb(InitPcl& _ip) : worker(_ip) {}
		InitPclTbb(InitPclTbb& other, tbb::split) : worker(other.worker) { reset(); }
		inline void operator() (const tbb::blocked_range<size_t>& range)
		{
			InitPclRes tmp;
			worker.work(range.begin(), tmp);
			res.join(tmp);
		}
		inline void join(const InitPclTbb& other) noexcept { res.join(other.res); }
	};

	struct MapPclToBgMeshTbb
	{
		MapPclToBgMeshRes res;
		MapPclToBgMesh& worker;
		inline void reset() { res.react_force.reset(); }
		MapPclToBgMeshTbb(MapPclToBgMesh& _mptm) : worker(_mptm) {}
		MapPclToBgMeshTbb(MapPclToBgMeshTbb& other, tbb::split) : worker(other.worker) { reset(); }
		inline void operator() (const tbb::blocked_range<size_t>& range)
		{
			MapPclToBgMeshRes tmp;
			worker.work(range.begin(), tmp);
			res.join(tmp);
		}
		inline void join(const MapPclToBgMeshTbb& other) noexcept { res.join(other.res); }
	};

	struct MapBgMeshToPclTbb
	{
		MapBgMeshToPclRes res;
		MapBgMeshToPcl& worker;
		inline void reset() { res.pcl_num = 0; }
		MapBgMeshToPclTbb(MapBgMeshToPcl& _mmtp) : worker(_mmtp) {}
		MapBgMeshToPclTbb(MapBgMeshToPclTbb& other, tbb::split) : worker(other.worker) { reset(); }
		inline void operator() (const tbb::blocked_range<size_t>& range)
		{
			MapBgMeshToPclRes tmp;
			worker.work(range.begin(), tmp);
			res.join(tmp);
		}
		inline void join(const MapBgMeshToPclTbb& other) noexcept { res.join(other.res); }
	};
}

#endif