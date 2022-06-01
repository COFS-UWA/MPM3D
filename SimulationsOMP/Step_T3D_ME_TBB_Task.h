#ifndef __Step_T3D_ME_Task_h__
#define __Step_T3D_ME_Task_h__

#define TBB_SUPPRESS_DEPRECATED_MESSAGES 1
#include <tbb/tbb.h>

#include "RigidObject/RigidCylinder.h"
#include "RigidObject/RigidCone.h"
#include "RigidObject/RigidCube.h"
#include "RigidObject/RigidObjectByT3DMesh.h"

class Step_T3D_ME_TBB;

namespace Step_T3D_ME_TBB_Task
{
	constexpr size_t init_pcl_task_num_per_thread = 100;
	constexpr size_t map_pcl_to_mesh_task_num_per_thread = 100;
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
		typedef Model_T3D_ME_mt::ShapeFunc ShapeFunc;
		typedef Model_T3D_ME_mt::Displacement Displacement;
		typedef Model_T3D_ME_mt::Position Position;
		
		Step_T3D_ME_TBB &stp;
		const double* pcl_m;
		size_t* in_pcl_in_elems, * in_prev_pcl_ids;
		size_t task_num;
	
	public:
		InitPcl(Step_T3D_ME_TBB &_stp) : stp(_stp) {}
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

	class MapPclToBgMesh
	{
	protected:
		typedef Model_T3D_ME_mt::Force Force;
		typedef Model_T3D_ME_mt::Acceleration Acceleration;
		typedef Model_T3D_ME_mt::Velocity Velocity;
		typedef Model_T3D_ME_mt::Displacement Displacement;
		typedef Model_T3D_ME_mt::Position Position;
		typedef Model_T3D_ME_mt::Stress Stress;
		typedef Model_T3D_ME_mt::ShapeFunc ShapeFunc;
		typedef Model_T3D_ME_mt::DShapeFuncABC DShapeFuncABC;
		typedef Model_T3D_ME_mt::DShapeFuncD DShapeFuncD;
		typedef Model_T3D_ME_mt::ElemNodeVM ElemNodeVM;

		Step_T3D_ME_TBB& stp;

		// pcl data
		const double* pcl_m;
		const Force* pcl_bf;
		const Force* pcl_t;
		double* pcl_vol;
		// bg mesh data
		const DShapeFuncABC *elem_dN_abc;
		const double *elem_vol;
		double* elem_pcl_m;
		double* elem_density;
		ElemNodeVM* elem_node_vm;
		Force* elem_node_force;
		// pcl range
		const size_t* pcl_in_elems;
		const size_t* prev_pcl_ids;

		// pcl_vars0
		size_t *pcl_index0;
		double *pcl_density0;
		Velocity *pcl_v0;
		Displacement *pcl_disp0;
		Stress *pcl_stress0;
		ShapeFunc *pcl_N0;
		// pcl_vars1
		const size_t *pcl_index1;
		const double *pcl_density1;
		const Velocity *pcl_v1;
		const Displacement *pcl_disp1;
		const Stress *pcl_stress1;
		const ShapeFunc *pcl_N1;
		
		size_t task_num;

	public:
		MapPclToBgMesh(Step_T3D_ME_TBB& _stp) : stp(_stp) {}
		void init() noexcept;
		void update(size_t tsk_num) noexcept;
		void work(size_t tsk_id, MapPclToBgMeshRes &res) const;
		inline tbb::task* operator() (tbb::task& parent, size_t tsk_id, MapPclToBgMeshRes &res) const
		{ work(tsk_id, res); return nullptr; }
	};
	
	class UpdateAccelerationAndVelocity
	{
	protected:
		typedef Model_T3D_ME_mt::Force Force;
		typedef Model_T3D_ME_mt::ElemNodeVM ElemNodeVM;
		typedef Model_T3D_ME_mt::Acceleration Acceleration;
		typedef Model_T3D_ME_mt::Velocity Velocity;
		typedef Model_T3D_ME_mt::NodeHasVBC NodeHasVBC;
		typedef Model_T3D_ME_mt::NodeVBCVec NodeVBCVec;

		Step_T3D_ME_TBB& stp;

		const double* elem_pcl_m;
		const Force* elem_node_force;
		const ElemNodeVM *elem_node_vm;
		Acceleration *node_a;
		double *node_am;
		Velocity *node_v;
		const NodeHasVBC *node_has_vbc;
		const NodeVBCVec *node_vbc_vec;
		// node ranges
		const size_t* node_ids;
		const size_t* node_elem_offs;
		
		size_t four_elem_num, task_num;

	public:
		UpdateAccelerationAndVelocity(Step_T3D_ME_TBB &_stp) : stp(_stp) {}
		void init() noexcept;
		void update(size_t tsk_num) noexcept;
		void work(size_t tsk_id) const;
		inline void operator() (const tbb::blocked_range<size_t>& range) const { work(range.begin()); }
		inline tbb::task* operator() (tbb::task& parent, size_t tsk_id) const { work(tsk_id); return nullptr; }
	};

	class CalElemDeAndMapToNode
	{
	protected:
		typedef Model_T3D_ME_mt::ElemNodeIndex ElemNodeIndex;
		typedef Model_T3D_ME_mt::Velocity Velocity;
		typedef Model_T3D_ME_mt::DShapeFuncABC DShapeFuncABC;
		typedef Model_T3D_ME_mt::StrainInc StrainInc;

		Step_T3D_ME_TBB &stp;

		const ElemNodeIndex* elem_node_id;
		const double* elem_pcl_m;
		const DShapeFuncABC *elem_dN_abc; 
		const Velocity* node_v;
		StrainInc* elem_de;
		double* elem_m_de_vol;
		// elem ranges
		const size_t *elem_ids;

		size_t elem_num, task_num;

	public:
		CalElemDeAndMapToNode(Step_T3D_ME_TBB &_stp) : stp(_stp) {}
		void init() noexcept;
		void update(size_t tsk_num) noexcept;
		void work(size_t tsk_id) const;
		inline void operator() (const tbb::blocked_range<size_t>& range) const { work(range.begin()); }
		inline tbb::task* operator() (tbb::task& parent, size_t tsk_id) const { work(tsk_id); return nullptr; }
	};

	class CalNodeDe
	{
	protected:
		Step_T3D_ME_TBB &stp;

		const double* elem_m_de_vol;
		const double* node_am;
		double* node_de_vol;
		// node ranges
		const size_t* node_ids;
		const size_t* node_elem_offs;
		
		size_t four_elem_num, task_num;

	public:
		CalNodeDe(Step_T3D_ME_TBB &_stp) : stp(_stp) {}
		void init() noexcept;
		void update(size_t tsk_num) noexcept;
		void work(size_t tsk_id) const;
		inline void operator() (const tbb::blocked_range<size_t>& range) const { work(range.begin()); }
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
		typedef Model_T3D_ME_mt::ElemNodeIndex ElemNodeIndex;
		typedef Model_T3D_ME_mt::Acceleration Acceleration;
		typedef Model_T3D_ME_mt::Velocity Velocity;
		typedef Model_T3D_ME_mt::Displacement Displacement;
		typedef Model_T3D_ME_mt::Position Position;
		typedef Model_T3D_ME_mt::Stress Stress;
		typedef Model_T3D_ME_mt::Strain Strain;
		typedef Model_T3D_ME_mt::StrainInc StrainInc;
		typedef Model_T3D_ME_mt::ShapeFunc ShapeFunc;
		
		Step_T3D_ME_TBB& stp;

		// pcls
		const Position* pcl_pos;
		MatModel::MaterialModel** pcl_mat_model;
		// bg mesh
		const ElemNodeIndex *elem_node_id;
		const Acceleration *node_a;
		const Velocity *node_v;
		double *elem_density;
		StrainInc *elem_de;
		const double *node_de_vol;
		// pcl vars
		const size_t *pcl_index0;
		double* pcl_density0;
		Velocity* pcl_v0;
		Displacement* pcl_disp0;
		ShapeFunc*pcl_N0;
		Stress* pcl_stress0;
		Strain* pcl_strain0;
		Strain* pcl_estrain0;
		Strain* pcl_pstrain0;
		//
		const Strain* pcl_strain1;
		const Strain* pcl_estrain1;
		const Strain* pcl_pstrain1;
		// pcl ranges
		const size_t* pcl_in_elems, *prev_pcl_ids;
		size_t* in_pcl_in_elems, *in_prev_pcl_ids;

		size_t pcl_num, task_num;
		
	public:
		MapBgMeshToPcl(Step_T3D_ME_TBB &_stp) : stp(_stp) {}
		void init() noexcept;
		void update(size_t tsk_num) noexcept;
		void work(size_t tsk_id, MapBgMeshToPclRes &res) const;
		inline tbb::task* operator() (tbb::task& parent, size_t tsk_id, MapBgMeshToPclRes &res) const
		{ work(tsk_id, res); return nullptr; }
	};
	
	class ContactRigidBody
	{
	protected:
		typedef Model_T3D_ME_mt::Force Force;
		typedef Model_T3D_ME_mt::Displacement Displacement;
		typedef Model_T3D_ME_mt::Position Position;
		typedef Model_T3D_ME_mt::ShapeFunc ShapeFunc;
		typedef Model_T3D_ME_mt::SortedPclVarArrays SortedPclVarArrays;

		Step_T3D_ME_TBB& stp;

		RigidCone* prco;
		RigidCube* prcu;
		RigidCylinder* prcy;
		RigidObjectByT3DMesh* prmesh;

		ContactModel3D* pcm;

		const size_t* pcl_in_elems;

		const Position* pcl_pos;
		const double* pcl_vol;
		Force* elem_node_force;

		const size_t* pcl_index;
		const Displacement* pcl_disp;
		const ShapeFunc* pcl_N;

	public:
		ContactRigidBody(Step_T3D_ME_TBB& _stp) : stp(_stp) {}
		void init(double max_pcl_vol) noexcept;
		void update() noexcept;
		inline bool has_rigid_cone() const noexcept { return prco != nullptr; }
		inline bool has_rigid_cube() const noexcept { return prcu != nullptr; }
		inline bool has_rigid_cylinder() const noexcept { return prcy != nullptr; }
		inline bool has_rigid_mesh() const noexcept { return prmesh != nullptr; }
		void apply_rigid_cone(size_t p_id0, size_t p_id1, Force3D& rc_cf) const noexcept;
		void apply_rigid_cube(size_t p_id0, size_t p_id1, Force3D& rc_cf) const noexcept;
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
		MapPclToBgMesh &worker;
		inline void reset() { res.react_force.reset(); }
		MapPclToBgMeshTbb(MapPclToBgMesh &_mptm) : worker(_mptm) {}
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
		MapBgMeshToPcl &worker;
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

	class MapBgMeshToPcl0
	{
	protected:
		typedef Model_T3D_ME_mt::ElemNodeIndex ElemNodeIndex;
		typedef Model_T3D_ME_mt::Acceleration Acceleration;
		typedef Model_T3D_ME_mt::Velocity Velocity;
		typedef Model_T3D_ME_mt::Displacement Displacement;
		typedef Model_T3D_ME_mt::Position Position;
		typedef Model_T3D_ME_mt::ShapeFunc ShapeFunc;

		Step_T3D_ME_TBB& stp;

		// pcls
		const Position* pcl_pos;
		// bg mesh
		const ElemNodeIndex* elem_node_id;
		const Acceleration* node_a;
		const Velocity* node_v;
		double* elem_density;
		const double* node_de_vol;
		// pcl vars
		const size_t* pcl_index0;
		double* pcl_density0;
		Velocity* pcl_v0;
		Displacement* pcl_disp0;
		ShapeFunc* pcl_N0;

		// pcl ranges
		const size_t* pcl_in_elems, * prev_pcl_ids;
		size_t* in_pcl_in_elems, * in_prev_pcl_ids;

		size_t pcl_num, task_num;

	public:
		MapBgMeshToPcl0(Step_T3D_ME_TBB& _stp) : stp(_stp) {}
		void init() noexcept;
		void update(size_t tsk_num) noexcept;
		void work(size_t tsk_id, MapBgMeshToPclRes& res) const;
		inline tbb::task* operator() (tbb::task& parent, size_t tsk_id, MapBgMeshToPclRes &res) const
		{ work(tsk_id, res); return nullptr; }
	};

	class MapBgMeshToPcl1
	{
	protected:
		typedef Model_T3D_ME_mt::ElemNodeIndex ElemNodeIndex;
		typedef Model_T3D_ME_mt::Acceleration Acceleration;
		typedef Model_T3D_ME_mt::Velocity Velocity;
		typedef Model_T3D_ME_mt::Displacement Displacement;
		typedef Model_T3D_ME_mt::Position Position;
		typedef Model_T3D_ME_mt::Stress Stress;
		typedef Model_T3D_ME_mt::Strain Strain;
		typedef Model_T3D_ME_mt::StrainInc StrainInc;
		typedef Model_T3D_ME_mt::ShapeFunc ShapeFunc;

		Step_T3D_ME_TBB& stp;

		// pcls
		MatModel::MaterialModel** pcl_mat_model;
		// bg mesh
		const ElemNodeIndex* elem_node_id;
		StrainInc* elem_de;
		const double* node_de_vol;
		// pcl vars
		const size_t* pcl_index0;
		Stress* pcl_stress0;
		Strain* pcl_strain0;
		Strain* pcl_estrain0;
		Strain* pcl_pstrain0;
		//
		const Strain* pcl_strain1;
		const Strain* pcl_estrain1;
		const Strain* pcl_pstrain1;
		// pcl ranges
		const size_t* pcl_in_elems, * prev_pcl_ids;

		size_t pcl_num, task_num;

	public:
		MapBgMeshToPcl1(Step_T3D_ME_TBB& _stp) : stp(_stp) {}
		void init() noexcept;
		void update(size_t tsk_num) noexcept;
		void work(size_t tsk_id) const;
		inline tbb::task* operator() (tbb::task& parent, size_t tsk_id) const { work(tsk_id); return nullptr; }
	};
}

#endif