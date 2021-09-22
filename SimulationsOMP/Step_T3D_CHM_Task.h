#ifndef __Step_T3D_CHM_Task_h__
#define __Step_T3D_CHM_Task_h__

#define TBB_SUPPRESS_DEPRECATED_MESSAGES 1
#include <tbb/tbb.h>

#include "ParallelUtils.h"
#include "MSDRadixSortUtils.h"
#include "SortParticleTask.h"
#include "SortTehMeshNodeTask.h"
#include "Model_T3D_CHM_mt.h"

namespace Step_T3D_CHM_Task
{
	using MSDRadixSortUtils::block_low;

	constexpr size_t task_num_per_thread = 4;
	constexpr size_t min_pcl_num_per_init_pcl_task = 100;
	constexpr size_t min_pcl_num_per_map_pcl_to_mesh_task = 100;
	constexpr size_t min_pcl_num_per_contact_rigid_rect_task = 100;
	constexpr size_t min_elem_num_per_update_a_and_v_task = 20;
	constexpr size_t min_elem_num_per_cal_elem_de_task = 20;
	constexpr size_t min_elem_num_per_cal_node_de_task = 20;
	constexpr size_t min_pcl_num_per_map_mesh_to_pcl_task = 100;

	struct CalData
	{
		using Force = Model_T3D_CHM_mt::Force;
		using Position = Model_T3D_CHM_mt::Position;
		using SortedPclVarArrays = Model_T3D_CHM_mt::SortedPclVarArrays;
		using ElemNodeIndex = Model_T3D_CHM_mt::ElemNodeIndex;
		using DShapeFuncABC = Model_T3D_CHM_mt::DShapeFuncABC;
		using DShapeFuncD = Model_T3D_CHM_mt::DShapeFuncD;
		using StrainInc = Model_T3D_CHM_mt::StrainInc;
		using ElemNodeVM = Model_T3D_CHM_mt::ElemNodeVM;
		using Acceleration = Model_T3D_CHM_mt::Acceleration;
		using Velocity = Model_T3D_CHM_mt::Velocity;
		using NodeHasVBC = Model_T3D_CHM_mt::NodeHasVBC;

		Model_T3D_CHM_mt *pmodel;

		// pcl data
		const double* pcl_m_s;
		const double* pcl_density_s;
		const double* pcl_vol_s;
		const Force* pcl_bf_s;
		const Force* pcl_bf_f;
		const Force* pcl_t;
		const Position* pcl_pos;
		double *pcl_vol;
		MatModel::MaterialModel** pcl_mat_model;

		SortedPclVarArrays spvas[2];

		// elem data
		const ElemNodeIndex* elem_node_id;
		const DShapeFuncABC* elem_N_abc;
		const DShapeFuncD* elem_N_d;
		const double* elem_vol;

		// element calculation data
		double* elem_density_f;
		double* elem_pcl_n;
		double* elem_pcl_m_s;
		double* elem_pcl_m_f;
		StrainInc* elem_de;
		double* elem_p;
		double* elem_n2_miu_div_k_vol;
		Force* elem_seep_force;
		double* elem_m_de_vol_s;
		double* elem_m_de_vol_f;

		// element-node data
		ElemNodeVM* elem_node_vm_s;
		ElemNodeVM* elem_node_vm_f;
		Force* elem_node_force_s;
		Force* elem_node_force_f;

		// node data
		Position* node_pos;
		Acceleration* node_a_s;
		Acceleration* node_a_f;
		Velocity* node_v_s;
		Velocity* node_v_f;
		NodeHasVBC* node_has_vbc_s;
		NodeHasVBC* node_has_vbc_f;
		double* node_am_s;
		double* node_am_f;
		double* node_de_vol_s;
		double* node_de_vol_f;

		double Kf, miu, k;

#ifdef _DEBUG
		size_t elem_num;
		size_t node_num;
		size_t ori_pcl_num;
#endif
		
		// cal data
		size_t thread_num;
		MSDRadixSortUtils::RadixBinBlockMemArray thread_bin_blocks_mem;
		SortParticleMem pcl_sort_mem;
		SortTehMeshNodeMem node_sort_mem;
		
		// data changed during computation
		double dt;
		size_t sorted_pcl_var_id;
		size_t prev_valid_pcl_num;
		size_t valid_pcl_num;
		size_t valid_elem_num;

		void set_model(Model_T3D_CHM_mt& md) noexcept;
	};
	
	class InitPcl;
	struct InitPclRes
	{
		size_t pcl_num;
		InitPcl* init_pcl;
		InitPclRes() : pcl_num(0), init_pcl(nullptr) {}
		InitPclRes(InitPclRes &other, tbb::split) :
			pcl_num(0), init_pcl(other.init_pcl) {}
		void operator() (tbb::blocked_range<size_t>& range);
		inline void join(const InitPclRes &other)
		{ pcl_num += other.pcl_num; }
	};

	class InitPcl
	{
	protected:
		typedef Model_T3D_CHM_mt::ShapeFunc ShapeFunc;
		typedef Model_T3D_CHM_mt::Displacement Displacement;
		typedef Model_T3D_CHM_mt::Position Position;
		CalData& cd;
		size_t task_num;

	public:
		InitPcl(CalData& _cd) : cd(_cd) {}
		inline size_t get_task_num() const noexcept { return task_num; }
		size_t work(size_t wk_id) const;
		inline tbb::task* operator() (tbb::task& parent, size_t wk_id, InitPclRes &res) const
		{
			res.pcl_num = work(wk_id);
			return nullptr;
		}

		void init(size_t thread_num) noexcept
		{
			task_num = ParallelUtils::cal_task_num<
				min_pcl_num_per_init_pcl_task, task_num_per_thread>(
					thread_num, cd.prev_valid_pcl_num);
		}
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

		CalData& cd;
		const double* pcl_m_s;
		const double* pcl_vol_s;
		const Force* pcl_bf_s;
		const Force* pcl_bf_f;
		const Force* pcl_t;
		double* pcl_vol;
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

		const size_t *pcl_in_elem;
		const size_t *cur_to_prev_pcl;
		
		size_t pcl_num, task_num;

	public:
		MapPclToBgMesh(CalData &_cd) : cd(_cd) {}
		void init() noexcept
		{
			// pcl data
			pcl_m_s = cd.pcl_m_s;
			pcl_vol_s = cd.pcl_vol_s;
			pcl_bf_s = cd.pcl_bf_s;
			pcl_bf_f = cd.pcl_bf_f;
			pcl_t = cd.pcl_t;
			pcl_vol = cd.pcl_vol;
			// bg mesh data
			elem_N_abc = cd.elem_N_abc;
			elem_vol = cd.elem_vol;
			elem_pcl_m_s = cd.elem_pcl_m_s;
			elem_pcl_m_f = cd.elem_pcl_m_f;
			elem_pcl_n = cd.elem_pcl_n;
			elem_density_f = cd.elem_density_f;
			elem_p = cd.elem_p;
			elem_node_vm_s = cd.elem_node_vm_s;
			elem_node_vm_f = cd.elem_node_vm_f;
			elem_node_force_s = cd.elem_node_force_s;
			elem_node_force_f = cd.elem_node_force_f;
		}
		void update(size_t thread_num) noexcept
		{
			const auto& spva0 = cd.spvas[cd.sorted_pcl_var_id];
			const auto& spva1 = cd.spvas[cd.sorted_pcl_var_id ^ 1];
			pcl_index0 = spva0.pcl_index;
			pcl_n0 = spva0.pcl_n;
			pcl_density_f0 = spva0.pcl_density_f;
			pcl_v_s0 = spva0.pcl_v_s;
			pcl_v_f0 = spva0.pcl_v_f;
			pcl_u_s0 = spva0.pcl_u_s;
			pcl_u_f0 = spva0.pcl_u_f;
			pcl_stress0 = spva0.pcl_stress;
			pcl_N0 = spva0.pcl_N;
			pcl_index1 = spva1.pcl_index;
			pcl_n1 = spva1.pcl_n;
			pcl_density_f1 = spva1.pcl_density_f;
			pcl_v_s1 = spva1.pcl_v_s;
			pcl_v_f1 = spva1.pcl_v_f;
			pcl_u_s1 = spva1.pcl_u_s;
			pcl_u_f1 = spva1.pcl_u_f;
			pcl_stress1 = spva1.pcl_stress;
			pcl_p1 = spva1.pcl_p;
			pcl_N1 = spva1.pcl_N;

			SortParticleMem& pcl_sort_mem = cd.pcl_sort_mem;
			pcl_in_elem = pcl_sort_mem.res_keys;
			cur_to_prev_pcl = pcl_sort_mem.res_vals;

			pcl_num = cd.valid_pcl_num;
			task_num = ParallelUtils::cal_task_num<min_pcl_num_per_init_pcl_task, task_num_per_thread>(thread_num, pcl_num);
		}
		inline size_t get_task_num() const noexcept { return task_num; }
		
		void work(size_t wk_id) const;
		void operator() (const tbb::blocked_range<size_t>& range) const { work(range.begin()); }
		tbb::task* operator() (tbb::task& parent, size_t wk_id) const { work(wk_id); return nullptr; }
	};

	//class ContactRigidRect
	//{
	//protected:
	//	typedef Model_T3D_CHM_mt::Force Force;
	//	typedef Model_T3D_CHM_mt::Position Position;
	//	typedef Model_T3D_CHM_mt::Displacement Displacement;
	//	typedef Model_T3D_CHM_mt::ShapeFunc ShapeFunc;

	//	CalData &cd;
	//	//RigidRect *prr;
	//	double K_cont;
	//	const Position *pcl_pos;
	//	const double *pcl_vol;
	//	Force* elem_node_force;

	//	const size_t* pcl_in_elem;
	//	const size_t* pcl_index;
	//	const ShapeFunc* pcl_N;
	//	const Displacement* pcl_disp;
	//	size_t valid_pcl_num;
	//	size_t task_num;

	//public:
	//	ContactRigidRect(CalData& _cd) : cd(_cd) {}
	//	inline void init(Model_T3D_CHM_mt &md) noexcept
	//	{
	//		//prr = &md.get_rigid_rect();
	//		//K_cont = md.get_Kn_cont();
	//		pcl_pos = cd.pcl_pos;
	//		pcl_vol = cd.pcl_vol;
	//		elem_node_force = cd.elem_node_force;
	//	}
	//	inline void update(size_t thread_num) noexcept
	//	{
	//		pcl_in_elem = cd.pcl_sort_mem.res_keys;
	//		const auto& spva0 = cd.spvas[cd.sorted_pcl_var_id];
	//		pcl_index = spva0.pcl_index;
	//		pcl_N = spva0.pcl_N;
	//		pcl_disp = spva0.pcl_disp;
	//		valid_pcl_num = cd.valid_pcl_num;
	//		task_num = ParallelUtils::cal_task_num<
	//			pcl_num_per_contact_rigid_rect_task, task_num_per_thread>(
	//				valid_pcl_num, thread_num);
	//	}
	//	inline size_t get_task_num() const noexcept { return task_num; }
	//	void operator() (size_t wk_id, Force3D &rr_cf) const;
	//};
	
	class UpdateAccelerationAndVelocity
	{
	protected:
		typedef Model_T3D_CHM_mt::Force Force;
		typedef Model_T3D_CHM_mt::ElemNodeVM ElemNodeVM;
		typedef Model_T3D_CHM_mt::Acceleration Acceleration;
		typedef Model_T3D_CHM_mt::Velocity Velocity;
		typedef Model_T3D_CHM_mt::NodeHasVBC NodeHasVBC;

		CalData& cd;

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
		const size_t* node_has_elem;
		const size_t* node_elem_pair;

		size_t four_elem_num, task_num;

	public:
		UpdateAccelerationAndVelocity(CalData& _cd) : cd(_cd) {}
		
		inline size_t get_task_num() const noexcept { return task_num; }
		inline void operator() (const tbb::blocked_range<size_t>& range) const
		{ work(range.begin()); }
		inline tbb::task* operator() (tbb::task& parent, size_t wk_id) const
		{ work(wk_id); return nullptr; }

		void init() noexcept
		{
			elem_pcl_m_s = cd.elem_pcl_m_s;
			elem_pcl_m_f = cd.elem_pcl_m_f;
			elem_node_force_s = cd.elem_node_force_s;
			elem_node_force_f = cd.elem_node_force_f;
			elem_node_vm_s = cd.elem_node_vm_s;
			elem_node_vm_f = cd.elem_node_vm_f;
			node_am_s = cd.node_am_s;
			node_am_f = cd.node_am_f;
			node_a_s = cd.node_a_s;
			node_a_f = cd.node_a_f;
			node_v_s = cd.node_v_s;
			node_v_f = cd.node_v_f;
			node_has_vbc_s = cd.node_has_vbc_s;
			node_has_vbc_f = cd.node_has_vbc_f;
			node_has_elem = cd.node_sort_mem.res_keys;
			node_elem_pair = cd.node_sort_mem.res_vals;
		}
		void update(size_t thread_num) noexcept
		{
			four_elem_num = cd.valid_elem_num * 4;
			task_num = ParallelUtils::cal_task_num<min_elem_num_per_update_a_and_v_task, task_num_per_thread>(thread_num, four_elem_num);
		}
		void work(size_t wk_id) const;
	};

	class CalElemDeAndMapToNode
	{
	protected:
		typedef Model_T3D_CHM_mt::ElemNodeIndex ElemNodeIndex;
		typedef Model_T3D_CHM_mt::Velocity Velocity;
		typedef Model_T3D_CHM_mt::DShapeFuncABC DShapeFuncABC;
		typedef Model_T3D_CHM_mt::StrainInc StrainInc;

		CalData& cd;

		const ElemNodeIndex* elem_node_id;
		const DShapeFuncABC *elem_N_abc;
		const double* elem_pcl_m_s, *elem_pcl_m_f;
		const double* elem_pcl_n;
		const Velocity *node_v_s, *node_v_f;
		StrainInc* elem_de;
		double* elem_m_de_vol_s, *elem_m_de_vol_f;
		const size_t *valid_elems;

		size_t elem_num, task_num;

	public:
		CalElemDeAndMapToNode(CalData &_cd) : cd(_cd) {}
		inline size_t get_task_num() const noexcept { return task_num; }
		inline void operator() (const tbb::blocked_range<size_t>& range) const
		{ work(range.begin()); }
		inline tbb::task* operator() (tbb::task& parent, size_t wk_id) const
		{ work(wk_id); return nullptr; }

		void init() noexcept
		{
			elem_node_id = cd.elem_node_id;
			elem_N_abc = cd.elem_N_abc;
			elem_pcl_m_s = cd.elem_pcl_m_s;
			elem_pcl_m_f = cd.elem_pcl_m_f;
			elem_pcl_n = cd.elem_pcl_n;
			node_v_s = cd.node_v_s;
			node_v_f = cd.node_v_f;
			elem_de = cd.elem_de;
			elem_m_de_vol_s = cd.elem_m_de_vol_s;
			elem_m_de_vol_f = cd.elem_m_de_vol_f;
			valid_elems = cd.node_sort_mem.res_elems;
		}
		void update(size_t thread_num) noexcept
		{
			elem_num = cd.valid_elem_num;
			task_num = ParallelUtils::cal_task_num<
				min_elem_num_per_cal_elem_de_task, task_num_per_thread>(thread_num, elem_num);
		}
		void work(size_t wk_id) const;
	};

	class CalNodeDe
	{
	protected:
		CalData& cd;

		const size_t* node_has_elem;
		const size_t* node_elem_pair;
		const double* elem_m_de_vol_s, *elem_m_de_vol_f;
		const double* node_am_s, *node_am_f;
		double* node_de_vol_s, *node_de_vol_f;

		size_t four_elem_num, task_num;

	public:
		CalNodeDe(CalData &_cd) : cd(_cd) {}
		inline size_t get_task_num() const noexcept { return task_num; }
		inline void operator() (const tbb::blocked_range<size_t>& range) const
		{ work(range.begin()); }
		inline tbb::task* operator() (tbb::task& parent, size_t wk_id) const
		{ work(wk_id); return nullptr; }

		void init() noexcept
		{
			node_has_elem = cd.node_sort_mem.res_keys;
			node_elem_pair = cd.node_sort_mem.res_vals;
			elem_m_de_vol_s = cd.elem_m_de_vol_s;
			elem_m_de_vol_f = cd.elem_m_de_vol_f;
			node_am_s = cd.node_am_s;
			node_am_f = cd.node_am_f;
			node_de_vol_s = cd.node_de_vol_s;
			node_de_vol_f = cd.node_de_vol_f;
		}
		void update(size_t thread_num) noexcept
		{
			four_elem_num = cd.valid_elem_num * 4;
			task_num = ParallelUtils::cal_task_num<
				min_elem_num_per_cal_node_de_task, task_num_per_thread>(
					thread_num, four_elem_num);
		}
		void work(size_t wk_id) const;
	};

	class MapBgMeshToPcl;
	struct MapBgMeshToPclRes
	{
		size_t pcl_num;
		MapBgMeshToPcl* map_bg_mesh_to_pcl;
		MapBgMeshToPclRes() : pcl_num(0), map_bg_mesh_to_pcl(nullptr) {}
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
		
		CalData& cd;

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

		const size_t* pcl_in_elem;
		const size_t* cur_to_prev_pcl;
		size_t* new_pcl_in_elem;
		size_t* new_cur_to_prev_pcl;

		size_t pcl_num, task_num;
		
	public:
		MapBgMeshToPcl(CalData& _cd) : cd(_cd) {}

		inline size_t get_task_num() const noexcept { return task_num; }
		inline tbb::task *operator() (tbb::task &parent, size_t wk_id, MapBgMeshToPclRes &res) const
		{ res.pcl_num = work(wk_id); return nullptr; }
		size_t work(size_t wk_id) const;

		void init() noexcept
		{
			elem_node_id = cd.elem_node_id;
			node_a_s = cd.node_a_s;
			node_a_f = cd.node_a_f;
			node_v_s = cd.node_v_s;
			node_v_f = cd.node_v_f;
			elem_density_f = cd.elem_density_f;
			elem_pcl_n = cd.elem_pcl_n;
			elem_p = cd.elem_p;
			elem_de = cd.elem_de;
			node_de_vol_s = cd.node_de_vol_s;
			node_de_vol_f = cd.node_de_vol_f;
			pcl_pos = cd.pcl_pos;
			pcl_mat_model = cd.pcl_mat_model;
		}
		void update(size_t thread_num) noexcept
		{
			const auto& spva0 = cd.spvas[cd.sorted_pcl_var_id];
			const auto& spva1 = cd.spvas[cd.sorted_pcl_var_id ^ 1];
			pcl_index0 = spva0.pcl_index;
			pcl_density_f0 = spva0.pcl_density_f;
			pcl_n0 = spva0.pcl_n;
			pcl_p0 = spva0.pcl_p;
			pcl_N0 = spva0.pcl_N;
			pcl_v_s0 = spva0.pcl_v_s;
			pcl_v_f0 = spva0.pcl_v_f;
			pcl_u_s0 = spva0.pcl_u_s;
			pcl_u_f0 = spva0.pcl_u_f;
			pcl_stress0 = spva0.pcl_stress;
			pcl_strain0 = spva0.pcl_strain;
			pcl_estrain0 = spva0.pcl_estrain;
			pcl_pstrain0 = spva0.pcl_pstrain;
			pcl_strain1 = spva1.pcl_strain;
			pcl_estrain1 = spva1.pcl_estrain;
			pcl_pstrain1 = spva1.pcl_pstrain;
			
			SortParticleMem& pcl_sort_mem = cd.pcl_sort_mem;
			pcl_in_elem = pcl_sort_mem.res_keys;
			cur_to_prev_pcl = pcl_sort_mem.res_vals;
			pcl_sort_mem.update_key_and_val();
			new_pcl_in_elem = pcl_sort_mem.ori_keys;
			new_cur_to_prev_pcl = pcl_sort_mem.ori_vals;

			// update task_Num
			pcl_num = cd.valid_pcl_num;
			task_num = ParallelUtils::cal_task_num<
				min_pcl_num_per_map_mesh_to_pcl_task, task_num_per_thread>(
					thread_num, pcl_num);
		}
	};
}

#endif