#ifndef __Step_T3D_ME_Task_h__
#define __Step_T3D_ME_Task_h__

#include "tbb/task.h"
#include "tbb/task_scheduler_init.h"

#include "ParallelUtils.h"
#include "MSDRadixSortUtils.h"
#include "SortParticleTask.h"
#include "SortTehMeshNodeTask.h"
#include "Model_T3D_ME_mt.h"

namespace Step_T3D_ME_Task
{
	using MSDRadixSortUtils::block_low;

	constexpr size_t task_num_per_thread = 3;
	constexpr size_t pcl_num_per_init_pcl_task = 100;
	constexpr size_t pcl_num_per_map_pcl_to_mesh_task = 100;
	constexpr size_t pcl_num_per_contact_rigid_rect_task = 100;
	constexpr size_t elem_num_per_update_a_and_v_task = 20;
	constexpr size_t elem_num_per_cal_elem_de_task = 20;
	constexpr size_t elem_num_per_cal_node_de_task = 20;
	constexpr size_t pcl_num_per_map_mesh_to_pcl_task = 100;

	struct CalData
	{
		Model_T3D_ME_mt *pmodel;

		// pcl data
		const double* pcl_m;
		const Model_T3D_ME_mt::Force* pcl_bf;
		const Model_T3D_ME_mt::Force* pcl_t;
		const Model_T3D_ME_mt::Position* pcl_pos;
		double* pcl_vol;
		MatModel::MaterialModel** pcl_mat_model;

		Model_T3D_ME_mt::SortedPclVarArrays spvas[2];

		// elem data
		const Model_T3D_ME_mt::ElemNodeIndex* elem_node_id;
		const Model_T3D_ME_mt::DShapeFuncABC* elem_dN_abc;
		const Model_T3D_ME_mt::DShapeFuncD* elem_dN_d;
		const double* elem_vol;

		double* elem_pcl_m;
		double* elem_density;
		Model_T3D_ME_mt::StrainInc* elem_de;
		double* elem_m_de_vol;
		
		// elem node data
		MSDRadixSortUtils::RadixBinBlockMemArray thread_bin_blocks_mem;
		Model_T3D_ME_mt::ElemNodeVM* elem_node_vm;
		Model_T3D_ME_mt::Force* elem_node_force;

		// node data
		Model_T3D_ME_mt::Acceleration* node_a;
		Model_T3D_ME_mt::Velocity* node_v;
		Model_T3D_ME_mt::NodeHasVBC* node_has_vbc;
		double* node_am;
		double* node_de_vol;
		
#ifdef _DEBUG
		size_t elem_num;
		size_t node_num;
		size_t ori_pcl_num;
#endif
		
		// cal data
		size_t thread_num;
		SortParticleMem pcl_sort_mem;
		SortTehMeshNodeMem node_sort_mem;
		
		// data changed during computation
		double dt;
		size_t sorted_pcl_var_id;
		size_t prev_valid_pcl_num;
		size_t valid_pcl_num;
		size_t valid_elem_num;

		void set_model(Model_T3D_ME_mt& md) noexcept;
	};
	
	class InitPcl
	{
	protected:
		typedef Model_T3D_ME_mt::ShapeFunc ShapeFunc;
		typedef Model_T3D_ME_mt::Displacement Displacement;
		typedef Model_T3D_ME_mt::Position Position;
		CalData& cd;
		size_t task_num;
	public:
		InitPcl(CalData& _cd) : cd(_cd) {}
		inline void init(size_t thread_num) noexcept
		{
			task_num = ParallelUtils::cal_task_num<
				pcl_num_per_init_pcl_task, task_num_per_thread>(
					cd.prev_valid_pcl_num, thread_num);
		}
		inline size_t get_task_num() const noexcept { return task_num; }
		void operator() (size_t wk_id, size_t &pcl_in_mesh_num) const;
	};
	
	class MapPclToBgMesh
	{
	protected:
		typedef Model_T3D_ME_mt::Force Force;
		typedef Model_T3D_ME_mt::Acceleration Acceleration;
		typedef Model_T3D_ME_mt::Velocity Velocity;
		typedef Model_T3D_ME_mt::Displacement Displacement;
		typedef Model_T3D_ME_mt::Stress Stress;
		typedef Model_T3D_ME_mt::ShapeFunc ShapeFunc;
		typedef Model_T3D_ME_mt::DShapeFuncABC DShapeFuncABC;
		typedef Model_T3D_ME_mt::DShapeFuncD DShapeFuncD;
		typedef Model_T3D_ME_mt::ElemNodeVM ElemNodeVM;

		CalData& cd;
		const double* pcl_m;
		const Force* pcl_bf;
		const Force* pcl_t;
		double* pcl_vol;
		const DShapeFuncABC *elem_dN_abc;
		const double *elem_vol;
		double* elem_pcl_m;
		double* elem_density;
		ElemNodeVM* elem_node_vm;
		Force* elem_node_force;

		size_t *pcl_index0;
		double *pcl_density0;
		Velocity *pcl_v0;
		Displacement *pcl_disp0;
		Stress *pcl_stress0;
		ShapeFunc *pcl_N0;
		const size_t *pcl_index1;
		const double *pcl_density1;
		const Velocity *pcl_v1;
		const Displacement *pcl_disp1;
		const Stress *pcl_stress1;
		const ShapeFunc *pcl_N1;

		size_t valid_pcl_num;
		size_t task_num;
		const size_t* pcl_in_elem;
		const size_t* cur_to_prev_pcl;
		
	public:
		MapPclToBgMesh(CalData &_cd) : cd(_cd) {}
		inline void init() noexcept
		{
			pcl_m = cd.pcl_m;
			pcl_bf = cd.pcl_bf;
			pcl_t = cd.pcl_t;
			pcl_vol = cd.pcl_vol;
			elem_dN_abc = cd.elem_dN_abc;
			elem_vol = cd.elem_vol;
			elem_pcl_m = cd.elem_pcl_m;
			elem_density = cd.elem_density;
			elem_node_vm = cd.elem_node_vm;
			elem_node_force = cd.elem_node_force;
		}
		inline void update(size_t thread_num) noexcept
		{
			const auto& spva0 = cd.spvas[cd.sorted_pcl_var_id];
			const auto& spva1 = cd.spvas[cd.sorted_pcl_var_id ^ 1];
			pcl_index0 = spva0.pcl_index;
			pcl_density0 = spva0.pcl_density;
			pcl_v0 = spva0.pcl_v;
			pcl_disp0 = spva0.pcl_disp;
			pcl_stress0 = spva0.pcl_stress;
			pcl_N0 = spva0.pcl_N;
			pcl_index1 = spva1.pcl_index;
			pcl_density1 = spva1.pcl_density;
			pcl_v1 = spva1.pcl_v;
			pcl_disp1 = spva1.pcl_disp;
			pcl_stress1 = spva1.pcl_stress;
			pcl_N1 = spva1.pcl_N;

			valid_pcl_num = cd.valid_pcl_num;
			task_num = ParallelUtils::cal_task_num<
				pcl_num_per_map_pcl_to_mesh_task, task_num_per_thread>(
					valid_pcl_num, thread_num);

			SortParticleMem& pcl_sort_mem = cd.pcl_sort_mem;
			pcl_in_elem = pcl_sort_mem.res_keys;
			cur_to_prev_pcl = pcl_sort_mem.res_vals;
		}
		inline size_t get_task_num() const noexcept { return task_num; }
		void operator() (size_t wk_id) const;
	};

	class ContactRigidRect
	{
	protected:
		typedef Model_T3D_ME_mt::Force Force;
		typedef Model_T3D_ME_mt::Position Position;
		typedef Model_T3D_ME_mt::Displacement Displacement;
		typedef Model_T3D_ME_mt::ShapeFunc ShapeFunc;

		CalData &cd;
		//RigidRect *prr;
		double K_cont;
		const Position *pcl_pos;
		const double *pcl_vol;
		Force* elem_node_force;

		const size_t* pcl_in_elem;
		const size_t* pcl_index;
		const ShapeFunc* pcl_N;
		const Displacement* pcl_disp;
		size_t valid_pcl_num;
		size_t task_num;

	public:
		ContactRigidRect(CalData& _cd) : cd(_cd) {}
		inline void init(Model_T3D_ME_mt &md) noexcept
		{
			//prr = &md.get_rigid_rect();
			//K_cont = md.get_Kn_cont();
			pcl_pos = cd.pcl_pos;
			pcl_vol = cd.pcl_vol;
			elem_node_force = cd.elem_node_force;
		}
		inline void update(size_t thread_num) noexcept
		{
			pcl_in_elem = cd.pcl_sort_mem.res_keys;
			const auto& spva0 = cd.spvas[cd.sorted_pcl_var_id];
			pcl_index = spva0.pcl_index;
			pcl_N = spva0.pcl_N;
			pcl_disp = spva0.pcl_disp;
			valid_pcl_num = cd.valid_pcl_num;
			task_num = ParallelUtils::cal_task_num<
				pcl_num_per_contact_rigid_rect_task, task_num_per_thread>(
					valid_pcl_num, thread_num);
		}
		inline size_t get_task_num() const noexcept { return task_num; }
		void operator() (size_t wk_id, Force3D &rr_cf) const;
	};
	
	class UpdateAccelerationAndVelocity
	{
	protected:
		typedef Model_T3D_ME_mt::Force Force;
		typedef Model_T3D_ME_mt::ElemNodeVM ElemNodeVM;
		typedef Model_T3D_ME_mt::Acceleration Acceleration;
		typedef Model_T3D_ME_mt::Velocity Velocity;
		typedef Model_T3D_ME_mt::NodeHasVBC NodeHasVBC;

		CalData& cd;
		const double* elem_pcl_m;
		const Force* elem_node_force;
		const ElemNodeVM *elem_node_vm;
		Acceleration *node_a;
		double *node_am;
		Velocity *node_v;
		NodeHasVBC *node_has_vbc;
		const size_t* node_has_elem;
		const size_t* node_elem_pair;
		
		size_t four_valid_elem_num;
		size_t task_num;

	public:
		UpdateAccelerationAndVelocity(CalData& _cd) : cd(_cd) {}
		inline void init() noexcept
		{
			elem_pcl_m = cd.elem_pcl_m;
			elem_node_force = cd.elem_node_force;
			elem_node_vm = cd.elem_node_vm;
			node_a = cd.node_a;
			node_am = cd.node_am;
			node_v = cd.node_v;
			node_has_vbc = cd.node_has_vbc;
			node_has_elem = cd.node_sort_mem.res_keys;
			node_elem_pair = cd.node_sort_mem.res_vals;
		}
		inline void update(size_t thread_num) noexcept
		{
			four_valid_elem_num = cd.valid_elem_num * 4;
			task_num = ParallelUtils::cal_task_num<
				elem_num_per_update_a_and_v_task, task_num_per_thread>(
					four_valid_elem_num, thread_num);
		}
		inline size_t get_task_num() const noexcept { return task_num; }
		void operator() (size_t wk_id) const;
	};

	class CalElemDeAndMapToNode
	{
	protected:
		typedef Model_T3D_ME_mt::ElemNodeIndex ElemNodeIndex;
		typedef Model_T3D_ME_mt::Velocity Velocity;
		typedef Model_T3D_ME_mt::DShapeFuncABC DShapeFuncABC;
		typedef Model_T3D_ME_mt::StrainInc StrainInc;

		CalData& cd;
		const ElemNodeIndex* elem_node_id;
		const double* elem_pcl_m;
		const DShapeFuncABC *elem_dN_abc; 
		const Velocity* node_v;
		StrainInc* elem_de;
		double* elem_m_de_vol;
		const size_t *valid_elems;

		size_t valid_elem_num;
		size_t task_num;

	public:
		CalElemDeAndMapToNode(CalData &_cd) : cd(_cd) {}
		inline void init() noexcept
		{
			elem_node_id = cd.elem_node_id;
			elem_pcl_m = cd.elem_pcl_m;
			elem_dN_abc = cd.elem_dN_abc;
			node_v = cd.node_v;
			elem_de = cd.elem_de;
			elem_m_de_vol = cd.elem_m_de_vol;
			valid_elems = cd.node_sort_mem.res_elems;
		}
		inline void update(size_t thread_num) noexcept
		{
			valid_elem_num = cd.valid_elem_num;
			task_num = ParallelUtils::cal_task_num<
				elem_num_per_cal_elem_de_task, task_num_per_thread>(
					valid_elem_num, thread_num);
		}
		inline size_t get_task_num() const noexcept { return task_num; }
		void operator() (size_t wk_id) const;
	};

	class CalNodeDe
	{
	protected:
		CalData& cd;
		const size_t* node_has_elem;
		const size_t* node_elem_pair;
		const double* elem_m_de_vol;
		const double* node_am;
		double* node_de_vol;

		size_t four_valid_elem_num;
		size_t task_num;

	public:
		CalNodeDe(CalData &_cd) : cd(_cd) {}
		inline void init() noexcept
		{
			node_has_elem = cd.node_sort_mem.res_keys;
			node_elem_pair = cd.node_sort_mem.res_vals;
			elem_m_de_vol = cd.elem_m_de_vol;
			node_am = cd.node_am;
			node_de_vol = cd.node_de_vol;
		}
		inline void update(size_t thread_num) noexcept
		{
			four_valid_elem_num = cd.valid_elem_num * 4;
			task_num = ParallelUtils::cal_task_num<
				elem_num_per_cal_node_de_task, task_num_per_thread>(
					four_valid_elem_num, thread_num);
		}
		inline size_t get_task_num() const noexcept { return task_num; }
		void operator() (size_t wk_id) const;
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
		
		CalData& cd;
		const ElemNodeIndex *elem_node_id;
		const Acceleration *node_a;
		const Velocity *node_v;
		double *elem_density;
		StrainInc *elem_de;
		const double *node_de_vol;
		const Position *pcl_pos;
		MatModel::MaterialModel **pcl_mat_model;

		const size_t *pcl_index0;
		double* pcl_density0;
		Velocity* pcl_v0;
		Displacement* pcl_disp0 ;
		ShapeFunc*pcl_N0;
		Stress* pcl_stress0;
		Strain* pcl_strain0;
		Strain* pcl_estrain0;
		Strain* pcl_pstrain0;
		const Strain* pcl_strain1;
		const Strain* pcl_estrain1;
		const Strain* pcl_pstrain1;

		size_t valid_pcl_num;
		size_t task_num;
		const size_t* pcl_in_elem;
		const size_t* cur_to_prev_pcl;
		size_t* new_pcl_in_elem;
		size_t* new_cur_to_prev_pcl;

	public:
		MapBgMeshToPcl(CalData& _cd) : cd(_cd) {}
		inline void init() noexcept
		{
			elem_node_id = cd.elem_node_id;
			node_a = cd.node_a;
			node_v = cd.node_v;
			elem_density = cd.elem_density;
			elem_de = cd.elem_de;
			node_de_vol = cd.node_de_vol;
			pcl_pos = cd.pcl_pos;
			pcl_mat_model = cd.pcl_mat_model;
		}
		inline void update(size_t thread_num) noexcept
		{
			const auto& spva0 = cd.spvas[cd.sorted_pcl_var_id];
			const auto& spva1 = cd.spvas[cd.sorted_pcl_var_id ^ 1];
			pcl_index0 = spva0.pcl_index;
			pcl_density0 = spva0.pcl_density;
			pcl_v0 = spva0.pcl_v;
			pcl_disp0 = spva0.pcl_disp;
			pcl_N0 = spva0.pcl_N;
			pcl_stress0 = spva0.pcl_stress;
			pcl_strain0 = spva0.pcl_strain;
			pcl_estrain0 = spva0.pcl_estrain;
			pcl_pstrain0 = spva0.pcl_pstrain;
			pcl_strain1 = spva1.pcl_strain;
			pcl_estrain1 = spva1.pcl_estrain;
			pcl_pstrain1 = spva1.pcl_pstrain;
			
			valid_pcl_num = cd.valid_pcl_num;
			task_num = ParallelUtils::cal_task_num<
				pcl_num_per_map_mesh_to_pcl_task, task_num_per_thread>(
					valid_pcl_num, thread_num);

			SortParticleMem& pcl_sort_mem = cd.pcl_sort_mem;
			pcl_in_elem = pcl_sort_mem.res_keys;
			cur_to_prev_pcl = pcl_sort_mem.res_vals;
			pcl_sort_mem.update_key_and_val();
			new_pcl_in_elem = pcl_sort_mem.ori_keys;
			new_cur_to_prev_pcl = pcl_sort_mem.ori_vals;
		}
		inline size_t get_task_num() const noexcept { return task_num; }
		void operator() (size_t wk_id, size_t &pcl_in_mesh_num) const;
	};
}

#endif