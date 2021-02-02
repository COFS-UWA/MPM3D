#include "SimulationsOMP_pcp.h"

#include <fstream>

#include "Step_T2D_ME_TBB_task.h"
#include "Step_T2D_ME_TBB.h"

#define one_third (1.0/3.0)
#define N_min (1.0e-8)
#define Block_Low(th_id, th_num, data_num) ((th_id)*(data_num)/(th_num))

#ifdef _DEBUG
static std::fstream res_file_t2d_me_tbb;
#endif

Step_T2D_ME_TBB::Step_T2D_ME_TBB(const char* _name) :
	Step_TBB(_name, "Step_T2D_ME_TBB", &cal_substep_func_T2D_ME_TBB) {}

Step_T2D_ME_TBB::~Step_T2D_ME_TBB() {}

int Step_T2D_ME_TBB::init_calculation()
{
#ifdef _DEBUG
	res_file_t2d_me_tbb.open("step_t2d_me_tbb.txt", std::ios::out | std::ios::binary);
#endif


	return 0;
}

int Step_T2D_ME_TBB::finalize_calculation()
{
	return 0;
}

int cal_substep_func_T2D_ME_TBB(void* _self)
{
	//Step_T2D_ME_TBB& self = *(Step_T2D_ME_TBB *)(_self);
	//Model_T2D_ME_mt& md = *(Model_T2D_ME_mt*)(self.model);

	//const double* const pcl_m = self.pcl_m;
	//const Force* const pcl_bf = self.pcl_bf;
	//const Force* const pcl_t = self.pcl_t;
	//const Position* const pcl_pos = self.pcl_pos;
	//double* const pcl_vol = self.pcl_vol;
	//MatModel::MaterialModel** const pcl_mat_model = self.pcl_mat_model;

	//const size_t thread_num = self.thread_num;
	//ThreadData& thd = self.thread_datas[my_th_id];
	//SortedPclVarArrays& spva0 = self.sorted_pcl_var_arrays[thd.sorted_pcl_var_id];
	//thd.sorted_pcl_var_id ^= 1;
	//SortedPclVarArrays& spva1 = self.sorted_pcl_var_arrays[thd.sorted_pcl_var_id];

	//size_t* const pcl_index0 = spva0.pcl_index;
	//double* const pcl_density0 = spva0.pcl_density;
	//Displacement* const pcl_disp0 = spva0.pcl_disp;
	//Velocity* const pcl_v0 = spva0.pcl_v;
	//Stress* const pcl_stress0 = spva0.pcl_stress;
	//Strain* const pcl_strain0 = spva0.pcl_strain;
	//Strain* const pcl_estrain0 = spva0.pcl_estrain;
	//Strain* const pcl_pstrain0 = spva0.pcl_pstrain;
	//ShapeFunc* const pcl_N0 = spva0.pcl_N;

	//size_t* const pcl_index1 = spva1.pcl_index;
	//double* const pcl_density1 = spva1.pcl_density;
	//Displacement* const pcl_disp1 = spva1.pcl_disp;
	//Velocity* const pcl_v1 = spva1.pcl_v;
	//Stress* const pcl_stress1 = spva1.pcl_stress;
	//Strain* const pcl_strain1 = spva1.pcl_strain;
	//Strain* const pcl_estrain1 = spva1.pcl_estrain;
	//Strain* const pcl_pstrain1 = spva1.pcl_pstrain;
	//ShapeFunc* const pcl_N1 = spva1.pcl_N;

	//ElemNodeIndex* const elem_node_id = self.elem_node_id;
	//DShapeFuncABC* const elem_dN_abc = self.elem_dN_abc;
	//DShapeFuncD* const elem_dN_d = self.elem_dN_d;
	//double* const elem_vol = self.elem_vol;

	//double* const elem_density = self.elem_density;
	//double* const elem_pcl_m = self.elem_pcl_m;
	//StrainInc* const elem_de = self.elem_de;
	//double* const elem_m_de_vol = self.elem_m_de_vol;

	//ElemNodeVM* const elem_node_vm = self.elem_node_vm;
	//ElemNodeForce* const elem_node_force = self.elem_node_force;

	//Acceleration* const node_a = self.node_a;
	//Velocity* const node_v = self.node_v;
	//NodeHasVBC* const node_has_vbc = self.node_has_vbc;
	//double* const node_am = self.node_am;
	//double* const node_de_vol = self.node_de_vol;

	//union
	//{
	//	struct
	//	{
	//		size_t* prev_pcl_id0;
	//		size_t* prev_pcl_id1;
	//		size_t* pcl_in_elem0;
	//		size_t* pcl_in_elem1;
	//		size_t* node_has_elem0;
	//		size_t* node_has_elem1;
	//		size_t* node_elem_pair0;
	//		size_t* node_elem_pair1;
	//	};
	//	struct
	//	{
	//		size_t prev_pcl_id_ui0;
	//		size_t prev_pcl_id_ui1;
	//		size_t pcl_in_elem_ui0;
	//		size_t pcl_in_elem_ui1;
	//		size_t node_has_elem_ui0;
	//		size_t node_has_elem_ui1;
	//		size_t node_elem_pair_ui0;
	//		size_t node_elem_pair_ui1;
	//	};
	//};

	//prev_pcl_id0 = self.prev_pcl_ids[thd.sorted_pcl_in_elem_id];
	//prev_pcl_id1 = self.prev_pcl_ids[thd.sorted_pcl_in_elem_id ^ 1];
	//pcl_in_elem0 = self.pcl_in_elems[thd.sorted_pcl_in_elem_id];
	//pcl_in_elem1 = self.pcl_in_elems[thd.sorted_pcl_in_elem_id ^ 1];
	//node_has_elem0 = self.node_has_elems[0];
	//node_has_elem1 = self.node_has_elems[1];
	//node_elem_pair0 = self.node_elem_pairs[0];
	//node_elem_pair1 = self.node_elem_pairs[1];
	//
	//using Step_T2D_ME_TBB_task::MapPclToBgMeshTask;
	return 0;
}
