#include "TestsParallel_pcp.h"

#include "test_simulations_omp.h"
#include "test_model_view_omp.h"

int main(int argc, char *argv[])
{
	//test_von_mises(argc, argv);

	//test_Step_OMP(argc, argv);

	//test_t2d_me_mt_test1(argc, argv);
	
	//test_t2d_me_mt_test2(argc, argv);
	//test_t2d_me_mt_test2_result(argc, argv);
	
	//test_t2d_me_s_test2(argc, argv);
	//test_t2d_me_s_test2_result(argc, argv);

	//test_t2d_me_mt_cap_compression(argc, argv);
	//test_t2d_me_mt_cap_compression_result(argc, argv);

	//test_t2d_me_mt_strip_footing_smaller(argc, argv);
	//test_t2d_me_mt_strip_footing_smaller_result(argc, argv);

	//test_t2d_me_mt_strip_footing(argc, argv);
	//test_t2d_me_mt_strip_footing_result(argc, argv);

	//test_t3d_me_mt_1d_compression(argc, argv);
	//test_t3d_me_mt_1d_compression_result(argc, argv);

	//test_t3d_me_s_1d_compression(argc, argv);
	//test_t3d_me_s_1d_compression_result(argc, argv);

	//test_omp_barrier_time();
	//test_rigid_cylinder(argc, argv);
	//test_rigid_cone(argc, argv);

	//test_contact_model_3d(argc, argv);
	
	//test_t3d_me_mt_cap_compression(argc, argv);
	//test_t3d_me_mt_cap_compression_result(argc, argv);
	//test_t3d_me_mt_cap_compression_restart(argc, argv);

	//test_t3d_me_mt_block_sliding(argc, argv);
	//test_t3d_me_mt_block_sliding_result(argc, argv);

	test_t3d_me_mt_cylinder_foundation_create_model(argc, argv);
	test_t3d_me_mt_cylinder_foundation(argc, argv);
	
	//system("pause");
	return 0;
}
