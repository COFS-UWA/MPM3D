#include "TestsParallel_pcp.h"

#include "test_simulations_omp.h"
#include "test_model_view_omp.h"

int main(int argc, char *argv[])
{
	//test_Step_OMP(argc, argv);

	//test_t2d_me_mt_test1(argc, argv);
	
	//test_t2d_me_mt_test2(argc, argv);
	//test_t2d_me_mt_test2_result(argc, argv);
	
	//test_t2d_me_s_test2(argc, argv);
	//test_t2d_me_s_test2_result(argc, argv);

	test_t2d_me_mt_strip_footing_smaller(argc, argv);
	//test_t2d_me_mt_strip_footing_smaller_result(argc, argv);

	//system("pause");
	return 0;
}