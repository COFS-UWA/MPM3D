#include "Tests_pcp.h"

#include "test_utilities.h"
#include "test_geometry.h"
#include "test_material_models.h"
#include "test_simulations.h"
#include "test_model_view.h"

int main(int argc, char *argv[])
{
	//test_stack_and_link_list();

	//test_model_container();

	//test_solve_functions();

	//test_tetrahedron_mesh();
	//test_searching_grid3d1();
	//test_searching_grid3d2();
	//test_searching_grid3d3();

	//test_Model_T3D_ME_s_display(argc, argv);

	test_t3d_me_s_1d_compression(argc, argv);
	
	//test_PospSingleFrame_display(argc, argv);
	//test_PospMPM3DApp(argc, argv);
	
	//system("pause");
	return 0;
}