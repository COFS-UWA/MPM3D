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
	//test_mcc_get_Su();
	//test_mcc_compression();
	//test_undrained_mcc();

	//test_solve_functions();

	//test_searching_grid2d1();

	//test_tetrahedron_mesh();
	//test_searching_grid3d1();
	//test_searching_grid3d2();
	//test_searching_grid3d3();

	//test_random_point_for_Model_T2D_ME_s(argc, argv);
	//test_random_point_for_Model_T2D_ME_s_result(argc, argv);

	//test_Model_T2D_ME_s_display(argc, argv);
	//test_Model_T3D_ME_s_display(argc, argv);

	//test_t3d_me_s_1d_compression(argc, argv);
	//test_t3d_me_s_1d_compression_result(argc, argv);

	//test_t3d_chm_s_1d_consolidation(argc, argv);
	//test_t3d_chm_s_1d_consolidation_result(argc, argv);

	//test_fem_t3d_me_s_1d_compression(argc, argv);
	//test_fem_t3d_me_s_1d_compression2(argc, argv);

	//test_PospSingleFrame_display(argc, argv);
	//test_PospMPM3DApp(argc, argv);

	//test_t2d_me_s_geostatic(argc, argv);
	//test_t2d_me_s_geostatic_result(argc, argv);
	//test_t2d_me_s_geostatic_restart(argc, argv);
	//test_t2d_me_s_geostatic_restart_result(argc, argv);

	//test_t2d_chm_s_geostatic_mcc(argc, argv);
	//test_t2d_chm_s_geostatic_mcc_result(argc, argv);

	//test_t2d_me_s_1d_compression(argc, argv);
	//test_t2d_me_s_1d_compression_static_result(argc, argv);
	//test_t2d_me_s_1d_compression_ani_result(argc, argv);

	//test_t2d_me_s_1d_compression_horizontal(argc, argv);
	//test_t2d_me_s_1d_compression_horizontal_result(argc, argv);

	//test_t2d_chm_s_1d_consolidation(argc, argv);
	//test_t2d_chm_s_1d_compression_static_result(argc, argv);
	//test_t2d_chm_s_1d_compression_ani_result(argc, argv);

	//test_t2d_chm_s_ud_compression(argc, argv);
	//test_t2d_chm_s_ud_compression_result(argc, argv);

	// test me restart
	//test_t2d_me_s_test_restart_geo_step(argc, argv);
	//test_t2d_me_s_test_restart_geo_step_result(argc, argv);
	//test_t2d_me_s_test_restart_penetration_step(argc, argv);
	//test_t2d_me_s_test_restart_penetration_step_result(argc, argv);

	// test chm restart
	//test_t2d_chm_s_test_restart_geo_step(argc, argv);
	//test_t2d_chm_s_test_restart_geo_step_result(argc, argv);
	//test_t2d_chm_s_test_restart_consolidation_step(argc, argv);
	//test_t2d_chm_s_test_restart_consolidation_step_result(argc, argv);

	// rigid cicle
	//test_t2d_me_s_test_rigid_circle(argc, argv);
	//test_t2d_me_s_test_rigid_circle_result(argc, argv);
	//test_t2d_chm_s_test_rigid_circle(argc, argv);
	//test_t2d_chm_s_test_rigid_circle_result(argc, argv);

	// parallelism
	//test_t2d_me_p_geostatic(argc, argv);
	//test_t2d_me_p_geostatic_result(argc, argv);

	//test_t2d_me_p_test(argc, argv);
	//test_t2d_me_p_test_result(argc, argv);
	
	//test_t2d_me_p_1d_compression(argc, argv);
	//test_t2d_me_p_1d_compression_result(argc, argv);

	//test_t2d_me_p_rigid_circle_penetration(argc, argv);
	//test_t2d_me_p_rigid_circle_penetration_result(argc, argv);

	// pipe embedment simulation for conference
	// coupled hydro-mechanics
	// geostatic step
	//test_t2d_chm_s_pipe_conference_geo(argc, argv);
	//test_t2d_chm_s_pipe_conference_geo_result(argc, argv);
	// penetration step 1
	//test_t2d_chm_s_pipe_conference_restart1(argc, argv);
	//test_t2d_chm_s_pipe_conference_restart1_result(argc, argv);
	// penetration step 2
	//test_t2d_chm_s_pipe_conference_restart2(argc, argv);
	//test_t2d_chm_s_pipe_conference_restart2_result(argc, argv);

	// completely drained
	//test_t2d_me_s_pipe_conference_geo(argc, argv);
	//test_t2d_me_s_pipe_conference_geo_result(argc, argv);
	
	//test_t2d_me_s_pipe_conference_geo_undrained(argc, argv);
	//test_t2d_me_s_pipe_conference_geo_undrained_result(argc, argv);

	//test_t2d_me_s_pipe_conference_restart(argc, argv);
	//test_t2d_me_s_pipe_conference_restart_result(argc, argv);
	//test_t2d_me_s_pipe_conference_restart2(argc, argv);
	//test_t2d_me_s_pipe_conference_restart_result2(argc, argv);

	//test_t2d_me_s_pipe_conference_drained_no_geostress(argc, argv);
	//test_t2d_me_s_pipe_conference_drained_no_geostress_result(argc, argv);
	
	//test_t2d_chm_s_pipe_conference_undrained_no_geostress(argc, argv);
	//test_t2d_chm_s_pipe_conference_undrained_no_geostress_result(argc, argv);

	//test_t2d_me_s_plate_with_hole(argc, argv);
	test_t2d_me_s_plate_with_hole_result(argc, argv);

	//test_t2d_me_p_pipe_conference_geo(argc, argv);
	//test_t2d_me_p_pipe_conference_geo_result(argc, argv);
	//test_t2d_me_p_pipe_conference_restart(argc, argv);
	//test_t2d_me_p_pipe_conference_restart_result(argc, argv);
	
	//test_RigidTetrahedronMesh_display(argc, argv);
	//test_RigidTetrahedronMesh_intersection(argc, argv);
	//test_RigidTetrahedronMesh_bg_grid(argc, argv);
	//test_RigidTetrahedronMesh_bg_grid2(argc, argv);
	//test_RigidTetrahedronMesh_close_to_boundary(argc, argv);
	//test_RigidTetrahedronMesh_search_dist(argc, argv);

	//test_t3d_me_s_cap_compression(argc, argv);
	//test_t3d_me_s_cap_compression_result(argc, argv);

	//system("pause");
	return 0;
}