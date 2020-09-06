#ifndef __Test_Simulations_h__
#define __Test_Simulations_h__

void test_solve_functions();

void test_random_point_for_Model_T2D_ME_s(int argc, char** argv);

void test_Model_T2D_ME_s_display(int argc, char** argv);
void test_Model_T3D_ME_s_display(int argc, char **argv);

void test_t2d_me_s_geostatic(int argc, char** argv);
void test_t2d_me_s_geostatic_restart(int argc, char** argv);

void test_t2d_chm_s_geostatic_mcc(int argc, char** argv);

void test_t2d_me_s_1d_compression(int argc, char **argv);
void test_t2d_chm_s_1d_consolidation(int argc, char** argv);

void test_t2d_me_s_test_restart_geo_step(int argc, char** argv);
void test_t2d_me_s_test_restart_penetration_step(int argc, char** argv);

void test_t2d_me_s_test_rigid_circle(int argc, char** argv);
void test_t2d_chm_s_test_rigid_circle(int argc, char** argv);

void test_t2d_chm_s_test_restart_geo_step(int argc, char** argv);
void test_t2d_chm_s_test_restart_consolidation_step(int argc, char** argv);

void test_t3d_me_s_1d_compression(int argc, char **argv);
void test_t3d_chm_s_1d_consolidation(int argc, char** argv);

// test on a simple cube
void test_fem_t3d_me_s_1d_compression(int argc, char** argv);

// test on more complex column
void test_fem_t3d_me_s_1d_compression2(int argc, char** argv);

// parallelsim
void test_t2d_me_p_geostatic(int argc, char** argv);
void test_t2d_me_p_test(int argc, char** argv);
void test_t2d_me_p_1d_compression(int argc, char** argv);
void test_t2d_me_p_rigid_circle_penetration(int argc, char** argv);

// shallow pipe embedment
// chm simulation
void test_t2d_chm_s_pipe_conference_geo(int argc, char** argv);
void test_t2d_chm_s_pipe_conference_restart1(int argc, char** argv);
void test_t2d_chm_s_pipe_conference_restart2(int argc, char** argv);

void test_t2d_me_s_pipe_conference_no_geostress(int argc, char** argv);

// completely drained or undrained
void test_t2d_me_s_pipe_conference_geo(int argc, char** argv);
void test_t2d_me_s_pipe_conference_restart(int argc, char** argv);
void test_t2d_me_s_pipe_conference_restart2(int argc, char** argv);
void test_t2d_me_s_pipe_conference_geo_undrained(int argc, char** argv);

void test_t2d_me_p_pipe_conference_geo(int argc, char** argv);
void test_t2d_me_p_pipe_conference_restart(int argc, char** argv);

void test_RigidTetrahedronMesh_intersection(int argc, char** argv);
void test_RigidTetrahedronMesh_bg_grid(int argc, char** argv);
void test_RigidTetrahedronMesh_bg_grid2(int argc, char** argv);
void test_RigidTetrahedronMesh_close_to_boundary(int argc, char** argv);
void test_RigidTetrahedronMesh_search_dist(int argc, char** argv);

void test_t3d_me_s_cap_compression(int argc, char** argv);

#endif