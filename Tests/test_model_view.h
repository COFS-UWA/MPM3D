#ifndef __Test_Model_View_h__
#define __Test_Model_View_h__

void test_PospMPM3DApp(int argc, char **argv);
void test_PospSingleFrame_display(int argc, char **argv);

void test_random_point_for_Model_T2D_ME_s_result(int argc, char** argv);

void test_t2d_me_s_geostatic_result(int argc, char** argv);
void test_t2d_me_s_geostatic_restart_result(int argc, char** argv);

void test_t2d_chm_s_geostatic_mcc_result(int argc, char** argv);

void test_t2d_me_s_1d_compression_static_result(int argc, char** argv);
void test_t2d_me_s_1d_compression_ani_result(int argc, char** argv);
void test_t2d_me_s_1d_compression_horizontal_result(int argc, char** argv);
void test_t2d_me_s_plate_with_hole_result(int argc, char** argv);
void test_t2d_me_s_cap_compression_result(int argc, char** argv);
void test_t2d_me_s_t_bar_smaller_soil_result(int argc, char** argv);
void test_t2d_me_s_shallow_foundation_result(int argc, char** argv);

void test_t2d_chm_s_1d_consolidation_static_result(int argc, char** argv);
void test_t2d_chm_s_1d_consolidation_ani_result(int argc, char** argv);
void test_t2d_chm_s_ud_compression_result(int argc, char** argv);

void test_t2d_me_s_test_restart_geo_step_result(int argc, char** argv);
void test_t2d_me_s_test_restart_penetration_step_result(int argc, char** argv);

void test_t2d_chm_s_test_restart_geo_step_result(int argc, char** argv);
void test_t2d_chm_s_test_restart_consolidation_step_result(int argc, char** argv);

void test_t2d_me_s_test_rigid_circle_result(int argc, char** argv);
void test_t2d_chm_s_test_rigid_circle_result(int argc, char** argv);

void test_t2d_chm_s_t_bar_smaller_soil_result(int argc, char** argv);
void test_t2d_chm_s_ud_t_bar_smaller_soil_result(int argc, char** argv);

void test_t3d_me_s_1d_compression_result(int argc, char **argv);
void test_t3d_chm_s_1d_consolidation_result(int argc, char **argv);

// parallelism
void test_t2d_me_p_geostatic_result(int argc, char** argv);
void test_t2d_me_p_test_result(int argc, char** argv);
void test_t2d_me_p_1d_compression_result(int argc, char** argv);
void test_t2d_me_p_rigid_circle_penetration_result(int argc, char** argv);

// shallow pipe embedment
// chm
void test_t2d_chm_s_pipe_conference_geo_result(int argc, char** argv);
void test_t2d_chm_s_pipe_conference_restart1_result(int argc, char** argv);
void test_t2d_chm_s_pipe_conference_restart2_result(int argc, char** argv);

// completely drained or undrained
void test_t2d_me_s_pipe_conference_geo_result(int argc, char** argv);
void test_t2d_me_s_pipe_conference_restart_result(int argc, char** argv);
void test_t2d_me_s_pipe_conference_restart_result2(int argc, char** argv);
void test_t2d_me_s_pipe_conference_geo_undrained_result(int argc, char** argv);

void test_t2d_me_s_pipe_conference_drained_no_geostress_result(int argc, char** argv);
void test_t2d_chm_s_pipe_conference_undrained_no_geostress_result(int argc, char** argv);

void test_t2d_me_p_pipe_conference_geo_result(int argc, char** argv);
void test_t2d_me_p_pipe_conference_restart_result(int argc, char** argv);

void test_t3d_me_s_cap_compression_result(int argc, char** argv);

#endif