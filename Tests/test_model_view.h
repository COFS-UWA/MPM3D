#ifndef __Test_Model_View_h__
#define __Test_Model_View_h__

void test_PospMPM3DApp(int argc, char **argv);
void test_PospSingleFrame_display(int argc, char **argv);

void test_t2d_me_s_geostatic_result(int argc, char** argv);
void test_t2d_chm_s_geostatic_mcc_result(int argc, char** argv);

void test_t2d_me_s_1d_compression_static_result(int argc, char** argv);
void test_t2d_me_s_1d_compression_ani_result(int argc, char** argv);

void test_t2d_chm_s_1d_compression_static_result(int argc, char** argv);
void test_t2d_chm_s_1d_compression_ani_result(int argc, char** argv);

void test_t2d_me_s_test_restart_geo_step_result(int argc, char** argv);
void test_t2d_me_s_test_restart_penetration_step_result(int argc, char** argv);

void test_t2d_chm_s_test_restart_geo_step_result(int argc, char** argv);
void test_t2d_chm_s_test_restart_consolidation_step_result(int argc, char** argv);

void test_t2d_chm_s_t_bar_conference_geo_result(int argc, char** argv);
void test_t2d_chm_s_t_bar_conference_restart1_result(int argc, char** argv);
void test_t2d_chm_s_t_bar_conference_restart2_result(int argc, char** argv);

void test_t3d_me_s_1d_compression_result(int argc, char **argv);
void test_t3d_chm_s_1d_consolidation_result(int argc, char **argv);

#endif