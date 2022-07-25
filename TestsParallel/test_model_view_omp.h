#ifndef __Model_View_OMP_h__
#define __Model_View_OMP_h__

void test_t2d_me_mt_test1_result(int argc, char** argv);
void test_t2d_me_s_test2_result(int argc, char** argv);
void test_t2d_me_mt_1d_compression_result(int argc, char** argv);
void test_t2d_me_mt_cap_compression_result(int argc, char** argv);
void test_t2d_me_mt_block_collision_result(int argc, char** argv);
void test_t2d_me_mt_block_sliding_result(int argc, char** argv);
void test_t2d_me_mt_strip_footing_smaller_result(int argc, char** argv);
void test_t2d_me_mt_strip_footing_result(int argc, char** argv);
void test_t2d_me_mt_geostatic_result(int argc, char** argv);

void test_t3d_me_mt_1d_compression_result(int argc, char** argv);
void test_t3d_me_s_1d_compression_result(int argc, char** argv);
void test_t3d_me_mt_1d_geostatic_result(int argc, char** argv);

void test_t3d_me_mt_cap_compression_result(int argc, char** argv);
void test_t3d_me_mt_cap_compression_result_div(int argc, char** argv);
void test_t3d_me_mt_triaxial_compression_result(int argc, char** argv);

void test_t3d_me_mt_block_sliding_result(int argc, char** argv);

void test_t3d_me_mt_cylinder_foundation_result(int argc, char** argv);
void test_t3d_me_mt_cylinder_foundation_result_den(int argc, char** argv);

void test_t3d_me_mt_test_rigid_mesh_result(int argc, char** argv);

void test_t2d_chm_mt_1d_consolidation_static_result(int argc, char** argv);
void test_t2d_chm_mt_1d_consolidation_ani_result(int argc, char** argv);

void test_t2d_chm_mt_geostatic_result(int argc, char** argv);

void test_t2d_chm_mt_block_sliding_result(int argc, char** argv);

void test_t2d_chm_mt_test_rigid_circle_result(int argc, char** argv);

void test_t2d_chm_mt_pipe_conference_result(int argc, char** argv);
void test_t2d_chm_mt_pipe_conference_restart_result(int argc, char** argv);

void test_t2d_chm_mt_pipe_conference_den_result(int argc, char** argv);
void test_t2d_chm_mt_pipe_embedment_result(int argc, char** argv);

void test_t3d_chm_mt_triaxial_compression_result(int argc, char** argv);

void test_t3d_me_mt_cylinder_foundation_result_ch_den2(int argc, char** argv);

void test_t3d_me_mt_spudcan_coarse_result(int argc, char** argv);

void test_t3d_chm_mt_1d_consolidation_result(int argc, char** argv);
void test_t3d_chm_mt_1d_geostatic_result(int argc, char** argv);

void test_t2d_me_tbb_1d_compression_result(int argc, char** argv);
void test_t2d_me_tbb_cap_compression_result(int argc, char** argv);
void test_t2d_chm_tbb_1d_consolidation_result(int argc, char** argv);

void test_t3d_me_tbb_1d_compression_result(int argc, char** argv);
void test_t3d_me_tbb_cap_compression_result(int argc, char** argv);

void test_t3d_chm_tbb_1d_consolidation_result(int argc, char** argv);
void test_t3d_chm_tbb_cap_compression_result(int argc, char** argv);

void test_t3d_me_tbb_piezofoundation_result(int argc, char** argv);
void test_t3d_chm_tbb_piezofoundation_result(int argc, char** argv);

void test_t3d_me_mt_spudcan_geo_result(int argc, char** argv);
void test_t3d_me_mt_spudcan_result(int argc, char** argv);

void test_t3d_chm_mt_spudcan_geo_result(int argc, char** argv);
void test_t3d_chm_mt_spudcan_result(int argc, char** argv);

void test_t3d_me_mt_spudcan_cy_geo_result(int argc, char** argv);
void test_t3d_me_mt_spudcan_cy_result(int argc, char** argv);

void test_t3d_chm_mt_spudcan_cy_geo_result(int argc, char** argv);
void test_t3d_chm_mt_spudcan_cy_result(int argc, char** argv);

void test_t3d_me_mt_spudcan_cy_Hossain_2006_result(int argc, char** argv);

#endif