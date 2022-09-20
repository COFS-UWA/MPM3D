#ifndef __Test_Simulations_OMP_h__
#define __Test_Simulations_OMP_h__

void test_von_mises(int argc, char** argv);

void test_rigid_cylinder(int argc, char** argv);
void test_rigid_cone(int argc, char** argv);
void test_rigid_mesh_contact(int argc, char** argv);
void test_rigid_mesh_contact2(int argc, char** argv);
void test_rigid_mesh_contact_2d(int argc, char** argv);
void test_contact_model_3d(int argc, char** argv);

void test_mohr_coulomb();
void test_sand_hypoplasticity_integration();
void test_sand_hypoplasticity_Herleand_Gudehus_1999();
void test_sand_hypoplasticity_Wichtmann_2019();
void test_sand_hypoplasticity_triaxial();
void test_triaxial_secant();
void test_sand_hypoplasticity_wrapper();

void test_t3d_me_mt_triaxial_compression(int argc, char** argv);

void test_Step_OMP(int argc, char** argv);

void test_t2d_me_mt_test1(int argc, char** argv);
void test_t2d_me_s_test2(int argc, char** argv);
void test_t2d_me_mt_1d_compression(int argc, char** argv);
void test_t2d_me_mt_cap_compression(int argc, char** argv);
void test_t2d_me_mt_strip_footing_smaller(int argc, char** argv);
void test_t2d_me_mt_strip_footing(int argc, char** argv);
void test_t2d_me_mt_block_sliding(int argc, char** argv);
void test_t3d_me_mt_weird_block_sliding(int argc, char** argv);
void test_t2d_me_mt_block_collision(int argc, char** argv);
void test_t2d_me_mt_geostatic(int argc, char** argv);
void test_t2d_me_mt_column_collapse(int argc, char** argv);

void test_t3d_me_mt_1d_compression(int argc, char** argv);
void test_t3d_me_s_1d_compression(int argc, char** argv);
void test_t3d_me_mt_1d_geostatic(int argc, char** argv);
void test_t3d_me_mt_cylinder_bcs(int argc, char** argv);
void test_t3d_me_mt_cylinder_bcs_result(int argc, char** argv);
void test_t3d_chm_mt_cylinder_bcs(int argc, char** argv);
void test_t3d_chm_mt_cylinder_bcs2(int argc, char** argv);
void test_t3d_chm_mt_cylinder_bcs_result(int argc, char** argv);

void test_t3d_me_mt_cap_compression(int argc, char** argv);
void test_t3d_me_mt_cap_compression_restart(int argc, char** argv);
void test_t3d_me_mt_cap_compression_restart2(int argc, char **argv);

void test_omp_barrier_time();
void test_t3d_me_mt_block_sliding(int argc, char** argv);

void test_t3d_me_mt_cylinder_foundation_create_model(int argc, char** argv);
void test_t3d_me_mt_cylinder_foundation(int argc, char** argv);
void test_t3d_me_mt_cylinder_foundation_restart(int argc, char** argv);
void test_t3d_me_mt_cylinder_foundation_restart2(int argc, char** argv);

void test_t3d_me_mt_cylinder_foundation_create_model_den(int argc, char** argv);
void test_t3d_me_mt_cylinder_foundation_den(int argc, char** argv);
void test_t3d_me_mt_cylinder_foundation_restart_den(int argc, char** argv);
//void test_t3d_me_mt_cylinder_foundation_restart2_den(int argc, char** argv);

void test_t3d_rigid_mesh(int argc, char** argv);
void test_t3d_me_mt_test_rigid_mesh(int argc, char** argv);

void test_t2d_chm_mt_1d_consolidation(int argc, char** argv);
void test_t2d_chm_mt_1d_consolidation_restart(int argc, char** argv);

void test_t2d_chm_mt_geostatic(int argc, char** argv);

void test_t2d_chm_mt_block_sliding(int argc, char** argv);

void test_t2d_chm_mt_test_rigid_circle(int argc, char** argv);

void test_t2d_chm_mt_pipe_conference(int argc, char** argv);
void test_t2d_chm_mt_pipe_conference_restart(int argc, char** argv);
void test_t2d_chm_mt_pipe_conference_restart1(int argc, char** argv);

void test_t2d_chm_mt_pipe_conference_den(int argc, char** argv);
void test_t2d_chm_mt_pipe_conference_den_restart(int argc, char** argv);

void test_t2d_chm_mt_pipe_embedment(int argc, char** argv);
void test_t2d_chm_mt_pipe_embedment_restart(int argc, char** argv);

void test_t3d_me_mt_cylinder_foundation_create_model_ch_den(int argc, char** argv);
void test_t3d_me_mt_cylinder_foundation_ch_den(int argc, char** argv);
void test_t3d_me_mt_cylinder_foundation_restart_ch_den(int argc, char** argv);
void test_t3d_me_mt_cylinder_foundation_restart_ch_den2(int argc, char** argv);

void test_t3d_me_mt_test_spudcan_model(int argc, char** argv);
void test_t3d_me_mt_spudcan_coarse_model(int argc, char** argv);
void test_t3d_me_mt_spudcan_coarse(int argc, char** argv);

void test_t3d_chm_mt_1d_consolidation(int argc, char** argv);
void test_t3d_chm_mt_1d_consolidation_restart(int argc, char** argv);
void test_t3d_chm_mt_1d_geostatic(int argc, char** argv);
void test_t3d_chm_mt_1d_geostatic_restart(int argc, char** argv);

void test_sort_pcl_task();
void test_sort_pcl_task2();
void test_sort_tri_mesh_node_task();
void test_sort_tri_mesh_node_task2();
void test_sort_teh_mesh_node_task();

void test_t2d_me_tbb_1d_compression(int argc, char** argv);
void test_t2d_me_tbb_cap_compression(int argc, char** argv);
void test_t2d_chm_tbb_1d_consolidation(int argc, char** argv);

void test_t3d_me_tbb_1d_compression(int argc, char** argv);
void test_t3d_me_tbb_cap_compression(int argc, char** argv);
void test_t3d_chm_tbb_1d_consolidation(int argc, char** argv);
void test_t3d_chm_tbb_cap_compression(int argc, char** argv);

void test_t3d_me_mt_spudcan_sand_hypo_model(int argc, char** argv);
void test_t3d_me_mt_spudcan_sand_hypo(int argc, char** argv);
void test_t3d_me_mt_spudcan_sand_hypo_result(int argc, char** argv);

void test_t3d_chm_mt_spudcan_sand_hypo_model(int argc, char** argv);
void test_t3d_chm_mt_spudcan_sand_hypo(int argc, char** argv);
void test_t3d_chm_mt_spudcan_sand_hypo_result(int argc, char** argv);

void test_t3d_me_mt_piezofoundation_model(int argc, char** argv);
void test_t3d_me_mt_piezofoundation_geo(int argc, char** argv);
void test_t3d_me_mt_piezofoundation(int argc, char** argv);
void test_t3d_me_mt_piezofoundation2(int argc, char** argv);
void test_t3d_me_mt_piezofoundation_geo_result(int argc, char** argv);
void test_t3d_me_mt_piezofoundation_result(int argc, char** argv);

void test_t3d_me_tbb_piezofoundation_sim_mat_model(int argc, char** argv);
void test_t3d_me_tbb_piezofoundation_sim_mat(int argc, char** argv);
void test_t3d_me_tbb_piezofoundation_sim_mat_restart(int argc, char **argv);
void test_t3d_me_mt_piezofoundation_sim_mat(int argc, char** argv);
void test_t3d_me_mt_piezofoundation_sim_mat_restart(int argc, char** argv);

void test_t3d_chm_mt_cap_compression(int argc, char** argv);
void test_t3d_chm_mt_cap_compression_restart(int argc, char** argv);
void test_t3d_chm_mt_cap_compression_result(int argc, char** argv);
void test_t3d_chm_mt_triaxial_compression(int argc, char** argv);

void test_t3d_chm_tbb_piezofoundation_sim_mat_model(int argc, char** argv);
void test_t3d_chm_tbb_piezofoundation_sim_mat(int argc, char** argv);

void test_t3d_me_mt_spudcan_model(int argc, char** argv);
void test_t3d_me_mt_spudcan_geostatic(int argc, char** argv);
void test_t3d_me_mt_spudcan(int argc, char** argv);

void test_t3d_chm_mt_spudcan_model(int argc, char** argv);
void test_t3d_chm_mt_spudcan_geostatic(int argc, char** argv);
void test_t3d_chm_mt_spudcan(int argc, char** argv);

void test_t3d_me_mt_spudcan_cy_model(int argc, char** argv);
void test_t3d_me_mt_spudcan_cy_geostatic(int argc, char** argv);
void test_t3d_me_mt_spudcan_cy(int argc, char** argv);

void test_t3d_chm_mt_spudcan_cy_model(int argc, char** argv);
void test_t3d_chm_mt_spudcan_cy_geostatic(int argc, char** argv);
void test_t3d_chm_mt_spudcan_cy(int argc, char** argv);

void test_t3d_chm_mt_spudcan_cy_HV_model(int argc, char** argv);
void test_t3d_chm_mt_spudcan_cy_HV_geostatic(int argc, char** argv);
void test_t3d_chm_mt_spudcan_cy_HV(int argc, char** argv);

void test_t3d_me_mt_spudcan_cy_Hossain_2006_model(int argc, char** argv);
void test_t3d_me_mt_spudcan_cy_Hossain_2006(int argc, char** argv);
void test_t3d_me_mt_spudcan_cy_Hossain_2006_restart(int argc, char** argv);

#endif