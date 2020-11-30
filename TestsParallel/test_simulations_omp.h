#ifndef __Test_Simulations_OMP_h__
#define __Test_Simulations_OMP_h__

void test_von_mises(int argc, char** argv);

void test_rigid_cylinder(int argc, char** argv);
void test_rigid_cone(int argc, char** argv);
void test_contact_model_3d(int argc, char** argv);

void test_Step_OMP(int argc, char** argv);

void test_t2d_me_mt_test1(int argc, char** argv);
void test_t2d_me_mt_test2(int argc, char** argv);
void test_t2d_me_s_test2(int argc, char** argv);
void test_t2d_me_mt_cap_compression(int argc, char** argv);
void test_t2d_me_mt_strip_footing_smaller(int argc, char** argv);
void test_t2d_me_mt_strip_footing(int argc, char** argv);

void test_t3d_me_mt_1d_compression(int argc, char** argv);
void test_t3d_me_s_1d_compression(int argc, char** argv);

void test_t3d_me_mt_cap_compression(int argc, char** argv);
void test_t3d_me_mt_cap_compression_restart(int argc, char** argv);

void test_omp_barrier_time();
void test_t3d_me_mt_cylinder_foundation_create_model(int argc, char** argv);
void test_t3d_me_mt_cylinder_foundation(int argc, char** argv);

void test_t3d_me_mt_block_sliding(int argc, char** argv);

#endif