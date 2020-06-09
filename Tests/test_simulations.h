#ifndef __Test_Simulations_h__
#define __Test_Simulations_h__

void test_solve_functions();
void test_Model_T3D_ME_s_display(int argc, char **argv);

void test_t3d_me_s_1d_compression(int argc, char **argv);
void test_t3d_chm_s_1d_consolidation(int argc, char** argv);

// test on a simple cube
void test_fem_t3d_me_s_1d_compression(int argc, char** argv);

// test on more complex column
void test_fem_t3d_me_s_1d_compression2(int argc, char** argv);

#endif