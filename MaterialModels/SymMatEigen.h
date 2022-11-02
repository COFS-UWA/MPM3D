#ifndef __Symmetric_Matrix_Eigen_Value_h__
#define __Symmetric_Matrix_Eigen_Value_h__

#include "MaterialModelPrecision.h"

void cal_sym_mat_eigen(const __Float_Type__ mat_va[6], __Float_Type__ ev[3]);

// evecs^T * mat_va * evecs = ev
// evecs * ev * evecs^T = mat_va
void cal_sym_mat_eigen(const __Float_Type__ mat_va[6], __Float_Type__ ev[3], __Float_Type__ evecs[3][3]);

void rotate_sym_mat_to_eigen_mat(const __Float_Type__ mat[6], const __Float_Type__ evecs[3][3], __Float_Type__ evs[3]);

void rotate_eigen_mat_to_sym_mat(const __Float_Type__ evs[3], const __Float_Type__ evecs[3][3], __Float_Type__ mat[6]);

#endif