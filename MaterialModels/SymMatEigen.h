#ifndef __Symmetric_Matrix_Eigen_Value_h__
#define __Symmetric_Matrix_Eigen_Value_h__

#include "MaterialModelPrecision.h"

void cal_sym_mat_eigen(const __Float_Type__ mat_va[6], __Float_Type__ ev[3]);

void cal_sym_mat_eigen(const __Float_Type__ mat_va[6], __Float_Type__ ev[3], __Float_Type__ evecs[3][3]);

#endif