#ifndef __Symmetric_Matrix_Eigen_Value_h__
#define __Symmetric_Matrix_Eigen_Value_h__

#include "MaterialModelPrecision.h"

void cal_sym_mat_eigen(const __Float_Type__ mat_va[6], __Float_Type__ ev[3]);

// evecs[3][3]^T * mat_va[3][3] * evecs[3][3] = ev[3][3]
// mat_va[3][3] = evecs[3][3] * ev[3][3] * evecs[3][3]^T
// [v1, v2, v3]
void cal_sym_mat_eigen(const __Float_Type__ mat_va[6], __Float_Type__ ev[3], __Float_Type__ evecs[3][3]);

#endif