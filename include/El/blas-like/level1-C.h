/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_BLAS_LEVEL1_C_H
#define EL_BLAS_LEVEL1_C_H

#include "El/core/DistMatrix-C.h"

#ifdef __cplusplus
extern "C" {
#endif

/* B = A^H 
   ------- */
ElError ElAdjointMatrix_i( ElConstMatrix_i A, ElMatrix_i B );
ElError ElAdjointMatrix_s( ElConstMatrix_s A, ElMatrix_s B );
ElError ElAdjointMatrix_d( ElConstMatrix_d A, ElMatrix_d B );
ElError ElAdjointMatrix_c( ElConstMatrix_c A, ElMatrix_c B );
ElError ElAdjointMatrix_z( ElConstMatrix_z A, ElMatrix_z B );

ElError ElAdjointDistMatrix_i( ElConstDistMatrix_i A, ElDistMatrix_i B );
ElError ElAdjointDistMatrix_s( ElConstDistMatrix_s A, ElDistMatrix_s B );
ElError ElAdjointDistMatrix_d( ElConstDistMatrix_d A, ElDistMatrix_d B );
ElError ElAdjointDistMatrix_c( ElConstDistMatrix_c A, ElDistMatrix_c B );
ElError ElAdjointDistMatrix_z( ElConstDistMatrix_z A, ElDistMatrix_z B );

/* Y := alpha A X + Y 
   ------------------ */
ElError ElAxpyMatrix_i( ElInt alpha, ElConstMatrix_i X, ElMatrix_i Y );
ElError ElAxpyMatrix_s( float alpha, ElConstMatrix_s X, ElMatrix_s Y );
ElError ElAxpyMatrix_d( double alpha, ElConstMatrix_d X, ElMatrix_d Y );
ElError ElAxpyMatrix_c( complex_float alpha, ElConstMatrix_c X, ElMatrix_c Y );
ElError ElAxpyMatrix_z( complex_double alpha, ElConstMatrix_z X, ElMatrix_z Y );

ElError ElAxpyDistMatrix_i
( ElInt alpha, ElConstDistMatrix_i X, ElDistMatrix_i Y );
ElError ElAxpyDistMatrix_s
( float alpha, ElConstDistMatrix_s X, ElDistMatrix_s Y );
ElError ElAxpyDistMatrix_d
( double alpha, ElConstDistMatrix_d X, ElDistMatrix_d Y );
ElError ElAxpyDistMatrix_c
( complex_float alpha, ElConstDistMatrix_c X, ElDistMatrix_c Y );
ElError ElAxpyDistMatrix_z
( complex_double alpha, ElConstDistMatrix_z X, ElDistMatrix_z Y );

/* tri(Y) := tri(alpha A X + Y)
   ---------------------------- */
ElError ElAxpyTriangleMatrix_i
( ElUpperOrLower uplo, ElInt alpha, ElConstMatrix_i X, ElMatrix_i Y );
ElError ElAxpyTriangleMatrix_s
( ElUpperOrLower uplo, float alpha, ElConstMatrix_s X, ElMatrix_s Y );
ElError ElAxpyTriangleMatrix_d
( ElUpperOrLower uplo, double alpha, ElConstMatrix_d X, ElMatrix_d Y );
ElError ElAxpyTriangleMatrix_c
( ElUpperOrLower uplo, complex_float alpha, ElConstMatrix_c X, ElMatrix_c Y );
ElError ElAxpyTriangleMatrix_z
( ElUpperOrLower uplo, complex_double alpha, ElConstMatrix_z X, ElMatrix_z Y );

ElError ElAxpyTriangleDistMatrix_i
( ElUpperOrLower uplo, ElInt alpha, 
  ElConstDistMatrix_i X, ElDistMatrix_i Y );
ElError ElAxpyTriangleDistMatrix_s
( ElUpperOrLower uplo, float alpha, 
  ElConstDistMatrix_s X, ElDistMatrix_s Y );
ElError ElAxpyTriangleDistMatrix_d
( ElUpperOrLower uplo, double alpha, 
  ElConstDistMatrix_d X, ElDistMatrix_d Y );
ElError ElAxpyTriangleDistMatrix_c
( ElUpperOrLower uplo, complex_float alpha, 
  ElConstDistMatrix_c X, ElDistMatrix_c Y );
ElError ElAxpyTriangleDistMatrix_z
( ElUpperOrLower uplo, complex_double alpha, 
  ElConstDistMatrix_z X, ElDistMatrix_z Y );

/* B = A 
   ----- */
ElError ElCopyMatrix_i( ElConstMatrix_i A, ElMatrix_i B );
ElError ElCopyMatrix_s( ElConstMatrix_s A, ElMatrix_s B );
ElError ElCopyMatrix_d( ElConstMatrix_d A, ElMatrix_d B );
ElError ElCopyMatrix_c( ElConstMatrix_c A, ElMatrix_c B );
ElError ElCopyMatrix_z( ElConstMatrix_z A, ElMatrix_z B );

ElError ElCopyDistMatrix_i( ElConstDistMatrix_i A, ElDistMatrix_i B );
ElError ElCopyDistMatrix_s( ElConstDistMatrix_s A, ElDistMatrix_s B );
ElError ElCopyDistMatrix_d( ElConstDistMatrix_d A, ElDistMatrix_d B );
ElError ElCopyDistMatrix_c( ElConstDistMatrix_c A, ElDistMatrix_c B );
ElError ElCopyDistMatrix_z( ElConstDistMatrix_z A, ElDistMatrix_z B );

/* DiagonalScale 
   ------------- */
ElError ElDiagonalScaleMatrix_i
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_i d, ElMatrix_i X );
ElError ElDiagonalScaleMatrix_s
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_s d, ElMatrix_s X );
ElError ElDiagonalScaleMatrix_d
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_d d, ElMatrix_d X );
ElError ElDiagonalScaleMatrix_c
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_c d, ElMatrix_c X );
ElError ElDiagonalScaleMatrix_z
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstMatrix_z d, ElMatrix_z X );

ElError ElDiagonalScaleDistMatrix_i
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_i d, ElDistMatrix_i X );
ElError ElDiagonalScaleDistMatrix_s
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_s d, ElDistMatrix_s X );
ElError ElDiagonalScaleDistMatrix_d
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_d d, ElDistMatrix_d X );
ElError ElDiagonalScaleDistMatrix_c
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_c d, ElDistMatrix_c X );
ElError ElDiagonalScaleDistMatrix_z
( ElLeftOrRight side, ElOrientation orientation, 
  ElConstDistMatrix_z d, ElDistMatrix_z X );

/* B = A^T 
   ------- */
ElError ElTransposeMatrix_i( ElConstMatrix_i A, ElMatrix_i B );
ElError ElTransposeMatrix_s( ElConstMatrix_s A, ElMatrix_s B );
ElError ElTransposeMatrix_d( ElConstMatrix_d A, ElMatrix_d B );
ElError ElTransposeMatrix_c( ElConstMatrix_c A, ElMatrix_c B );
ElError ElTransposeMatrix_z( ElConstMatrix_z A, ElMatrix_z B );

ElError ElTransposeDistMatrix_i( ElConstDistMatrix_i A, ElDistMatrix_i B );
ElError ElTransposeDistMatrix_s( ElConstDistMatrix_s A, ElDistMatrix_s B );
ElError ElTransposeDistMatrix_d( ElConstDistMatrix_d A, ElDistMatrix_d B );
ElError ElTransposeDistMatrix_c( ElConstDistMatrix_c A, ElDistMatrix_c B );
ElError ElTransposeDistMatrix_z( ElConstDistMatrix_z A, ElDistMatrix_z B );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_BLAS_LEVEL1_C_H */
