/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_LEVEL3_C_H
#define EL_BLAS_LEVEL3_C_H

#include "El/core/DistMatrix.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Gemm
   ==== */
typedef enum {
  EL_GEMM_DEFAULT, 
  EL_GEMM_SUMMA_A,
  EL_GEMM_SUMMA_B,
  EL_GEMM_SUMMA_C,
  EL_GEMM_SUMMA_DOT,
  EL_GEMM_CANNON
} ElGemmAlgorithm;

EL_EXPORT ElError ElGemm_i
( ElOrientation orientationOfA, ElOrientation orientationOfB,
  ElInt alpha, ElConstMatrix_i A, ElConstMatrix_i B,
  ElInt beta,  ElMatrix_i C );
EL_EXPORT ElError ElGemm_s
( ElOrientation orientationOfA, ElOrientation orientationOfB,
  float alpha, ElConstMatrix_s A, ElConstMatrix_s B,
  float beta,  ElMatrix_s C );
EL_EXPORT ElError ElGemm_d
( ElOrientation orientationOfA, ElOrientation orientationOfB,
  double alpha, ElConstMatrix_d A, ElConstMatrix_d B,
  double beta,  ElMatrix_d C );
EL_EXPORT ElError ElGemm_c
( ElOrientation orientationOfA, ElOrientation orientationOfB,
  complex_float alpha, ElConstMatrix_c A, ElConstMatrix_c B,
  complex_float beta,  ElMatrix_c C );
EL_EXPORT ElError ElGemm_z
( ElOrientation orientationOfA, ElOrientation orientationOfB,
  complex_double alpha, ElConstMatrix_z A, ElConstMatrix_z B,
  complex_double beta,  ElMatrix_z C );

EL_EXPORT ElError ElGemmDist_i
( ElOrientation orientationOfA, ElOrientation orientationOfB,
  ElInt alpha, ElConstDistMatrix_i A, ElConstDistMatrix_i B,
  ElInt beta,  ElDistMatrix_i C );
EL_EXPORT ElError ElGemmDist_s
( ElOrientation orientationOfA, ElOrientation orientationOfB,
  float alpha, ElConstDistMatrix_s A, ElConstDistMatrix_s B,
  float beta,  ElDistMatrix_s C );
EL_EXPORT ElError ElGemmDist_d
( ElOrientation orientationOfA, ElOrientation orientationOfB,
  double alpha, ElConstDistMatrix_d A, ElConstDistMatrix_d B,
  double beta,  ElDistMatrix_d C );
EL_EXPORT ElError ElGemmDist_c
( ElOrientation orientationOfA, ElOrientation orientationOfB,
  complex_float alpha, ElConstDistMatrix_c A, ElConstDistMatrix_c B,
  complex_float beta,  ElDistMatrix_c C );
EL_EXPORT ElError ElGemmDist_z
( ElOrientation orientationOfA, ElOrientation orientationOfB,
  complex_double alpha, ElConstDistMatrix_z A, ElConstDistMatrix_z B,
  complex_double beta,  ElDistMatrix_z C );

/* Expert version
   ^^^^^^^^^^^^^^ */
EL_EXPORT ElError ElGemmXDist_i
( ElOrientation orientationOfA, ElOrientation orientationOfB,
  ElInt alpha, ElConstDistMatrix_i A, ElConstDistMatrix_i B,
  ElInt beta,  ElDistMatrix_i C, ElGemmAlgorithm alg );
EL_EXPORT ElError ElGemmXDist_s
( ElOrientation orientationOfA, ElOrientation orientationOfB,
  float alpha, ElConstDistMatrix_s A, ElConstDistMatrix_s B,
  float beta,  ElDistMatrix_s C, ElGemmAlgorithm alg );
EL_EXPORT ElError ElGemmXDist_d
( ElOrientation orientationOfA, ElOrientation orientationOfB,
  double alpha, ElConstDistMatrix_d A, ElConstDistMatrix_d B,
  double beta,  ElDistMatrix_d C, ElGemmAlgorithm alg );
EL_EXPORT ElError ElGemmXDist_c
( ElOrientation orientationOfA, ElOrientation orientationOfB,
  complex_float alpha, ElConstDistMatrix_c A, ElConstDistMatrix_c B,
  complex_float beta,  ElDistMatrix_c C, ElGemmAlgorithm alg );
EL_EXPORT ElError ElGemmXDist_z
( ElOrientation orientationOfA, ElOrientation orientationOfB,
  complex_double alpha, ElConstDistMatrix_z A, ElConstDistMatrix_z B,
  complex_double beta,  ElDistMatrix_z C, ElGemmAlgorithm alg );

/* Hemm
   ==== */
EL_EXPORT ElError ElHemm_c
( ElLeftOrRight side, ElUpperOrLower uplo, 
  complex_float alpha, ElConstMatrix_c A, ElConstMatrix_c B,
  complex_float beta,  ElMatrix_c C );
EL_EXPORT ElError ElHemm_z
( ElLeftOrRight side, ElUpperOrLower uplo, 
  complex_double alpha, ElConstMatrix_z A, ElConstMatrix_z B,
  complex_double beta,  ElMatrix_z C );

EL_EXPORT ElError ElHemmDist_c
( ElLeftOrRight side, ElUpperOrLower uplo, 
  complex_float alpha, ElConstDistMatrix_c A, ElConstDistMatrix_c B,
  complex_float beta,  ElDistMatrix_c C );
EL_EXPORT ElError ElHemmDist_z
( ElLeftOrRight side, ElUpperOrLower uplo, 
  complex_double alpha, ElConstDistMatrix_z A, ElConstDistMatrix_z B,
  complex_double beta,  ElDistMatrix_z C );

/* Herk
   ==== */
EL_EXPORT ElError ElHerk_c
( ElUpperOrLower uplo, ElOrientation orientation,
  float alpha, ElConstMatrix_c A, float beta, ElMatrix_c C );
EL_EXPORT ElError ElHerk_z
( ElUpperOrLower uplo, ElOrientation orientation,
  double alpha, ElConstMatrix_z A, double beta, ElMatrix_z C );

EL_EXPORT ElError ElHerkDist_c
( ElUpperOrLower uplo, ElOrientation orientation,
  float alpha, ElConstDistMatrix_c A, float beta, ElDistMatrix_c C );
EL_EXPORT ElError ElHerkDist_z
( ElUpperOrLower uplo, ElOrientation orientation,
  double alpha, ElConstDistMatrix_z A, double beta, ElDistMatrix_z C );

EL_EXPORT ElError ElHerkSparse_c
( ElUpperOrLower uplo, ElOrientation orientation,
  float alpha, ElConstSparseMatrix_c A,
  float beta,  ElSparseMatrix_c C );
EL_EXPORT ElError ElHerkSparse_z
( ElUpperOrLower uplo, ElOrientation orientation,
  double alpha, ElConstSparseMatrix_z A,
  double beta,  ElSparseMatrix_z C );

EL_EXPORT ElError ElHerkDistSparse_c
( ElUpperOrLower uplo, ElOrientation orientation,
  float alpha, ElConstDistSparseMatrix_c A,
  float beta,  ElDistSparseMatrix_c C );
EL_EXPORT ElError ElHerkDistSparse_z
( ElUpperOrLower uplo, ElOrientation orientation,
  double alpha, ElConstDistSparseMatrix_z A,
  double beta,  ElDistSparseMatrix_z C );

/* Her2k
   ===== */
EL_EXPORT ElError ElHer2k_c
( ElUpperOrLower uplo, ElOrientation orientation,
  complex_float alpha, ElConstMatrix_c A, ElConstMatrix_c B,
  float beta,          ElMatrix_c C );
EL_EXPORT ElError ElHer2k_z
( ElUpperOrLower uplo, ElOrientation orientation,
  complex_double alpha, ElConstMatrix_z A, ElConstMatrix_z B,
  double beta,          ElMatrix_z C );

EL_EXPORT ElError ElHer2kDist_c
( ElUpperOrLower uplo, ElOrientation orientation,
  complex_float alpha, ElConstDistMatrix_c A, ElConstDistMatrix_c B,
  float beta,          ElDistMatrix_c C );
EL_EXPORT ElError ElHer2kDist_z
( ElUpperOrLower uplo, ElOrientation orientation,
  complex_double alpha, ElConstDistMatrix_z A, ElConstDistMatrix_z B,
  double beta,          ElDistMatrix_z C );

/* Multiply
   ======== */
EL_EXPORT ElError ElMultiply_i
( ElOrientation orientation,
  ElInt alpha, ElConstSparseMatrix_i A, ElConstMatrix_i X, 
  ElInt beta, ElMatrix_i Y );
EL_EXPORT ElError ElMultiply_s
( ElOrientation orientation,
  float alpha, ElConstSparseMatrix_s A, ElConstMatrix_s X, 
  float beta, ElMatrix_s Y );
EL_EXPORT ElError ElMultiply_d
( ElOrientation orientation,
  double alpha, ElConstSparseMatrix_d A, ElConstMatrix_d X, 
  double beta, ElMatrix_d Y );
EL_EXPORT ElError ElMultiply_c
( ElOrientation orientation,
  complex_float alpha, ElConstSparseMatrix_c A, ElConstMatrix_c X, 
  complex_float beta, ElMatrix_c Y );
EL_EXPORT ElError ElMultiply_z
( ElOrientation orientation,
  complex_double alpha, ElConstSparseMatrix_z A, ElConstMatrix_z X, 
  complex_double beta, ElMatrix_z Y );

EL_EXPORT ElError ElMultiplyDist_i
( ElOrientation orientation,
  ElInt alpha, ElConstDistSparseMatrix_i A, ElConstDistMultiVec_i X, 
  ElInt beta, ElDistMultiVec_i Y );
EL_EXPORT ElError ElMultiplyDist_s
( ElOrientation orientation,
  float alpha, ElConstDistSparseMatrix_s A, ElConstDistMultiVec_s X, 
  float beta, ElDistMultiVec_s Y );
EL_EXPORT ElError ElMultiplyDist_d
( ElOrientation orientation,
  double alpha, ElConstDistSparseMatrix_d A, ElConstDistMultiVec_d X, 
  double beta, ElDistMultiVec_d Y );
EL_EXPORT ElError ElMultiplyDist_c
( ElOrientation orientation,
  complex_float alpha, ElConstDistSparseMatrix_c A, ElConstDistMultiVec_c X, 
  complex_float beta, ElDistMultiVec_c Y );
EL_EXPORT ElError ElMultiplyDist_z
( ElOrientation orientation,
  complex_double alpha, ElConstDistSparseMatrix_z A, ElConstDistMultiVec_z X, 
  complex_double beta, ElDistMultiVec_z Y );

/* MultiShiftQuasiTrsm
   =================== */
EL_EXPORT ElError ElMultiShiftQuasiTrsm_s
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  float alpha, ElConstMatrix_s A, ElConstMatrix_s shifts, 
               ElMatrix_s B );
EL_EXPORT ElError ElMultiShiftQuasiTrsm_d
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  double alpha, ElConstMatrix_d A, ElConstMatrix_d shifts, 
                ElMatrix_d B );
EL_EXPORT ElError ElMultiShiftQuasiTrsm_c
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  complex_float alpha, ElConstMatrix_c A, ElConstMatrix_c shifts, 
                       ElMatrix_c B );
EL_EXPORT ElError ElMultiShiftQuasiTrsm_z
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  complex_double alpha, ElConstMatrix_z A, ElConstMatrix_z shifts, 
                        ElMatrix_z B );

EL_EXPORT ElError ElMultiShiftQuasiTrsmDist_s
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  float alpha, ElConstDistMatrix_s A, ElConstDistMatrix_s shifts, 
               ElDistMatrix_s B );
EL_EXPORT ElError ElMultiShiftQuasiTrsmDist_d
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  double alpha, ElConstDistMatrix_d A, ElConstDistMatrix_d shifts, 
                ElDistMatrix_d B );
EL_EXPORT ElError ElMultiShiftQuasiTrsmDist_c
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  complex_float alpha, ElConstDistMatrix_c A, ElConstDistMatrix_c shifts, 
                       ElDistMatrix_c B );
EL_EXPORT ElError ElMultiShiftQuasiTrsmDist_z
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  complex_double alpha, ElConstDistMatrix_z A, ElConstDistMatrix_z shifts, 
                        ElDistMatrix_z B );

/* MultiShiftTrsm
   ============== */
EL_EXPORT ElError ElMultiShiftTrsm_s
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  float alpha, ElMatrix_s A, ElConstMatrix_s shifts, 
               ElMatrix_s B );
EL_EXPORT ElError ElMultiShiftTrsm_d
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  double alpha, ElMatrix_d A, ElConstMatrix_d shifts, 
                ElMatrix_d B );
EL_EXPORT ElError ElMultiShiftTrsm_c
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  complex_float alpha, ElMatrix_c A, ElConstMatrix_c shifts, 
                       ElMatrix_c B );
EL_EXPORT ElError ElMultiShiftTrsm_z
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  complex_double alpha, ElMatrix_z A, ElConstMatrix_z shifts, 
                        ElMatrix_z B );

EL_EXPORT ElError ElMultiShiftTrsmDist_s
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  float alpha, ElConstDistMatrix_s A, ElConstDistMatrix_s shifts, 
               ElDistMatrix_s B );
EL_EXPORT ElError ElMultiShiftTrsmDist_d
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  double alpha, ElConstDistMatrix_d A, ElConstDistMatrix_d shifts, 
                ElDistMatrix_d B );
EL_EXPORT ElError ElMultiShiftTrsmDist_c
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  complex_float alpha, ElConstDistMatrix_c A, ElConstDistMatrix_c shifts, 
                       ElDistMatrix_c B );
EL_EXPORT ElError ElMultiShiftTrsmDist_z
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  complex_double alpha, ElConstDistMatrix_z A, ElConstDistMatrix_z shifts, 
                        ElDistMatrix_z B );

/* SafeMultiShiftTrsm
   ================== */
EL_EXPORT ElError ElSafeMultiShiftTrsm_s
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  float alpha, ElMatrix_s A, ElConstMatrix_s shifts, 
               ElMatrix_s B, ElMatrix_s scales );
EL_EXPORT ElError ElSafeMultiShiftTrsm_d
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  double alpha, ElMatrix_d A, ElConstMatrix_d shifts, 
                ElMatrix_d B, ElMatrix_d scales );
EL_EXPORT ElError ElSafeMultiShiftTrsm_c
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  complex_float alpha, ElMatrix_c A, ElConstMatrix_c shifts, 
                       ElMatrix_c B, ElMatrix_c scales );
EL_EXPORT ElError ElSafeMultiShiftTrsm_z
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  complex_double alpha, ElMatrix_z A, ElConstMatrix_z shifts, 
                        ElMatrix_z B, ElMatrix_z scales );

EL_EXPORT ElError ElSafeMultiShiftTrsmDist_s
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  float alpha, ElConstDistMatrix_s A, ElConstDistMatrix_s shifts, 
               ElDistMatrix_s B, ElDistMatrix_s scales );
EL_EXPORT ElError ElSafeMultiShiftTrsmDist_d
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  double alpha, ElConstDistMatrix_d A, ElConstDistMatrix_d shifts, 
                ElDistMatrix_d B, ElDistMatrix_d scales );
EL_EXPORT ElError ElSafeMultiShiftTrsmDist_c
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  complex_float alpha, ElConstDistMatrix_c A, ElConstDistMatrix_c shifts, 
                       ElDistMatrix_c B, ElDistMatrix_c scales );
EL_EXPORT ElError ElSafeMultiShiftTrsmDist_z
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  complex_double alpha, ElConstDistMatrix_z A, ElConstDistMatrix_z shifts, 
                        ElDistMatrix_z B, ElDistMatrix_z scales );

/* QuasiTrsm
   ========= */
EL_EXPORT ElError ElQuasiTrsm_s
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  float alpha, ElConstMatrix_s A, ElMatrix_s B );
EL_EXPORT ElError ElQuasiTrsm_d
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  double alpha, ElConstMatrix_d A, ElMatrix_d B );
EL_EXPORT ElError ElQuasiTrsm_c
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  complex_float alpha, ElConstMatrix_c A, ElMatrix_c B );
EL_EXPORT ElError ElQuasiTrsm_z
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  complex_double alpha, ElConstMatrix_z A, ElMatrix_z B );

EL_EXPORT ElError ElQuasiTrsmDist_s
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  float alpha, ElConstDistMatrix_s A, ElDistMatrix_s B );
EL_EXPORT ElError ElQuasiTrsmDist_d
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  double alpha, ElConstDistMatrix_d A, ElDistMatrix_d B );
EL_EXPORT ElError ElQuasiTrsmDist_c
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  complex_float alpha, ElConstDistMatrix_c A, ElDistMatrix_c B );
EL_EXPORT ElError ElQuasiTrsmDist_z
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  complex_double alpha, ElConstDistMatrix_z A, ElDistMatrix_z B );

/* Symm
   ==== */
EL_EXPORT ElError ElSymm_s
( ElLeftOrRight side, ElUpperOrLower uplo, 
  float alpha, ElConstMatrix_s A, ElConstMatrix_s B,
  float beta,  ElMatrix_s C );
EL_EXPORT ElError ElSymm_d
( ElLeftOrRight side, ElUpperOrLower uplo, 
  double alpha, ElConstMatrix_d A, ElConstMatrix_d B,
  double beta,  ElMatrix_d C );
EL_EXPORT ElError ElSymm_c
( ElLeftOrRight side, ElUpperOrLower uplo, 
  complex_float alpha, ElConstMatrix_c A, ElConstMatrix_c B,
  complex_float beta,  ElMatrix_c C );
EL_EXPORT ElError ElSymm_z
( ElLeftOrRight side, ElUpperOrLower uplo, 
  complex_double alpha, ElConstMatrix_z A, ElConstMatrix_z B,
  complex_double beta,  ElMatrix_z C );

EL_EXPORT ElError ElSymmDist_s
( ElLeftOrRight side, ElUpperOrLower uplo, 
  float alpha, ElConstDistMatrix_s A, ElConstDistMatrix_s B,
  float beta,  ElDistMatrix_s C );
EL_EXPORT ElError ElSymmDist_d
( ElLeftOrRight side, ElUpperOrLower uplo, 
  double alpha, ElConstDistMatrix_d A, ElConstDistMatrix_d B,
  double beta,  ElDistMatrix_d C );
EL_EXPORT ElError ElSymmDist_c
( ElLeftOrRight side, ElUpperOrLower uplo, 
  complex_float alpha, ElConstDistMatrix_c A, ElConstDistMatrix_c B,
  complex_float beta,  ElDistMatrix_c C );
EL_EXPORT ElError ElSymmDist_z
( ElLeftOrRight side, ElUpperOrLower uplo, 
  complex_double alpha, ElConstDistMatrix_z A, ElConstDistMatrix_z B,
  complex_double beta,  ElDistMatrix_z C );

/* Syrk
   ==== */
EL_EXPORT ElError ElSyrk_s
( ElUpperOrLower uplo, ElOrientation orientation,
  float alpha, ElConstMatrix_s A, 
  float beta,  ElMatrix_s C );
EL_EXPORT ElError ElSyrk_d
( ElUpperOrLower uplo, ElOrientation orientation,
  double alpha, ElConstMatrix_d A, 
  double beta,  ElMatrix_d C );
EL_EXPORT ElError ElSyrk_c
( ElUpperOrLower uplo, ElOrientation orientation,
  complex_float alpha, ElConstMatrix_c A, 
  complex_float beta,  ElMatrix_c C );
EL_EXPORT ElError ElSyrk_z
( ElUpperOrLower uplo, ElOrientation orientation,
  complex_double alpha, ElConstMatrix_z A, 
  complex_double beta,  ElMatrix_z C );

EL_EXPORT ElError ElSyrkDist_s
( ElUpperOrLower uplo, ElOrientation orientation,
  float alpha, ElConstDistMatrix_s A, 
  float beta,  ElDistMatrix_s C );
EL_EXPORT ElError ElSyrkDist_d
( ElUpperOrLower uplo, ElOrientation orientation,
  double alpha, ElConstDistMatrix_d A, 
  double beta,  ElDistMatrix_d C );
EL_EXPORT ElError ElSyrkDist_c
( ElUpperOrLower uplo, ElOrientation orientation,
  complex_float alpha, ElConstDistMatrix_c A, 
  complex_float beta,  ElDistMatrix_c C );
EL_EXPORT ElError ElSyrkDist_z
( ElUpperOrLower uplo, ElOrientation orientation,
  complex_double alpha, ElConstDistMatrix_z A, 
  complex_double beta,  ElDistMatrix_z C );

EL_EXPORT ElError ElSyrkSparse_s
( ElUpperOrLower uplo, ElOrientation orientation,
  float alpha, ElConstSparseMatrix_s A,
  float beta,  ElSparseMatrix_s C );
EL_EXPORT ElError ElSyrkSparse_d
( ElUpperOrLower uplo, ElOrientation orientation,
  double alpha, ElConstSparseMatrix_d A,
  double beta,  ElSparseMatrix_d C );
EL_EXPORT ElError ElSyrkSparse_c
( ElUpperOrLower uplo, ElOrientation orientation,
  complex_float alpha, ElConstSparseMatrix_c A,
  complex_float beta,  ElSparseMatrix_c C );
EL_EXPORT ElError ElSyrkSparse_z
( ElUpperOrLower uplo, ElOrientation orientation,
  complex_double alpha, ElConstSparseMatrix_z A,
  complex_double beta,  ElSparseMatrix_z C );

EL_EXPORT ElError ElSyrkDistSparse_s
( ElUpperOrLower uplo, ElOrientation orientation,
  float alpha, ElConstDistSparseMatrix_s A,
  float beta,  ElDistSparseMatrix_s C );
EL_EXPORT ElError ElSyrkDistSparse_d
( ElUpperOrLower uplo, ElOrientation orientation,
  double alpha, ElConstDistSparseMatrix_d A,
  double beta,  ElDistSparseMatrix_d C );
EL_EXPORT ElError ElSyrkDistSparse_c
( ElUpperOrLower uplo, ElOrientation orientation,
  complex_float alpha, ElConstDistSparseMatrix_c A,
  complex_float beta,  ElDistSparseMatrix_c C );
EL_EXPORT ElError ElSyrkDistSparse_z
( ElUpperOrLower uplo, ElOrientation orientation,
  complex_double alpha, ElConstDistSparseMatrix_z A,
  complex_double beta,  ElDistSparseMatrix_z C );

/* Syr2k
   ===== */
EL_EXPORT ElError ElSyr2k_s
( ElUpperOrLower uplo, ElOrientation orientation,
  float alpha, ElConstMatrix_s A, ElConstMatrix_s B,
  float beta,  ElMatrix_s C );
EL_EXPORT ElError ElSyr2k_d
( ElUpperOrLower uplo, ElOrientation orientation,
  double alpha, ElConstMatrix_d A, ElConstMatrix_d B,
  double beta,  ElMatrix_d C );
EL_EXPORT ElError ElSyr2k_c
( ElUpperOrLower uplo, ElOrientation orientation,
  complex_float alpha, ElConstMatrix_c A, ElConstMatrix_c B,
  complex_float beta,  ElMatrix_c C );
EL_EXPORT ElError ElSyr2k_z
( ElUpperOrLower uplo, ElOrientation orientation,
  complex_double alpha, ElConstMatrix_z A, ElConstMatrix_z B,
  complex_double beta,  ElMatrix_z C );

EL_EXPORT ElError ElSyr2kDist_s
( ElUpperOrLower uplo, ElOrientation orientation,
  float alpha, ElConstDistMatrix_s A, ElConstDistMatrix_s B,
  float beta,  ElDistMatrix_s C );
EL_EXPORT ElError ElSyr2kDist_d
( ElUpperOrLower uplo, ElOrientation orientation,
  double alpha, ElConstDistMatrix_d A, ElConstDistMatrix_d B,
  double beta,  ElDistMatrix_d C );
EL_EXPORT ElError ElSyr2kDist_c
( ElUpperOrLower uplo, ElOrientation orientation,
  complex_float alpha, ElConstDistMatrix_c A, ElConstDistMatrix_c B,
  complex_float beta,  ElDistMatrix_c C );
EL_EXPORT ElError ElSyr2kDist_z
( ElUpperOrLower uplo, ElOrientation orientation,
  complex_double alpha, ElConstDistMatrix_z A, ElConstDistMatrix_z B,
  complex_double beta,  ElDistMatrix_z C );

/* Trdtrmm
   ======= */
EL_EXPORT ElError ElTrdtrmm_s
( ElUpperOrLower uplo, ElMatrix_s A );
EL_EXPORT ElError ElTrdtrmm_d
( ElUpperOrLower uplo, ElMatrix_d A );
EL_EXPORT ElError ElTrdtrmm_c
( ElUpperOrLower uplo, ElMatrix_c A, bool conjugate );
EL_EXPORT ElError ElTrdtrmm_z
( ElUpperOrLower uplo, ElMatrix_z A, bool conjugate );

EL_EXPORT ElError ElTrdtrmmDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A );
EL_EXPORT ElError ElTrdtrmmDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A );
EL_EXPORT ElError ElTrdtrmmDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A, bool conjugate );
EL_EXPORT ElError ElTrdtrmmDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A, bool conjugate );

/* TrdtrmmQuasi
   ============ */
/* TODO: Come up with a better name for this and Trdtrmm and make consistent
         with the C++ API */
EL_EXPORT ElError ElTrdtrmmQuasi_s
( ElUpperOrLower uplo, ElMatrix_s A, ElConstMatrix_s dOff );
EL_EXPORT ElError ElTrdtrmmQuasi_d
( ElUpperOrLower uplo, ElMatrix_d A, ElConstMatrix_d dOff );
EL_EXPORT ElError ElTrdtrmmQuasi_c
( ElUpperOrLower uplo, ElMatrix_c A, ElConstMatrix_c dOff, bool conjugate );
EL_EXPORT ElError ElTrdtrmmQuasi_z
( ElUpperOrLower uplo, ElMatrix_z A, ElConstMatrix_z dOff, bool conjugate );

EL_EXPORT ElError ElTrdtrmmQuasiDist_s
( ElUpperOrLower uplo, 
  ElDistMatrix_s A, ElConstDistMatrix_s dOff );
EL_EXPORT ElError ElTrdtrmmQuasiDist_d
( ElUpperOrLower uplo, 
  ElDistMatrix_d A, ElConstDistMatrix_d dOff );
EL_EXPORT ElError ElTrdtrmmQuasiDist_c
( ElUpperOrLower uplo, 
  ElDistMatrix_c A, ElConstDistMatrix_c dOff, bool conjugate );
EL_EXPORT ElError ElTrdtrmmQuasiDist_z
( ElUpperOrLower uplo, 
  ElDistMatrix_z A, ElConstDistMatrix_z dOff, bool conjugate );

/* Trmm
   ==== */
EL_EXPORT ElError ElTrmm_s
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  float alpha, ElConstMatrix_s A, ElMatrix_s B );
EL_EXPORT ElError ElTrmm_d
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  double alpha, ElConstMatrix_d A, ElMatrix_d B );
EL_EXPORT ElError ElTrmm_c
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  complex_float alpha, ElConstMatrix_c A, ElMatrix_c B );
EL_EXPORT ElError ElTrmm_z
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  complex_double alpha, ElConstMatrix_z A, ElMatrix_z B );

EL_EXPORT ElError ElTrmmDist_s
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  float alpha, ElConstDistMatrix_s A, ElDistMatrix_s B );
EL_EXPORT ElError ElTrmmDist_d
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  double alpha, ElConstDistMatrix_d A, ElDistMatrix_d B );
EL_EXPORT ElError ElTrmmDist_c
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  complex_float alpha, ElConstDistMatrix_c A, ElDistMatrix_c B );
EL_EXPORT ElError ElTrmmDist_z
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  complex_double alpha, ElConstDistMatrix_z A, ElDistMatrix_z B );

/* Trrk
   ==== */
/* TRiangular Rank-K update */
EL_EXPORT ElError ElTrrk_s
( ElUpperOrLower uplo, ElOrientation orientA, ElOrientation orientB, 
  float alpha, ElConstMatrix_s A, ElConstMatrix_s B, 
  float beta,                     ElMatrix_s C );
EL_EXPORT ElError ElTrrk_d
( ElUpperOrLower uplo, ElOrientation orientA, ElOrientation orientB, 
  double alpha, ElConstMatrix_d A, ElConstMatrix_d B, 
  double beta,                     ElMatrix_d C );
EL_EXPORT ElError ElTrrk_c
( ElUpperOrLower uplo, ElOrientation orientA, ElOrientation orientB, 
  complex_float alpha, ElConstMatrix_c A, ElConstMatrix_c B, 
  complex_float beta,                     ElMatrix_c C );
EL_EXPORT ElError ElTrrk_z
( ElUpperOrLower uplo, ElOrientation orientA, ElOrientation orientB, 
  complex_double alpha, ElConstMatrix_z A, ElConstMatrix_z B, 
  complex_double beta,                     ElMatrix_z C );

EL_EXPORT ElError ElTrrkDist_s
( ElUpperOrLower uplo, ElOrientation orientA, ElOrientation orientB, 
  float alpha, ElConstDistMatrix_s A, ElConstDistMatrix_s B, 
  float beta,                         ElDistMatrix_s C );
EL_EXPORT ElError ElTrrkDist_d
( ElUpperOrLower uplo, ElOrientation orientA, ElOrientation orientB, 
  double alpha, ElConstDistMatrix_d A, ElConstDistMatrix_d B, 
  double beta,                         ElDistMatrix_d C );
EL_EXPORT ElError ElTrrkDist_c
( ElUpperOrLower uplo, ElOrientation orientA, ElOrientation orientB, 
  complex_float alpha, ElConstDistMatrix_c A, ElConstDistMatrix_c B, 
  complex_float beta,                         ElDistMatrix_c C );
EL_EXPORT ElError ElTrrkDist_z
( ElUpperOrLower uplo, ElOrientation orientA, ElOrientation orientB, 
  complex_double alpha, ElConstDistMatrix_z A, ElConstDistMatrix_z B, 
  complex_double beta,                         ElDistMatrix_z C );

/* Trr2k
   ===== */
/* TRiangular Rank-2K update */
/* An optimized sequential implementation does not yet exist */
/*
EL_EXPORT ElError ElTrr2k_s
( ElUpperOrLower uplo, 
  ElOrientation orientA, ElOrientation orientB, 
  ElOrientation orientC, ElOrientation orientD,
  float alpha, ElConstMatrix_s A, ElConstMatrix_s B, 
               ElConstMatrix_s C, ElConstMatrix_s D,
  float beta,                     ElMatrix_s E );
EL_EXPORT ElError ElTrr2k_d
( ElUpperOrLower uplo, 
  ElOrientation orientA, ElOrientation orientB, 
  ElOrientation orientC, ElOrientation orientD, 
  double alpha, ElConstMatrix_d A, ElConstMatrix_d B, 
                ElConstMatrix_d C, ElConstMatrix_d D,
  double beta,                     ElMatrix_d E );
EL_EXPORT ElError ElTrr2k_c
( ElUpperOrLower uplo, 
  ElOrientation orientA, ElOrientation orientB, 
  ElOrientation orientC, ElOrientation orientD, 
  complex_float alpha, ElConstMatrix_c A, ElConstMatrix_c B, 
                       ElConstMatrix_c C, ElConstMatrix_c D,
  complex_float beta,                     ElMatrix_c E );
EL_EXPORT ElError ElTrr2k_z
( ElUpperOrLower uplo,
  ElOrientation orientA, ElOrientation orientB, 
  ElOrientation orientC, ElOrientation orientD, 
  complex_double alpha, ElConstMatrix_z A, ElConstMatrix_z B, 
                        ElConstMatrix_z C, ElConstMatrix_z D,
  complex_double beta,                     ElMatrix_z E );
*/

EL_EXPORT ElError ElTrr2kDist_s
( ElUpperOrLower uplo, 
  ElOrientation orientA, ElOrientation orientB, 
  ElOrientation orientC, ElOrientation orientD,
  float alpha, ElConstDistMatrix_s A, ElConstDistMatrix_s B, 
  float beta,  ElConstDistMatrix_s C, ElConstDistMatrix_s D,
  float gamma,                        ElDistMatrix_s E );
EL_EXPORT ElError ElTrr2kDist_d
( ElUpperOrLower uplo, 
  ElOrientation orientA, ElOrientation orientB, 
  ElOrientation orientC, ElOrientation orientD, 
  double alpha, ElConstDistMatrix_d A, ElConstDistMatrix_d B, 
  double beta,  ElConstDistMatrix_d C, ElConstDistMatrix_d D,
  double gamma,                        ElDistMatrix_d E );
EL_EXPORT ElError ElTrr2kDist_c
( ElUpperOrLower uplo, 
  ElOrientation orientA, ElOrientation orientB, 
  ElOrientation orientC, ElOrientation orientD, 
  complex_float alpha, ElConstDistMatrix_c A, ElConstDistMatrix_c B, 
  complex_float beta,  ElConstDistMatrix_c C, ElConstDistMatrix_c D,
  complex_float gamma,                        ElDistMatrix_c E );
EL_EXPORT ElError ElTrr2kDist_z
( ElUpperOrLower uplo,
  ElOrientation orientA, ElOrientation orientB, 
  ElOrientation orientC, ElOrientation orientD, 
  complex_double alpha, ElConstDistMatrix_z A, ElConstDistMatrix_z B, 
  complex_double beta,  ElConstDistMatrix_z C, ElConstDistMatrix_z D,
  complex_double gamma,                        ElDistMatrix_z E );

/* Trsm
   ==== */
EL_EXPORT ElError ElTrsm_s
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  float alpha, ElConstMatrix_s A, ElMatrix_s B );
EL_EXPORT ElError ElTrsm_d
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  double alpha, ElConstMatrix_d A, ElMatrix_d B );
EL_EXPORT ElError ElTrsm_c
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  complex_float alpha, ElConstMatrix_c A, ElMatrix_c B );
EL_EXPORT ElError ElTrsm_z
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  complex_double alpha, ElConstMatrix_z A, ElMatrix_z B );

EL_EXPORT ElError ElTrsmDist_s
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  float alpha, ElConstDistMatrix_s A, ElDistMatrix_s B );
EL_EXPORT ElError ElTrsmDist_d
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  double alpha, ElConstDistMatrix_d A, ElDistMatrix_d B );
EL_EXPORT ElError ElTrsmDist_c
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  complex_float alpha, ElConstDistMatrix_c A, ElDistMatrix_c B );
EL_EXPORT ElError ElTrsmDist_z
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  complex_double alpha, ElConstDistMatrix_z A, ElDistMatrix_z B );

/* Trstrm
   ====== */
EL_EXPORT ElError ElTrstrm_s
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  float alpha, ElConstMatrix_s A, ElMatrix_s B );
EL_EXPORT ElError ElTrstrm_d
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  double alpha, ElConstMatrix_d A, ElMatrix_d B );
EL_EXPORT ElError ElTrstrm_c
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  complex_float alpha, ElConstMatrix_c A, ElMatrix_c B );
EL_EXPORT ElError ElTrstrm_z
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  complex_double alpha, ElConstMatrix_z A, ElMatrix_z B );

EL_EXPORT ElError ElTrstrmDist_s
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  float alpha, ElConstDistMatrix_s A, ElDistMatrix_s B );
EL_EXPORT ElError ElTrstrmDist_d
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  double alpha, ElConstDistMatrix_d A, ElDistMatrix_d B );
EL_EXPORT ElError ElTrstrmDist_c
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  complex_float alpha, ElConstDistMatrix_c A, ElDistMatrix_c B );
EL_EXPORT ElError ElTrstrmDist_z
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  complex_double alpha, ElConstDistMatrix_z A, ElDistMatrix_z B );

/* Trtrmm
   ====== */
EL_EXPORT ElError ElTrtrmm_s
( ElUpperOrLower uplo, ElMatrix_s A );
EL_EXPORT ElError ElTrtrmm_d
( ElUpperOrLower uplo, ElMatrix_d A );
EL_EXPORT ElError ElTrtrmm_c
( ElUpperOrLower uplo, ElMatrix_c A, bool conjugate );
EL_EXPORT ElError ElTrtrmm_z
( ElUpperOrLower uplo, ElMatrix_z A, bool conjugate );

EL_EXPORT ElError ElTrtrmmDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A );
EL_EXPORT ElError ElTrtrmmDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A );
EL_EXPORT ElError ElTrtrmmDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A, bool conjugate );
EL_EXPORT ElError ElTrtrmmDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A, bool conjugate );

/* TwoSidedTrmm
   ============ */
EL_EXPORT ElError ElTwoSidedTrmm_s
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElMatrix_s A, ElConstMatrix_s B );
EL_EXPORT ElError ElTwoSidedTrmm_d
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElMatrix_d A, ElConstMatrix_d B );
EL_EXPORT ElError ElTwoSidedTrmm_c
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElMatrix_c A, ElConstMatrix_c B );
EL_EXPORT ElError ElTwoSidedTrmm_z
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElMatrix_z A, ElConstMatrix_z B );

EL_EXPORT ElError ElTwoSidedTrmmDist_s
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, 
  ElDistMatrix_s A, ElConstDistMatrix_s B );
EL_EXPORT ElError ElTwoSidedTrmmDist_d
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, 
  ElDistMatrix_d A, ElConstDistMatrix_d B );
EL_EXPORT ElError ElTwoSidedTrmmDist_c
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, 
  ElDistMatrix_c A, ElConstDistMatrix_c B );
EL_EXPORT ElError ElTwoSidedTrmmDist_z
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, 
  ElDistMatrix_z A, ElConstDistMatrix_z B );

/* TwoSidedTrsm
   ============ */
EL_EXPORT ElError ElTwoSidedTrsm_s
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElMatrix_s A, ElConstMatrix_s B );
EL_EXPORT ElError ElTwoSidedTrsm_d
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElMatrix_d A, ElConstMatrix_d B );
EL_EXPORT ElError ElTwoSidedTrsm_c
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElMatrix_c A, ElConstMatrix_c B );
EL_EXPORT ElError ElTwoSidedTrsm_z
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElMatrix_z A, ElConstMatrix_z B );

EL_EXPORT ElError ElTwoSidedTrsmDist_s
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, 
  ElDistMatrix_s A, ElConstDistMatrix_s B );
EL_EXPORT ElError ElTwoSidedTrsmDist_d
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, 
  ElDistMatrix_d A, ElConstDistMatrix_d B );
EL_EXPORT ElError ElTwoSidedTrsmDist_c
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, 
  ElDistMatrix_c A, ElConstDistMatrix_c B );
EL_EXPORT ElError ElTwoSidedTrsmDist_z
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, 
  ElDistMatrix_z A, ElConstDistMatrix_z B );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_BLAS_LEVEL3_C_H */
