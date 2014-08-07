/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_BLAS_LEVEL3_C_H
#define EL_BLAS_LEVEL3_C_H

#include "El/core/DistMatrix-C.h"

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

ElError ElGemm_i
( ElOrientation orientationOfA, ElOrientation orientationOfB,
  ElInt alpha, ElConstMatrix_i A, ElConstMatrix_i B,
  ElInt beta,  ElMatrix_i C );
ElError ElGemm_s
( ElOrientation orientationOfA, ElOrientation orientationOfB,
  float alpha, ElConstMatrix_s A, ElConstMatrix_s B,
  float beta,  ElMatrix_s C );
ElError ElGemm_d
( ElOrientation orientationOfA, ElOrientation orientationOfB,
  double alpha, ElConstMatrix_d A, ElConstMatrix_d B,
  double beta,  ElMatrix_d C );
ElError ElGemm_c
( ElOrientation orientationOfA, ElOrientation orientationOfB,
  complex_float alpha, ElConstMatrix_c A, ElConstMatrix_c B,
  complex_float beta,  ElMatrix_c C );
ElError ElGemm_z
( ElOrientation orientationOfA, ElOrientation orientationOfB,
  complex_double alpha, ElConstMatrix_z A, ElConstMatrix_z B,
  complex_double beta,  ElMatrix_z C );

ElError ElGemmDist_i
( ElOrientation orientationOfA, ElOrientation orientationOfB,
  ElInt alpha, ElConstDistMatrix_i A, ElConstDistMatrix_i B,
  ElInt beta,  ElDistMatrix_i C );
ElError ElGemmDist_s
( ElOrientation orientationOfA, ElOrientation orientationOfB,
  float alpha, ElConstDistMatrix_s A, ElConstDistMatrix_s B,
  float beta,  ElDistMatrix_s C );
ElError ElGemmDist_d
( ElOrientation orientationOfA, ElOrientation orientationOfB,
  double alpha, ElConstDistMatrix_d A, ElConstDistMatrix_d B,
  double beta,  ElDistMatrix_d C );
ElError ElGemmDist_c
( ElOrientation orientationOfA, ElOrientation orientationOfB,
  complex_float alpha, ElConstDistMatrix_c A, ElConstDistMatrix_c B,
  complex_float beta,  ElDistMatrix_c C );
ElError ElGemmDist_z
( ElOrientation orientationOfA, ElOrientation orientationOfB,
  complex_double alpha, ElConstDistMatrix_z A, ElConstDistMatrix_z B,
  complex_double beta,  ElDistMatrix_z C );

/* Expert version
   ^^^^^^^^^^^^^^ */
ElError ElGemmXDist_i
( ElOrientation orientationOfA, ElOrientation orientationOfB,
  ElInt alpha, ElConstDistMatrix_i A, ElConstDistMatrix_i B,
  ElInt beta,  ElDistMatrix_i C, ElGemmAlgorithm alg );
ElError ElGemmXDist_s
( ElOrientation orientationOfA, ElOrientation orientationOfB,
  float alpha, ElConstDistMatrix_s A, ElConstDistMatrix_s B,
  float beta,  ElDistMatrix_s C, ElGemmAlgorithm alg );
ElError ElGemmXDist_d
( ElOrientation orientationOfA, ElOrientation orientationOfB,
  double alpha, ElConstDistMatrix_d A, ElConstDistMatrix_d B,
  double beta,  ElDistMatrix_d C, ElGemmAlgorithm alg );
ElError ElGemmXDist_c
( ElOrientation orientationOfA, ElOrientation orientationOfB,
  complex_float alpha, ElConstDistMatrix_c A, ElConstDistMatrix_c B,
  complex_float beta,  ElDistMatrix_c C, ElGemmAlgorithm alg );
ElError ElGemmXDist_z
( ElOrientation orientationOfA, ElOrientation orientationOfB,
  complex_double alpha, ElConstDistMatrix_z A, ElConstDistMatrix_z B,
  complex_double beta,  ElDistMatrix_z C, ElGemmAlgorithm alg );

/* Hemm
   ==== */
ElError ElHemm_c
( ElLeftOrRight side, ElUpperOrLower uplo, 
  complex_float alpha, ElConstMatrix_c A, ElConstMatrix_c B,
  complex_float beta,  ElMatrix_c C );
ElError ElHemm_z
( ElLeftOrRight side, ElUpperOrLower uplo, 
  complex_double alpha, ElConstMatrix_z A, ElConstMatrix_z B,
  complex_double beta,  ElMatrix_z C );

ElError ElHemmDist_c
( ElLeftOrRight side, ElUpperOrLower uplo, 
  complex_float alpha, ElConstDistMatrix_c A, ElConstDistMatrix_c B,
  complex_float beta,  ElDistMatrix_c C );
ElError ElHemmDist_z
( ElLeftOrRight side, ElUpperOrLower uplo, 
  complex_double alpha, ElConstDistMatrix_z A, ElConstDistMatrix_z B,
  complex_double beta,  ElDistMatrix_z C );

/* Herk
   ==== */
ElError ElHerk_c
( ElUpperOrLower uplo, ElOrientation orientation,
  complex_float alpha, ElConstMatrix_c A, 
  complex_float beta,  ElMatrix_c C );
ElError ElHerk_z
( ElUpperOrLower uplo, ElOrientation orientation,
  complex_double alpha, ElConstMatrix_z A, 
  complex_double beta,  ElMatrix_z C );

ElError ElHerkDist_c
( ElUpperOrLower uplo, ElOrientation orientation,
  complex_float alpha, ElConstDistMatrix_c A, 
  complex_float beta,  ElDistMatrix_c C );
ElError ElHerkDist_z
( ElUpperOrLower uplo, ElOrientation orientation,
  complex_double alpha, ElConstDistMatrix_z A, 
  complex_double beta,  ElDistMatrix_z C );

/* Her2k
   ===== */
ElError ElHer2k_c
( ElUpperOrLower uplo, ElOrientation orientation,
  complex_float alpha, ElConstMatrix_c A, ElConstMatrix_c B,
  complex_float beta,  ElMatrix_c C );
ElError ElHer2k_z
( ElUpperOrLower uplo, ElOrientation orientation,
  complex_double alpha, ElConstMatrix_z A, ElConstMatrix_z B,
  complex_double beta,  ElMatrix_z C );

ElError ElHer2kDist_c
( ElUpperOrLower uplo, ElOrientation orientation,
  complex_float alpha, ElConstDistMatrix_c A, ElConstDistMatrix_c B,
  complex_float beta,  ElDistMatrix_c C );
ElError ElHer2kDist_z
( ElUpperOrLower uplo, ElOrientation orientation,
  complex_double alpha, ElConstDistMatrix_z A, ElConstDistMatrix_z B,
  complex_double beta,  ElDistMatrix_z C );

/* MultiShiftQuasiTrsm
   =================== */
ElError ElMultiShiftQuasiTrsm_s
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  float alpha, ElConstMatrix_s A, ElConstMatrix_s shifts, 
               ElMatrix_s B );
ElError ElMultiShiftQuasiTrsm_d
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  double alpha, ElConstMatrix_d A, ElConstMatrix_d shifts, 
                ElMatrix_d B );
ElError ElMultiShiftQuasiTrsm_c
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  complex_float alpha, ElConstMatrix_c A, ElConstMatrix_c shifts, 
                       ElMatrix_c B );
ElError ElMultiShiftQuasiTrsm_z
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  complex_double alpha, ElConstMatrix_z A, ElConstMatrix_z shifts, 
                        ElMatrix_z B );


/* NOTE: 'A' and 'B' must be in [MC,MR] distributions, while
         'shifts' must be in a [VR,STAR] distribution */
ElError ElMultiShiftQuasiTrsmDist_s
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  float alpha, ElConstDistMatrix_s A, ElConstDistMatrix_s shifts, 
               ElDistMatrix_s B );
ElError ElMultiShiftQuasiTrsmDist_d
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  double alpha, ElConstDistMatrix_d A, ElConstDistMatrix_d shifts, 
                ElDistMatrix_d B );
ElError ElMultiShiftQuasiTrsmDist_c
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  complex_float alpha, ElConstDistMatrix_c A, ElConstDistMatrix_c shifts, 
                       ElDistMatrix_c B );
ElError ElMultiShiftQuasiTrsmDist_z
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  complex_double alpha, ElConstDistMatrix_z A, ElConstDistMatrix_z shifts, 
                        ElDistMatrix_z B );

/* MultiShiftTrsm
   ============== */
ElError ElMultiShiftTrsm_s
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  float alpha, ElMatrix_s A, ElConstMatrix_s shifts, 
               ElMatrix_s B );
ElError ElMultiShiftTrsm_d
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  double alpha, ElMatrix_d A, ElConstMatrix_d shifts, 
                ElMatrix_d B );
ElError ElMultiShiftTrsm_c
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  complex_float alpha, ElMatrix_c A, ElConstMatrix_c shifts, 
                       ElMatrix_c B );
ElError ElMultiShiftTrsm_z
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  complex_double alpha, ElMatrix_z A, ElConstMatrix_z shifts, 
                        ElMatrix_z B );

/* NOTE: 'A' and 'B' must be in [MC,MR] distributions, while
         'shifts' must be in a [VR,STAR] distribution */
ElError ElMultiShiftTrsmDist_s
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  float alpha, ElConstDistMatrix_s A, ElConstDistMatrix_s shifts, 
               ElDistMatrix_s B );
ElError ElMultiShiftTrsmDist_d
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  double alpha, ElConstDistMatrix_d A, ElConstDistMatrix_d shifts, 
                ElDistMatrix_d B );
ElError ElMultiShiftTrsmDist_c
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  complex_float alpha, ElConstDistMatrix_c A, ElConstDistMatrix_c shifts, 
                       ElDistMatrix_c B );
ElError ElMultiShiftTrsmDist_z
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  complex_double alpha, ElConstDistMatrix_z A, ElConstDistMatrix_z shifts, 
                        ElDistMatrix_z B );

/* QuasiTrsm
   ========= */
ElError ElQuasiTrsm_s
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  float alpha, ElConstMatrix_s A, ElMatrix_s B );
ElError ElQuasiTrsm_d
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  double alpha, ElConstMatrix_d A, ElMatrix_d B );
ElError ElQuasiTrsm_c
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  complex_float alpha, ElConstMatrix_c A, ElMatrix_c B );
ElError ElQuasiTrsm_z
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  complex_double alpha, ElConstMatrix_z A, ElMatrix_z B );

/* NOTE: 'A' and 'B' must be in [MC,MR] distributions */
ElError ElQuasiTrsmDist_s
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  float alpha, ElConstDistMatrix_s A, ElDistMatrix_s B );
ElError ElQuasiTrsmDist_d
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  double alpha, ElConstDistMatrix_d A, ElDistMatrix_d B );
ElError ElQuasiTrsmDist_c
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  complex_float alpha, ElConstDistMatrix_c A, ElDistMatrix_c B );
ElError ElQuasiTrsmDist_z
( ElLeftOrRight side, ElUpperOrLower uplo, ElOrientation orientation,
  complex_double alpha, ElConstDistMatrix_z A, ElDistMatrix_z B );

/* Symm
   ==== */
ElError ElSymm_s
( ElLeftOrRight side, ElUpperOrLower uplo, 
  float alpha, ElConstMatrix_s A, ElConstMatrix_s B,
  float beta,  ElMatrix_s C );
ElError ElSymm_d
( ElLeftOrRight side, ElUpperOrLower uplo, 
  double alpha, ElConstMatrix_d A, ElConstMatrix_d B,
  double beta,  ElMatrix_d C );
ElError ElSymm_c
( ElLeftOrRight side, ElUpperOrLower uplo, 
  complex_float alpha, ElConstMatrix_c A, ElConstMatrix_c B,
  complex_float beta,  ElMatrix_c C );
ElError ElSymm_z
( ElLeftOrRight side, ElUpperOrLower uplo, 
  complex_double alpha, ElConstMatrix_z A, ElConstMatrix_z B,
  complex_double beta,  ElMatrix_z C );

ElError ElSymmDist_s
( ElLeftOrRight side, ElUpperOrLower uplo, 
  float alpha, ElConstDistMatrix_s A, ElConstDistMatrix_s B,
  float beta,  ElDistMatrix_s C );
ElError ElSymmDist_d
( ElLeftOrRight side, ElUpperOrLower uplo, 
  double alpha, ElConstDistMatrix_d A, ElConstDistMatrix_d B,
  double beta,  ElDistMatrix_d C );
ElError ElSymmDist_c
( ElLeftOrRight side, ElUpperOrLower uplo, 
  complex_float alpha, ElConstDistMatrix_c A, ElConstDistMatrix_c B,
  complex_float beta,  ElDistMatrix_c C );
ElError ElSymmDist_z
( ElLeftOrRight side, ElUpperOrLower uplo, 
  complex_double alpha, ElConstDistMatrix_z A, ElConstDistMatrix_z B,
  complex_double beta,  ElDistMatrix_z C );

/* Syrk
   ==== */
ElError ElSyrk_s
( ElUpperOrLower uplo, ElOrientation orientation,
  float alpha, ElConstMatrix_s A, 
  float beta,  ElMatrix_s C );
ElError ElSyrk_d
( ElUpperOrLower uplo, ElOrientation orientation,
  double alpha, ElConstMatrix_d A, 
  double beta,  ElMatrix_d C );
ElError ElSyrk_c
( ElUpperOrLower uplo, ElOrientation orientation,
  complex_float alpha, ElConstMatrix_c A, 
  complex_float beta,  ElMatrix_c C );
ElError ElSyrk_z
( ElUpperOrLower uplo, ElOrientation orientation,
  complex_double alpha, ElConstMatrix_z A, 
  complex_double beta,  ElMatrix_z C );

ElError ElSyrkDist_s
( ElUpperOrLower uplo, ElOrientation orientation,
  float alpha, ElConstDistMatrix_s A, 
  float beta,  ElDistMatrix_s C );
ElError ElSyrkDist_d
( ElUpperOrLower uplo, ElOrientation orientation,
  double alpha, ElConstDistMatrix_d A, 
  double beta,  ElDistMatrix_d C );
ElError ElSyrkDist_c
( ElUpperOrLower uplo, ElOrientation orientation,
  complex_float alpha, ElConstDistMatrix_c A, 
  complex_float beta,  ElDistMatrix_c C );
ElError ElSyrkDist_z
( ElUpperOrLower uplo, ElOrientation orientation,
  complex_double alpha, ElConstDistMatrix_z A, 
  complex_double beta,  ElDistMatrix_z C );

/* Syr2k
   ===== */
ElError ElSyr2k_s
( ElUpperOrLower uplo, ElOrientation orientation,
  float alpha, ElConstMatrix_s A, ElConstMatrix_s B,
  float beta,  ElMatrix_s C );
ElError ElSyr2k_d
( ElUpperOrLower uplo, ElOrientation orientation,
  double alpha, ElConstMatrix_d A, ElConstMatrix_d B,
  double beta,  ElMatrix_d C );
ElError ElSyr2k_c
( ElUpperOrLower uplo, ElOrientation orientation,
  complex_float alpha, ElConstMatrix_c A, ElConstMatrix_c B,
  complex_float beta,  ElMatrix_c C );
ElError ElSyr2k_z
( ElUpperOrLower uplo, ElOrientation orientation,
  complex_double alpha, ElConstMatrix_z A, ElConstMatrix_z B,
  complex_double beta,  ElMatrix_z C );

ElError ElSyr2kDist_s
( ElUpperOrLower uplo, ElOrientation orientation,
  float alpha, ElConstDistMatrix_s A, ElConstDistMatrix_s B,
  float beta,  ElDistMatrix_s C );
ElError ElSyr2kDist_d
( ElUpperOrLower uplo, ElOrientation orientation,
  double alpha, ElConstDistMatrix_d A, ElConstDistMatrix_d B,
  double beta,  ElDistMatrix_d C );
ElError ElSyr2kDist_c
( ElUpperOrLower uplo, ElOrientation orientation,
  complex_float alpha, ElConstDistMatrix_c A, ElConstDistMatrix_c B,
  complex_float beta,  ElDistMatrix_c C );
ElError ElSyr2kDist_z
( ElUpperOrLower uplo, ElOrientation orientation,
  complex_double alpha, ElConstDistMatrix_z A, ElConstDistMatrix_z B,
  complex_double beta,  ElDistMatrix_z C );

/* Trdtrmm
   ======= */
ElError ElTrdtrmm_s( ElUpperOrLower uplo, ElMatrix_s A, bool conjugate );
ElError ElTrdtrmm_d( ElUpperOrLower uplo, ElMatrix_d A, bool conjugate );
ElError ElTrdtrmm_c( ElUpperOrLower uplo, ElMatrix_c A, bool conjugate );
ElError ElTrdtrmm_z( ElUpperOrLower uplo, ElMatrix_z A, bool conjugate );

/* NOTE: 'A' must be in a [MC,MR] distribution */
ElError ElTrdtrmmDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A, bool conjugate );
ElError ElTrdtrmmDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A, bool conjugate );
ElError ElTrdtrmmDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A, bool conjugate );
ElError ElTrdtrmmDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A, bool conjugate );

/* TrdtrmmQuasi
   ============ */
/* TODO: Come up with a better name for this and Trdtrmm and make consistent
         with the C++ API */
ElError ElTrdtrmmQuasi_s
( ElUpperOrLower uplo, ElMatrix_s A, ElConstMatrix_s dOff, bool conjugate );
ElError ElTrdtrmmQuasi_d
( ElUpperOrLower uplo, ElMatrix_d A, ElConstMatrix_d dOff, bool conjugate );
ElError ElTrdtrmmQuasi_c
( ElUpperOrLower uplo, ElMatrix_c A, ElConstMatrix_c dOff, bool conjugate );
ElError ElTrdtrmmQuasi_z
( ElUpperOrLower uplo, ElMatrix_z A, ElConstMatrix_z dOff, bool conjugate );

/* NOTE: 'A' must be in a [MC,MR] distribution, while 
         'dOff' must be in a [MD,STAR] distribution */
ElError ElTrdtrmmQuasiDist_s
( ElUpperOrLower uplo, 
  ElDistMatrix_s A, ElConstDistMatrix_s dOff, bool conjugate );
ElError ElTrdtrmmQuasiDist_d
( ElUpperOrLower uplo, 
  ElDistMatrix_d A, ElConstDistMatrix_d dOff, bool conjugate );
ElError ElTrdtrmmQuasiDist_c
( ElUpperOrLower uplo, 
  ElDistMatrix_c A, ElConstDistMatrix_c dOff, bool conjugate );
ElError ElTrdtrmmQuasiDist_z
( ElUpperOrLower uplo, 
  ElDistMatrix_z A, ElConstDistMatrix_z dOff, bool conjugate );

/* Trmm
   ==== */
ElError ElTrmm_s
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  float alpha, ElConstMatrix_s A, ElMatrix_s B );
ElError ElTrmm_d
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  double alpha, ElConstMatrix_d A, ElMatrix_d B );
ElError ElTrmm_c
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  complex_float alpha, ElConstMatrix_c A, ElMatrix_c B );
ElError ElTrmm_z
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  complex_double alpha, ElConstMatrix_z A, ElMatrix_z B );

ElError ElTrmmDist_s
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  float alpha, ElConstDistMatrix_s A, ElDistMatrix_s B );
ElError ElTrmmDist_d
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  double alpha, ElConstDistMatrix_d A, ElDistMatrix_d B );
ElError ElTrmmDist_c
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  complex_float alpha, ElConstDistMatrix_c A, ElDistMatrix_c B );
ElError ElTrmmDist_z
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  complex_double alpha, ElConstDistMatrix_z A, ElDistMatrix_z B );

/* Trrk
   ==== */
/* TRiangular Rank-K update */
ElError ElTrrk_s
( ElUpperOrLower uplo, ElOrientation orientA, ElOrientation orientB, 
  float alpha, ElConstMatrix_s A, ElConstMatrix_s B, 
  float beta,                     ElMatrix_s C );
ElError ElTrrk_d
( ElUpperOrLower uplo, ElOrientation orientA, ElOrientation orientB, 
  double alpha, ElConstMatrix_d A, ElConstMatrix_d B, 
  double beta,                     ElMatrix_d C );
ElError ElTrrk_c
( ElUpperOrLower uplo, ElOrientation orientA, ElOrientation orientB, 
  complex_float alpha, ElConstMatrix_c A, ElConstMatrix_c B, 
  complex_float beta,                     ElMatrix_c C );
ElError ElTrrk_z
( ElUpperOrLower uplo, ElOrientation orientA, ElOrientation orientB, 
  complex_double alpha, ElConstMatrix_z A, ElConstMatrix_z B, 
  complex_double beta,                     ElMatrix_z C );

ElError ElTrrkDist_s
( ElUpperOrLower uplo, ElOrientation orientA, ElOrientation orientB, 
  float alpha, ElConstDistMatrix_s A, ElConstDistMatrix_s B, 
  float beta,                         ElDistMatrix_s C );
ElError ElTrrkDist_d
( ElUpperOrLower uplo, ElOrientation orientA, ElOrientation orientB, 
  double alpha, ElConstDistMatrix_d A, ElConstDistMatrix_d B, 
  double beta,                         ElDistMatrix_d C );
ElError ElTrrkDist_c
( ElUpperOrLower uplo, ElOrientation orientA, ElOrientation orientB, 
  complex_float alpha, ElConstDistMatrix_c A, ElConstDistMatrix_c B, 
  complex_float beta,                         ElDistMatrix_c C );
ElError ElTrrkDist_z
( ElUpperOrLower uplo, ElOrientation orientA, ElOrientation orientB, 
  complex_double alpha, ElConstDistMatrix_z A, ElConstDistMatrix_z B, 
  complex_double beta,                         ElDistMatrix_z C );

/* Trr2k
   ===== */
/* TRiangular Rank-2K update */
ElError ElTrr2k_s
( ElUpperOrLower uplo, 
  ElOrientation orientA, ElOrientation orientB, 
  ElOrientation orientC, ElOrientation orientD,
  float alpha, ElConstMatrix_s A, ElConstMatrix_s B, 
               ElConstMatrix_s C, ElConstMatrix_s D,
  float beta,                     ElMatrix_s E );
ElError ElTrr2k_d
( ElUpperOrLower uplo, 
  ElOrientation orientA, ElOrientation orientB, 
  ElOrientation orientC, ElOrientation orientD, 
  double alpha, ElConstMatrix_d A, ElConstMatrix_d B, 
                ElConstMatrix_d C, ElConstMatrix_d D,
  double beta,                     ElMatrix_d E );
ElError ElTrr2k_c
( ElUpperOrLower uplo, 
  ElOrientation orientA, ElOrientation orientB, 
  ElOrientation orientC, ElOrientation orientD, 
  complex_float alpha, ElConstMatrix_c A, ElConstMatrix_c B, 
                       ElConstMatrix_c C, ElConstMatrix_c D,
  complex_float beta,                     ElMatrix_c E );
ElError ElTrr2k_z
( ElUpperOrLower uplo,
  ElOrientation orientA, ElOrientation orientB, 
  ElOrientation orientC, ElOrientation orientD, 
  complex_double alpha, ElConstMatrix_z A, ElConstMatrix_z B, 
                        ElConstMatrix_z C, ElConstMatrix_z D,
  complex_double beta,                     ElMatrix_z E );

ElError ElTrr2kDist_s
( ElUpperOrLower uplo, 
  ElOrientation orientA, ElOrientation orientB, 
  ElOrientation orientC, ElOrientation orientD,
  float alpha, ElConstDistMatrix_s A, ElConstDistMatrix_s B, 
               ElConstDistMatrix_s C, ElConstDistMatrix_s D,
  float beta,                         ElDistMatrix_s E );
ElError ElTrr2kDist_d
( ElUpperOrLower uplo, 
  ElOrientation orientA, ElOrientation orientB, 
  ElOrientation orientC, ElOrientation orientD, 
  double alpha, ElConstDistMatrix_d A, ElConstDistMatrix_d B, 
                ElConstDistMatrix_d C, ElConstDistMatrix_d D,
  double beta,                         ElDistMatrix_d E );
ElError ElTrr2kDist_c
( ElUpperOrLower uplo, 
  ElOrientation orientA, ElOrientation orientB, 
  ElOrientation orientC, ElOrientation orientD, 
  complex_float alpha, ElConstDistMatrix_c A, ElConstDistMatrix_c B, 
                       ElConstDistMatrix_c C, ElConstDistMatrix_c D,
  complex_float beta,                         ElDistMatrix_c E );
ElError ElTrr2kDist_z
( ElUpperOrLower uplo,
  ElOrientation orientA, ElOrientation orientB, 
  ElOrientation orientC, ElOrientation orientD, 
  complex_double alpha, ElConstDistMatrix_z A, ElConstDistMatrix_z B, 
                        ElConstDistMatrix_z C, ElConstDistMatrix_z D,
  complex_double beta,                         ElDistMatrix_z E );

/* Trsm
   ==== */
ElError ElTrsm_s
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  float alpha, ElConstMatrix_s A, ElMatrix_s B );
ElError ElTrsm_d
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  double alpha, ElConstMatrix_d A, ElMatrix_d B );
ElError ElTrsm_c
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  complex_float alpha, ElConstMatrix_c A, ElMatrix_c B );
ElError ElTrsm_z
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  complex_double alpha, ElConstMatrix_z A, ElMatrix_z B );

ElError ElTrsmDist_s
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  float alpha, ElConstDistMatrix_s A, ElDistMatrix_s B );
ElError ElTrsmDist_d
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  double alpha, ElConstDistMatrix_d A, ElDistMatrix_d B );
ElError ElTrsmDist_c
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  complex_float alpha, ElConstDistMatrix_c A, ElDistMatrix_c B );
ElError ElTrsmDist_z
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  complex_double alpha, ElConstDistMatrix_z A, ElDistMatrix_z B );

/* Trstrm
   ====== */
ElError ElTrstrm_s
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  float alpha, ElConstMatrix_s A, ElMatrix_s B );
ElError ElTrstrm_d
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  double alpha, ElConstMatrix_d A, ElMatrix_d B );
ElError ElTrstrm_c
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  complex_float alpha, ElConstMatrix_c A, ElMatrix_c B );
ElError ElTrstrm_z
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  complex_double alpha, ElConstMatrix_z A, ElMatrix_z B );

ElError ElTrstrmDist_s
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  float alpha, ElConstDistMatrix_s A, ElDistMatrix_s B );
ElError ElTrstrmDist_d
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  double alpha, ElConstDistMatrix_d A, ElDistMatrix_d B );
ElError ElTrstrmDist_c
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  complex_float alpha, ElConstDistMatrix_c A, ElDistMatrix_c B );
ElError ElTrstrmDist_z
( ElLeftOrRight side, ElUpperOrLower uplo, 
  ElOrientation orientation, ElUnitOrNonUnit diag,
  complex_double alpha, ElConstDistMatrix_z A, ElDistMatrix_z B );

/* Trtrmm
   ====== */
ElError ElTrtrmm_s( ElUpperOrLower uplo, ElMatrix_s A, bool conjugate );
ElError ElTrtrmm_d( ElUpperOrLower uplo, ElMatrix_d A, bool conjugate );
ElError ElTrtrmm_c( ElUpperOrLower uplo, ElMatrix_c A, bool conjugate );
ElError ElTrtrmm_z( ElUpperOrLower uplo, ElMatrix_z A, bool conjugate );

ElError ElTrtrmmDist_s( ElUpperOrLower uplo, ElDistMatrix_s A, bool conjugate );
ElError ElTrtrmmDist_d( ElUpperOrLower uplo, ElDistMatrix_d A, bool conjugate );
ElError ElTrtrmmDist_c( ElUpperOrLower uplo, ElDistMatrix_c A, bool conjugate );
ElError ElTrtrmmDist_z( ElUpperOrLower uplo, ElDistMatrix_z A, bool conjugate );

/* TwoSidedTrmm
   ============ */
ElError ElTwoSidedTrmm_s
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElMatrix_s A, ElConstMatrix_s B );
ElError ElTwoSidedTrmm_d
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElMatrix_d A, ElConstMatrix_d B );
ElError ElTwoSidedTrmm_c
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElMatrix_c A, ElConstMatrix_c B );
ElError ElTwoSidedTrmm_z
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElMatrix_z A, ElConstMatrix_z B );

ElError ElTwoSidedTrmmDist_s
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, 
  ElDistMatrix_s A, ElConstDistMatrix_s B );
ElError ElTwoSidedTrmmDist_d
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, 
  ElDistMatrix_d A, ElConstDistMatrix_d B );
ElError ElTwoSidedTrmmDist_c
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, 
  ElDistMatrix_c A, ElConstDistMatrix_c B );
ElError ElTwoSidedTrmmDist_z
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, 
  ElDistMatrix_z A, ElConstDistMatrix_z B );

/* TwoSidedTrsm
   ============ */
ElError ElTwoSidedTrsm_s
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElMatrix_s A, ElConstMatrix_s B );
ElError ElTwoSidedTrsm_d
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElMatrix_d A, ElConstMatrix_d B );
ElError ElTwoSidedTrsm_c
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElMatrix_c A, ElConstMatrix_c B );
ElError ElTwoSidedTrsm_z
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, ElMatrix_z A, ElConstMatrix_z B );

ElError ElTwoSidedTrsmDist_s
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, 
  ElDistMatrix_s A, ElConstDistMatrix_s B );
ElError ElTwoSidedTrsmDist_d
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, 
  ElDistMatrix_d A, ElConstDistMatrix_d B );
ElError ElTwoSidedTrsmDist_c
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, 
  ElDistMatrix_c A, ElConstDistMatrix_c B );
ElError ElTwoSidedTrsmDist_z
( ElUpperOrLower uplo, ElUnitOrNonUnit diag, 
  ElDistMatrix_z A, ElConstDistMatrix_z B );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_BLAS_LEVEL3_C_H */
