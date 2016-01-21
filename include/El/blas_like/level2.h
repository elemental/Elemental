/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_LEVEL2_C_H
#define EL_BLAS_LEVEL2_C_H

#include "El/core/DistMatrix.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Gemv
   ==== */
EL_EXPORT ElError ElGemv_i
( ElOrientation orientation, 
  ElInt alpha, ElConstMatrix_i A, ElConstMatrix_i x, 
  ElInt beta,                     ElMatrix_i y );
EL_EXPORT ElError ElGemv_s
( ElOrientation orientation, 
  float alpha, ElConstMatrix_s A, ElConstMatrix_s x, 
  float beta,                     ElMatrix_s y );
EL_EXPORT ElError ElGemv_d
( ElOrientation orientation, 
  double alpha, ElConstMatrix_d A, ElConstMatrix_d x, 
  double beta,                     ElMatrix_d y );
EL_EXPORT ElError ElGemv_c
( ElOrientation orientation, 
  complex_float alpha, ElConstMatrix_c A, ElConstMatrix_c x, 
  complex_float beta,                     ElMatrix_c y );
EL_EXPORT ElError ElGemv_z
( ElOrientation orientation, 
  complex_double alpha, ElConstMatrix_z A, ElConstMatrix_z x, 
  complex_double beta,                     ElMatrix_z y );

EL_EXPORT ElError ElGemvDist_i
( ElOrientation orientation, 
  ElInt alpha, ElConstDistMatrix_i A, ElConstDistMatrix_i x, 
  ElInt beta,                         ElDistMatrix_i y );
EL_EXPORT ElError ElGemvDist_s
( ElOrientation orientation, 
  float alpha, ElConstDistMatrix_s A, ElConstDistMatrix_s x, 
  float beta,                         ElDistMatrix_s y );
EL_EXPORT ElError ElGemvDist_d
( ElOrientation orientation, 
  double alpha, ElConstDistMatrix_d A, ElConstDistMatrix_d x, 
  double beta,                         ElDistMatrix_d y );
EL_EXPORT ElError ElGemvDist_c
( ElOrientation orientation, 
  complex_float alpha, ElConstDistMatrix_c A, ElConstDistMatrix_c x, 
  complex_float beta,                         ElDistMatrix_c y );
EL_EXPORT ElError ElGemvDist_z
( ElOrientation orientation, 
  complex_double alpha, ElConstDistMatrix_z A, ElConstDistMatrix_z x, 
  complex_double beta,                         ElDistMatrix_z y );

/* Ger
   === */
EL_EXPORT ElError ElGer_i
( ElInt alpha, ElConstMatrix_i x, ElConstMatrix_i y, ElMatrix_i A );
EL_EXPORT ElError ElGer_s
( float alpha, ElConstMatrix_s x, ElConstMatrix_s y, ElMatrix_s A );
EL_EXPORT ElError ElGer_d
( double alpha, ElConstMatrix_d x, ElConstMatrix_d y, ElMatrix_d A );
EL_EXPORT ElError ElGer_c
( complex_float alpha, ElConstMatrix_c x, ElConstMatrix_c y, ElMatrix_c A );
EL_EXPORT ElError ElGer_z
( complex_double alpha, ElConstMatrix_z x, ElConstMatrix_z y, ElMatrix_z A );

EL_EXPORT ElError ElGerDist_i
( ElInt alpha, ElConstDistMatrix_i x, ElConstDistMatrix_i y, 
               ElDistMatrix_i A );
EL_EXPORT ElError ElGerDist_s
( float alpha, ElConstDistMatrix_s x, ElConstDistMatrix_s y, 
               ElDistMatrix_s A );
EL_EXPORT ElError ElGerDist_d
( double alpha, ElConstDistMatrix_d x, ElConstDistMatrix_d y, 
                ElDistMatrix_d A );
EL_EXPORT ElError ElGerDist_c
( complex_float alpha, ElConstDistMatrix_c x, ElConstDistMatrix_c y, 
                       ElDistMatrix_c A );
EL_EXPORT ElError ElGerDist_z
( complex_double alpha, ElConstDistMatrix_z x, ElConstDistMatrix_z y, 
                        ElDistMatrix_z A );

/* Geru
   ==== */
EL_EXPORT ElError ElGeru_c
( complex_float alpha, ElConstMatrix_c x, ElConstMatrix_c y, ElMatrix_c A );
EL_EXPORT ElError ElGeru_z
( complex_double alpha, ElConstMatrix_z x, ElConstMatrix_z y, ElMatrix_z A );

EL_EXPORT ElError ElGeruDist_c
( complex_float alpha, ElConstDistMatrix_c x, ElConstDistMatrix_c y, 
                       ElDistMatrix_c A );
EL_EXPORT ElError ElGeruDist_z
( complex_double alpha, ElConstDistMatrix_z x, ElConstDistMatrix_z y, 
                        ElDistMatrix_z A );

/* Hemv
   ==== */
EL_EXPORT ElError ElHemv_c
( ElUpperOrLower uplo, 
  complex_float alpha, ElConstMatrix_c A, ElConstMatrix_c x, 
  complex_float beta,                     ElMatrix_c y );
EL_EXPORT ElError ElHemv_z
( ElUpperOrLower uplo, 
  complex_double alpha, ElConstMatrix_z A, ElConstMatrix_z x, 
  complex_double beta,                     ElMatrix_z y );

/* NOTE: 'A', 'x', and 'y' must be in [MC,MR] distributions */
EL_EXPORT ElError ElHemvDist_c
( ElUpperOrLower uplo, 
  complex_float alpha, ElConstDistMatrix_c A, ElConstDistMatrix_c x, 
  complex_float beta,                         ElDistMatrix_c y );
EL_EXPORT ElError ElHemvDist_z
( ElUpperOrLower uplo, 
  complex_double alpha, ElConstDistMatrix_z A, ElConstDistMatrix_z x, 
  complex_double beta,                         ElDistMatrix_z y );

/* Her 
   === */
EL_EXPORT ElError ElHer_c
( ElUpperOrLower uplo, float alpha, ElConstMatrix_c x, ElMatrix_c A );
EL_EXPORT ElError ElHer_z
( ElUpperOrLower uplo, double alpha, ElConstMatrix_z x, ElMatrix_z A );

/* NOTE: 'A' and 'x' must be in [MC,MR] distributions */
EL_EXPORT ElError ElHerDist_c
( ElUpperOrLower uplo, float alpha, 
  ElConstDistMatrix_c x, ElDistMatrix_c A );
EL_EXPORT ElError ElHerDist_z
( ElUpperOrLower uplo, double alpha, 
  ElConstDistMatrix_z x, ElDistMatrix_z A );

/* Her2 
   ==== */
EL_EXPORT ElError ElHer2_c
( ElUpperOrLower uplo, complex_float alpha, 
  ElConstMatrix_c x, ElConstMatrix_c y, ElMatrix_c A );
EL_EXPORT ElError ElHer2_z
( ElUpperOrLower uplo, complex_double alpha, 
  ElConstMatrix_z x, ElConstMatrix_z y, ElMatrix_z A );

/* NOTE: 'A', 'x', and 'y' must be in [MC,MR] distributions */
EL_EXPORT ElError ElHer2Dist_c
( ElUpperOrLower uplo, complex_float alpha, 
  ElConstDistMatrix_c x, ElConstDistMatrix_c y, ElDistMatrix_c A );
EL_EXPORT ElError ElHer2Dist_z
( ElUpperOrLower uplo, complex_double alpha, 
  ElConstDistMatrix_z x, ElConstDistMatrix_z y, ElDistMatrix_z A );

/* QuasiTrsv
   ========= */
EL_EXPORT ElError ElQuasiTrsv_s
( ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_s A, ElMatrix_s x );
EL_EXPORT ElError ElQuasiTrsv_d
( ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_d A, ElMatrix_d x );
EL_EXPORT ElError ElQuasiTrsv_c
( ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstMatrix_c A, ElMatrix_c x );
EL_EXPORT ElError ElQuasiTrsv_z
( ElUpperOrLower uplo, ElOrientation orientation,
  ElConstMatrix_z A, ElMatrix_z x );

EL_EXPORT ElError ElQuasiTrsvDist_s
( ElUpperOrLower uplo, ElOrientation orientation,
  ElConstDistMatrix_s A, ElDistMatrix_s x );
EL_EXPORT ElError ElQuasiTrsvDist_d
( ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_d A, ElDistMatrix_d x );
EL_EXPORT ElError ElQuasiTrsvDist_c
( ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_c A, ElDistMatrix_c x );
EL_EXPORT ElError ElQuasiTrsvDist_z
( ElUpperOrLower uplo, ElOrientation orientation, 
  ElConstDistMatrix_z A, ElDistMatrix_z x );

/* Symv
   ==== */
typedef struct {
  ElInt bsize;
  bool avoidTrmvBasedLocalSymv;
} ElSymvCtrl;
EL_EXPORT ElError ElSymvCtrlDefault_s( ElSymvCtrl* ctrl );
EL_EXPORT ElError ElSymvCtrlDefault_d( ElSymvCtrl* ctrl );
EL_EXPORT ElError ElSymvCtrlDefault_c( ElSymvCtrl* ctrl );
EL_EXPORT ElError ElSymvCtrlDefault_z( ElSymvCtrl* ctrl );

EL_EXPORT ElError ElSymv_i
( ElUpperOrLower uplo, 
  ElInt alpha, ElConstMatrix_i A, ElConstMatrix_i x, 
  ElInt beta,                     ElMatrix_i y );
EL_EXPORT ElError ElSymv_s
( ElUpperOrLower uplo, 
  float alpha, ElConstMatrix_s A, ElConstMatrix_s x, 
  float beta,                     ElMatrix_s y );
EL_EXPORT ElError ElSymv_d
( ElUpperOrLower uplo, 
  double alpha, ElConstMatrix_d A, ElConstMatrix_d x, 
  double beta,                     ElMatrix_d y );
EL_EXPORT ElError ElSymv_c
( ElUpperOrLower uplo, 
  complex_float alpha, ElConstMatrix_c A, ElConstMatrix_c x, 
  complex_float beta,                     ElMatrix_c y );
EL_EXPORT ElError ElSymv_z
( ElUpperOrLower uplo, 
  complex_double alpha, ElConstMatrix_z A, ElConstMatrix_z x, 
  complex_double beta,                     ElMatrix_z y );

EL_EXPORT ElError ElSymvDist_i
( ElUpperOrLower uplo, 
  ElInt alpha, ElConstDistMatrix_i A, ElConstDistMatrix_i x, 
  ElInt beta,                         ElDistMatrix_i y );
EL_EXPORT ElError ElSymvDist_s
( ElUpperOrLower uplo, 
  float alpha, ElConstDistMatrix_s A, ElConstDistMatrix_s x, 
  float beta,                         ElDistMatrix_s y );
EL_EXPORT ElError ElSymvDist_d
( ElUpperOrLower uplo, 
  double alpha, ElConstDistMatrix_d A, ElConstDistMatrix_d x, 
  double beta,                         ElDistMatrix_d y );
EL_EXPORT ElError ElSymvDist_c
( ElUpperOrLower uplo, 
  complex_float alpha, ElConstDistMatrix_c A, ElConstDistMatrix_c x, 
  complex_float beta,                         ElDistMatrix_c y );
EL_EXPORT ElError ElSymvDist_z
( ElUpperOrLower uplo, 
  complex_double alpha, ElConstDistMatrix_z A, ElConstDistMatrix_z x, 
  complex_double beta,                         ElDistMatrix_z y );

/* Syr 
   === */
EL_EXPORT ElError ElSyr_s
( ElUpperOrLower uplo, float alpha, ElConstMatrix_s x, ElMatrix_s A );
EL_EXPORT ElError ElSyr_d
( ElUpperOrLower uplo, double alpha, ElConstMatrix_d x, ElMatrix_d A );
EL_EXPORT ElError ElSyr_c
( ElUpperOrLower uplo, complex_float alpha, ElConstMatrix_c x, ElMatrix_c A );
EL_EXPORT ElError ElSyr_z
( ElUpperOrLower uplo, complex_double alpha, ElConstMatrix_z x, ElMatrix_z A );

EL_EXPORT ElError ElSyrDist_s
( ElUpperOrLower uplo, float alpha, 
  ElConstDistMatrix_s x, ElDistMatrix_s A );
EL_EXPORT ElError ElSyrDist_d
( ElUpperOrLower uplo, double alpha, 
  ElConstDistMatrix_d x, ElDistMatrix_d A );
EL_EXPORT ElError ElSyrDist_c
( ElUpperOrLower uplo, complex_float alpha, 
  ElConstDistMatrix_c x, ElDistMatrix_c A );
EL_EXPORT ElError ElSyrDist_z
( ElUpperOrLower uplo, complex_double alpha, 
  ElConstDistMatrix_z x, ElDistMatrix_z A );

/* Syr2 
   ==== */
EL_EXPORT ElError ElSyr2_s
( ElUpperOrLower uplo, float alpha, 
  ElConstMatrix_s x, ElConstMatrix_s y, ElMatrix_s A );
EL_EXPORT ElError ElSyr2_d
( ElUpperOrLower uplo, double alpha, 
  ElConstMatrix_d x, ElConstMatrix_d y, ElMatrix_d A );
EL_EXPORT ElError ElSyr2_c
( ElUpperOrLower uplo, complex_float alpha, 
  ElConstMatrix_c x, ElConstMatrix_c y, ElMatrix_c A );
EL_EXPORT ElError ElSyr2_z
( ElUpperOrLower uplo, complex_double alpha, 
  ElConstMatrix_z x, ElConstMatrix_z y, ElMatrix_z A );

EL_EXPORT ElError ElSyr2Dist_s
( ElUpperOrLower uplo, float alpha, 
  ElConstDistMatrix_s x, ElConstDistMatrix_s y, ElDistMatrix_s A );
EL_EXPORT ElError ElSyr2Dist_d
( ElUpperOrLower uplo, double alpha, 
  ElConstDistMatrix_d x, ElConstDistMatrix_d y, ElDistMatrix_d A );
EL_EXPORT ElError ElSyr2Dist_c
( ElUpperOrLower uplo, complex_float alpha, 
  ElConstDistMatrix_c x, ElConstDistMatrix_c y, ElDistMatrix_c A );
EL_EXPORT ElError ElSyr2Dist_z
( ElUpperOrLower uplo, complex_double alpha, 
  ElConstDistMatrix_z x, ElConstDistMatrix_z y, ElDistMatrix_z A );

/* Trmv
   ==== */
EL_EXPORT ElError ElTrmv_s
( ElUpperOrLower uplo, ElOrientation orientation, ElUnitOrNonUnit diag, 
  ElConstMatrix_s A, ElMatrix_s x );
EL_EXPORT ElError ElTrmv_d
( ElUpperOrLower uplo, ElOrientation orientation, ElUnitOrNonUnit diag, 
  ElConstMatrix_d A, ElMatrix_d x );
EL_EXPORT ElError ElTrmv_c
( ElUpperOrLower uplo, ElOrientation orientation, ElUnitOrNonUnit diag, 
  ElConstMatrix_c A, ElMatrix_c x );
EL_EXPORT ElError ElTrmv_z
( ElUpperOrLower uplo, ElOrientation orientation, ElUnitOrNonUnit diag, 
  ElConstMatrix_z A, ElMatrix_z x );
/* NOTE: Distributed Trmv not implemented, use Trmm */

/* Trr
   === */
EL_EXPORT ElError ElTrr_i
( ElUpperOrLower uplo, ElInt alpha, 
  ElConstMatrix_i x, ElConstMatrix_i y, ElMatrix_i A );
EL_EXPORT ElError ElTrr_s
( ElUpperOrLower uplo, float alpha, 
  ElConstMatrix_s x, ElConstMatrix_s y, ElMatrix_s A );
EL_EXPORT ElError ElTrr_d
( ElUpperOrLower uplo, double alpha, 
  ElConstMatrix_d x, ElConstMatrix_d y, ElMatrix_d A );
EL_EXPORT ElError ElTrr_c
( ElUpperOrLower uplo, complex_float alpha, 
  ElConstMatrix_c x, ElConstMatrix_c y, ElMatrix_c A, bool conjugate );
EL_EXPORT ElError ElTrr_z
( ElUpperOrLower uplo, complex_double alpha, 
  ElConstMatrix_z x, ElConstMatrix_z y, ElMatrix_z A, bool conjugate );

EL_EXPORT ElError ElTrrDist_i
( ElUpperOrLower uplo, ElInt alpha, 
  ElConstDistMatrix_i x, ElConstDistMatrix_i y, ElDistMatrix_i A );
EL_EXPORT ElError ElTrrDist_s
( ElUpperOrLower uplo, float alpha, 
  ElConstDistMatrix_s x, ElConstDistMatrix_s y, ElDistMatrix_s A );
EL_EXPORT ElError ElTrrDist_d
( ElUpperOrLower uplo, double alpha, 
  ElConstDistMatrix_d x, ElConstDistMatrix_d y, ElDistMatrix_d A );
EL_EXPORT ElError ElTrrDist_c
( ElUpperOrLower uplo, complex_float alpha, 
  ElConstDistMatrix_c x, ElConstDistMatrix_c y, ElDistMatrix_c A, 
  bool conjugate );
EL_EXPORT ElError ElTrrDist_z
( ElUpperOrLower uplo, complex_double alpha, 
  ElConstDistMatrix_z x, ElConstDistMatrix_z y, ElDistMatrix_z A, 
  bool conjugate );

/* Trr2
   ==== */
EL_EXPORT ElError ElTrr2_i
( ElUpperOrLower uplo, ElInt alpha, 
  ElConstMatrix_i X, ElConstMatrix_i Y, ElMatrix_i A );
EL_EXPORT ElError ElTrr2_s
( ElUpperOrLower uplo, float alpha, 
  ElConstMatrix_s X, ElConstMatrix_s Y, ElMatrix_s A );
EL_EXPORT ElError ElTrr2_d
( ElUpperOrLower uplo, double alpha, 
  ElConstMatrix_d X, ElConstMatrix_d Y, ElMatrix_d A );
EL_EXPORT ElError ElTrr2_c
( ElUpperOrLower uplo, complex_float alpha, 
  ElConstMatrix_c X, ElConstMatrix_c Y, ElMatrix_c A, bool conjugate );
EL_EXPORT ElError ElTrr2_z
( ElUpperOrLower uplo, complex_double alpha, 
  ElConstMatrix_z X, ElConstMatrix_z Y, ElMatrix_z A, bool conjugate );

EL_EXPORT ElError ElTrr2Dist_i
( ElUpperOrLower uplo, ElInt alpha, 
  ElConstDistMatrix_i X, ElConstDistMatrix_i Y, ElDistMatrix_i A );
EL_EXPORT ElError ElTrr2Dist_s
( ElUpperOrLower uplo, float alpha, 
  ElConstDistMatrix_s X, ElConstDistMatrix_s Y, ElDistMatrix_s A );
EL_EXPORT ElError ElTrr2Dist_d
( ElUpperOrLower uplo, double alpha, 
  ElConstDistMatrix_d X, ElConstDistMatrix_d Y, ElDistMatrix_d A );
EL_EXPORT ElError ElTrr2Dist_c
( ElUpperOrLower uplo, complex_float alpha, 
  ElConstDistMatrix_c X, ElConstDistMatrix_c Y, ElDistMatrix_c A, 
  bool conjugate );
EL_EXPORT ElError ElTrr2Dist_z
( ElUpperOrLower uplo, complex_double alpha, 
  ElConstDistMatrix_z X, ElConstDistMatrix_z Y, ElDistMatrix_z A, 
  bool conjugate );

/* Trsv
   ==== */
EL_EXPORT ElError ElTrsv_s
( ElUpperOrLower uplo, ElOrientation orientation, ElUnitOrNonUnit diag, 
  ElConstMatrix_s A, ElMatrix_s x );
EL_EXPORT ElError ElTrsv_d
( ElUpperOrLower uplo, ElOrientation orientation, ElUnitOrNonUnit diag, 
  ElConstMatrix_d A, ElMatrix_d x );
EL_EXPORT ElError ElTrsv_c
( ElUpperOrLower uplo, ElOrientation orientation, ElUnitOrNonUnit diag, 
  ElConstMatrix_c A, ElMatrix_c x );
EL_EXPORT ElError ElTrsv_z
( ElUpperOrLower uplo, ElOrientation orientation, ElUnitOrNonUnit diag, 
  ElConstMatrix_z A, ElMatrix_z x );

EL_EXPORT ElError ElTrsvDist_s
( ElUpperOrLower uplo, ElOrientation orientation, ElUnitOrNonUnit diag, 
  ElConstDistMatrix_s A, ElDistMatrix_s x );
EL_EXPORT ElError ElTrsvDist_d
( ElUpperOrLower uplo, ElOrientation orientation, ElUnitOrNonUnit diag, 
  ElConstDistMatrix_d A, ElDistMatrix_d x );
EL_EXPORT ElError ElTrsvDist_c
( ElUpperOrLower uplo, ElOrientation orientation, ElUnitOrNonUnit diag, 
  ElConstDistMatrix_c A, ElDistMatrix_c x );
EL_EXPORT ElError ElTrsvDist_z
( ElUpperOrLower uplo, ElOrientation orientation, ElUnitOrNonUnit diag, 
  ElConstDistMatrix_z A, ElDistMatrix_z x );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_BLAS_LEVEL2_C_H */
