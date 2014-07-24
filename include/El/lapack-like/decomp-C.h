/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_LAPACK_DECOMP_C_H
#define EL_LAPACK_DECOMP_C_H

#include "El/lapack-like/funcs-C.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  EL_AXBX=1,
  EL_ABX=2, 
  EL_BAX=3
} ElHermitianGenDefiniteEigType;

typedef struct {
  ElInt cutoff;
  ElInt maxInnerIts, maxOuterIts;
  float tol;
  float spreadFactor;
  bool progress;
} ElHermitianSdcCtrl_s;
ElError ElHermitianSdcCtrlDefault_s( ElHermitianSdcCtrl_s* ctrl );

typedef struct {
  ElInt cutoff;
  ElInt maxInnerIts, maxOuterIts;
  double tol;
  double spreadFactor;
  bool progress;
} ElHermitianSdcCtrl_d;
ElError ElHermitianSdcCtrlDefault_d( ElHermitianSdcCtrl_d* ctrl );

typedef struct {
  ElHermitianTridiagCtrl tridiagCtrl;
  ElHermitianSdcCtrl_s sdcCtrl;
  bool useSdc;
} ElHermitianEigCtrl_s;
ElError ElHermitianEigCtrlDefault_s( ElHermitianEigCtrl_s* ctrl );

typedef struct {
  ElHermitianTridiagCtrl tridiagCtrl;
  ElHermitianSdcCtrl_d sdcCtrl;
  bool useSdc;
} ElHermitianEigCtrl_d;
ElError ElHermitianEigCtrlDefault_d( ElHermitianEigCtrl_d* ctrl );

typedef struct {
  bool qdwh;
  bool colPiv;
  ElInt maxIts;
  ElInt numIts;
} ElPolarCtrl;
ElError ElPolarCtrlDefault( ElPolarCtrl* ctrl );

typedef struct {
  bool aed;
  ElInt blockHeight, blockWidth;
} ElHessQrCtrl;
ElError ElHessQrCtrlDefault( ElHessQrCtrl* ctrl );

typedef struct {
  ElInt cutoff;
  ElInt maxInnerIts, maxOuterIts;
  float tol;
  float spreadFactor;
  bool random;
  bool progress;
  ElSignCtrl_s signCtrl;
} ElSdcCtrl_s;
ElError ElSdcCtrlDefault_s( ElSdcCtrl_s* ctrl );

typedef struct {
  ElInt cutoff;
  ElInt maxInnerIts, maxOuterIts;
  double tol;
  double spreadFactor;
  bool random;
  bool progress;
  ElSignCtrl_d signCtrl;
} ElSdcCtrl_d;
ElError ElSdcCtrlDefault_d( ElSdcCtrl_d* ctrl );

typedef struct {
  bool useSdc;
  ElHessQrCtrl qrCtrl;
  ElSdcCtrl_s sdcCtrl;
} ElSchurCtrl_s;
ElError ElSchurCtrlDefault_s( ElSchurCtrl_s* ctrl );

typedef struct {
  bool useSdc;
  ElHessQrCtrl qrCtrl;
  ElSdcCtrl_d sdcCtrl;
} ElSchurCtrl_d;
ElError ElSchurCtrlDefault_d( ElSchurCtrl_d* ctrl );

/* Hermitian eigenvalue solvers
   ============================ */

/* Compute all of the eigenvalues
   ------------------------------ */
ElError ElHermitianEig_s
( ElUpperOrLower uplo, ElMatrix_s A, ElMatrix_s w, ElSortType sortType );
ElError ElHermitianEig_d
( ElUpperOrLower uplo, ElMatrix_d A, ElMatrix_d w, ElSortType sortType );
ElError ElHermitianEig_c
( ElUpperOrLower uplo, ElMatrix_c A, ElMatrix_s w, ElSortType sortType );
ElError ElHermitianEig_z
( ElUpperOrLower uplo, ElMatrix_z A, ElMatrix_d w, ElSortType sortType );

/* NOTE: 'A' must be in a [MC,MR] distribution, while 
         'w' must be in a [VR,STAR] distribution */
/*
ElError ElHermitianEigDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A, ElDistMatrix_s w, 
  ElSortType sortType );
*/
ElError ElHermitianEigDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A, ElDistMatrix_d w,
  ElSortType sortType );
/*
ElError ElHermitianEigDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A, ElDistMatrix_s w,
  ElSortType sortType );
*/
ElError ElHermitianEigDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A, ElDistMatrix_d w,
  ElSortType sortType );

/* TODO: Expert version */

/* Compute the entire eigenvalue decomposition 
   ------------------------------------------- */
ElError ElHermitianEigPair_s
( ElUpperOrLower uplo, ElMatrix_s A, ElMatrix_s w, ElMatrix_s Z,
  ElSortType sortType );
ElError ElHermitianEigPair_d
( ElUpperOrLower uplo, ElMatrix_d A, ElMatrix_d w, ElMatrix_d Z,
  ElSortType sortType );
ElError ElHermitianEigPair_c
( ElUpperOrLower uplo, ElMatrix_c A, ElMatrix_s w, ElMatrix_c Z,
  ElSortType sortType );
ElError ElHermitianEigPair_z
( ElUpperOrLower uplo, ElMatrix_z A, ElMatrix_d w, ElMatrix_z Z,
  ElSortType sortType );

/* NOTE: 'A' and 'Z' must be in a [MC,MR] distribution, while 
         'w' must be in a [VR,STAR] distribution */
/*
ElError ElHermitianEigPairDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A, ElDistMatrix_s w, ElDistMatrix_s Z,
  ElSortType sortType );
*/
ElError ElHermitianEigPairDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A, ElDistMatrix_d w, ElDistMatrix_d Z,
  ElSortType sortType );
/*
ElError ElHermitianEigPairDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A, ElDistMatrix_s w, ElDistMatrix_c Z,
  ElSortType sortType );
*/
ElError ElHermitianEigPairDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A, ElDistMatrix_d w, ElDistMatrix_z Z,
  ElSortType sortType );

/* TODO: Expert version */

/* Compute the eigenvalues within a half-open interval
   --------------------------------------------------- */
ElError ElHermitianEigInterval_s
( ElUpperOrLower uplo, ElMatrix_s A, ElMatrix_s w, float lower, float upper,
  ElSortType sortType );
ElError ElHermitianEigInterval_d
( ElUpperOrLower uplo, ElMatrix_d A, ElMatrix_d w, double lower, double upper,
  ElSortType sortType );
ElError ElHermitianEigInterval_c
( ElUpperOrLower uplo, ElMatrix_c A, ElMatrix_s w, float lower, float upper,
  ElSortType sortType );
ElError ElHermitianEigInterval_z
( ElUpperOrLower uplo, ElMatrix_z A, ElMatrix_d w, double lower, double upper,
  ElSortType sortType );

/* NOTE: 'A' must be in a [MC,MR] distribution, while 
         'w' must be in a [VR,STAR] distribution */
/*
ElError ElHermitianEigIntervalDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A, ElDistMatrix_s w, 
  float lower, float upper );
*/
ElError ElHermitianEigIntervalDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A, ElDistMatrix_d w,
  double lower, double upper, ElSortType sortType );
/*
ElError ElHermitianEigIntervalDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A, ElDistMatrix_s w,
  float lower, float upper, ElSortType sortType );
*/
ElError ElHermitianEigIntervalDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A, ElDistMatrix_d w,
  double lower, double upper, ElSortType sortType );

/* TODO: Expert version */

/* Compute the eigenpairs within a half-open interval
   -------------------------------------------------- */
ElError ElHermitianEigPairInterval_s
( ElUpperOrLower uplo, ElMatrix_s A, ElMatrix_s w, ElMatrix_s Z, 
  float lower, float upper, ElSortType sortType );
ElError ElHermitianEigPairInterval_d
( ElUpperOrLower uplo, ElMatrix_d A, ElMatrix_d w, ElMatrix_d Z, 
  double lower, double upper, ElSortType sortType );
ElError ElHermitianEigPairInterval_c
( ElUpperOrLower uplo, ElMatrix_c A, ElMatrix_s w, ElMatrix_c Z,
  float lower, float upper, ElSortType sortType );
ElError ElHermitianEigPairInterval_z
( ElUpperOrLower uplo, ElMatrix_z A, ElMatrix_d w, ElMatrix_z Z,
  double lower, double upper, ElSortType sortType );

/* NOTE: 'A' and 'Z' must be in a [MC,MR] distribution, while 
         'w' must be in a [VR,STAR] distribution */
/*
ElError ElHermitianEigPairIntervalDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A, ElDistMatrix_s w, ElDistMatrix_s Z,
  float lower, float upper, ElSortType sortType );
*/
ElError ElHermitianEigPairIntervalDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A, ElDistMatrix_d w, ElDistMatrix_d Z,
  double lower, double upper, ElSortType sortType );
/*
ElError ElHermitianEigPairIntervalDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A, ElDistMatrix_s w, ElDistMatrix_c Z,
  float lower, float upper, ElSortType sortType );
*/
ElError ElHermitianEigPairIntervalDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A, ElDistMatrix_d w, ElDistMatrix_z Z,
  double lower, double upper, ElSortType sortType );

/* TODO: Expert version */

/* Compute the eigenvalues within an index range
   --------------------------------------------- */
ElError ElHermitianEigIndices_s
( ElUpperOrLower uplo, ElMatrix_s A, ElMatrix_s w, ElInt lower, ElInt upper,
  ElSortType sortType );
ElError ElHermitianEigIndices_d
( ElUpperOrLower uplo, ElMatrix_d A, ElMatrix_d w, ElInt lower, ElInt upper,
  ElSortType sortType );
ElError ElHermitianEigIndices_c
( ElUpperOrLower uplo, ElMatrix_c A, ElMatrix_s w, ElInt lower, ElInt upper,
  ElSortType sortType );
ElError ElHermitianEigIndices_z
( ElUpperOrLower uplo, ElMatrix_z A, ElMatrix_d w, ElInt lower, ElInt upper,
  ElSortType sortType );

/* NOTE: 'A' must be in a [MC,MR] distribution, while 
         'w' must be in a [VR,STAR] distribution */
/*
ElError ElHermitianEigIndicesDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A, ElDistMatrix_s w, 
  ElInt lower, ElInt upper, ElSortType sortType );
*/
ElError ElHermitianEigIndicesDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A, ElDistMatrix_d w,
  ElInt lower, ElInt upper, ElSortType sortType );
/*
ElError ElHermitianEigIndicesDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A, ElDistMatrix_s w,
  ElInt lower, ElInt upper, ElSortType sortType );
*/
ElError ElHermitianEigIndicesDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A, ElDistMatrix_d w,
  ElInt lower, ElInt upper, ElSortType sortType );

/* TODO: Expert version */

/* Compute the eigenpairs within an index range
   -------------------------------------------- */
ElError ElHermitianEigPairIndices_s
( ElUpperOrLower uplo, ElMatrix_s A, ElMatrix_s w, ElMatrix_s Z, 
  ElInt lower, ElInt upper, ElSortType sortType );
ElError ElHermitianEigPairIndices_d
( ElUpperOrLower uplo, ElMatrix_d A, ElMatrix_d w, ElMatrix_d Z, 
  ElInt lower, ElInt upper, ElSortType sortType );
ElError ElHermitianEigPairIndices_c
( ElUpperOrLower uplo, ElMatrix_c A, ElMatrix_s w, ElMatrix_c Z,
  ElInt lower, ElInt upper, ElSortType sortType );
ElError ElHermitianEigPairIndices_z
( ElUpperOrLower uplo, ElMatrix_z A, ElMatrix_d w, ElMatrix_z Z,
  ElInt lower, ElInt upper, ElSortType sortType );

/* NOTE: 'A' and 'Z' must be in a [MC,MR] distribution, while 
         'w' must be in a [VR,STAR] distribution */
/*
ElError ElHermitianEigPairIndicesDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A, ElDistMatrix_s w, ElDistMatrix_s Z,
  ElInt lower, ElInt upper, ElSortType sortType );
*/
ElError ElHermitianEigPairIndicesDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A, ElDistMatrix_d w, ElDistMatrix_d Z,
  ElInt lower, ElInt upper, ElSortType sortType );
/*
ElError ElHermitianEigPairIndicesDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A, ElDistMatrix_s w, ElDistMatrix_c Z,
  ElInt lower, ElInt upper, ElSortType sortType );
*/
ElError ElHermitianEigPairIndicesDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A, ElDistMatrix_d w, ElDistMatrix_z Z,
  ElInt lower, ElInt upper, ElSortType sortType );

/* TODO: Expert version */

/* Hermitian generalized-definite eigensolvers
   =========================================== */
/* TODO: Mirror approach from above */

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_LAPACK_DECOMP_C_H */
