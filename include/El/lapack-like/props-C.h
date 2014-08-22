/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_LAPACK_PROPS_C_H
#define EL_LAPACK_PROPS_C_H

#include "El/core/DistMatrix-C.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
  EL_PS_TWO_NORM,
  EL_PS_ONE_NORM
} ElPseudospecNorm;

typedef struct {
  ElInt realSize, imagSize;

  ElInt imgSaveFreq, numSaveFreq, imgDispFreq;
  ElInt imgSaveCount, numSaveCount, imgDispCount;
  const char *imgBase, *numBase;
  ElFileFormat imgFormat, numFormat; 
  bool itCounts;
} ElSnapshotCtrl;
ElError ElSnapshotCtrlDefault( ElSnapshotCtrl* ctrl );
/* NOTE: Since conversion from SnapshotCtrl involves deep copies of char* */
ElError ElSnapshotCtrlDestroy( const ElSnapshotCtrl* ctrl );

typedef struct {
  ElPseudospecNorm norm;
  ElInt blockWidth;

  bool schur;
  bool forceComplexSchur;
  bool forceComplexPs;
  ElSchurCtrl_s schurCtrl;

  ElInt maxIts;
  float tol;
  bool deflate;
  
  bool arnoldi;
  ElInt basisSize;
  bool reorthog;

  bool progress;

  ElSnapshotCtrl snapCtrl;
} ElPseudospecCtrl_s;
ElError ElPseudospecCtrlDefault_s( ElPseudospecCtrl_s* ctrl );
/* NOTE: Since conversion from SnapshotCtrl involves deep copies of char* */
ElError ElPseudospecCtrlDestroy_s( const ElPseudospecCtrl_s* ctrl );

typedef struct {
  ElPseudospecNorm norm;
  ElInt blockWidth;

  bool schur;
  bool forceComplexSchur;
  bool forceComplexPs;
  ElSchurCtrl_d schurCtrl;

  ElInt maxIts;
  double tol;
  bool deflate;
  
  bool arnoldi;
  ElInt basisSize;
  bool reorthog;

  bool progress;

  ElSnapshotCtrl snapCtrl;
} ElPseudospecCtrl_d;
ElError ElPseudospecCtrlDefault_d( ElPseudospecCtrl_d* ctrl );
/* NOTE: Since conversion from SnapshotCtrl involves deep copies of char* */
ElError ElPseudospecCtrlDestroy_d( const ElPseudospecCtrl_d* ctrl );

/* Condition number
   ================ */
ElError ElCondition_s( ElConstMatrix_s A, ElNormType type, float* cond );
ElError ElCondition_d( ElConstMatrix_d A, ElNormType type, double* cond );
ElError ElCondition_c( ElConstMatrix_c A, ElNormType type, float* cond );
ElError ElCondition_z( ElConstMatrix_z A, ElNormType type, double* cond );

ElError ElConditionDist_s
( ElConstDistMatrix_s A, ElNormType type, float* cond );
ElError ElConditionDist_d
( ElConstDistMatrix_d A, ElNormType type, double* cond );
ElError ElConditionDist_c
( ElConstDistMatrix_c A, ElNormType type, float* cond );
ElError ElConditionDist_z
( ElConstDistMatrix_z A, ElNormType type, double* cond );

/* Frobenius
   --------- */
ElError ElFrobeniusCondition_s( ElConstMatrix_s A, float* cond );
ElError ElFrobeniusCondition_d( ElConstMatrix_d A, double* cond );
ElError ElFrobeniusCondition_c( ElConstMatrix_c A, float* cond );
ElError ElFrobeniusCondition_z( ElConstMatrix_z A, double* cond );

ElError ElFrobeniusConditionDist_s( ElConstDistMatrix_s A, float* cond );
ElError ElFrobeniusConditionDist_d( ElConstDistMatrix_d A, double* cond );
ElError ElFrobeniusConditionDist_c( ElConstDistMatrix_c A, float* cond );
ElError ElFrobeniusConditionDist_z( ElConstDistMatrix_z A, double* cond );

/* Infinity
   -------- */
ElError ElInfinityCondition_s( ElConstMatrix_s A, float* cond );
ElError ElInfinityCondition_d( ElConstMatrix_d A, double* cond );
ElError ElInfinityCondition_c( ElConstMatrix_c A, float* cond );
ElError ElInfinityCondition_z( ElConstMatrix_z A, double* cond );

ElError ElInfinityConditionDist_s( ElConstDistMatrix_s A, float* cond );
ElError ElInfinityConditionDist_d( ElConstDistMatrix_d A, double* cond );
ElError ElInfinityConditionDist_c( ElConstDistMatrix_c A, float* cond );
ElError ElInfinityConditionDist_z( ElConstDistMatrix_z A, double* cond );

/* Max
   --- */
ElError ElMaxCondition_s( ElConstMatrix_s A, float* cond );
ElError ElMaxCondition_d( ElConstMatrix_d A, double* cond );
ElError ElMaxCondition_c( ElConstMatrix_c A, float* cond );
ElError ElMaxCondition_z( ElConstMatrix_z A, double* cond );

ElError ElMaxConditionDist_s( ElConstDistMatrix_s A, float* cond );
ElError ElMaxConditionDist_d( ElConstDistMatrix_d A, double* cond );
ElError ElMaxConditionDist_c( ElConstDistMatrix_c A, float* cond );
ElError ElMaxConditionDist_z( ElConstDistMatrix_z A, double* cond );

/* One
   --- */
ElError ElOneCondition_s( ElConstMatrix_s A, float* cond );
ElError ElOneCondition_d( ElConstMatrix_d A, double* cond );
ElError ElOneCondition_c( ElConstMatrix_c A, float* cond );
ElError ElOneCondition_z( ElConstMatrix_z A, double* cond );

ElError ElOneConditionDist_s( ElConstDistMatrix_s A, float* cond );
ElError ElOneConditionDist_d( ElConstDistMatrix_d A, double* cond );
ElError ElOneConditionDist_c( ElConstDistMatrix_c A, float* cond );
ElError ElOneConditionDist_z( ElConstDistMatrix_z A, double* cond );

/* Two
   --- */
ElError ElTwoCondition_s( ElConstMatrix_s A, float* cond );
ElError ElTwoCondition_d( ElConstMatrix_d A, double* cond );
ElError ElTwoCondition_c( ElConstMatrix_c A, float* cond );
ElError ElTwoCondition_z( ElConstMatrix_z A, double* cond );

ElError ElTwoConditionDist_s( ElConstDistMatrix_s A, float* cond );
ElError ElTwoConditionDist_d( ElConstDistMatrix_d A, double* cond );
ElError ElTwoConditionDist_c( ElConstDistMatrix_c A, float* cond );
ElError ElTwoConditionDist_z( ElConstDistMatrix_z A, double* cond );

/* Determinant
   =========== */
/* Return the result in a safer, expanded format
   --------------------------------------------- */
/* TODO: Version which allows for overwriting? How to name? */
ElError ElSafeDeterminant_s( ElConstMatrix_s A, ElSafeProduct_s* prod );
ElError ElSafeDeterminant_d( ElConstMatrix_d A, ElSafeProduct_d* prod );
ElError ElSafeDeterminant_c( ElConstMatrix_c A, ElSafeProduct_c* prod );
ElError ElSafeDeterminant_z( ElConstMatrix_z A, ElSafeProduct_z* prod );

ElError ElSafeDeterminantDist_s( ElConstDistMatrix_s A, ElSafeProduct_s* prod );
ElError ElSafeDeterminantDist_d( ElConstDistMatrix_d A, ElSafeProduct_d* prod );
ElError ElSafeDeterminantDist_c( ElConstDistMatrix_c A, ElSafeProduct_c* prod );
ElError ElSafeDeterminantDist_z( ElConstDistMatrix_z A, ElSafeProduct_z* prod );

/* TODO: Version which allows for overwriting? How to name? */
ElError ElSafeHPDDeterminant_s
( ElUpperOrLower uplo, ElConstMatrix_s A, ElSafeProduct_s* prod );
ElError ElSafeHPDDeterminant_d
( ElUpperOrLower uplo, ElConstMatrix_d A, ElSafeProduct_d* prod );
ElError ElSafeHPDDeterminant_c
( ElUpperOrLower uplo, ElConstMatrix_c A, ElSafeProduct_s* prod );
ElError ElSafeHPDDeterminant_z
( ElUpperOrLower uplo, ElConstMatrix_z A, ElSafeProduct_d* prod );

ElError ElSafeHPDDeterminantDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, ElSafeProduct_s* prod );
ElError ElSafeHPDDeterminantDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, ElSafeProduct_d* prod );
ElError ElSafeHPDDeterminantDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, ElSafeProduct_s* prod );
ElError ElSafeHPDDeterminantDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, ElSafeProduct_d* prod );

/* Return the direct result
   ------------------------ */
/* TODO: Version which allows for overwriting? How to name? */
ElError ElDeterminant_s( ElConstMatrix_s A, float* prod );
ElError ElDeterminant_d( ElConstMatrix_d A, double* prod );
ElError ElDeterminant_c( ElConstMatrix_c A, complex_float* prod );
ElError ElDeterminant_z( ElConstMatrix_z A, complex_double* prod );

ElError ElDeterminantDist_s( ElConstDistMatrix_s A, float* prod );
ElError ElDeterminantDist_d( ElConstDistMatrix_d A, double* prod );
ElError ElDeterminantDist_c( ElConstDistMatrix_c A, complex_float* prod );
ElError ElDeterminantDist_z( ElConstDistMatrix_z A, complex_double* prod );

/* TODO: Version which allows for overwriting? How to name? */
ElError ElHPDDeterminant_s
( ElUpperOrLower uplo, ElConstMatrix_s A, float* prod );
ElError ElHPDDeterminant_d
( ElUpperOrLower uplo, ElConstMatrix_d A, double* prod );
ElError ElHPDDeterminant_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float* prod );
ElError ElHPDDeterminant_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double* prod );

ElError ElHPDDeterminantDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, float* prod );
ElError ElHPDDeterminantDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, double* prod );
ElError ElHPDDeterminantDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float* prod );
ElError ElHPDDeterminantDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double* prod );

/* TODO: Determinant after Cholesky or LU */

/* Inertia
   ======= */
ElError ElInertia_s
( ElUpperOrLower uplo, ElMatrix_s A, ElLDLPivotType pivotType, 
  ElInertiaType* inertia );
ElError ElInertia_d
( ElUpperOrLower uplo, ElMatrix_d A, ElLDLPivotType pivotType,
  ElInertiaType* inertia );
ElError ElInertia_c
( ElUpperOrLower uplo, ElMatrix_c A, ElLDLPivotType pivotType,
  ElInertiaType* inertia );
ElError ElInertia_z
( ElUpperOrLower uplo, ElMatrix_z A, ElLDLPivotType pivotType,
  ElInertiaType* inertia );

ElError ElInertiaDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A, ElLDLPivotType pivotType,
  ElInertiaType* inertia );
ElError ElInertiaDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A, ElLDLPivotType pivotType,
  ElInertiaType* inertia );
ElError ElInertiaDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A, ElLDLPivotType pivotType,
  ElInertiaType* inertia );
ElError ElInertiaDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A, ElLDLPivotType pivotType,
  ElInertiaType* inertia );

/* Norm
   ==== */
ElError ElNorm_s( ElConstMatrix_s A, ElNormType normType, float* norm );
ElError ElNorm_d( ElConstMatrix_d A, ElNormType normType, double* norm );
ElError ElNorm_c( ElConstMatrix_c A, ElNormType normType, float* norm );
ElError ElNorm_z( ElConstMatrix_z A, ElNormType normType, double* norm );

ElError ElNormDist_s
( ElConstDistMatrix_s A, ElNormType normType, float* norm );
ElError ElNormDist_d
( ElConstDistMatrix_d A, ElNormType normType, double* norm );
ElError ElNormDist_c
( ElConstDistMatrix_c A, ElNormType normType, float* norm );
ElError ElNormDist_z
( ElConstDistMatrix_z A, ElNormType normType, double* norm );

ElError ElHermitianNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, ElNormType normType, float* norm );
ElError ElHermitianNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, ElNormType normType, double* norm );

ElError ElHermitianNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, ElNormType normType, 
  float* norm );
ElError ElHermitianNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, ElNormType normType, 
  double* norm );

ElError ElSymmetricNorm_s
( ElUpperOrLower uplo, ElConstMatrix_s A, ElNormType normType, float* norm );
ElError ElSymmetricNorm_d
( ElUpperOrLower uplo, ElConstMatrix_d A, ElNormType normType, double* norm );
ElError ElSymmetricNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, ElNormType normType, float* norm );
ElError ElSymmetricNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, ElNormType normType, double* norm );

ElError ElSymmetricNormDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, ElNormType normType, 
  float* norm );
ElError ElSymmetricNormDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, ElNormType normType, 
  double* norm );
ElError ElSymmetricNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, ElNormType normType, 
  float* norm );
ElError ElSymmetricNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, ElNormType normType, 
  double* norm );

/* Entrywise norm
   -------------- */
ElError ElEntrywiseNorm_s( ElConstMatrix_s A, float p, float* norm );
ElError ElEntrywiseNorm_d( ElConstMatrix_d A, double p, double* norm );
ElError ElEntrywiseNorm_c( ElConstMatrix_c A, float p, float* norm );
ElError ElEntrywiseNorm_z( ElConstMatrix_z A, double p, double* norm );

ElError ElEntrywiseNormDist_s( ElConstDistMatrix_s A, float p, float* norm );
ElError ElEntrywiseNormDist_d( ElConstDistMatrix_d A, double p, double* norm );
ElError ElEntrywiseNormDist_c( ElConstDistMatrix_c A, float p, float* norm );
ElError ElEntrywiseNormDist_z( ElConstDistMatrix_z A, double p, double* norm );

ElError ElHermitianEntrywiseNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float p, float* norm );
ElError ElHermitianEntrywiseNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double p, double* norm );

ElError ElHermitianEntrywiseNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float p, float* norm );
ElError ElHermitianEntrywiseNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double p, double* norm );

ElError ElSymmetricEntrywiseNorm_s
( ElUpperOrLower uplo, ElConstMatrix_s A, float p, float* norm );
ElError ElSymmetricEntrywiseNorm_d
( ElUpperOrLower uplo, ElConstMatrix_d A, double p, double* norm );
ElError ElSymmetricEntrywiseNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float p, float* norm );
ElError ElSymmetricEntrywiseNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double p, double* norm );

ElError ElSymmetricEntrywiseNormDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, float p, float* norm );
ElError ElSymmetricEntrywiseNormDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, double p, double* norm );
ElError ElSymmetricEntrywiseNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float p, float* norm );
ElError ElSymmetricEntrywiseNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double p, double* norm );

/* Entrywise one norm
   ------------------ */
ElError ElEntrywiseOneNorm_s( ElConstMatrix_s A, float* norm );
ElError ElEntrywiseOneNorm_d( ElConstMatrix_d A, double* norm );
ElError ElEntrywiseOneNorm_c( ElConstMatrix_c A, float* norm );
ElError ElEntrywiseOneNorm_z( ElConstMatrix_z A, double* norm );

ElError ElEntrywiseOneNormDist_s( ElConstDistMatrix_s A, float* norm );
ElError ElEntrywiseOneNormDist_d( ElConstDistMatrix_d A, double* norm );
ElError ElEntrywiseOneNormDist_c( ElConstDistMatrix_c A, float* norm );
ElError ElEntrywiseOneNormDist_z( ElConstDistMatrix_z A, double* norm );

ElError ElHermitianEntrywiseOneNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float* norm );
ElError ElHermitianEntrywiseOneNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double* norm );

ElError ElHermitianEntrywiseOneNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float* norm );
ElError ElHermitianEntrywiseOneNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double* norm );

ElError ElSymmetricEntrywiseOneNorm_s
( ElUpperOrLower uplo, ElConstMatrix_s A, float* norm );
ElError ElSymmetricEntrywiseOneNorm_d
( ElUpperOrLower uplo, ElConstMatrix_d A, double* norm );
ElError ElSymmetricEntrywiseOneNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float* norm );
ElError ElSymmetricEntrywiseOneNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double* norm );

ElError ElSymmetricEntrywiseOneNormDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, float* norm );
ElError ElSymmetricEntrywiseOneNormDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, double* norm );
ElError ElSymmetricEntrywiseOneNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float* norm );
ElError ElSymmetricEntrywiseOneNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double* norm );

/* Frobenius norm
   -------------- */
ElError ElFrobeniusNorm_s( ElConstMatrix_s A, float* norm );
ElError ElFrobeniusNorm_d( ElConstMatrix_d A, double* norm );
ElError ElFrobeniusNorm_c( ElConstMatrix_c A, float* norm );
ElError ElFrobeniusNorm_z( ElConstMatrix_z A, double* norm );

ElError ElFrobeniusNormDist_s( ElConstDistMatrix_s A, float* norm );
ElError ElFrobeniusNormDist_d( ElConstDistMatrix_d A, double* norm );
ElError ElFrobeniusNormDist_c( ElConstDistMatrix_c A, float* norm );
ElError ElFrobeniusNormDist_z( ElConstDistMatrix_z A, double* norm );

ElError ElHermitianFrobeniusNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float* norm );
ElError ElHermitianFrobeniusNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double* norm );

ElError ElHermitianFrobeniusNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float* norm );
ElError ElHermitianFrobeniusNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double* norm );

ElError ElSymmetricFrobeniusNorm_s
( ElUpperOrLower uplo, ElConstMatrix_s A, float* norm );
ElError ElSymmetricFrobeniusNorm_d
( ElUpperOrLower uplo, ElConstMatrix_d A, double* norm );
ElError ElSymmetricFrobeniusNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float* norm );
ElError ElSymmetricFrobeniusNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double* norm );

ElError ElSymmetricFrobeniusNormDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, float* norm );
ElError ElSymmetricFrobeniusNormDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, double* norm );
ElError ElSymmetricFrobeniusNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float* norm );
ElError ElSymmetricFrobeniusNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double* norm );

/* Infinity norm
   ------------- */
ElError ElInfinityNorm_s( ElConstMatrix_s A, float* norm );
ElError ElInfinityNorm_d( ElConstMatrix_d A, double* norm );
ElError ElInfinityNorm_c( ElConstMatrix_c A, float* norm );
ElError ElInfinityNorm_z( ElConstMatrix_z A, double* norm );

ElError ElInfinityNormDist_s( ElConstDistMatrix_s A, float* norm );
ElError ElInfinityNormDist_d( ElConstDistMatrix_d A, double* norm );
ElError ElInfinityNormDist_c( ElConstDistMatrix_c A, float* norm );
ElError ElInfinityNormDist_z( ElConstDistMatrix_z A, double* norm );

ElError ElHermitianInfinityNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float* norm );
ElError ElHermitianInfinityNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double* norm );

ElError ElHermitianInfinityNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float* norm );
ElError ElHermitianInfinityNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double* norm );

ElError ElSymmetricInfinityNorm_s
( ElUpperOrLower uplo, ElConstMatrix_s A, float* norm );
ElError ElSymmetricInfinityNorm_d
( ElUpperOrLower uplo, ElConstMatrix_d A, double* norm );
ElError ElSymmetricInfinityNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float* norm );
ElError ElSymmetricInfinityNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double* norm );

ElError ElSymmetricInfinityNormDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, float* norm );
ElError ElSymmetricInfinityNormDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, double* norm );
ElError ElSymmetricInfinityNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float* norm );
ElError ElSymmetricInfinityNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double* norm );

/* Ky-Fan norm
   ----------- */
ElError ElKyFanNorm_s( ElConstMatrix_s A, ElInt k, float* norm );
ElError ElKyFanNorm_d( ElConstMatrix_d A, ElInt k, double* norm );
ElError ElKyFanNorm_c( ElConstMatrix_c A, ElInt k, float* norm );
ElError ElKyFanNorm_z( ElConstMatrix_z A, ElInt k, double* norm );

ElError ElKyFanNormDist_s( ElConstDistMatrix_s A, ElInt k, float* norm );
ElError ElKyFanNormDist_d( ElConstDistMatrix_d A, ElInt k, double* norm );
ElError ElKyFanNormDist_c( ElConstDistMatrix_c A, ElInt k, float* norm );
ElError ElKyFanNormDist_z( ElConstDistMatrix_z A, ElInt k, double* norm );

ElError ElHermitianKyFanNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, ElInt k, float* norm );
ElError ElHermitianKyFanNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, ElInt k, double* norm );

ElError ElHermitianKyFanNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, ElInt k, float* norm );
ElError ElHermitianKyFanNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, ElInt k, double* norm );

ElError ElSymmetricKyFanNorm_s
( ElUpperOrLower uplo, ElConstMatrix_s A, ElInt k, float* norm );
ElError ElSymmetricKyFanNorm_d
( ElUpperOrLower uplo, ElConstMatrix_d A, ElInt k, double* norm );
ElError ElSymmetricKyFanNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, ElInt k, float* norm );
ElError ElSymmetricKyFanNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, ElInt k, double* norm );

ElError ElSymmetricKyFanNormDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, ElInt k, float* norm );
ElError ElSymmetricKyFanNormDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, ElInt k, double* norm );
ElError ElSymmetricKyFanNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, ElInt k, float* norm );
ElError ElSymmetricKyFanNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, ElInt k, double* norm );

/* Ky-Fan-Schatten norm
   -------------------- */
ElError ElKyFanSchattenNorm_s
( ElConstMatrix_s A, ElInt k, float p, float* norm );
ElError ElKyFanSchattenNorm_d
( ElConstMatrix_d A, ElInt k, double p, double* norm );
ElError ElKyFanSchattenNorm_c
( ElConstMatrix_c A, ElInt k, float p, float* norm );
ElError ElKyFanSchattenNorm_z
( ElConstMatrix_z A, ElInt k, double p, double* norm );

ElError ElKyFanSchattenNormDist_s
( ElConstDistMatrix_s A, ElInt k, float p, float* norm );
ElError ElKyFanSchattenNormDist_d
( ElConstDistMatrix_d A, ElInt k, double p, double* norm );
ElError ElKyFanSchattenNormDist_c
( ElConstDistMatrix_c A, ElInt k, float p, float* norm );
ElError ElKyFanSchattenNormDist_z
( ElConstDistMatrix_z A, ElInt k, double p, double* norm );

ElError ElHermitianKyFanSchattenNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, ElInt k, float p, float* norm );
ElError ElHermitianKyFanSchattenNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, ElInt k, double p, double* norm );

ElError ElHermitianKyFanSchattenNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, ElInt k, float p, float* norm );
ElError ElHermitianKyFanSchattenNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, ElInt k, double p, double* norm );

ElError ElSymmetricKyFanSchattenNorm_s
( ElUpperOrLower uplo, ElConstMatrix_s A, ElInt k, float p, float* norm );
ElError ElSymmetricKyFanSchattenNorm_d
( ElUpperOrLower uplo, ElConstMatrix_d A, ElInt k, double p, double* norm );
ElError ElSymmetricKyFanSchattenNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, ElInt k, float p, float* norm );
ElError ElSymmetricKyFanSchattenNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, ElInt k, double p, double* norm );

ElError ElSymmetricKyFanSchattenNormDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, ElInt k, float p, float* norm );
ElError ElSymmetricKyFanSchattenNormDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, ElInt k, double p, double* norm );
ElError ElSymmetricKyFanSchattenNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, ElInt k, float p, float* norm );
ElError ElSymmetricKyFanSchattenNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, ElInt k, double p, double* norm );

/* Max norm
   -------- */
ElError ElMaxNorm_s( ElConstMatrix_s A, float* norm );
ElError ElMaxNorm_d( ElConstMatrix_d A, double* norm );
ElError ElMaxNorm_c( ElConstMatrix_c A, float* norm );
ElError ElMaxNorm_z( ElConstMatrix_z A, double* norm );

ElError ElMaxNormDist_s( ElConstDistMatrix_s A, float* norm );
ElError ElMaxNormDist_d( ElConstDistMatrix_d A, double* norm );
ElError ElMaxNormDist_c( ElConstDistMatrix_c A, float* norm );
ElError ElMaxNormDist_z( ElConstDistMatrix_z A, double* norm );

ElError ElHermitianMaxNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float* norm );
ElError ElHermitianMaxNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double* norm );

ElError ElHermitianMaxNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float* norm );
ElError ElHermitianMaxNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double* norm );

ElError ElSymmetricMaxNorm_s
( ElUpperOrLower uplo, ElConstMatrix_s A, float* norm );
ElError ElSymmetricMaxNorm_d
( ElUpperOrLower uplo, ElConstMatrix_d A, double* norm );
ElError ElSymmetricMaxNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float* norm );
ElError ElSymmetricMaxNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double* norm );

ElError ElSymmetricMaxNormDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, float* norm );
ElError ElSymmetricMaxNormDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, double* norm );
ElError ElSymmetricMaxNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float* norm );
ElError ElSymmetricMaxNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double* norm );

/* Nuclear norm
   ------------ */
ElError ElNuclearNorm_s( ElConstMatrix_s A, float* norm );
ElError ElNuclearNorm_d( ElConstMatrix_d A, double* norm );
ElError ElNuclearNorm_c( ElConstMatrix_c A, float* norm );
ElError ElNuclearNorm_z( ElConstMatrix_z A, double* norm );

ElError ElNuclearNormDist_s( ElConstDistMatrix_s A, float* norm );
ElError ElNuclearNormDist_d( ElConstDistMatrix_d A, double* norm );
ElError ElNuclearNormDist_c( ElConstDistMatrix_c A, float* norm );
ElError ElNuclearNormDist_z( ElConstDistMatrix_z A, double* norm );

ElError ElHermitianNuclearNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float* norm );
ElError ElHermitianNuclearNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double* norm );

ElError ElHermitianNuclearNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float* norm );
ElError ElHermitianNuclearNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double* norm );

ElError ElSymmetricNuclearNorm_s
( ElUpperOrLower uplo, ElConstMatrix_s A, float* norm );
ElError ElSymmetricNuclearNorm_d
( ElUpperOrLower uplo, ElConstMatrix_d A, double* norm );
ElError ElSymmetricNuclearNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float* norm );
ElError ElSymmetricNuclearNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double* norm );

ElError ElSymmetricNuclearNormDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, float* norm );
ElError ElSymmetricNuclearNormDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, double* norm );
ElError ElSymmetricNuclearNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float* norm );
ElError ElSymmetricNuclearNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double* norm );

/* One norm
   -------- */
ElError ElOneNorm_s( ElConstMatrix_s A, float* norm );
ElError ElOneNorm_d( ElConstMatrix_d A, double* norm );
ElError ElOneNorm_c( ElConstMatrix_c A, float* norm );
ElError ElOneNorm_z( ElConstMatrix_z A, double* norm );

ElError ElOneNormDist_s( ElConstDistMatrix_s A, float* norm );
ElError ElOneNormDist_d( ElConstDistMatrix_d A, double* norm );
ElError ElOneNormDist_c( ElConstDistMatrix_c A, float* norm );
ElError ElOneNormDist_z( ElConstDistMatrix_z A, double* norm );

ElError ElHermitianOneNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float* norm );
ElError ElHermitianOneNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double* norm );

ElError ElHermitianOneNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float* norm );
ElError ElHermitianOneNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double* norm );

ElError ElSymmetricOneNorm_s
( ElUpperOrLower uplo, ElConstMatrix_s A, float* norm );
ElError ElSymmetricOneNorm_d
( ElUpperOrLower uplo, ElConstMatrix_d A, double* norm );
ElError ElSymmetricOneNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float* norm );
ElError ElSymmetricOneNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double* norm );

ElError ElSymmetricOneNormDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, float* norm );
ElError ElSymmetricOneNormDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, double* norm );
ElError ElSymmetricOneNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float* norm );
ElError ElSymmetricOneNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double* norm );

/* Schatten norm
   ------------- */
ElError ElSchattenNorm_s( ElConstMatrix_s A, float p, float* norm );
ElError ElSchattenNorm_d( ElConstMatrix_d A, double p, double* norm );
ElError ElSchattenNorm_c( ElConstMatrix_c A, float p, float* norm );
ElError ElSchattenNorm_z( ElConstMatrix_z A, double p, double* norm );

ElError ElSchattenNormDist_s( ElConstDistMatrix_s A, float p, float* norm );
ElError ElSchattenNormDist_d( ElConstDistMatrix_d A, double p, double* norm );
ElError ElSchattenNormDist_c( ElConstDistMatrix_c A, float p, float* norm );
ElError ElSchattenNormDist_z( ElConstDistMatrix_z A, double p, double* norm );

ElError ElHermitianSchattenNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float p, float* norm );
ElError ElHermitianSchattenNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double p, double* norm );

ElError ElHermitianSchattenNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float p, float* norm );
ElError ElHermitianSchattenNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double p, double* norm );

ElError ElSymmetricSchattenNorm_s
( ElUpperOrLower uplo, ElConstMatrix_s A, float p, float* norm );
ElError ElSymmetricSchattenNorm_d
( ElUpperOrLower uplo, ElConstMatrix_d A, double p, double* norm );
ElError ElSymmetricSchattenNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float p, float* norm );
ElError ElSymmetricSchattenNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double p, double* norm );

ElError ElSymmetricSchattenNormDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, float p, float* norm );
ElError ElSymmetricSchattenNormDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, double p, double* norm );
ElError ElSymmetricSchattenNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float p, float* norm );
ElError ElSymmetricSchattenNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double p, double* norm );

/* Two norm
   -------- */
ElError ElTwoNorm_s( ElConstMatrix_s A, float* norm );
ElError ElTwoNorm_d( ElConstMatrix_d A, double* norm );
ElError ElTwoNorm_c( ElConstMatrix_c A, float* norm );
ElError ElTwoNorm_z( ElConstMatrix_z A, double* norm );

ElError ElTwoNormDist_s( ElConstDistMatrix_s A, float* norm );
ElError ElTwoNormDist_d( ElConstDistMatrix_d A, double* norm );
ElError ElTwoNormDist_c( ElConstDistMatrix_c A, float* norm );
ElError ElTwoNormDist_z( ElConstDistMatrix_z A, double* norm );

ElError ElHermitianTwoNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float* norm );
ElError ElHermitianTwoNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double* norm );

ElError ElHermitianTwoNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float* norm );
ElError ElHermitianTwoNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double* norm );

ElError ElSymmetricTwoNorm_s
( ElUpperOrLower uplo, ElConstMatrix_s A, float* norm );
ElError ElSymmetricTwoNorm_d
( ElUpperOrLower uplo, ElConstMatrix_d A, double* norm );
ElError ElSymmetricTwoNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float* norm );
ElError ElSymmetricTwoNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double* norm );

ElError ElSymmetricTwoNormDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, float* norm );
ElError ElSymmetricTwoNormDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, double* norm );
ElError ElSymmetricTwoNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float* norm );
ElError ElSymmetricTwoNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double* norm );

/* Zero norm
   --------- */
ElError ElZeroNorm_i( ElConstMatrix_s A, ElInt tol, ElInt* norm );
ElError ElZeroNorm_s( ElConstMatrix_s A, float tol, ElInt* norm );
ElError ElZeroNorm_d( ElConstMatrix_d A, double tol, ElInt* norm );
ElError ElZeroNorm_c( ElConstMatrix_c A, float tol, ElInt* norm );
ElError ElZeroNorm_z( ElConstMatrix_z A, double tol, ElInt* norm );

ElError ElZeroNormDist_i( ElConstDistMatrix_s A, ElInt tol, ElInt* norm );
ElError ElZeroNormDist_s( ElConstDistMatrix_s A, float tol, ElInt* norm );
ElError ElZeroNormDist_d( ElConstDistMatrix_d A, double tol, ElInt* norm );
ElError ElZeroNormDist_c( ElConstDistMatrix_c A, float tol, ElInt* norm );
ElError ElZeroNormDist_z( ElConstDistMatrix_z A, double tol, ElInt* norm );

/* Two-norm estimates
   ------------------ */
ElError ElTwoNormEstimate_s
( ElConstMatrix_s A, float tol, ElInt maxIts, float* normEst );
ElError ElTwoNormEstimate_d
( ElConstMatrix_d A, double tol, ElInt maxIts, double* normEst );
ElError ElTwoNormEstimate_c
( ElConstMatrix_c A, float tol, ElInt maxIts, float* normEst );
ElError ElTwoNormEstimate_z
( ElConstMatrix_z A, double tol, ElInt maxIts, double* normEst );

ElError ElTwoNormEstimateDist_s
( ElConstDistMatrix_s A, float tol, ElInt maxIts, float* normEst );
ElError ElTwoNormEstimateDist_d
( ElConstDistMatrix_d A, double tol, ElInt maxIts, double* normEst );
ElError ElTwoNormEstimateDist_c
( ElConstDistMatrix_c A, float tol, ElInt maxIts, float* normEst );
ElError ElTwoNormEstimateDist_z
( ElConstDistMatrix_z A, double tol, ElInt maxIts, double* normEst );

ElError ElSymmetricTwoNormEstimate_s
( ElUpperOrLower uplo, ElConstMatrix_s A, float tol, ElInt maxIts, 
  float* normEst );
ElError ElSymmetricTwoNormEstimate_d
( ElUpperOrLower uplo, ElConstMatrix_d A, double tol, ElInt maxIts, 
  double* normEst );
ElError ElSymmetricTwoNormEstimate_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float tol, ElInt maxIts, 
  float* normEst );
ElError ElSymmetricTwoNormEstimate_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double tol, ElInt maxIts, 
  double* normEst );

ElError ElSymmetricTwoNormEstimateDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, float tol, ElInt maxIts, 
  float* normEst );
ElError ElSymmetricTwoNormEstimateDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, double tol, ElInt maxIts, 
  double* normEst );
ElError ElSymmetricTwoNormEstimateDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float tol, ElInt maxIts, 
  float* normEst );
ElError ElSymmetricTwoNormEstimateDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double tol, ElInt maxIts, 
  double* normEst );

ElError ElHermitianTwoNormEstimate_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float tol, ElInt maxIts, 
  float* normEst );
ElError ElHermitianTwoNormEstimate_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double tol, ElInt maxIts, 
  double* normEst );

ElError ElHermitianTwoNormEstimateDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float tol, ElInt maxIts, 
  float* normEst );
ElError ElHermitianTwoNormEstimateDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double tol, ElInt maxIts, 
  double* normEst );

/* Pseudospectra
   ============= */
/* Automatic window choice
   ----------------------- */
ElError ElPseudospectralWindowAuto_s
( ElConstMatrix_s A, ElMatrix_s invNormMap, ElInt realSize, ElInt imagSize );
ElError ElPseudospectralWindowAuto_d
( ElConstMatrix_d A, ElMatrix_d invNormMap, ElInt realSize, ElInt imagSize );
ElError ElPseudospectralWindowAuto_c
( ElConstMatrix_c A, ElMatrix_s invNormMap, ElInt realSize, ElInt imagSize );
ElError ElPseudospectralWindowAuto_z
( ElConstMatrix_z A, ElMatrix_d invNormMap, ElInt realSize, ElInt imagSize );

ElError ElPseudospectralWindowAutoDist_s
( ElConstDistMatrix_s A, ElDistMatrix_s invNormMap, 
  ElInt realSize, ElInt imagSize );
ElError ElPseudospectralWindowAutoDist_d
( ElConstDistMatrix_d A, ElDistMatrix_d invNormMap, 
  ElInt realSize, ElInt imagSize );
ElError ElPseudospectralWindowAutoDist_c
( ElConstDistMatrix_c A, ElDistMatrix_s invNormMap, 
  ElInt realSize, ElInt imagSize );
ElError ElPseudospectralWindowAutoDist_z
( ElConstDistMatrix_z A, ElDistMatrix_d invNormMap, 
  ElInt realSize, ElInt imagSize );

/* Expert version
   ^^^^^^^^^^^^^^ */
ElError ElPseudospectralWindowAutoX_s
( ElConstMatrix_s A, ElMatrix_s invNormMap, ElInt realSize, ElInt imagSize,
  ElPseudospecCtrl_s ctrl );
ElError ElPseudospectralWindowAutoX_d
( ElConstMatrix_d A, ElMatrix_d invNormMap, ElInt realSize, ElInt imagSize,
  ElPseudospecCtrl_d ctrl );
ElError ElPseudospectralWindowAutoX_c
( ElConstMatrix_c A, ElMatrix_s invNormMap, ElInt realSize, ElInt imagSize,
  ElPseudospecCtrl_s ctrl );
ElError ElPseudospectralWindowAutoX_z
( ElConstMatrix_z A, ElMatrix_d invNormMap, ElInt realSize, ElInt imagSize,
  ElPseudospecCtrl_d ctrl );

ElError ElPseudospectralWindowAutoXDist_s
( ElConstDistMatrix_s A, ElDistMatrix_s invNormMap, 
  ElInt realSize, ElInt imagSize, ElPseudospecCtrl_s ctrl );
ElError ElPseudospectralWindowAutoXDist_d
( ElConstDistMatrix_d A, ElDistMatrix_d invNormMap, 
  ElInt realSize, ElInt imagSize, ElPseudospecCtrl_d ctrl );
ElError ElPseudospectralWindowAutoXDist_c
( ElConstDistMatrix_c A, ElDistMatrix_s invNormMap, 
  ElInt realSize, ElInt imagSize, ElPseudospecCtrl_s ctrl );
ElError ElPseudospectralWindowAutoXDist_z
( ElConstDistMatrix_z A, ElDistMatrix_d invNormMap, 
  ElInt realSize, ElInt imagSize, ElPseudospecCtrl_d ctrl );

/* Manual window choice
   -------------------- */
ElError ElPseudospectralWindow_s
( ElConstMatrix_s A, ElMatrix_s invNormMap, 
  complex_float center, float realWidth, float imagWidth,
  ElInt realSize, ElInt imagSize );
ElError ElPseudospectralWindow_d
( ElConstMatrix_d A, ElMatrix_d invNormMap, 
  complex_double center, double realWidth, double imagWidth,
  ElInt realSize, ElInt imagSize );
ElError ElPseudospectralWindow_c
( ElConstMatrix_c A, ElMatrix_s invNormMap, 
  complex_float center, float realWidth, float imagWidth,
  ElInt realSize, ElInt imagSize );
ElError ElPseudospectralWindow_z
( ElConstMatrix_z A, ElMatrix_d invNormMap, 
  complex_double center, double realWidth, double imagWidth,
  ElInt realSize, ElInt imagSize );

ElError ElPseudospectralWindowDist_s
( ElConstDistMatrix_s A, ElDistMatrix_s invNormMap, 
  complex_float center, float realWidth, float imagWidth,
  ElInt realSize, ElInt imagSize );
ElError ElPseudospectralWindowDist_d
( ElConstDistMatrix_d A, ElDistMatrix_d invNormMap, 
  complex_double center, double realWidth, double imagWidth,
  ElInt realSize, ElInt imagSize );
ElError ElPseudospectralWindowDist_c
( ElConstDistMatrix_c A, ElDistMatrix_s invNormMap, 
  complex_float center, float realWidth, float imagWidth,
  ElInt realSize, ElInt imagSize );
ElError ElPseudospectralWindowDist_z
( ElConstDistMatrix_z A, ElDistMatrix_d invNormMap, 
  complex_double center, double realWidth, double imagWidth,
  ElInt realSize, ElInt imagSize );

/* Expert version
   ^^^^^^^^^^^^^^ */
ElError ElPseudospectralWindowX_s
( ElConstMatrix_s A, ElMatrix_s invNormMap, 
  complex_float center, float realWidth, float imagWidth,
  ElInt realSize, ElInt imagSize, ElPseudospecCtrl_s ctrl );
ElError ElPseudospectralWindowX_d
( ElConstMatrix_d A, ElMatrix_d invNormMap, 
  complex_double center, double realWidth, double imagWidth,
  ElInt realSize, ElInt imagSize, ElPseudospecCtrl_d ctrl );
ElError ElPseudospectralWindowX_c
( ElConstMatrix_c A, ElMatrix_s invNormMap, 
  complex_float center, float realWidth, float imagWidth,
  ElInt realSize, ElInt imagSize, ElPseudospecCtrl_s ctrl );
ElError ElPseudospectralWindowX_z
( ElConstMatrix_z A, ElMatrix_d invNormMap, 
  complex_double center, double realWidth, double imagWidth,
  ElInt realSize, ElInt imagSize, ElPseudospecCtrl_d ctrl );

ElError ElPseudospectralWindowXDist_s
( ElConstDistMatrix_s A, ElDistMatrix_s invNormMap, 
  complex_float center, float realWidth, float imagWidth,
  ElInt realSize, ElInt imagSize, ElPseudospecCtrl_s ctrl );
ElError ElPseudospectralWindowXDist_d
( ElConstDistMatrix_d A, ElDistMatrix_d invNormMap, 
  complex_double center, double realWidth, double imagWidth,
  ElInt realSize, ElInt imagSize, ElPseudospecCtrl_d ctrl );
ElError ElPseudospectralWindowXDist_c
( ElConstDistMatrix_c A, ElDistMatrix_s invNormMap, 
  complex_float center, float realWidth, float imagWidth,
  ElInt realSize, ElInt imagSize, ElPseudospecCtrl_s ctrl );
ElError ElPseudospectralWindowXDist_z
( ElConstDistMatrix_z A, ElDistMatrix_d invNormMap, 
  complex_double center, double realWidth, double imagWidth,
  ElInt realSize, ElInt imagSize, ElPseudospecCtrl_d ctrl );

/* Point cloud
   ----------- */
ElError ElPseudospectralCloud_s
( ElConstMatrix_s A, ElConstMatrix_c shifts, ElMatrix_s invNormMap );
ElError ElPseudospectralCloud_d
( ElConstMatrix_d A, ElConstMatrix_z shifts, ElMatrix_d invNormMap );
ElError ElPseudospectralCloud_c
( ElConstMatrix_c A, ElConstMatrix_c shifts, ElMatrix_s invNormMap );
ElError ElPseudospectralCloud_z
( ElConstMatrix_z A, ElConstMatrix_z shifts, ElMatrix_d invNormMap );

ElError ElPseudospectralCloudDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_c shifts, 
  ElDistMatrix_s invNormMap );
ElError ElPseudospectralCloudDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_z shifts, 
  ElDistMatrix_d invNormMap );
ElError ElPseudospectralCloudDist_c
( ElConstDistMatrix_c A, ElConstDistMatrix_c shifts, 
  ElDistMatrix_s invNormMap );
ElError ElPseudospectralCloudDist_z
( ElConstDistMatrix_z A, ElConstDistMatrix_z shifts, 
  ElDistMatrix_d invNormMap );

/* Expert version
   ^^^^^^^^^^^^^^ */
ElError ElPseudospectralCloudX_s
( ElConstMatrix_s A, ElConstMatrix_c shifts, ElMatrix_s invNormMap,
  ElPseudospecCtrl_s ctrl );
ElError ElPseudospectralCloudX_d
( ElConstMatrix_d A, ElConstMatrix_z shifts, ElMatrix_d invNormMap,
  ElPseudospecCtrl_d ctrl );
ElError ElPseudospectralCloudX_c
( ElConstMatrix_c A, ElConstMatrix_c shifts, ElMatrix_s invNormMap,
  ElPseudospecCtrl_s ctrl );
ElError ElPseudospectralCloudX_z
( ElConstMatrix_z A, ElConstMatrix_z shifts, ElMatrix_d invNormMap,
  ElPseudospecCtrl_d ctrl );

ElError ElPseudospectralCloudXDist_s
( ElConstDistMatrix_s A, ElConstDistMatrix_c shifts, 
  ElDistMatrix_s invNormMap, ElPseudospecCtrl_s ctrl );
ElError ElPseudospectralCloudXDist_d
( ElConstDistMatrix_d A, ElConstDistMatrix_z shifts, 
  ElDistMatrix_d invNormMap, ElPseudospecCtrl_d ctrl );
ElError ElPseudospectralCloudXDist_c
( ElConstDistMatrix_c A, ElConstDistMatrix_c shifts, 
  ElDistMatrix_s invNormMap, ElPseudospecCtrl_s ctrl );
ElError ElPseudospectralCloudXDist_z
( ElConstDistMatrix_z A, ElConstDistMatrix_z shifts, 
  ElDistMatrix_d invNormMap, ElPseudospecCtrl_d ctrl );

/* Trace
   ===== */
ElError ElTrace_i( ElConstMatrix_i A, ElInt* trace );
ElError ElTrace_s( ElConstMatrix_s A, float* trace );
ElError ElTrace_d( ElConstMatrix_d A, double* trace );
ElError ElTrace_c( ElConstMatrix_c A, complex_float* trace );
ElError ElTrace_z( ElConstMatrix_z A, complex_double* trace );

ElError ElTraceDist_i( ElConstDistMatrix_i A, ElInt* trace );
ElError ElTraceDist_s( ElConstDistMatrix_s A, float* trace );
ElError ElTraceDist_d( ElConstDistMatrix_d A, double* trace );
ElError ElTraceDist_c( ElConstDistMatrix_c A, complex_float* trace );
ElError ElTraceDist_z( ElConstDistMatrix_z A, complex_double* trace );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_LAPACK_PROPS_C_H */
