/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LAPACK_PROPS_C_H
#define EL_LAPACK_PROPS_C_H

#include "El/core/DistMatrix.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Condition number
   ================ */
EL_EXPORT ElError ElCondition_s
( ElConstMatrix_s A, ElNormType type, float* cond );
EL_EXPORT ElError ElCondition_d
( ElConstMatrix_d A, ElNormType type, double* cond );
EL_EXPORT ElError ElCondition_c
( ElConstMatrix_c A, ElNormType type, float* cond );
EL_EXPORT ElError ElCondition_z
( ElConstMatrix_z A, ElNormType type, double* cond );
EL_EXPORT ElError ElConditionDist_s
( ElConstDistMatrix_s A, ElNormType type, float* cond );
EL_EXPORT ElError ElConditionDist_d
( ElConstDistMatrix_d A, ElNormType type, double* cond );
EL_EXPORT ElError ElConditionDist_c
( ElConstDistMatrix_c A, ElNormType type, float* cond );
EL_EXPORT ElError ElConditionDist_z
( ElConstDistMatrix_z A, ElNormType type, double* cond );

/* Frobenius
   --------- */
EL_EXPORT ElError ElFrobeniusCondition_s( ElConstMatrix_s A, float* cond );
EL_EXPORT ElError ElFrobeniusCondition_d( ElConstMatrix_d A, double* cond );
EL_EXPORT ElError ElFrobeniusCondition_c( ElConstMatrix_c A, float* cond );
EL_EXPORT ElError ElFrobeniusCondition_z( ElConstMatrix_z A, double* cond );

EL_EXPORT ElError ElFrobeniusConditionDist_s
( ElConstDistMatrix_s A, float* cond );
EL_EXPORT ElError ElFrobeniusConditionDist_d
( ElConstDistMatrix_d A, double* cond );
EL_EXPORT ElError ElFrobeniusConditionDist_c
( ElConstDistMatrix_c A, float* cond );
EL_EXPORT ElError ElFrobeniusConditionDist_z
( ElConstDistMatrix_z A, double* cond );

/* Infinity
   -------- */
EL_EXPORT ElError ElInfinityCondition_s( ElConstMatrix_s A, float* cond );
EL_EXPORT ElError ElInfinityCondition_d( ElConstMatrix_d A, double* cond );
EL_EXPORT ElError ElInfinityCondition_c( ElConstMatrix_c A, float* cond );
EL_EXPORT ElError ElInfinityCondition_z( ElConstMatrix_z A, double* cond );
EL_EXPORT ElError ElInfinityConditionDist_s
( ElConstDistMatrix_s A, float* cond );
EL_EXPORT ElError ElInfinityConditionDist_d
( ElConstDistMatrix_d A, double* cond );
EL_EXPORT ElError ElInfinityConditionDist_c
( ElConstDistMatrix_c A, float* cond );
EL_EXPORT ElError ElInfinityConditionDist_z
( ElConstDistMatrix_z A, double* cond );

/* Max
   --- */
EL_EXPORT ElError ElMaxCondition_s( ElConstMatrix_s A, float* cond );
EL_EXPORT ElError ElMaxCondition_d( ElConstMatrix_d A, double* cond );
EL_EXPORT ElError ElMaxCondition_c( ElConstMatrix_c A, float* cond );
EL_EXPORT ElError ElMaxCondition_z( ElConstMatrix_z A, double* cond );
EL_EXPORT ElError ElMaxConditionDist_s( ElConstDistMatrix_s A, float* cond );
EL_EXPORT ElError ElMaxConditionDist_d( ElConstDistMatrix_d A, double* cond );
EL_EXPORT ElError ElMaxConditionDist_c( ElConstDistMatrix_c A, float* cond );
EL_EXPORT ElError ElMaxConditionDist_z( ElConstDistMatrix_z A, double* cond );

/* One
   --- */
EL_EXPORT ElError ElOneCondition_s( ElConstMatrix_s A, float* cond );
EL_EXPORT ElError ElOneCondition_d( ElConstMatrix_d A, double* cond );
EL_EXPORT ElError ElOneCondition_c( ElConstMatrix_c A, float* cond );
EL_EXPORT ElError ElOneCondition_z( ElConstMatrix_z A, double* cond );
EL_EXPORT ElError ElOneConditionDist_s( ElConstDistMatrix_s A, float* cond );
EL_EXPORT ElError ElOneConditionDist_d( ElConstDistMatrix_d A, double* cond );
EL_EXPORT ElError ElOneConditionDist_c( ElConstDistMatrix_c A, float* cond );
EL_EXPORT ElError ElOneConditionDist_z( ElConstDistMatrix_z A, double* cond );

/* Two
   --- */
EL_EXPORT ElError ElTwoCondition_s( ElConstMatrix_s A, float* cond );
EL_EXPORT ElError ElTwoCondition_d( ElConstMatrix_d A, double* cond );
EL_EXPORT ElError ElTwoCondition_c( ElConstMatrix_c A, float* cond );
EL_EXPORT ElError ElTwoCondition_z( ElConstMatrix_z A, double* cond );
EL_EXPORT ElError ElTwoConditionDist_s( ElConstDistMatrix_s A, float* cond );
EL_EXPORT ElError ElTwoConditionDist_d( ElConstDistMatrix_d A, double* cond );
EL_EXPORT ElError ElTwoConditionDist_c( ElConstDistMatrix_c A, float* cond );
EL_EXPORT ElError ElTwoConditionDist_z( ElConstDistMatrix_z A, double* cond );

/* Determinant
   =========== */
/* Return the result in a safer, expanded format
   --------------------------------------------- */
/* TODO: Version which allows for overwriting? How to name? */
EL_EXPORT ElError ElSafeDeterminant_s
( ElConstMatrix_s A, ElSafeProduct_s* prod );
EL_EXPORT ElError ElSafeDeterminant_d
( ElConstMatrix_d A, ElSafeProduct_d* prod );
EL_EXPORT ElError ElSafeDeterminant_c
( ElConstMatrix_c A, ElSafeProduct_c* prod );
EL_EXPORT ElError ElSafeDeterminant_z
( ElConstMatrix_z A, ElSafeProduct_z* prod );
EL_EXPORT ElError ElSafeDeterminantDist_s
( ElConstDistMatrix_s A, ElSafeProduct_s* prod );
EL_EXPORT ElError ElSafeDeterminantDist_d
( ElConstDistMatrix_d A, ElSafeProduct_d* prod );
EL_EXPORT ElError ElSafeDeterminantDist_c
( ElConstDistMatrix_c A, ElSafeProduct_c* prod );
EL_EXPORT ElError ElSafeDeterminantDist_z
( ElConstDistMatrix_z A, ElSafeProduct_z* prod );

/* TODO: Version which allows for overwriting? How to name? */
EL_EXPORT ElError ElSafeHPDDeterminant_s
( ElUpperOrLower uplo, ElConstMatrix_s A, ElSafeProduct_s* prod );
EL_EXPORT ElError ElSafeHPDDeterminant_d
( ElUpperOrLower uplo, ElConstMatrix_d A, ElSafeProduct_d* prod );
EL_EXPORT ElError ElSafeHPDDeterminant_c
( ElUpperOrLower uplo, ElConstMatrix_c A, ElSafeProduct_s* prod );
EL_EXPORT ElError ElSafeHPDDeterminant_z
( ElUpperOrLower uplo, ElConstMatrix_z A, ElSafeProduct_d* prod );
EL_EXPORT ElError ElSafeHPDDeterminantDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, ElSafeProduct_s* prod );
EL_EXPORT ElError ElSafeHPDDeterminantDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, ElSafeProduct_d* prod );
EL_EXPORT ElError ElSafeHPDDeterminantDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, ElSafeProduct_s* prod );
EL_EXPORT ElError ElSafeHPDDeterminantDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, ElSafeProduct_d* prod );

/* Return the direct result
   ------------------------ */
/* TODO: Version which allows for overwriting? How to name? */
EL_EXPORT ElError ElDeterminant_s( ElConstMatrix_s A, float* prod );
EL_EXPORT ElError ElDeterminant_d( ElConstMatrix_d A, double* prod );
EL_EXPORT ElError ElDeterminant_c( ElConstMatrix_c A, complex_float* prod );
EL_EXPORT ElError ElDeterminant_z( ElConstMatrix_z A, complex_double* prod );
EL_EXPORT ElError ElDeterminantDist_s
( ElConstDistMatrix_s A, float* prod );
EL_EXPORT ElError ElDeterminantDist_d
( ElConstDistMatrix_d A, double* prod );
EL_EXPORT ElError ElDeterminantDist_c
( ElConstDistMatrix_c A, complex_float* prod );
EL_EXPORT ElError ElDeterminantDist_z
( ElConstDistMatrix_z A, complex_double* prod );

/* TODO: Version which allows for overwriting? How to name? */
EL_EXPORT ElError ElHPDDeterminant_s
( ElUpperOrLower uplo, ElConstMatrix_s A, float* prod );
EL_EXPORT ElError ElHPDDeterminant_d
( ElUpperOrLower uplo, ElConstMatrix_d A, double* prod );
EL_EXPORT ElError ElHPDDeterminant_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float* prod );
EL_EXPORT ElError ElHPDDeterminant_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double* prod );
EL_EXPORT ElError ElHPDDeterminantDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, float* prod );
EL_EXPORT ElError ElHPDDeterminantDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, double* prod );
EL_EXPORT ElError ElHPDDeterminantDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float* prod );
EL_EXPORT ElError ElHPDDeterminantDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double* prod );

/* TODO: Determinant after Cholesky or LU */

/* Inertia
   ======= */
EL_EXPORT ElError ElInertia_s
( ElUpperOrLower uplo, ElMatrix_s A, ElInertiaType* inertia );
EL_EXPORT ElError ElInertia_d
( ElUpperOrLower uplo, ElMatrix_d A, ElInertiaType* inertia );
EL_EXPORT ElError ElInertia_c
( ElUpperOrLower uplo, ElMatrix_c A, ElInertiaType* inertia );
EL_EXPORT ElError ElInertia_z
( ElUpperOrLower uplo, ElMatrix_z A, ElInertiaType* inertia );
EL_EXPORT ElError ElInertiaDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A, ElInertiaType* inertia );
EL_EXPORT ElError ElInertiaDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A, ElInertiaType* inertia );
EL_EXPORT ElError ElInertiaDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A, ElInertiaType* inertia );
EL_EXPORT ElError ElInertiaDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A, ElInertiaType* inertia );

/* TODO: Expert versions */

/* Norm
   ==== */
EL_EXPORT ElError ElNorm_s
( ElConstMatrix_s A, ElNormType normType, float* norm );
EL_EXPORT ElError ElNorm_d
( ElConstMatrix_d A, ElNormType normType, double* norm );
EL_EXPORT ElError ElNorm_c
( ElConstMatrix_c A, ElNormType normType, float* norm );
EL_EXPORT ElError ElNorm_z
( ElConstMatrix_z A, ElNormType normType, double* norm );
EL_EXPORT ElError ElNormDist_s
( ElConstDistMatrix_s A, ElNormType normType, float* norm );
EL_EXPORT ElError ElNormDist_d
( ElConstDistMatrix_d A, ElNormType normType, double* norm );
EL_EXPORT ElError ElNormDist_c
( ElConstDistMatrix_c A, ElNormType normType, float* norm );
EL_EXPORT ElError ElNormDist_z
( ElConstDistMatrix_z A, ElNormType normType, double* norm );

EL_EXPORT ElError ElHermitianNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, ElNormType normType, float* norm );
EL_EXPORT ElError ElHermitianNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, ElNormType normType, double* norm );
EL_EXPORT ElError ElHermitianNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, ElNormType normType, 
  float* norm );
EL_EXPORT ElError ElHermitianNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, ElNormType normType, 
  double* norm );

EL_EXPORT ElError ElSymmetricNorm_s
( ElUpperOrLower uplo, ElConstMatrix_s A, ElNormType normType, float* norm );
EL_EXPORT ElError ElSymmetricNorm_d
( ElUpperOrLower uplo, ElConstMatrix_d A, ElNormType normType, double* norm );
EL_EXPORT ElError ElSymmetricNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, ElNormType normType, float* norm );
EL_EXPORT ElError ElSymmetricNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, ElNormType normType, double* norm );
EL_EXPORT ElError ElSymmetricNormDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, ElNormType normType, 
  float* norm );
EL_EXPORT ElError ElSymmetricNormDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, ElNormType normType, 
  double* norm );
EL_EXPORT ElError ElSymmetricNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, ElNormType normType, 
  float* norm );
EL_EXPORT ElError ElSymmetricNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, ElNormType normType, 
  double* norm );

/* Entrywise norm
   -------------- */
EL_EXPORT ElError ElEntrywiseNorm_s
( ElConstMatrix_s A, float p, float* norm );
EL_EXPORT ElError ElEntrywiseNorm_d
( ElConstMatrix_d A, double p, double* norm );
EL_EXPORT ElError ElEntrywiseNorm_c
( ElConstMatrix_c A, float p, float* norm );
EL_EXPORT ElError ElEntrywiseNorm_z
( ElConstMatrix_z A, double p, double* norm );

EL_EXPORT ElError ElEntrywiseNormDist_s
( ElConstDistMatrix_s A, float p, float* norm );
EL_EXPORT ElError ElEntrywiseNormDist_d
( ElConstDistMatrix_d A, double p, double* norm );
EL_EXPORT ElError ElEntrywiseNormDist_c
( ElConstDistMatrix_c A, float p, float* norm );
EL_EXPORT ElError ElEntrywiseNormDist_z
( ElConstDistMatrix_z A, double p, double* norm );

EL_EXPORT ElError ElEntrywiseNormSparse_s
( ElConstSparseMatrix_s A, float p, float* norm );
EL_EXPORT ElError ElEntrywiseNormSparse_d
( ElConstSparseMatrix_d A, double p, double* norm );
EL_EXPORT ElError ElEntrywiseNormSparse_c
( ElConstSparseMatrix_c A, float p, float* norm );
EL_EXPORT ElError ElEntrywiseNormSparse_z
( ElConstSparseMatrix_z A, double p, double* norm );

EL_EXPORT ElError ElEntrywiseNormDistSparse_s
( ElConstDistSparseMatrix_s A, float p, float* norm );
EL_EXPORT ElError ElEntrywiseNormDistSparse_d
( ElConstDistSparseMatrix_d A, double p, double* norm );
EL_EXPORT ElError ElEntrywiseNormDistSparse_c
( ElConstDistSparseMatrix_c A, float p, float* norm );
EL_EXPORT ElError ElEntrywiseNormDistSparse_z
( ElConstDistSparseMatrix_z A, double p, double* norm );

EL_EXPORT ElError ElEntrywiseNormDistMultiVec_s
( ElConstDistMultiVec_s A, float p, float* norm );
EL_EXPORT ElError ElEntrywiseNormDistMultiVec_d
( ElConstDistMultiVec_d A, double p, double* norm );
EL_EXPORT ElError ElEntrywiseNormDistMultiVec_c
( ElConstDistMultiVec_c A, float p, float* norm );
EL_EXPORT ElError ElEntrywiseNormDistMultiVec_z
( ElConstDistMultiVec_z A, double p, double* norm );

EL_EXPORT ElError ElHermitianEntrywiseNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float p, float* norm );
EL_EXPORT ElError ElHermitianEntrywiseNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double p, double* norm );

EL_EXPORT ElError ElHermitianEntrywiseNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float p, float* norm );
EL_EXPORT ElError ElHermitianEntrywiseNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double p, double* norm );

EL_EXPORT ElError ElHermitianEntrywiseNormSparse_c
( ElUpperOrLower uplo, ElConstSparseMatrix_c A, float p, float* norm );
EL_EXPORT ElError ElHermitianEntrywiseNormSparse_z
( ElUpperOrLower uplo, ElConstSparseMatrix_z A, double p, double* norm );

EL_EXPORT ElError ElHermitianEntrywiseNormDistSparse_c
( ElUpperOrLower uplo, ElConstDistSparseMatrix_c A, float p, float* norm );
EL_EXPORT ElError ElHermitianEntrywiseNormDistSparse_z
( ElUpperOrLower uplo, ElConstDistSparseMatrix_z A, double p, double* norm );

EL_EXPORT ElError ElSymmetricEntrywiseNorm_s
( ElUpperOrLower uplo, ElConstMatrix_s A, float p, float* norm );
EL_EXPORT ElError ElSymmetricEntrywiseNorm_d
( ElUpperOrLower uplo, ElConstMatrix_d A, double p, double* norm );
EL_EXPORT ElError ElSymmetricEntrywiseNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float p, float* norm );
EL_EXPORT ElError ElSymmetricEntrywiseNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double p, double* norm );

EL_EXPORT ElError ElSymmetricEntrywiseNormDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, float p, float* norm );
EL_EXPORT ElError ElSymmetricEntrywiseNormDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, double p, double* norm );
EL_EXPORT ElError ElSymmetricEntrywiseNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float p, float* norm );
EL_EXPORT ElError ElSymmetricEntrywiseNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double p, double* norm );

EL_EXPORT ElError ElSymmetricEntrywiseNormSparse_s
( ElUpperOrLower uplo, ElConstSparseMatrix_s A, float p, float* norm );
EL_EXPORT ElError ElSymmetricEntrywiseNormSparse_d
( ElUpperOrLower uplo, ElConstSparseMatrix_d A, double p, double* norm );
EL_EXPORT ElError ElSymmetricEntrywiseNormSparse_c
( ElUpperOrLower uplo, ElConstSparseMatrix_c A, float p, float* norm );
EL_EXPORT ElError ElSymmetricEntrywiseNormSparse_z
( ElUpperOrLower uplo, ElConstSparseMatrix_z A, double p, double* norm );

EL_EXPORT ElError ElSymmetricEntrywiseNormDistSparse_s
( ElUpperOrLower uplo, ElConstDistSparseMatrix_s A, float p, float* norm );
EL_EXPORT ElError ElSymmetricEntrywiseNormDistSparse_d
( ElUpperOrLower uplo, ElConstDistSparseMatrix_d A, double p, double* norm );
EL_EXPORT ElError ElSymmetricEntrywiseNormDistSparse_c
( ElUpperOrLower uplo, ElConstDistSparseMatrix_c A, float p, float* norm );
EL_EXPORT ElError ElSymmetricEntrywiseNormDistSparse_z
( ElUpperOrLower uplo, ElConstDistSparseMatrix_z A, double p, double* norm );

/* Frobenius norm
   -------------- */
EL_EXPORT ElError ElFrobeniusNorm_s( ElConstMatrix_s A, float* norm );
EL_EXPORT ElError ElFrobeniusNorm_d( ElConstMatrix_d A, double* norm );
EL_EXPORT ElError ElFrobeniusNorm_c( ElConstMatrix_c A, float* norm );
EL_EXPORT ElError ElFrobeniusNorm_z( ElConstMatrix_z A, double* norm );
EL_EXPORT ElError ElFrobeniusNormDist_s( ElConstDistMatrix_s A, float* norm );
EL_EXPORT ElError ElFrobeniusNormDist_d( ElConstDistMatrix_d A, double* norm );
EL_EXPORT ElError ElFrobeniusNormDist_c( ElConstDistMatrix_c A, float* norm );
EL_EXPORT ElError ElFrobeniusNormDist_z( ElConstDistMatrix_z A, double* norm );
EL_EXPORT ElError ElFrobeniusNormSparse_s
( ElConstSparseMatrix_s A, float* norm );
EL_EXPORT ElError ElFrobeniusNormSparse_d
( ElConstSparseMatrix_d A, double* norm );
EL_EXPORT ElError ElFrobeniusNormSparse_c
( ElConstSparseMatrix_c A, float* norm );
EL_EXPORT ElError ElFrobeniusNormSparse_z
( ElConstSparseMatrix_z A, double* norm );
EL_EXPORT ElError ElFrobeniusNormDistSparse_s
( ElConstDistSparseMatrix_s A, float* norm );
EL_EXPORT ElError ElFrobeniusNormDistSparse_d
( ElConstDistSparseMatrix_d A, double* norm );
EL_EXPORT ElError ElFrobeniusNormDistSparse_c
( ElConstDistSparseMatrix_c A, float* norm );
EL_EXPORT ElError ElFrobeniusNormDistSparse_z
( ElConstDistSparseMatrix_z A, double* norm );
EL_EXPORT ElError ElFrobeniusNormDistMultiVec_s
( ElConstDistMultiVec_s A, float* norm );
EL_EXPORT ElError ElFrobeniusNormDistMultiVec_d
( ElConstDistMultiVec_d A, double* norm );
EL_EXPORT ElError ElFrobeniusNormDistMultiVec_c
( ElConstDistMultiVec_c A, float* norm );
EL_EXPORT ElError ElFrobeniusNormDistMultiVec_z
( ElConstDistMultiVec_z A, double* norm );

EL_EXPORT ElError ElHermitianFrobeniusNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float* norm );
EL_EXPORT ElError ElHermitianFrobeniusNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double* norm );
EL_EXPORT ElError ElHermitianFrobeniusNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float* norm );
EL_EXPORT ElError ElHermitianFrobeniusNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double* norm );
EL_EXPORT ElError ElHermitianFrobeniusNormSparse_c
( ElUpperOrLower uplo, ElConstSparseMatrix_c A, float* norm );
EL_EXPORT ElError ElHermitianFrobeniusNormSparse_z
( ElUpperOrLower uplo, ElConstSparseMatrix_z A, double* norm );
EL_EXPORT ElError ElHermitianFrobeniusNormDistSparse_c
( ElUpperOrLower uplo, ElConstDistSparseMatrix_c A, float* norm );
EL_EXPORT ElError ElHermitianFrobeniusNormDistSparse_z
( ElUpperOrLower uplo, ElConstDistSparseMatrix_z A, double* norm );

EL_EXPORT ElError ElSymmetricFrobeniusNorm_s
( ElUpperOrLower uplo, ElConstMatrix_s A, float* norm );
EL_EXPORT ElError ElSymmetricFrobeniusNorm_d
( ElUpperOrLower uplo, ElConstMatrix_d A, double* norm );
EL_EXPORT ElError ElSymmetricFrobeniusNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float* norm );
EL_EXPORT ElError ElSymmetricFrobeniusNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double* norm );
EL_EXPORT ElError ElSymmetricFrobeniusNormDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, float* norm );
EL_EXPORT ElError ElSymmetricFrobeniusNormDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, double* norm );
EL_EXPORT ElError ElSymmetricFrobeniusNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float* norm );
EL_EXPORT ElError ElSymmetricFrobeniusNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double* norm );
EL_EXPORT ElError ElSymmetricFrobeniusNormSparse_s
( ElUpperOrLower uplo, ElConstSparseMatrix_s A, float* norm );
EL_EXPORT ElError ElSymmetricFrobeniusNormSparse_d
( ElUpperOrLower uplo, ElConstSparseMatrix_d A, double* norm );
EL_EXPORT ElError ElSymmetricFrobeniusNormSparse_c
( ElUpperOrLower uplo, ElConstSparseMatrix_c A, float* norm );
EL_EXPORT ElError ElSymmetricFrobeniusNormSparse_z
( ElUpperOrLower uplo, ElConstSparseMatrix_z A, double* norm );
EL_EXPORT ElError ElSymmetricFrobeniusNormDistSparse_s
( ElUpperOrLower uplo, ElConstDistSparseMatrix_s A, float* norm );
EL_EXPORT ElError ElSymmetricFrobeniusNormDistSparse_d
( ElUpperOrLower uplo, ElConstDistSparseMatrix_d A, double* norm );
EL_EXPORT ElError ElSymmetricFrobeniusNormDistSparse_c
( ElUpperOrLower uplo, ElConstDistSparseMatrix_c A, float* norm );
EL_EXPORT ElError ElSymmetricFrobeniusNormDistSparse_z
( ElUpperOrLower uplo, ElConstDistSparseMatrix_z A, double* norm );

/* Infinity norm
   ------------- */
EL_EXPORT ElError ElInfinityNorm_s( ElConstMatrix_s A, float* norm );
EL_EXPORT ElError ElInfinityNorm_d( ElConstMatrix_d A, double* norm );
EL_EXPORT ElError ElInfinityNorm_c( ElConstMatrix_c A, float* norm );
EL_EXPORT ElError ElInfinityNorm_z( ElConstMatrix_z A, double* norm );
EL_EXPORT ElError ElInfinityNormDist_s( ElConstDistMatrix_s A, float* norm );
EL_EXPORT ElError ElInfinityNormDist_d( ElConstDistMatrix_d A, double* norm );
EL_EXPORT ElError ElInfinityNormDist_c( ElConstDistMatrix_c A, float* norm );
EL_EXPORT ElError ElInfinityNormDist_z( ElConstDistMatrix_z A, double* norm );
EL_EXPORT ElError ElInfinityNormSparse_s
( ElConstSparseMatrix_s A, float* norm );
EL_EXPORT ElError ElInfinityNormSparse_d
( ElConstSparseMatrix_d A, double* norm );
EL_EXPORT ElError ElInfinityNormSparse_c
( ElConstSparseMatrix_c A, float* norm );
EL_EXPORT ElError ElInfinityNormSparse_z
( ElConstSparseMatrix_z A, double* norm );
EL_EXPORT ElError ElInfinityNormDistSparse_s
( ElConstDistSparseMatrix_s A, float* norm );
EL_EXPORT ElError ElInfinityNormDistSparse_d
( ElConstDistSparseMatrix_d A, double* norm );
EL_EXPORT ElError ElInfinityNormDistSparse_c
( ElConstDistSparseMatrix_c A, float* norm );
EL_EXPORT ElError ElInfinityNormDistSparse_z
( ElConstDistSparseMatrix_z A, double* norm );

EL_EXPORT ElError ElHermitianInfinityNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float* norm );
EL_EXPORT ElError ElHermitianInfinityNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double* norm );
EL_EXPORT ElError ElHermitianInfinityNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float* norm );
EL_EXPORT ElError ElHermitianInfinityNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double* norm );

EL_EXPORT ElError ElSymmetricInfinityNorm_s
( ElUpperOrLower uplo, ElConstMatrix_s A, float* norm );
EL_EXPORT ElError ElSymmetricInfinityNorm_d
( ElUpperOrLower uplo, ElConstMatrix_d A, double* norm );
EL_EXPORT ElError ElSymmetricInfinityNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float* norm );
EL_EXPORT ElError ElSymmetricInfinityNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double* norm );
EL_EXPORT ElError ElSymmetricInfinityNormDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, float* norm );
EL_EXPORT ElError ElSymmetricInfinityNormDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, double* norm );
EL_EXPORT ElError ElSymmetricInfinityNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float* norm );
EL_EXPORT ElError ElSymmetricInfinityNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double* norm );

/* Ky-Fan norm
   ----------- */
EL_EXPORT ElError ElKyFanNorm_s( ElConstMatrix_s A, ElInt k, float* norm );
EL_EXPORT ElError ElKyFanNorm_d( ElConstMatrix_d A, ElInt k, double* norm );
EL_EXPORT ElError ElKyFanNorm_c( ElConstMatrix_c A, ElInt k, float* norm );
EL_EXPORT ElError ElKyFanNorm_z( ElConstMatrix_z A, ElInt k, double* norm );
EL_EXPORT ElError ElKyFanNormDist_s
( ElConstDistMatrix_s A, ElInt k, float* norm );
EL_EXPORT ElError ElKyFanNormDist_d
( ElConstDistMatrix_d A, ElInt k, double* norm );
EL_EXPORT ElError ElKyFanNormDist_c
( ElConstDistMatrix_c A, ElInt k, float* norm );
EL_EXPORT ElError ElKyFanNormDist_z
( ElConstDistMatrix_z A, ElInt k, double* norm );

EL_EXPORT ElError ElHermitianKyFanNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, ElInt k, float* norm );
EL_EXPORT ElError ElHermitianKyFanNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, ElInt k, double* norm );
EL_EXPORT ElError ElHermitianKyFanNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, ElInt k, float* norm );
EL_EXPORT ElError ElHermitianKyFanNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, ElInt k, double* norm );

EL_EXPORT ElError ElSymmetricKyFanNorm_s
( ElUpperOrLower uplo, ElConstMatrix_s A, ElInt k, float* norm );
EL_EXPORT ElError ElSymmetricKyFanNorm_d
( ElUpperOrLower uplo, ElConstMatrix_d A, ElInt k, double* norm );
EL_EXPORT ElError ElSymmetricKyFanNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, ElInt k, float* norm );
EL_EXPORT ElError ElSymmetricKyFanNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, ElInt k, double* norm );
EL_EXPORT ElError ElSymmetricKyFanNormDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, ElInt k, float* norm );
EL_EXPORT ElError ElSymmetricKyFanNormDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, ElInt k, double* norm );
EL_EXPORT ElError ElSymmetricKyFanNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, ElInt k, float* norm );
EL_EXPORT ElError ElSymmetricKyFanNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, ElInt k, double* norm );

/* Ky-Fan-Schatten norm
   -------------------- */
EL_EXPORT ElError ElKyFanSchattenNorm_s
( ElConstMatrix_s A, ElInt k, float p, float* norm );
EL_EXPORT ElError ElKyFanSchattenNorm_d
( ElConstMatrix_d A, ElInt k, double p, double* norm );
EL_EXPORT ElError ElKyFanSchattenNorm_c
( ElConstMatrix_c A, ElInt k, float p, float* norm );
EL_EXPORT ElError ElKyFanSchattenNorm_z
( ElConstMatrix_z A, ElInt k, double p, double* norm );
EL_EXPORT ElError ElKyFanSchattenNormDist_s
( ElConstDistMatrix_s A, ElInt k, float p, float* norm );
EL_EXPORT ElError ElKyFanSchattenNormDist_d
( ElConstDistMatrix_d A, ElInt k, double p, double* norm );
EL_EXPORT ElError ElKyFanSchattenNormDist_c
( ElConstDistMatrix_c A, ElInt k, float p, float* norm );
EL_EXPORT ElError ElKyFanSchattenNormDist_z
( ElConstDistMatrix_z A, ElInt k, double p, double* norm );

EL_EXPORT ElError ElHermitianKyFanSchattenNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, ElInt k, float p, float* norm );
EL_EXPORT ElError ElHermitianKyFanSchattenNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, ElInt k, double p, double* norm );
EL_EXPORT ElError ElHermitianKyFanSchattenNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, ElInt k, float p, float* norm );
EL_EXPORT ElError ElHermitianKyFanSchattenNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, ElInt k, double p, double* norm );

EL_EXPORT ElError ElSymmetricKyFanSchattenNorm_s
( ElUpperOrLower uplo, ElConstMatrix_s A, ElInt k, float p, float* norm );
EL_EXPORT ElError ElSymmetricKyFanSchattenNorm_d
( ElUpperOrLower uplo, ElConstMatrix_d A, ElInt k, double p, double* norm );
EL_EXPORT ElError ElSymmetricKyFanSchattenNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, ElInt k, float p, float* norm );
EL_EXPORT ElError ElSymmetricKyFanSchattenNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, ElInt k, double p, double* norm );
EL_EXPORT ElError ElSymmetricKyFanSchattenNormDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, ElInt k, float p, float* norm );
EL_EXPORT ElError ElSymmetricKyFanSchattenNormDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, ElInt k, double p, double* norm );
EL_EXPORT ElError ElSymmetricKyFanSchattenNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, ElInt k, float p, float* norm );
EL_EXPORT ElError ElSymmetricKyFanSchattenNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, ElInt k, double p, double* norm );

/* Max norm
   -------- */
EL_EXPORT ElError ElMaxNorm_i( ElConstMatrix_i A, ElInt* norm );
EL_EXPORT ElError ElMaxNorm_s( ElConstMatrix_s A, float* norm );
EL_EXPORT ElError ElMaxNorm_d( ElConstMatrix_d A, double* norm );
EL_EXPORT ElError ElMaxNorm_c( ElConstMatrix_c A, float* norm );
EL_EXPORT ElError ElMaxNorm_z( ElConstMatrix_z A, double* norm );
EL_EXPORT ElError ElMaxNormDist_i( ElConstDistMatrix_i A, ElInt* norm );
EL_EXPORT ElError ElMaxNormDist_s( ElConstDistMatrix_s A, float* norm );
EL_EXPORT ElError ElMaxNormDist_d( ElConstDistMatrix_d A, double* norm );
EL_EXPORT ElError ElMaxNormDist_c( ElConstDistMatrix_c A, float* norm );
EL_EXPORT ElError ElMaxNormDist_z( ElConstDistMatrix_z A, double* norm );
EL_EXPORT ElError ElMaxNormSparse_i( ElConstSparseMatrix_i A, ElInt* norm );
EL_EXPORT ElError ElMaxNormSparse_s( ElConstSparseMatrix_s A, float* norm );
EL_EXPORT ElError ElMaxNormSparse_d( ElConstSparseMatrix_d A, double* norm );
EL_EXPORT ElError ElMaxNormSparse_c( ElConstSparseMatrix_c A, float* norm );
EL_EXPORT ElError ElMaxNormSparse_z( ElConstSparseMatrix_z A, double* norm );
EL_EXPORT ElError ElMaxNormDistSparse_i
( ElConstDistSparseMatrix_i A, ElInt* norm );
EL_EXPORT ElError ElMaxNormDistSparse_s
( ElConstDistSparseMatrix_s A, float* norm );
EL_EXPORT ElError ElMaxNormDistSparse_d
( ElConstDistSparseMatrix_d A, double* norm );
EL_EXPORT ElError ElMaxNormDistSparse_c
( ElConstDistSparseMatrix_c A, float* norm );
EL_EXPORT ElError ElMaxNormDistSparse_z
( ElConstDistSparseMatrix_z A, double* norm );
EL_EXPORT ElError ElMaxNormDistMultiVec_i
( ElConstDistMultiVec_i A, ElInt* norm );
EL_EXPORT ElError ElMaxNormDistMultiVec_s
( ElConstDistMultiVec_s A, float* norm );
EL_EXPORT ElError ElMaxNormDistMultiVec_d
( ElConstDistMultiVec_d A, double* norm );
EL_EXPORT ElError ElMaxNormDistMultiVec_c
( ElConstDistMultiVec_c A, float* norm );
EL_EXPORT ElError ElMaxNormDistMultiVec_z
( ElConstDistMultiVec_z A, double* norm );

EL_EXPORT ElError ElHermitianMaxNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float* norm );
EL_EXPORT ElError ElHermitianMaxNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double* norm );
EL_EXPORT ElError ElHermitianMaxNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float* norm );
EL_EXPORT ElError ElHermitianMaxNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double* norm );
EL_EXPORT ElError ElHermitianMaxNormSparse_c
( ElUpperOrLower uplo, ElConstSparseMatrix_c A, float* norm );
EL_EXPORT ElError ElHermitianMaxNormSparse_z
( ElUpperOrLower uplo, ElConstSparseMatrix_z A, double* norm );
EL_EXPORT ElError ElHermitianMaxNormDistSparse_c
( ElUpperOrLower uplo, ElConstDistSparseMatrix_c A, float* norm );
EL_EXPORT ElError ElHermitianMaxNormDistSparse_z
( ElUpperOrLower uplo, ElConstDistSparseMatrix_z A, double* norm );

EL_EXPORT ElError ElSymmetricMaxNorm_i
( ElUpperOrLower uplo, ElConstMatrix_i A, ElInt* norm );
EL_EXPORT ElError ElSymmetricMaxNorm_s
( ElUpperOrLower uplo, ElConstMatrix_s A, float* norm );
EL_EXPORT ElError ElSymmetricMaxNorm_d
( ElUpperOrLower uplo, ElConstMatrix_d A, double* norm );
EL_EXPORT ElError ElSymmetricMaxNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float* norm );
EL_EXPORT ElError ElSymmetricMaxNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double* norm );
EL_EXPORT ElError ElSymmetricMaxNormDist_i
( ElUpperOrLower uplo, ElConstDistMatrix_i A, ElInt* norm );
EL_EXPORT ElError ElSymmetricMaxNormDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, float* norm );
EL_EXPORT ElError ElSymmetricMaxNormDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, double* norm );
EL_EXPORT ElError ElSymmetricMaxNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float* norm );
EL_EXPORT ElError ElSymmetricMaxNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double* norm );
EL_EXPORT ElError ElSymmetricMaxNormSparse_i
( ElUpperOrLower uplo, ElConstSparseMatrix_i A, ElInt* norm );
EL_EXPORT ElError ElSymmetricMaxNormSparse_s
( ElUpperOrLower uplo, ElConstSparseMatrix_s A, float* norm );
EL_EXPORT ElError ElSymmetricMaxNormSparse_d
( ElUpperOrLower uplo, ElConstSparseMatrix_d A, double* norm );
EL_EXPORT ElError ElSymmetricMaxNormSparse_c
( ElUpperOrLower uplo, ElConstSparseMatrix_c A, float* norm );
EL_EXPORT ElError ElSymmetricMaxNormSparse_z
( ElUpperOrLower uplo, ElConstSparseMatrix_z A, double* norm );
EL_EXPORT ElError ElSymmetricMaxNormDistSparse_i
( ElUpperOrLower uplo, ElConstDistSparseMatrix_i A, ElInt* norm );
EL_EXPORT ElError ElSymmetricMaxNormDistSparse_s
( ElUpperOrLower uplo, ElConstDistSparseMatrix_s A, float* norm );
EL_EXPORT ElError ElSymmetricMaxNormDistSparse_d
( ElUpperOrLower uplo, ElConstDistSparseMatrix_d A, double* norm );
EL_EXPORT ElError ElSymmetricMaxNormDistSparse_c
( ElUpperOrLower uplo, ElConstDistSparseMatrix_c A, float* norm );
EL_EXPORT ElError ElSymmetricMaxNormDistSparse_z
( ElUpperOrLower uplo, ElConstDistSparseMatrix_z A, double* norm );

/* Nuclear norm
   ------------ */
EL_EXPORT ElError ElNuclearNorm_s( ElConstMatrix_s A, float* norm );
EL_EXPORT ElError ElNuclearNorm_d( ElConstMatrix_d A, double* norm );
EL_EXPORT ElError ElNuclearNorm_c( ElConstMatrix_c A, float* norm );
EL_EXPORT ElError ElNuclearNorm_z( ElConstMatrix_z A, double* norm );
EL_EXPORT ElError ElNuclearNormDist_s( ElConstDistMatrix_s A, float* norm );
EL_EXPORT ElError ElNuclearNormDist_d( ElConstDistMatrix_d A, double* norm );
EL_EXPORT ElError ElNuclearNormDist_c( ElConstDistMatrix_c A, float* norm );
EL_EXPORT ElError ElNuclearNormDist_z( ElConstDistMatrix_z A, double* norm );

EL_EXPORT ElError ElHermitianNuclearNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float* norm );
EL_EXPORT ElError ElHermitianNuclearNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double* norm );
EL_EXPORT ElError ElHermitianNuclearNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float* norm );
EL_EXPORT ElError ElHermitianNuclearNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double* norm );

EL_EXPORT ElError ElSymmetricNuclearNorm_s
( ElUpperOrLower uplo, ElConstMatrix_s A, float* norm );
EL_EXPORT ElError ElSymmetricNuclearNorm_d
( ElUpperOrLower uplo, ElConstMatrix_d A, double* norm );
EL_EXPORT ElError ElSymmetricNuclearNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float* norm );
EL_EXPORT ElError ElSymmetricNuclearNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double* norm );
EL_EXPORT ElError ElSymmetricNuclearNormDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, float* norm );
EL_EXPORT ElError ElSymmetricNuclearNormDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, double* norm );
EL_EXPORT ElError ElSymmetricNuclearNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float* norm );
EL_EXPORT ElError ElSymmetricNuclearNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double* norm );

/* One norm
   -------- */
EL_EXPORT ElError ElOneNorm_s( ElConstMatrix_s A, float* norm );
EL_EXPORT ElError ElOneNorm_d( ElConstMatrix_d A, double* norm );
EL_EXPORT ElError ElOneNorm_c( ElConstMatrix_c A, float* norm );
EL_EXPORT ElError ElOneNorm_z( ElConstMatrix_z A, double* norm );
EL_EXPORT ElError ElOneNormDist_s( ElConstDistMatrix_s A, float* norm );
EL_EXPORT ElError ElOneNormDist_d( ElConstDistMatrix_d A, double* norm );
EL_EXPORT ElError ElOneNormDist_c( ElConstDistMatrix_c A, float* norm );
EL_EXPORT ElError ElOneNormDist_z( ElConstDistMatrix_z A, double* norm );
EL_EXPORT ElError ElOneNormSparse_s( ElConstSparseMatrix_s A, float* norm );
EL_EXPORT ElError ElOneNormSparse_d( ElConstSparseMatrix_d A, double* norm );
EL_EXPORT ElError ElOneNormSparse_c( ElConstSparseMatrix_c A, float* norm );
EL_EXPORT ElError ElOneNormSparse_z( ElConstSparseMatrix_z A, double* norm );
EL_EXPORT ElError ElOneNormDistSparse_s
( ElConstDistSparseMatrix_s A, float* norm );
EL_EXPORT ElError ElOneNormDistSparse_d
( ElConstDistSparseMatrix_d A, double* norm );
EL_EXPORT ElError ElOneNormDistSparse_c
( ElConstDistSparseMatrix_c A, float* norm );
EL_EXPORT ElError ElOneNormDistSparse_z
( ElConstDistSparseMatrix_z A, double* norm );

EL_EXPORT ElError ElHermitianOneNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float* norm );
EL_EXPORT ElError ElHermitianOneNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double* norm );
EL_EXPORT ElError ElHermitianOneNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float* norm );
EL_EXPORT ElError ElHermitianOneNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double* norm );

EL_EXPORT ElError ElSymmetricOneNorm_s
( ElUpperOrLower uplo, ElConstMatrix_s A, float* norm );
EL_EXPORT ElError ElSymmetricOneNorm_d
( ElUpperOrLower uplo, ElConstMatrix_d A, double* norm );
EL_EXPORT ElError ElSymmetricOneNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float* norm );
EL_EXPORT ElError ElSymmetricOneNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double* norm );
EL_EXPORT ElError ElSymmetricOneNormDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, float* norm );
EL_EXPORT ElError ElSymmetricOneNormDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, double* norm );
EL_EXPORT ElError ElSymmetricOneNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float* norm );
EL_EXPORT ElError ElSymmetricOneNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double* norm );

/* Schatten norm
   ------------- */
EL_EXPORT ElError ElSchattenNorm_s
( ElConstMatrix_s A, float p, float* norm );
EL_EXPORT ElError ElSchattenNorm_d
( ElConstMatrix_d A, double p, double* norm );
EL_EXPORT ElError ElSchattenNorm_c
( ElConstMatrix_c A, float p, float* norm );
EL_EXPORT ElError ElSchattenNorm_z
( ElConstMatrix_z A, double p, double* norm );
EL_EXPORT ElError ElSchattenNormDist_s
( ElConstDistMatrix_s A, float p, float* norm );
EL_EXPORT ElError ElSchattenNormDist_d
( ElConstDistMatrix_d A, double p, double* norm );
EL_EXPORT ElError ElSchattenNormDist_c
( ElConstDistMatrix_c A, float p, float* norm );
EL_EXPORT ElError ElSchattenNormDist_z
( ElConstDistMatrix_z A, double p, double* norm );

EL_EXPORT ElError ElHermitianSchattenNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float p, float* norm );
EL_EXPORT ElError ElHermitianSchattenNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double p, double* norm );
EL_EXPORT ElError ElHermitianSchattenNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float p, float* norm );
EL_EXPORT ElError ElHermitianSchattenNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double p, double* norm );

EL_EXPORT ElError ElSymmetricSchattenNorm_s
( ElUpperOrLower uplo, ElConstMatrix_s A, float p, float* norm );
EL_EXPORT ElError ElSymmetricSchattenNorm_d
( ElUpperOrLower uplo, ElConstMatrix_d A, double p, double* norm );
EL_EXPORT ElError ElSymmetricSchattenNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float p, float* norm );
EL_EXPORT ElError ElSymmetricSchattenNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double p, double* norm );
EL_EXPORT ElError ElSymmetricSchattenNormDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, float p, float* norm );
EL_EXPORT ElError ElSymmetricSchattenNormDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, double p, double* norm );
EL_EXPORT ElError ElSymmetricSchattenNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float p, float* norm );
EL_EXPORT ElError ElSymmetricSchattenNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double p, double* norm );

/* Two norm
   -------- */
EL_EXPORT ElError ElTwoNorm_s( ElConstMatrix_s A, float* norm );
EL_EXPORT ElError ElTwoNorm_d( ElConstMatrix_d A, double* norm );
EL_EXPORT ElError ElTwoNorm_c( ElConstMatrix_c A, float* norm );
EL_EXPORT ElError ElTwoNorm_z( ElConstMatrix_z A, double* norm );
EL_EXPORT ElError ElTwoNormDist_s( ElConstDistMatrix_s A, float* norm );
EL_EXPORT ElError ElTwoNormDist_d( ElConstDistMatrix_d A, double* norm );
EL_EXPORT ElError ElTwoNormDist_c( ElConstDistMatrix_c A, float* norm );
EL_EXPORT ElError ElTwoNormDist_z( ElConstDistMatrix_z A, double* norm );

EL_EXPORT ElError ElHermitianTwoNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float* norm );
EL_EXPORT ElError ElHermitianTwoNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double* norm );
EL_EXPORT ElError ElHermitianTwoNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float* norm );
EL_EXPORT ElError ElHermitianTwoNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double* norm );

EL_EXPORT ElError ElSymmetricTwoNorm_s
( ElUpperOrLower uplo, ElConstMatrix_s A, float* norm );
EL_EXPORT ElError ElSymmetricTwoNorm_d
( ElUpperOrLower uplo, ElConstMatrix_d A, double* norm );
EL_EXPORT ElError ElSymmetricTwoNorm_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float* norm );
EL_EXPORT ElError ElSymmetricTwoNorm_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double* norm );
EL_EXPORT ElError ElSymmetricTwoNormDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, float* norm );
EL_EXPORT ElError ElSymmetricTwoNormDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, double* norm );
EL_EXPORT ElError ElSymmetricTwoNormDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float* norm );
EL_EXPORT ElError ElSymmetricTwoNormDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double* norm );

/* Zero norm
   --------- */
EL_EXPORT ElError ElZeroNorm_i
( ElConstMatrix_i A, ElInt tol, ElInt* numNonzero );
EL_EXPORT ElError ElZeroNorm_s
( ElConstMatrix_s A, float tol, ElInt* numNonzero );
EL_EXPORT ElError ElZeroNorm_d
( ElConstMatrix_d A, double tol, ElInt* numNonzero );
EL_EXPORT ElError ElZeroNorm_c
( ElConstMatrix_c A, float tol, ElInt* numNonzero );
EL_EXPORT ElError ElZeroNorm_z
( ElConstMatrix_z A, double tol, ElInt* numNonzero );

EL_EXPORT ElError ElZeroNormDist_i
( ElConstDistMatrix_i A, ElInt tol, ElInt* numNonzero );
EL_EXPORT ElError ElZeroNormDist_s
( ElConstDistMatrix_s A, float tol, ElInt* numNonzero );
EL_EXPORT ElError ElZeroNormDist_d
( ElConstDistMatrix_d A, double tol, ElInt* numNonzero );
EL_EXPORT ElError ElZeroNormDist_c
( ElConstDistMatrix_c A, float tol, ElInt* numNonzero );
EL_EXPORT ElError ElZeroNormDist_z
( ElConstDistMatrix_z A, double tol, ElInt* numNonzero );

EL_EXPORT ElError ElZeroNormSparse_i
( ElConstSparseMatrix_i A, ElInt tol, ElInt* numNonzero );
EL_EXPORT ElError ElZeroNormSparse_s
( ElConstSparseMatrix_s A, float tol, ElInt* numNonzero );
EL_EXPORT ElError ElZeroNormSparse_d
( ElConstSparseMatrix_d A, double tol, ElInt* numNonzero );
EL_EXPORT ElError ElZeroNormSparse_c
( ElConstSparseMatrix_c A, float tol, ElInt* numNonzero );
EL_EXPORT ElError ElZeroNormSparse_z
( ElConstSparseMatrix_z A, double tol, ElInt* numNonzero );

EL_EXPORT ElError ElZeroNormDistSparse_i
( ElConstDistSparseMatrix_i A, ElInt tol, ElInt* numNonzero );
EL_EXPORT ElError ElZeroNormDistSparse_s
( ElConstDistSparseMatrix_s A, float tol, ElInt* numNonzero );
EL_EXPORT ElError ElZeroNormDistSparse_d
( ElConstDistSparseMatrix_d A, double tol, ElInt* numNonzero );
EL_EXPORT ElError ElZeroNormDistSparse_c
( ElConstDistSparseMatrix_c A, float tol, ElInt* numNonzero );
EL_EXPORT ElError ElZeroNormDistSparse_z
( ElConstDistSparseMatrix_z A, double tol, ElInt* numNonzero );

/* Two-norm estimates
   ------------------ */
EL_EXPORT ElError ElTwoNormEstimate_s
( ElConstMatrix_s A, float tol, ElInt maxIts, float* normEst );
EL_EXPORT ElError ElTwoNormEstimate_d
( ElConstMatrix_d A, double tol, ElInt maxIts, double* normEst );
EL_EXPORT ElError ElTwoNormEstimate_c
( ElConstMatrix_c A, float tol, ElInt maxIts, float* normEst );
EL_EXPORT ElError ElTwoNormEstimate_z
( ElConstMatrix_z A, double tol, ElInt maxIts, double* normEst );

EL_EXPORT ElError ElTwoNormEstimateDist_s
( ElConstDistMatrix_s A, float tol, ElInt maxIts, float* normEst );
EL_EXPORT ElError ElTwoNormEstimateDist_d
( ElConstDistMatrix_d A, double tol, ElInt maxIts, double* normEst );
EL_EXPORT ElError ElTwoNormEstimateDist_c
( ElConstDistMatrix_c A, float tol, ElInt maxIts, float* normEst );
EL_EXPORT ElError ElTwoNormEstimateDist_z
( ElConstDistMatrix_z A, double tol, ElInt maxIts, double* normEst );

EL_EXPORT ElError ElTwoNormEstimateSparse_s
( ElConstSparseMatrix_s A, ElInt basisSize, float* normEst );
EL_EXPORT ElError ElTwoNormEstimateSparse_d
( ElConstSparseMatrix_d A, ElInt basisSize, double* normEst );
EL_EXPORT ElError ElTwoNormEstimateSparse_c
( ElConstSparseMatrix_c A, ElInt basisSize, float* normEst );
EL_EXPORT ElError ElTwoNormEstimateSparse_z
( ElConstSparseMatrix_z A, ElInt basisSize, double* normEst );

EL_EXPORT ElError ElTwoNormEstimateDistSparse_s
( ElConstDistSparseMatrix_s A, ElInt basisSize, float* normEst );
EL_EXPORT ElError ElTwoNormEstimateDistSparse_d
( ElConstDistSparseMatrix_d A, ElInt basisSize, double* normEst );
EL_EXPORT ElError ElTwoNormEstimateDistSparse_c
( ElConstDistSparseMatrix_c A, ElInt basisSize, float* normEst );
EL_EXPORT ElError ElTwoNormEstimateDistSparse_z
( ElConstDistSparseMatrix_z A, ElInt basisSize, double* normEst );

EL_EXPORT ElError ElSymmetricTwoNormEstimate_s
( ElUpperOrLower uplo, ElConstMatrix_s A, float tol, ElInt maxIts, 
  float* normEst );
EL_EXPORT ElError ElSymmetricTwoNormEstimate_d
( ElUpperOrLower uplo, ElConstMatrix_d A, double tol, ElInt maxIts, 
  double* normEst );
EL_EXPORT ElError ElSymmetricTwoNormEstimate_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float tol, ElInt maxIts, 
  float* normEst );
EL_EXPORT ElError ElSymmetricTwoNormEstimate_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double tol, ElInt maxIts, 
  double* normEst );
EL_EXPORT ElError ElSymmetricTwoNormEstimateDist_s
( ElUpperOrLower uplo, ElConstDistMatrix_s A, float tol, ElInt maxIts, 
  float* normEst );
EL_EXPORT ElError ElSymmetricTwoNormEstimateDist_d
( ElUpperOrLower uplo, ElConstDistMatrix_d A, double tol, ElInt maxIts, 
  double* normEst );
EL_EXPORT ElError ElSymmetricTwoNormEstimateDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float tol, ElInt maxIts, 
  float* normEst );
EL_EXPORT ElError ElSymmetricTwoNormEstimateDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double tol, ElInt maxIts, 
  double* normEst );

EL_EXPORT ElError ElHermitianTwoNormEstimate_c
( ElUpperOrLower uplo, ElConstMatrix_c A, float tol, ElInt maxIts, 
  float* normEst );
EL_EXPORT ElError ElHermitianTwoNormEstimate_z
( ElUpperOrLower uplo, ElConstMatrix_z A, double tol, ElInt maxIts, 
  double* normEst );
EL_EXPORT ElError ElHermitianTwoNormEstimateDist_c
( ElUpperOrLower uplo, ElConstDistMatrix_c A, float tol, ElInt maxIts, 
  float* normEst );
EL_EXPORT ElError ElHermitianTwoNormEstimateDist_z
( ElUpperOrLower uplo, ElConstDistMatrix_z A, double tol, ElInt maxIts, 
  double* normEst );

/* Trace
   ===== */
EL_EXPORT ElError ElTrace_i( ElConstMatrix_i A, ElInt* trace );
EL_EXPORT ElError ElTrace_s( ElConstMatrix_s A, float* trace );
EL_EXPORT ElError ElTrace_d( ElConstMatrix_d A, double* trace );
EL_EXPORT ElError ElTrace_c( ElConstMatrix_c A, complex_float* trace );
EL_EXPORT ElError ElTrace_z( ElConstMatrix_z A, complex_double* trace );

EL_EXPORT ElError ElTraceDist_i( ElConstDistMatrix_i A, ElInt* trace );
EL_EXPORT ElError ElTraceDist_s( ElConstDistMatrix_s A, float* trace );
EL_EXPORT ElError ElTraceDist_d( ElConstDistMatrix_d A, double* trace );
EL_EXPORT ElError ElTraceDist_c( ElConstDistMatrix_c A, complex_float* trace );
EL_EXPORT ElError ElTraceDist_z( ElConstDistMatrix_z A, complex_double* trace );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_LAPACK_PROPS_C_H */
