/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_PERM_C_H
#define EL_PERM_C_H

#include "El/core/DistMatrix.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  ElInt align;
  MPI_Comm comm;

  int *sendCounts, *sendDispls,
      *recvCounts, *recvDispls;

  int numSendIdx, numRecvIdx;
  int *sendIdx, *sendRanks,
      *recvIdx, *recvRanks;
} ElPermutationMeta;

EL_EXPORT ElError ElPermutationMetaSet
( const ElDistMatrix_i p, const ElDistMatrix_i pInv, ElPermutationMeta* meta );
EL_EXPORT ElError ElPermutationMetaClear( const ElPermutationMeta* meta );

EL_EXPORT ElError ElPermutationMetaTotalSend
( const ElPermutationMeta* meta, int* total );
EL_EXPORT ElError ElPermutationMetaTotalRecv
( const ElPermutationMeta* meta, int* total );

EL_EXPORT ElError ElPermutationMetaScaleUp
( ElPermutationMeta* meta, ElInt length );
EL_EXPORT ElError ElPermutationMetaScaleDown
( ElPermutationMeta* meta, ElInt length );

/* Apply column pivots
   =================== */
EL_EXPORT ElError ElApplyColPivots_i
( ElMatrix_i A, ElConstMatrix_i pivots, ElInt offset );
EL_EXPORT ElError ElApplyColPivots_s
( ElMatrix_s A, ElConstMatrix_i pivots, ElInt offset );
EL_EXPORT ElError ElApplyColPivots_d
( ElMatrix_d A, ElConstMatrix_i pivots, ElInt offset );
EL_EXPORT ElError ElApplyColPivots_c
( ElMatrix_c A, ElConstMatrix_i pivots, ElInt offset );
EL_EXPORT ElError ElApplyColPivots_z
( ElMatrix_z A, ElConstMatrix_i pivots, ElInt offset );

EL_EXPORT ElError ElApplyColPivotsDist_i
( ElDistMatrix_i A, ElConstDistMatrix_i pivots, ElInt offset );
EL_EXPORT ElError ElApplyColPivotsDist_s
( ElDistMatrix_s A, ElConstDistMatrix_i pivots, ElInt offset );
EL_EXPORT ElError ElApplyColPivotsDist_d
( ElDistMatrix_d A, ElConstDistMatrix_i pivots, ElInt offset );
EL_EXPORT ElError ElApplyColPivotsDist_c
( ElDistMatrix_c A, ElConstDistMatrix_i pivots, ElInt offset );
EL_EXPORT ElError ElApplyColPivotsDist_z
( ElDistMatrix_z A, ElConstDistMatrix_i pivots, ElInt offset );

/* Apply row pivots
   ================ */
EL_EXPORT ElError ElApplyRowPivots_i
( ElMatrix_i A, ElConstMatrix_i pivots, ElInt offset );
EL_EXPORT ElError ElApplyRowPivots_s
( ElMatrix_s A, ElConstMatrix_i pivots, ElInt offset );
EL_EXPORT ElError ElApplyRowPivots_d
( ElMatrix_d A, ElConstMatrix_i pivots, ElInt offset );
EL_EXPORT ElError ElApplyRowPivots_c
( ElMatrix_c A, ElConstMatrix_i pivots, ElInt offset );
EL_EXPORT ElError ElApplyRowPivots_z
( ElMatrix_z A, ElConstMatrix_i pivots, ElInt offset );

EL_EXPORT ElError ElApplyRowPivotsDist_i
( ElDistMatrix_i A, ElConstDistMatrix_i pivots, ElInt offset );
EL_EXPORT ElError ElApplyRowPivotsDist_s
( ElDistMatrix_s A, ElConstDistMatrix_i pivots, ElInt offset );
EL_EXPORT ElError ElApplyRowPivotsDist_d
( ElDistMatrix_d A, ElConstDistMatrix_i pivots, ElInt offset );
EL_EXPORT ElError ElApplyRowPivotsDist_c
( ElDistMatrix_c A, ElConstDistMatrix_i pivots, ElInt offset );
EL_EXPORT ElError ElApplyRowPivotsDist_z
( ElDistMatrix_z A, ElConstDistMatrix_i pivots, ElInt offset );

/* Apply symmetric pivots
   ====================== */
EL_EXPORT ElError ElApplySymmetricPivots_i
( ElUpperOrLower uplo, ElMatrix_i A, ElConstMatrix_i p, ElInt offset );
EL_EXPORT ElError ElApplySymmetricPivots_s
( ElUpperOrLower uplo, ElMatrix_s A, ElConstMatrix_i p, ElInt offset );
EL_EXPORT ElError ElApplySymmetricPivots_d
( ElUpperOrLower uplo, ElMatrix_d A, ElConstMatrix_i p, ElInt offset );
EL_EXPORT ElError ElApplySymmetricPivots_c
( ElUpperOrLower uplo, ElMatrix_c A, ElConstMatrix_i p, 
  bool conjugate, ElInt offset );
EL_EXPORT ElError ElApplySymmetricPivots_z
( ElUpperOrLower uplo, ElMatrix_z A, ElConstMatrix_i p, 
  bool conjugate, ElInt offset );

EL_EXPORT ElError ElApplySymmetricPivotsDist_i
( ElUpperOrLower uplo, ElDistMatrix_i A, ElConstDistMatrix_i p, ElInt offset );
EL_EXPORT ElError ElApplySymmetricPivotsDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A, ElConstDistMatrix_i p, ElInt offset );
EL_EXPORT ElError ElApplySymmetricPivotsDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A, ElConstDistMatrix_i p, ElInt offset );
EL_EXPORT ElError ElApplySymmetricPivotsDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A, ElConstDistMatrix_i p, 
  bool conjugate, ElInt offset );
EL_EXPORT ElError ElApplySymmetricPivotsDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A, ElConstDistMatrix_i p, 
  bool conjugate, ElInt offset );

EL_EXPORT ElError ElApplyInverseSymmetricPivots_i
( ElUpperOrLower uplo, ElMatrix_i A, ElConstMatrix_i p, ElInt offset );
EL_EXPORT ElError ElApplyInverseSymmetricPivots_s
( ElUpperOrLower uplo, ElMatrix_s A, ElConstMatrix_i p, ElInt offset );
EL_EXPORT ElError ElApplyInverseSymmetricPivots_d
( ElUpperOrLower uplo, ElMatrix_d A, ElConstMatrix_i p, ElInt offset );
EL_EXPORT ElError ElApplyInverseSymmetricPivots_c
( ElUpperOrLower uplo, ElMatrix_c A, ElConstMatrix_i p, 
  bool conjugate, ElInt offset );
EL_EXPORT ElError ElApplyInverseSymmetricPivots_z
( ElUpperOrLower uplo, ElMatrix_z A, ElConstMatrix_i p, 
  bool conjugate, ElInt offset );

EL_EXPORT ElError ElApplyInverseSymmetricPivotsDist_i
( ElUpperOrLower uplo, ElDistMatrix_i A, ElConstDistMatrix_i p, ElInt offset );
EL_EXPORT ElError ElApplyInverseSymmetricPivotsDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A, ElConstDistMatrix_i p, ElInt offset );
EL_EXPORT ElError ElApplyInverseSymmetricPivotsDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A, ElConstDistMatrix_i p, ElInt offset );
EL_EXPORT ElError ElApplyInverseSymmetricPivotsDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A, ElConstDistMatrix_i p, 
  bool conjugate, ElInt offset );
EL_EXPORT ElError ElApplyInverseSymmetricPivotsDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A, ElConstDistMatrix_i p, 
  bool conjugate, ElInt offset );

/* Explicit permutation
   ==================== */
EL_EXPORT ElError ElExplicitPermutation( ElConstMatrix_i p, ElMatrix_i P );
EL_EXPORT ElError ElExplicitPermutationDist
( ElConstDistMatrix_i p, ElDistMatrix_i P );

/* Invert permutation
   ================== */
EL_EXPORT ElError ElInvertPermutation( ElConstMatrix_i p, ElMatrix_i pInv );
EL_EXPORT ElError ElInvertPermutationDist
( ElConstDistMatrix_i p, ElDistMatrix_i pInv );

/* Parity of a permutation
   ======================= */
EL_EXPORT ElError ElPermutationParity( ElConstMatrix_i p, bool* parity );
EL_EXPORT ElError ElPermutationParityDist
( ElConstDistMatrix_i p, bool* parity );

/* Permute columns
   =============== */

/* With only the forward permutation vector provided
   ------------------------------------------------- */
EL_EXPORT ElError ElPermuteCols_i( ElMatrix_i A, ElConstMatrix_i p );
EL_EXPORT ElError ElPermuteCols_s( ElMatrix_s A, ElConstMatrix_i p );
EL_EXPORT ElError ElPermuteCols_d( ElMatrix_d A, ElConstMatrix_i p );
EL_EXPORT ElError ElPermuteCols_c( ElMatrix_c A, ElConstMatrix_i p );
EL_EXPORT ElError ElPermuteCols_z( ElMatrix_z A, ElConstMatrix_i p );

EL_EXPORT ElError ElPermuteColsDist_i
( ElDistMatrix_i A, ElConstDistMatrix_i p );
EL_EXPORT ElError ElPermuteColsDist_s
( ElDistMatrix_s A, ElConstDistMatrix_i p );
EL_EXPORT ElError ElPermuteColsDist_d
( ElDistMatrix_d A, ElConstDistMatrix_i p );
EL_EXPORT ElError ElPermuteColsDist_c
( ElDistMatrix_c A, ElConstDistMatrix_i p );
EL_EXPORT ElError ElPermuteColsDist_z
( ElDistMatrix_z A, ElConstDistMatrix_i p );

EL_EXPORT ElError ElInversePermuteCols_i( ElMatrix_i A, ElConstMatrix_i p );
EL_EXPORT ElError ElInversePermuteCols_s( ElMatrix_s A, ElConstMatrix_i p );
EL_EXPORT ElError ElInversePermuteCols_d( ElMatrix_d A, ElConstMatrix_i p );
EL_EXPORT ElError ElInversePermuteCols_c( ElMatrix_c A, ElConstMatrix_i p );
EL_EXPORT ElError ElInversePermuteCols_z( ElMatrix_z A, ElConstMatrix_i p );

EL_EXPORT ElError ElInversePermuteColsDist_i
( ElDistMatrix_i A, ElConstDistMatrix_i p );
EL_EXPORT ElError ElInversePermuteColsDist_s
( ElDistMatrix_s A, ElConstDistMatrix_i p );
EL_EXPORT ElError ElInversePermuteColsDist_d
( ElDistMatrix_d A, ElConstDistMatrix_i p );
EL_EXPORT ElError ElInversePermuteColsDist_c
( ElDistMatrix_c A, ElConstDistMatrix_i p );
EL_EXPORT ElError ElInversePermuteColsDist_z
( ElDistMatrix_z A, ElConstDistMatrix_i p );

/* With the forward and inverse permutation vectors provided
   --------------------------------------------------------- */
EL_EXPORT ElError ElPermuteColsBoth_i
( ElMatrix_i A, ElConstMatrix_i p, ElConstMatrix_i pInv );
EL_EXPORT ElError ElPermuteColsBoth_s
( ElMatrix_s A, ElConstMatrix_i p, ElConstMatrix_i pInv );
EL_EXPORT ElError ElPermuteColsBoth_d
( ElMatrix_d A, ElConstMatrix_i p, ElConstMatrix_i pInv );
EL_EXPORT ElError ElPermuteColsBoth_c
( ElMatrix_c A, ElConstMatrix_i p, ElConstMatrix_i pInv );
EL_EXPORT ElError ElPermuteColsBoth_z
( ElMatrix_z A, ElConstMatrix_i p, ElConstMatrix_i pInv );

EL_EXPORT ElError ElPermuteColsBothDist_i
( ElDistMatrix_i A, ElConstDistMatrix_i p, ElConstDistMatrix_i pInv );
EL_EXPORT ElError ElPermuteColsBothDist_s
( ElDistMatrix_s A, ElConstDistMatrix_i p, ElConstDistMatrix_i pInv );
EL_EXPORT ElError ElPermuteColsBothDist_d
( ElDistMatrix_d A, ElConstDistMatrix_i p, ElConstDistMatrix_i pInv );
EL_EXPORT ElError ElPermuteColsBothDist_c
( ElDistMatrix_c A, ElConstDistMatrix_i p, ElConstDistMatrix_i pInv );
EL_EXPORT ElError ElPermuteColsBothDist_z
( ElDistMatrix_z A, ElConstDistMatrix_i p, ElConstDistMatrix_i pInv );

/* With the pivot plan provided
   ---------------------------- */
EL_EXPORT ElError ElPermuteColsMetaDist_i
( ElDistMatrix_i A, const ElPermutationMeta* meta );
EL_EXPORT ElError ElPermuteColsMetaDist_s
( ElDistMatrix_s A, const ElPermutationMeta* meta );
EL_EXPORT ElError ElPermuteColsMetaDist_d
( ElDistMatrix_d A, const ElPermutationMeta* meta );
EL_EXPORT ElError ElPermuteColsMetaDist_c
( ElDistMatrix_c A, const ElPermutationMeta* meta );
EL_EXPORT ElError ElPermuteColsMetaDist_z
( ElDistMatrix_z A, const ElPermutationMeta* meta );

/* Permute rows
   ============ */

/* With only the forward permutation vector provided
   ------------------------------------------------- */
EL_EXPORT ElError ElPermuteRows_i( ElMatrix_i A, ElConstMatrix_i p );
EL_EXPORT ElError ElPermuteRows_s( ElMatrix_s A, ElConstMatrix_i p );
EL_EXPORT ElError ElPermuteRows_d( ElMatrix_d A, ElConstMatrix_i p );
EL_EXPORT ElError ElPermuteRows_c( ElMatrix_c A, ElConstMatrix_i p );
EL_EXPORT ElError ElPermuteRows_z( ElMatrix_z A, ElConstMatrix_i p );

EL_EXPORT ElError ElPermuteRowsDist_i
( ElDistMatrix_i A, ElConstDistMatrix_i p );
EL_EXPORT ElError ElPermuteRowsDist_s
( ElDistMatrix_s A, ElConstDistMatrix_i p );
EL_EXPORT ElError ElPermuteRowsDist_d
( ElDistMatrix_d A, ElConstDistMatrix_i p );
EL_EXPORT ElError ElPermuteRowsDist_c
( ElDistMatrix_c A, ElConstDistMatrix_i p );
EL_EXPORT ElError ElPermuteRowsDist_z
( ElDistMatrix_z A, ElConstDistMatrix_i p );

EL_EXPORT ElError ElInversePermuteRows_i( ElMatrix_i A, ElConstMatrix_i p );
EL_EXPORT ElError ElInversePermuteRows_s( ElMatrix_s A, ElConstMatrix_i p );
EL_EXPORT ElError ElInversePermuteRows_d( ElMatrix_d A, ElConstMatrix_i p );
EL_EXPORT ElError ElInversePermuteRows_c( ElMatrix_c A, ElConstMatrix_i p );
EL_EXPORT ElError ElInversePermuteRows_z( ElMatrix_z A, ElConstMatrix_i p );

EL_EXPORT ElError ElInversePermuteRowsDist_i
( ElDistMatrix_i A, ElConstDistMatrix_i p );
EL_EXPORT ElError ElInversePermuteRowsDist_s
( ElDistMatrix_s A, ElConstDistMatrix_i p );
EL_EXPORT ElError ElInversePermuteRowsDist_d
( ElDistMatrix_d A, ElConstDistMatrix_i p );
EL_EXPORT ElError ElInversePermuteRowsDist_c
( ElDistMatrix_c A, ElConstDistMatrix_i p );
EL_EXPORT ElError ElInversePermuteRowsDist_z
( ElDistMatrix_z A, ElConstDistMatrix_i p );

/* With the forward and inverse permutation vectors provided
   --------------------------------------------------------- */
EL_EXPORT ElError ElPermuteRowsBoth_i
( ElMatrix_i A, ElConstMatrix_i p, ElConstMatrix_i pInv );
EL_EXPORT ElError ElPermuteRowsBoth_s
( ElMatrix_s A, ElConstMatrix_i p, ElConstMatrix_i pInv );
EL_EXPORT ElError ElPermuteRowsBoth_d
( ElMatrix_d A, ElConstMatrix_i p, ElConstMatrix_i pInv );
EL_EXPORT ElError ElPermuteRowsBoth_c
( ElMatrix_c A, ElConstMatrix_i p, ElConstMatrix_i pInv );
EL_EXPORT ElError ElPermuteRowsBoth_z
( ElMatrix_z A, ElConstMatrix_i p, ElConstMatrix_i pInv );

EL_EXPORT ElError ElPermuteRowsBothDist_i
( ElDistMatrix_i A, ElConstDistMatrix_i p, ElConstDistMatrix_i pInv );
EL_EXPORT ElError ElPermuteRowsBothDist_s
( ElDistMatrix_s A, ElConstDistMatrix_i p, ElConstDistMatrix_i pInv );
EL_EXPORT ElError ElPermuteRowsBothDist_d
( ElDistMatrix_d A, ElConstDistMatrix_i p, ElConstDistMatrix_i pInv );
EL_EXPORT ElError ElPermuteRowsBothDist_c
( ElDistMatrix_c A, ElConstDistMatrix_i p, ElConstDistMatrix_i pInv );
EL_EXPORT ElError ElPermuteRowsBothDist_z
( ElDistMatrix_z A, ElConstDistMatrix_i p, ElConstDistMatrix_i pInv );

/* With the pivot plan provided
   ---------------------------- */
EL_EXPORT ElError ElPermuteRowsMetaDist_i
( ElDistMatrix_i A, const ElPermutationMeta* meta );
EL_EXPORT ElError ElPermuteRowsMetaDist_s
( ElDistMatrix_s A, const ElPermutationMeta* meta );
EL_EXPORT ElError ElPermuteRowsMetaDist_d
( ElDistMatrix_d A, const ElPermutationMeta* meta );
EL_EXPORT ElError ElPermuteRowsMetaDist_c
( ElDistMatrix_c A, const ElPermutationMeta* meta );
EL_EXPORT ElError ElPermuteRowsMetaDist_z
( ElDistMatrix_z A, const ElPermutationMeta* meta );

/* Pivot parity
   ============ */
EL_EXPORT ElError ElPivotParity
( ElConstMatrix_i p, ElInt pivotOffset, bool* parity );
EL_EXPORT ElError ElPivotParityDist
( ElConstDistMatrix_i p, ElInt pivotOffset, bool* parity );

/* Convert a pivot sequence to partial permutation vector
   ====================================================== */
EL_EXPORT ElError ElPivotsToPartialPermutation
( ElConstMatrix_i pivots, ElMatrix_i p, ElMatrix_i pInv, 
  ElInt offset );
EL_EXPORT ElError ElPivotsToPartialPermutationDist
( ElConstDistMatrix_i pivots, ElDistMatrix_i p, ElDistMatrix_i pInv, 
  ElInt offset );

/* Convert a pivot sequence to a permutation vector
   ================================================ */
EL_EXPORT ElError ElPivotsToPermutation
( ElConstMatrix_i pivots, ElMatrix_i p, ElInt offset );
EL_EXPORT ElError ElPivotsToPermutationDist
( ElConstDistMatrix_i pivots, ElDistMatrix_i p, ElInt offset );

EL_EXPORT ElError ElPivotsToInversePermutation
( ElConstMatrix_i pivots, ElMatrix_i pInv, ElInt offset );
EL_EXPORT ElError ElPivotsToInversePermutationDist
( ElConstDistMatrix_i pivots, ElDistMatrix_i pInv, ElInt offset );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_PERM_C_H */
