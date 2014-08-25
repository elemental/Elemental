/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_PERM_C_H
#define EL_PERM_C_H

#include "El/core/DistMatrix-C.h"

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

ElError ElPermutationMetaSet
( const ElDistMatrix_i p, const ElDistMatrix_i pInv, ElPermutationMeta* meta );
ElError ElPermutationMetaClear( const ElPermutationMeta* meta );

ElError ElPermutationMetaTotalSend( const ElPermutationMeta* meta, int* total );
ElError ElPermutationMetaTotalRecv( const ElPermutationMeta* meta, int* total );

ElError ElPermutationMetaScaleUp( ElPermutationMeta* meta, ElInt length );
ElError ElPermutationMetaScaleDown( ElPermutationMeta* meta, ElInt length );

/* Apply column pivots
   =================== */
ElError ElApplyColPivots_i
( ElMatrix_i A, ElConstMatrix_i pivots, ElInt offset );
ElError ElApplyColPivots_s
( ElMatrix_s A, ElConstMatrix_i pivots, ElInt offset );
ElError ElApplyColPivots_d
( ElMatrix_d A, ElConstMatrix_i pivots, ElInt offset );
ElError ElApplyColPivots_c
( ElMatrix_c A, ElConstMatrix_i pivots, ElInt offset );
ElError ElApplyColPivots_z
( ElMatrix_z A, ElConstMatrix_i pivots, ElInt offset );

ElError ElApplyColPivotsDist_i
( ElDistMatrix_i A, ElConstDistMatrix_i pivots, ElInt offset );
ElError ElApplyColPivotsDist_s
( ElDistMatrix_s A, ElConstDistMatrix_i pivots, ElInt offset );
ElError ElApplyColPivotsDist_d
( ElDistMatrix_d A, ElConstDistMatrix_i pivots, ElInt offset );
ElError ElApplyColPivotsDist_c
( ElDistMatrix_c A, ElConstDistMatrix_i pivots, ElInt offset );
ElError ElApplyColPivotsDist_z
( ElDistMatrix_z A, ElConstDistMatrix_i pivots, ElInt offset );

/* Apply row pivots
   ================ */
ElError ElApplyRowPivots_i
( ElMatrix_i A, ElConstMatrix_i pivots, ElInt offset );
ElError ElApplyRowPivots_s
( ElMatrix_s A, ElConstMatrix_i pivots, ElInt offset );
ElError ElApplyRowPivots_d
( ElMatrix_d A, ElConstMatrix_i pivots, ElInt offset );
ElError ElApplyRowPivots_c
( ElMatrix_c A, ElConstMatrix_i pivots, ElInt offset );
ElError ElApplyRowPivots_z
( ElMatrix_z A, ElConstMatrix_i pivots, ElInt offset );

ElError ElApplyRowPivotsDist_i
( ElDistMatrix_i A, ElConstDistMatrix_i pivots, ElInt offset );
ElError ElApplyRowPivotsDist_s
( ElDistMatrix_s A, ElConstDistMatrix_i pivots, ElInt offset );
ElError ElApplyRowPivotsDist_d
( ElDistMatrix_d A, ElConstDistMatrix_i pivots, ElInt offset );
ElError ElApplyRowPivotsDist_c
( ElDistMatrix_c A, ElConstDistMatrix_i pivots, ElInt offset );
ElError ElApplyRowPivotsDist_z
( ElDistMatrix_z A, ElConstDistMatrix_i pivots, ElInt offset );

/* Apply symmetric pivots
   ====================== */
ElError ElApplySymmetricPivots_i
( ElUpperOrLower uplo, ElMatrix_i A, ElConstMatrix_i p, ElInt offset );
ElError ElApplySymmetricPivots_s
( ElUpperOrLower uplo, ElMatrix_s A, ElConstMatrix_i p, ElInt offset );
ElError ElApplySymmetricPivots_d
( ElUpperOrLower uplo, ElMatrix_d A, ElConstMatrix_i p, ElInt offset );
ElError ElApplySymmetricPivots_c
( ElUpperOrLower uplo, ElMatrix_c A, ElConstMatrix_i p, 
  bool conjugate, ElInt offset );
ElError ElApplySymmetricPivots_z
( ElUpperOrLower uplo, ElMatrix_z A, ElConstMatrix_i p, 
  bool conjugate, ElInt offset );

ElError ElApplySymmetricPivotsDist_i
( ElUpperOrLower uplo, ElDistMatrix_i A, ElConstDistMatrix_i p, ElInt offset );
ElError ElApplySymmetricPivotsDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A, ElConstDistMatrix_i p, ElInt offset );
ElError ElApplySymmetricPivotsDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A, ElConstDistMatrix_i p, ElInt offset );
ElError ElApplySymmetricPivotsDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A, ElConstDistMatrix_i p, 
  bool conjugate, ElInt offset );
ElError ElApplySymmetricPivotsDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A, ElConstDistMatrix_i p, 
  bool conjugate, ElInt offset );

ElError ElApplyInverseSymmetricPivots_i
( ElUpperOrLower uplo, ElMatrix_i A, ElConstMatrix_i p, ElInt offset );
ElError ElApplyInverseSymmetricPivots_s
( ElUpperOrLower uplo, ElMatrix_s A, ElConstMatrix_i p, ElInt offset );
ElError ElApplyInverseSymmetricPivots_d
( ElUpperOrLower uplo, ElMatrix_d A, ElConstMatrix_i p, ElInt offset );
ElError ElApplyInverseSymmetricPivots_c
( ElUpperOrLower uplo, ElMatrix_c A, ElConstMatrix_i p, 
  bool conjugate, ElInt offset );
ElError ElApplyInverseSymmetricPivots_z
( ElUpperOrLower uplo, ElMatrix_z A, ElConstMatrix_i p, 
  bool conjugate, ElInt offset );

ElError ElApplyInverseSymmetricPivotsDist_i
( ElUpperOrLower uplo, ElDistMatrix_i A, ElConstDistMatrix_i p, ElInt offset );
ElError ElApplyInverseSymmetricPivotsDist_s
( ElUpperOrLower uplo, ElDistMatrix_s A, ElConstDistMatrix_i p, ElInt offset );
ElError ElApplyInverseSymmetricPivotsDist_d
( ElUpperOrLower uplo, ElDistMatrix_d A, ElConstDistMatrix_i p, ElInt offset );
ElError ElApplyInverseSymmetricPivotsDist_c
( ElUpperOrLower uplo, ElDistMatrix_c A, ElConstDistMatrix_i p, 
  bool conjugate, ElInt offset );
ElError ElApplyInverseSymmetricPivotsDist_z
( ElUpperOrLower uplo, ElDistMatrix_z A, ElConstDistMatrix_i p, 
  bool conjugate, ElInt offset );

/* Explicit permutation
   ==================== */
ElError ElExplicitPermutation( ElConstMatrix_i p, ElMatrix_i P );
ElError ElExplicitPermutationDist( ElConstDistMatrix_i p, ElDistMatrix_i P );

/* Invert permutation
   ================== */
ElError ElInvertPermutation( ElConstMatrix_i p, ElMatrix_i pInv );
ElError ElInvertPermutationDist( ElConstDistMatrix_i p, ElDistMatrix_i pInv );

/* Permute columns
   =============== */

/* With only the forward permutation vector provided
   ------------------------------------------------- */
ElError ElPermuteCols_i( ElMatrix_i A, ElConstMatrix_i p );
ElError ElPermuteCols_s( ElMatrix_s A, ElConstMatrix_i p );
ElError ElPermuteCols_d( ElMatrix_d A, ElConstMatrix_i p );
ElError ElPermuteCols_c( ElMatrix_c A, ElConstMatrix_i p );
ElError ElPermuteCols_z( ElMatrix_z A, ElConstMatrix_i p );

ElError ElPermuteColsDist_i( ElDistMatrix_i A, ElConstDistMatrix_i p );
ElError ElPermuteColsDist_s( ElDistMatrix_s A, ElConstDistMatrix_i p );
ElError ElPermuteColsDist_d( ElDistMatrix_d A, ElConstDistMatrix_i p );
ElError ElPermuteColsDist_c( ElDistMatrix_c A, ElConstDistMatrix_i p );
ElError ElPermuteColsDist_z( ElDistMatrix_z A, ElConstDistMatrix_i p );

ElError ElInversePermuteCols_i( ElMatrix_i A, ElConstMatrix_i p );
ElError ElInversePermuteCols_s( ElMatrix_s A, ElConstMatrix_i p );
ElError ElInversePermuteCols_d( ElMatrix_d A, ElConstMatrix_i p );
ElError ElInversePermuteCols_c( ElMatrix_c A, ElConstMatrix_i p );
ElError ElInversePermuteCols_z( ElMatrix_z A, ElConstMatrix_i p );

ElError ElInversePermuteColsDist_i( ElDistMatrix_i A, ElConstDistMatrix_i p );
ElError ElInversePermuteColsDist_s( ElDistMatrix_s A, ElConstDistMatrix_i p );
ElError ElInversePermuteColsDist_d( ElDistMatrix_d A, ElConstDistMatrix_i p );
ElError ElInversePermuteColsDist_c( ElDistMatrix_c A, ElConstDistMatrix_i p );
ElError ElInversePermuteColsDist_z( ElDistMatrix_z A, ElConstDistMatrix_i p );

/* With the forward and inverse permutation vectors provided
   --------------------------------------------------------- */
ElError ElPermuteColsBoth_i
( ElMatrix_i A, ElConstMatrix_i p, ElConstMatrix_i pInv );
ElError ElPermuteColsBoth_s
( ElMatrix_s A, ElConstMatrix_i p, ElConstMatrix_i pInv );
ElError ElPermuteColsBoth_d
( ElMatrix_d A, ElConstMatrix_i p, ElConstMatrix_i pInv );
ElError ElPermuteColsBoth_c
( ElMatrix_c A, ElConstMatrix_i p, ElConstMatrix_i pInv );
ElError ElPermuteColsBoth_z
( ElMatrix_z A, ElConstMatrix_i p, ElConstMatrix_i pInv );

ElError ElPermuteColsBothDist_i
( ElDistMatrix_i A, ElConstDistMatrix_i p, ElConstDistMatrix_i pInv );
ElError ElPermuteColsBothDist_s
( ElDistMatrix_s A, ElConstDistMatrix_i p, ElConstDistMatrix_i pInv );
ElError ElPermuteColsBothDist_d
( ElDistMatrix_d A, ElConstDistMatrix_i p, ElConstDistMatrix_i pInv );
ElError ElPermuteColsBothDist_c
( ElDistMatrix_c A, ElConstDistMatrix_i p, ElConstDistMatrix_i pInv );
ElError ElPermuteColsBothDist_z
( ElDistMatrix_z A, ElConstDistMatrix_i p, ElConstDistMatrix_i pInv );

/* With the pivot plan provided
   ---------------------------- */
ElError ElPermuteColsMetaDist_i
( ElDistMatrix_i A, const ElPermutationMeta* meta );
ElError ElPermuteColsMetaDist_s
( ElDistMatrix_s A, const ElPermutationMeta* meta );
ElError ElPermuteColsMetaDist_d
( ElDistMatrix_d A, const ElPermutationMeta* meta );
ElError ElPermuteColsMetaDist_c
( ElDistMatrix_c A, const ElPermutationMeta* meta );
ElError ElPermuteColsMetaDist_z
( ElDistMatrix_z A, const ElPermutationMeta* meta );

/* Permute rows
   ============ */

/* With only the forward permutation vector provided
   ------------------------------------------------- */
ElError ElPermuteRows_i( ElMatrix_i A, ElConstMatrix_i p );
ElError ElPermuteRows_s( ElMatrix_s A, ElConstMatrix_i p );
ElError ElPermuteRows_d( ElMatrix_d A, ElConstMatrix_i p );
ElError ElPermuteRows_c( ElMatrix_c A, ElConstMatrix_i p );
ElError ElPermuteRows_z( ElMatrix_z A, ElConstMatrix_i p );

ElError ElPermuteRowsDist_i( ElDistMatrix_i A, ElConstDistMatrix_i p );
ElError ElPermuteRowsDist_s( ElDistMatrix_s A, ElConstDistMatrix_i p );
ElError ElPermuteRowsDist_d( ElDistMatrix_d A, ElConstDistMatrix_i p );
ElError ElPermuteRowsDist_c( ElDistMatrix_c A, ElConstDistMatrix_i p );
ElError ElPermuteRowsDist_z( ElDistMatrix_z A, ElConstDistMatrix_i p );

ElError ElInversePermuteRows_i( ElMatrix_i A, ElConstMatrix_i p );
ElError ElInversePermuteRows_s( ElMatrix_s A, ElConstMatrix_i p );
ElError ElInversePermuteRows_d( ElMatrix_d A, ElConstMatrix_i p );
ElError ElInversePermuteRows_c( ElMatrix_c A, ElConstMatrix_i p );
ElError ElInversePermuteRows_z( ElMatrix_z A, ElConstMatrix_i p );

ElError ElInversePermuteRowsDist_i( ElDistMatrix_i A, ElConstDistMatrix_i p );
ElError ElInversePermuteRowsDist_s( ElDistMatrix_s A, ElConstDistMatrix_i p );
ElError ElInversePermuteRowsDist_d( ElDistMatrix_d A, ElConstDistMatrix_i p );
ElError ElInversePermuteRowsDist_c( ElDistMatrix_c A, ElConstDistMatrix_i p );
ElError ElInversePermuteRowsDist_z( ElDistMatrix_z A, ElConstDistMatrix_i p );

/* With the forward and inverse permutation vectors provided
   --------------------------------------------------------- */
ElError ElPermuteRowsBoth_i
( ElMatrix_i A, ElConstMatrix_i p, ElConstMatrix_i pInv );
ElError ElPermuteRowsBoth_s
( ElMatrix_s A, ElConstMatrix_i p, ElConstMatrix_i pInv );
ElError ElPermuteRowsBoth_d
( ElMatrix_d A, ElConstMatrix_i p, ElConstMatrix_i pInv );
ElError ElPermuteRowsBoth_c
( ElMatrix_c A, ElConstMatrix_i p, ElConstMatrix_i pInv );
ElError ElPermuteRowsBoth_z
( ElMatrix_z A, ElConstMatrix_i p, ElConstMatrix_i pInv );

ElError ElPermuteRowsBothDist_i
( ElDistMatrix_i A, ElConstDistMatrix_i p, ElConstDistMatrix_i pInv );
ElError ElPermuteRowsBothDist_s
( ElDistMatrix_s A, ElConstDistMatrix_i p, ElConstDistMatrix_i pInv );
ElError ElPermuteRowsBothDist_d
( ElDistMatrix_d A, ElConstDistMatrix_i p, ElConstDistMatrix_i pInv );
ElError ElPermuteRowsBothDist_c
( ElDistMatrix_c A, ElConstDistMatrix_i p, ElConstDistMatrix_i pInv );
ElError ElPermuteRowsBothDist_z
( ElDistMatrix_z A, ElConstDistMatrix_i p, ElConstDistMatrix_i pInv );

/* With the pivot plan provided
   ---------------------------- */
ElError ElPermuteRowsMetaDist_i
( ElDistMatrix_i A, const ElPermutationMeta* meta );
ElError ElPermuteRowsMetaDist_s
( ElDistMatrix_s A, const ElPermutationMeta* meta );
ElError ElPermuteRowsMetaDist_d
( ElDistMatrix_d A, const ElPermutationMeta* meta );
ElError ElPermuteRowsMetaDist_c
( ElDistMatrix_c A, const ElPermutationMeta* meta );
ElError ElPermuteRowsMetaDist_z
( ElDistMatrix_z A, const ElPermutationMeta* meta );

/* Pivot parity
   ============ */
ElError ElPivotParity
( ElConstMatrix_i p, ElInt pivotOffset, bool* parity );
ElError ElPivotParityDist
( ElConstDistMatrix_i p, ElInt pivotOffset, bool* parity );

/* Convert a pivot sequence to partial permutation vectors
   ======================================================= */
ElError ElPivotsToPartialPermutation
( ElConstMatrix_i pivots, ElMatrix_i p, ElMatrix_i pInv, 
  ElInt offset );
ElError ElPivotsToPartialPermutationDist
( ElConstDistMatrix_i pivots, ElDistMatrix_i p, ElDistMatrix_i pInv, 
  ElInt offset );

/* Convert a pivot sequence to a permutation vector
   ================================================ */
ElError ElPivotsToPermutation
( ElConstMatrix_i pivots, ElMatrix_i p, ElInt offset );
ElError ElPivotsToPermutationDist
( ElConstDistMatrix_i pivots, ElDistMatrix_i p, ElInt offset );

ElError ElPivotsToInversePermutation
( ElConstMatrix_i pivots, ElMatrix_i pInv, ElInt offset );
ElError ElPivotsToInversePermutationDist
( ElConstDistMatrix_i pivots, ElDistMatrix_i pInv, ElInt offset );

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_PERM_C_H */
