/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "El.h"
using namespace El;

extern "C" {

ElError ElPermutationMetaSet
( const ElDistMatrix_i p, const ElDistMatrix_i pInv, ElPermutationMeta* meta )
{
    EL_TRY( 
      PermutationMeta metaCpp( *CReflect(p), *CReflect(pInv) );
      *meta = CReflect(metaCpp);
    )
}

ElError ElPermutationMetaClear( const ElPermutationMeta* meta )
{
    EL_TRY( 
      delete[] meta->sendCounts;
      delete[] meta->sendDispls;
      delete[] meta->recvCounts;
      delete[] meta->recvDispls;
      delete[] meta->sendIdx;
      delete[] meta->sendRanks;
      delete[] meta->recvIdx;
      delete[] meta->recvRanks;
    )
}

ElError ElPermutationMetaTotalSend( const ElPermutationMeta* meta, int* total )
{
    EL_TRY(
      int commSize;
      MPI_Comm_size( meta->comm, &commSize );
      *total = meta->sendDispls[commSize-1] + meta->sendCounts[commSize-1];
    )
}

ElError ElPermutationMetaTotalRecv( const ElPermutationMeta* meta, int* total )
{
    EL_TRY(
      int commSize;
      MPI_Comm_size( meta->comm, &commSize );
      *total = meta->recvDispls[commSize-1] + meta->recvCounts[commSize-1];
    )
}

ElError ElPermutationMetaScaleUp( ElPermutationMeta* meta, ElInt length )
{
    EL_TRY(
      int commSize;
      MPI_Comm_size( meta->comm, &commSize );
      for( int q=0; q<commSize; ++q )
      {
          meta->sendCounts[q] *= length;
          meta->sendDispls[q] *= length;
          meta->recvCounts[q] *= length;
          meta->recvDispls[q] *= length;
      }
    )
}

ElError ElPermutationMetaScaleDown( ElPermutationMeta* meta, ElInt length )
{
    EL_TRY(
      int commSize;
      MPI_Comm_size( meta->comm, &commSize );
      for( int q=0; q<commSize; ++q )
      {
          meta->sendCounts[q] /= length;
          meta->sendDispls[q] /= length;
          meta->recvCounts[q] /= length;
          meta->recvDispls[q] /= length;
      }
    )
}

// Explicit permutation
// ====================
ElError ElExplicitPermutation( ElConstMatrix_i p, ElMatrix_i P )
{ EL_TRY( ExplicitPermutation( *CReflect(p), *CReflect(P) ) ) }

ElError ElExplicitPermutationDist( ElConstDistMatrix_i p, ElDistMatrix_i P )
{ EL_TRY( ExplicitPermutation( *CReflect(p), *CReflect(P) ) ) }

// Invert permutation
// ==================
ElError ElInvertPermutation( ElConstMatrix_i p, ElMatrix_i pInv )
{ EL_TRY( InvertPermutation( *CReflect(p), *CReflect(pInv) ) ) }

ElError ElInvertPermutationDist( ElConstDistMatrix_i p, ElDistMatrix_i pInv )
{ EL_TRY( InvertPermutation( *CReflect(p), *CReflect(pInv) ) ) }

// Parity of a permutation
// =======================
ElError ElPermutationParity( ElConstMatrix_i p, bool* parity )
{ EL_TRY( *parity = PermutationParity( *CReflect(p) ) ) }

ElError ElPermutationParityDist( ElConstDistMatrix_i p, bool* parity )
{ EL_TRY( *parity = PermutationParity( *CReflect(p) ) ) }

// Pivot parity
// ============
ElError ElPivotParity
( ElConstMatrix_i p, ElInt pivotOffset, bool* parity )
{ EL_TRY( *parity = PivotParity( *CReflect(p), pivotOffset ) ) }

ElError ElPivotParityDist
( ElConstDistMatrix_i p, ElInt pivotOffset, bool* parity )
{ EL_TRY( *parity = PivotParity( *CReflect(p), pivotOffset ) ) }

// Convert a pivot sequence to partial permutation vectors
// =======================================================
ElError ElPivotsToPartialPermutation
( ElConstMatrix_i pivots, ElMatrix_i p, ElMatrix_i pInv, ElInt offset )
{ EL_TRY( PivotsToPartialPermutation
          ( *CReflect(pivots), *CReflect(p), *CReflect(pInv), offset ) ) }

ElError ElPivotsToPartialPermutationDist
( ElConstDistMatrix_i pivots, ElDistMatrix_i p, ElDistMatrix_i pInv, 
  ElInt offset )
{ EL_TRY( PivotsToPartialPermutation
          ( *CReflect(pivots), *CReflect(p), *CReflect(pInv), offset ) ) }

// Convert a pivot sequence to a permutation vector
// ================================================  
ElError ElPivotsToPermutation
( ElConstMatrix_i pivots, ElMatrix_i p, ElInt offset )
{ EL_TRY( PivotsToPermutation( *CReflect(pivots), *CReflect(p), offset ) ) }

ElError ElPivotsToPermutationDist
( ElConstDistMatrix_i pivots, ElDistMatrix_i p, ElInt offset )
{ EL_TRY( PivotsToPermutation( *CReflect(pivots), *CReflect(p), offset ) ) }

ElError ElPivotsToInversePermutation
( ElConstMatrix_i pivots, ElMatrix_i p, ElInt offset )
{ EL_TRY( PivotsToInversePermutation
          ( *CReflect(pivots), *CReflect(p), offset ) ) }

ElError ElPivotsToInversePermutationDist
( ElConstDistMatrix_i pivots, ElDistMatrix_i p, ElInt offset )
{ EL_TRY( PivotsToInversePermutation
          ( *CReflect(pivots), *CReflect(p), offset ) ) }

#define C_PROTO_BASE(SIG,SIGBASE,T) \
  /* Apply column pivots
     =================== */ \
  ElError ElApplyColPivots_ ## SIG \
  ( ElMatrix_ ## SIG A, ElConstMatrix_i pivots, ElInt offset ) \
  { EL_TRY( ApplyColPivots( *CReflect(A), *CReflect(pivots), offset ) ) } \
  ElError ElApplyColPivotsDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElConstDistMatrix_i pivots, ElInt offset ) \
  { EL_TRY( ApplyColPivots( *CReflect(A), *CReflect(pivots), offset ) ) } \
  ElError ElApplyInverseColPivots_ ## SIG \
  ( ElMatrix_ ## SIG A, ElConstMatrix_i pivots, ElInt offset ) \
  { EL_TRY( ApplyInverseColPivots \
           ( *CReflect(A), *CReflect(pivots), offset ) ) } \
  ElError ElApplyInverseColPivotsDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElConstDistMatrix_i pivots, ElInt offset ) \
  { EL_TRY( ApplyInverseColPivots \
            ( *CReflect(A), *CReflect(pivots), offset ) ) } \
  /* Apply row pivots
     ================ */ \
  ElError ElApplyRowPivots_ ## SIG \
  ( ElMatrix_ ## SIG A, ElConstMatrix_i pivots, ElInt offset ) \
  { EL_TRY( ApplyRowPivots( *CReflect(A), *CReflect(pivots), offset ) ) } \
  ElError ElApplyRowPivotsDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElConstDistMatrix_i pivots, ElInt offset ) \
  { EL_TRY( ApplyRowPivots( *CReflect(A), *CReflect(pivots), offset ) ) } \
  ElError ElApplyInverseRowPivots_ ## SIG \
  ( ElMatrix_ ## SIG A, ElConstMatrix_i pivots, ElInt offset ) \
  { EL_TRY( ApplyInverseRowPivots \
           ( *CReflect(A), *CReflect(pivots), offset ) ) } \
  ElError ElApplyInverseRowPivotsDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElConstDistMatrix_i pivots, ElInt offset ) \
  { EL_TRY( ApplyInverseRowPivots \
            ( *CReflect(A), *CReflect(pivots), offset ) ) } \
  /* Permute columns
     =============== */ \
  ElError ElPermuteCols_ ## SIG \
  ( ElMatrix_ ## SIG A, ElConstMatrix_i p ) \
  { EL_TRY( PermuteCols( *CReflect(A), *CReflect(p) ) ) } \
  ElError ElPermuteColsDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElConstDistMatrix_i p ) \
  { EL_TRY( PermuteCols( *CReflect(A), *CReflect(p) ) ) } \
  ElError ElInversePermuteCols_ ## SIG \
  ( ElMatrix_ ## SIG A, ElConstMatrix_i p ) \
  { EL_TRY( InversePermuteCols( *CReflect(A), *CReflect(p) ) ) } \
  ElError ElInversePermuteColsDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElConstDistMatrix_i p ) \
  { EL_TRY( InversePermuteCols( *CReflect(A), *CReflect(p) ) ) } \
  ElError ElPermuteColsBoth_ ## SIG \
  ( ElMatrix_ ## SIG A, ElConstMatrix_i p, ElConstMatrix_i pInv ) \
  { EL_TRY( PermuteCols( *CReflect(A), *CReflect(p), *CReflect(pInv) ) ) } \
  ElError ElPermuteColsBothDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElConstDistMatrix_i p, ElConstDistMatrix_i pInv ) \
  { EL_TRY( PermuteCols( *CReflect(A), *CReflect(p), *CReflect(pInv) ) ) } \
  ElError ElPermuteColsMetaDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, const ElPermutationMeta* meta ) \
  { EL_TRY( PermuteCols( *CReflect(A), CReflect(*meta) ) ) } \
  /* Permute rows
     ============ */ \
  ElError ElPermuteRows_ ## SIG \
  ( ElMatrix_ ## SIG A, ElConstMatrix_i p ) \
  { EL_TRY( PermuteRows( *CReflect(A), *CReflect(p) ) ) } \
  ElError ElPermuteRowsDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElConstDistMatrix_i p ) \
  { EL_TRY( PermuteRows( *CReflect(A), *CReflect(p) ) ) } \
  ElError ElInversePermuteRows_ ## SIG \
  ( ElMatrix_ ## SIG A, ElConstMatrix_i p ) \
  { EL_TRY( InversePermuteRows( *CReflect(A), *CReflect(p) ) ) } \
  ElError ElInversePermuteRowsDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElConstDistMatrix_i p ) \
  { EL_TRY( InversePermuteRows( *CReflect(A), *CReflect(p) ) ) } \
  ElError ElPermuteRowsBoth_ ## SIG \
  ( ElMatrix_ ## SIG A, ElConstMatrix_i p, ElConstMatrix_i pInv ) \
  { EL_TRY( PermuteRows( *CReflect(A), *CReflect(p), *CReflect(pInv) ) ) } \
  ElError ElPermuteRowsBothDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, ElConstDistMatrix_i p, ElConstDistMatrix_i pInv ) \
  { EL_TRY( PermuteRows( *CReflect(A), *CReflect(p), *CReflect(pInv) ) ) } \
  ElError ElPermuteRowsMetaDist_ ## SIG \
  ( ElDistMatrix_ ## SIG A, const ElPermutationMeta* meta ) \
  { EL_TRY( PermuteRows( *CReflect(A), CReflect(*meta) ) ) }

#define C_PROTO_NOCOMPLEX(SIG,Real) \
  C_PROTO_BASE(SIG,SIG,Real) \
  /* Apply symmetric pivots
     ====================== */ \
  ElError ElApplySymmetricPivots_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElConstMatrix_i p, \
    ElInt offset ) \
  { EL_TRY( ApplySymmetricPivots \
           ( CReflect(uplo), *CReflect(A), *CReflect(p), false, offset ) ) } \
  ElError ElApplySymmetricPivotsDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElConstDistMatrix_i p, \
    ElInt offset ) \
  { EL_TRY( ApplySymmetricPivots \
           ( CReflect(uplo), *CReflect(A), *CReflect(p), false, offset ) ) } \
  ElError ElApplyInverseSymmetricPivots_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElConstMatrix_i p, \
    ElInt offset ) \
  { EL_TRY( ApplyInverseSymmetricPivots \
           ( CReflect(uplo), *CReflect(A), *CReflect(p), false, offset ) ) } \
  ElError ElApplyInverseSymmetricPivotsDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElConstDistMatrix_i p, \
    ElInt offset ) \
  { EL_TRY( ApplyInverseSymmetricPivots \
           ( CReflect(uplo), *CReflect(A), *CReflect(p), false, offset ) ) }

#define C_PROTO_INT(SIG,T) C_PROTO_NOCOMPLEX(SIG,T)
#define C_PROTO_REAL(SIG,Real) C_PROTO_NOCOMPLEX(SIG,Real)

#define C_PROTO_COMPLEX(SIG,SIGBASE,F) \
  C_PROTO_BASE(SIG,SIGBASE,F) \
  /* Apply symmetric pivots
     ====================== */ \
  ElError ElApplySymmetricPivots_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElConstMatrix_i p, \
    bool conjugate, ElInt offset ) \
  { EL_TRY( ApplySymmetricPivots \
      ( CReflect(uplo), *CReflect(A), *CReflect(p), conjugate, offset ) ) } \
  ElError ElApplySymmetricPivotsDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElConstDistMatrix_i p, \
    bool conjugate, ElInt offset ) \
  { EL_TRY( ApplySymmetricPivots \
      ( CReflect(uplo), *CReflect(A), *CReflect(p), conjugate, offset ) ) } \
  ElError ElApplyInverseSymmetricPivots_ ## SIG \
  ( ElUpperOrLower uplo, ElMatrix_ ## SIG A, ElConstMatrix_i p, \
    bool conjugate, ElInt offset ) \
  { EL_TRY( ApplyInverseSymmetricPivots \
      ( CReflect(uplo), *CReflect(A), *CReflect(p), conjugate, offset ) ) } \
  ElError ElApplyInverseSymmetricPivotsDist_ ## SIG \
  ( ElUpperOrLower uplo, ElDistMatrix_ ## SIG A, ElConstDistMatrix_i p, \
    bool conjugate, ElInt offset ) \
  { EL_TRY( ApplyInverseSymmetricPivots \
      ( CReflect(uplo), *CReflect(A), *CReflect(p), conjugate, offset ) ) }

#include "El/macros/CInstantiate.h"

} // extern "C"
