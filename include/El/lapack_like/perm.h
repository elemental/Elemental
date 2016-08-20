/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_PERM_C_H
#define EL_PERM_C_H

#include <El/core/DistMatrix.h>

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
( const ElDistMatrix_i p,
  const ElDistMatrix_i pInv,
        ElInt permAlign,
        MPI_Comm comm,
        ElPermutationMeta* meta );
EL_EXPORT ElError ElPermutationMetaClear( const ElPermutationMeta* meta );

EL_EXPORT ElError ElPermutationMetaTotalSend
( const ElPermutationMeta* meta, int* total );
EL_EXPORT ElError ElPermutationMetaTotalRecv
( const ElPermutationMeta* meta, int* total );

EL_EXPORT ElError ElPermutationMetaScaleUp
( ElPermutationMeta* meta, ElInt length );
EL_EXPORT ElError ElPermutationMetaScaleDown
( ElPermutationMeta* meta, ElInt length );

/* Anonymous placeholders for Permutation and DistPermutation
   ---------------------------------------------------------- */
typedef struct ElPermutationDummy* ElPermutation;
typedef struct ElDistPermutationDummy* ElDistPermutation;

typedef const struct ElPermutationDummy* ElConstPermutation;
typedef const struct ElDistPermutationDummy* ElConstDistPermutation;

EL_EXPORT ElError ElPermutationCreate( ElPermutation* P );
EL_EXPORT ElError ElDistPermutationCreate( ElDistPermutation* P, ElGrid g );

EL_EXPORT ElError ElPermutationDestroy( ElConstPermutation P );
EL_EXPORT ElError ElDistPermutationDestroy( ElConstDistPermutation P );

EL_EXPORT ElError ElPermutationEmpty( ElPermutation P );
EL_EXPORT ElError ElDistPermutationEmpty( ElDistPermutation P );

EL_EXPORT ElError ElPermutationMakeIdentity
( ElPermutation P, ElInt size );
EL_EXPORT ElError ElDistPermutationMakeIdentity
( ElDistPermutation P, ElInt size );

EL_EXPORT ElError ElPermutationReserveSwaps
( ElPermutation P, ElInt maxSwaps );
EL_EXPORT ElError ElDistPermutationReserveSwaps
( ElDistPermutation P, ElInt maxSwaps );

EL_EXPORT ElError ElPermutationMakeArbitrary( ElPermutation P );
EL_EXPORT ElError ElDistPermutationMakeArbitrary( ElDistPermutation P );

EL_EXPORT ElError ElPermutationSwap
( ElPermutation P, ElInt origin, ElInt dest );
EL_EXPORT ElError ElDistPermutationSwap
( ElDistPermutation P, ElInt origin, ElInt dest );

EL_EXPORT ElError ElPermutationSwapSequence
( ElPermutation P, ElConstPermutation PAppend, ElInt offset );
EL_EXPORT ElError ElDistPermutationSwapSequence
( ElDistPermutation P, ElConstDistPermutation PAppend, ElInt offset );

/* TODO: Support for ElMatrix_i and ElDistMatrix_i swap vectors */

EL_EXPORT ElError ElPermutationSetImage
( ElPermutation P, ElInt origin, ElInt dest );
EL_EXPORT ElError ElDistPermutationSetImage
( ElDistPermutation P, ElInt origin, ElInt dest );

EL_EXPORT ElError ElPermutationHeight
( ElConstPermutation P, ElInt* height );
EL_EXPORT ElError ElDistPermutationHeight
( ElConstDistPermutation P, ElInt* height );
EL_EXPORT ElError ElPermutationWidth
( ElConstPermutation P, ElInt* width );
EL_EXPORT ElError ElDistPermutationWidth
( ElConstDistPermutation P, ElInt* width );

EL_EXPORT ElError ElPermutationParity
( ElConstPermutation P, bool* parity );
EL_EXPORT ElError ElDistPermutationParity
( ElConstDistPermutation P, bool* parity );

EL_EXPORT ElError ElPermutationIsSwapSequence
( ElConstPermutation P, bool* isSwap );
EL_EXPORT ElError ElDistPermutationIsSwapSequence
( ElConstDistPermutation P, bool* isSwap );

EL_EXPORT ElError ElPermutationIsImplicitSwapSequence
( ElConstPermutation P, bool* isImplicit );
EL_EXPORT ElError ElDistPermutationIsImplicitSwapSequence
( ElConstDistPermutation P, bool* isImplicit );

EL_EXPORT ElError ElPermutationPermuteCols_i
( ElConstPermutation P, ElMatrix_i A, ElInt offset );
EL_EXPORT ElError ElPermutationPermuteCols_s
( ElConstPermutation P, ElMatrix_s A, ElInt offset );
EL_EXPORT ElError ElPermutationPermuteCols_d
( ElConstPermutation P, ElMatrix_d A, ElInt offset );
EL_EXPORT ElError ElPermutationPermuteCols_c
( ElConstPermutation P, ElMatrix_c A, ElInt offset );
EL_EXPORT ElError ElPermutationPermuteCols_z
( ElConstPermutation P, ElMatrix_z A, ElInt offset );

EL_EXPORT ElError ElDistPermutationPermuteCols_i
( ElConstDistPermutation P, ElDistMatrix_i A, ElInt offset );
EL_EXPORT ElError ElDistPermutationPermuteCols_s
( ElConstDistPermutation P, ElDistMatrix_s A, ElInt offset );
EL_EXPORT ElError ElDistPermutationPermuteCols_d
( ElConstDistPermutation P, ElDistMatrix_d A, ElInt offset );
EL_EXPORT ElError ElDistPermutationPermuteCols_c
( ElConstDistPermutation P, ElDistMatrix_c A, ElInt offset );
EL_EXPORT ElError ElDistPermutationPermuteCols_z
( ElConstDistPermutation P, ElDistMatrix_z A, ElInt offset );

EL_EXPORT ElError ElPermutationPermuteRows_i
( ElConstPermutation P, ElMatrix_i A, ElInt offset );
EL_EXPORT ElError ElPermutationPermuteRows_s
( ElConstPermutation P, ElMatrix_s A, ElInt offset );
EL_EXPORT ElError ElPermutationPermuteRows_d
( ElConstPermutation P, ElMatrix_d A, ElInt offset );
EL_EXPORT ElError ElPermutationPermuteRows_c
( ElConstPermutation P, ElMatrix_c A, ElInt offset );
EL_EXPORT ElError ElPermutationPermuteRows_z
( ElConstPermutation P, ElMatrix_z A, ElInt offset );

EL_EXPORT ElError ElDistPermutationPermuteRows_i
( ElConstDistPermutation P, ElDistMatrix_i A, ElInt offset );
EL_EXPORT ElError ElDistPermutationPermuteRows_s
( ElConstDistPermutation P, ElDistMatrix_s A, ElInt offset );
EL_EXPORT ElError ElDistPermutationPermuteRows_d
( ElConstDistPermutation P, ElDistMatrix_d A, ElInt offset );
EL_EXPORT ElError ElDistPermutationPermuteRows_c
( ElConstDistPermutation P, ElDistMatrix_c A, ElInt offset );
EL_EXPORT ElError ElDistPermutationPermuteRows_z
( ElConstDistPermutation P, ElDistMatrix_z A, ElInt offset );

EL_EXPORT ElError ElPermutationPermuteSymmetrically_i
( ElConstPermutation P, ElUpperOrLower uplo, ElMatrix_i A, ElInt offset );
EL_EXPORT ElError ElPermutationPermuteSymmetrically_s
( ElConstPermutation P, ElUpperOrLower uplo, ElMatrix_s A, ElInt offset );
EL_EXPORT ElError ElPermutationPermuteSymmetrically_d
( ElConstPermutation P, ElUpperOrLower uplo, ElMatrix_d A, ElInt offset );
EL_EXPORT ElError ElPermutationPermuteSymmetrically_c
( ElConstPermutation P, ElUpperOrLower uplo, ElMatrix_c A,
  bool conjugate, ElInt offset );
EL_EXPORT ElError ElPermutationPermuteSymmetrically_z
( ElConstPermutation P, ElUpperOrLower uplo, ElMatrix_z A,
  bool conjugate, ElInt offset );

EL_EXPORT ElError ElDistPermutationPermuteSymmetrically_i
( ElConstDistPermutation P,
  ElUpperOrLower uplo, ElDistMatrix_i A, ElInt offset );
EL_EXPORT ElError ElDistPermutationPermuteSymmetrically_s
( ElConstDistPermutation P,
  ElUpperOrLower uplo, ElDistMatrix_s A, ElInt offset );
EL_EXPORT ElError ElDistPermutationPermuteSymmetrically_d
( ElConstDistPermutation P,
  ElUpperOrLower uplo, ElDistMatrix_d A, ElInt offset );
EL_EXPORT ElError ElDistPermutationPermuteSymmetrically_c
( ElConstDistPermutation P,
  ElUpperOrLower uplo, ElDistMatrix_c A,
  bool conjugate, ElInt offset );
EL_EXPORT ElError ElDistPermutationPermuteSymmetrically_z
( ElConstDistPermutation P,
  ElUpperOrLower uplo, ElDistMatrix_z A,
  bool conjugate, ElInt offset );

EL_EXPORT ElError ElPermutationExplicitVector
( ElConstPermutation P, ElMatrix_i PVec );
EL_EXPORT ElError ElDistPermutationExplicitVector
( ElConstDistPermutation P, ElDistMatrix_i PVec );

EL_EXPORT ElError ElPermutationExplicitMatrix
( ElConstPermutation P, ElMatrix_i PMat );
EL_EXPORT ElError ElDistPermutationExplicitMatrix
( ElConstDistPermutation P, ElDistMatrix_i PMat );

#ifdef __cplusplus
} // extern "C"
#endif

#ifdef __cplusplus
#include <El/lapack_like/perm/CReflect.hpp>
#endif

#endif /* ifndef EL_PERM_C_H */
