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

EL_EXPORT ElError ElPermutationCreate( ElPermutation* p );
EL_EXPORT ElError ElDistPermutationCreate( ElDistPermutation* p, ElGrid g );

EL_EXPORT ElError ElPermutationDestroy( ElConstPermutation p );
EL_EXPORT ElError ElDistPermutationDestroy( ElConstDistPermutation p );

/* TODO: Extend interfaces to Permutation and DistPermutation */

#ifdef __cplusplus
} // extern "C"
#endif

#endif /* ifndef EL_PERM_C_H */
