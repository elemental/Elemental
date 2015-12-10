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
( const ElDistMatrix_i p,
  const ElDistMatrix_i pInv,
        ElInt permAlign,
        MPI_Comm comm,
        ElPermutationMeta* meta )
{
    EL_TRY( 
      PermutationMeta metaCpp( *CReflect(p), *CReflect(pInv), permAlign, comm );
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

ElError ElPermutationCreate( ElPermutation* p )
{ EL_TRY( *p = CReflect( new Permutation ) ) }

ElError ElPermutationDestroy( ElConstPermutation p )
{ EL_TRY( delete CReflect(p) ) }

ElError ElDistPermutationCreate( ElDistPermutation* p, ElGrid grid )
{ EL_TRY( *p = CReflect( new DistPermutation(*CReflect(grid)) ) ) }

ElError ElDistPermutationDestroy( ElConstDistPermutation p )
{ EL_TRY( delete CReflect(p) ) }

/* TODO: Extend interfaces to Permutation and DistPermutation */

} // extern "C"
