/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
#include <El.h>
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
ElError ElDistPermutationCreate( ElDistPermutation* p, ElGrid grid )
{ EL_TRY( *p = CReflect( new DistPermutation(*CReflect(grid)) ) ) }

ElError ElPermutationDestroy( ElConstPermutation p )
{ EL_TRY( delete CReflect(p) ) }
ElError ElDistPermutationDestroy( ElConstDistPermutation p )
{ EL_TRY( delete CReflect(p) ) }

ElError ElPermutationEmpty( ElPermutation P )
{ EL_TRY( CReflect(P)->Empty() ) }
ElError ElDistPermutationEmpty( ElDistPermutation P )
{ EL_TRY( CReflect(P)->Empty() ) }

ElError ElPermutationMakeIdentity( ElPermutation P, ElInt size )
{ EL_TRY( CReflect(P)->MakeIdentity(size) ) }
ElError ElDistPermutationMakeIdentity( ElDistPermutation P, ElInt size )
{ EL_TRY( CReflect(P)->MakeIdentity(size) ) }

ElError ElPermutationReserveSwaps( ElPermutation P, ElInt maxSwaps )
{ EL_TRY( CReflect(P)->ReserveSwaps(maxSwaps) ) }
ElError ElDistPermutationReserveSwaps( ElDistPermutation P, ElInt maxSwaps )
{ EL_TRY( CReflect(P)->ReserveSwaps(maxSwaps) ) }

ElError ElPermutationMakeArbitrary( ElPermutation P )
{ EL_TRY( CReflect(P)->MakeArbitrary() ) }
ElError ElDistPermutationMakeArbitrary( ElDistPermutation P )
{ EL_TRY( CReflect(P)->MakeArbitrary() ) }

ElError ElPermutationSwap
( ElPermutation P, ElInt origin, ElInt dest )
{ EL_TRY( CReflect(P)->Swap(origin,dest) ) }
ElError ElDistPermutationSwap
( ElDistPermutation P, ElInt origin, ElInt dest )
{ EL_TRY( CReflect(P)->Swap(origin,dest) ) }

ElError ElPermutationSwapSequence
( ElPermutation P, ElConstPermutation PAppend, ElInt offset )
{ EL_TRY( CReflect(P)->SwapSequence(*CReflect(PAppend),offset) ) }
ElError ElDistPermutationSwapSequence
( ElDistPermutation P, ElConstDistPermutation PAppend, ElInt offset )
{ EL_TRY( CReflect(P)->SwapSequence(*CReflect(PAppend),offset) ) }

ElError ElPermutationSetImage
( ElPermutation P, ElInt origin, ElInt dest )
{ EL_TRY( CReflect(P)->SetImage(origin,dest) ) }
ElError ElDistPermutationSetImage
( ElDistPermutation P, ElInt origin, ElInt dest )
{ EL_TRY( CReflect(P)->SetImage(origin,dest) ) }

ElError ElPermutationHeight( ElConstPermutation P, ElInt* height )
{ EL_TRY( *height = CReflect(P)->Height() ) }
ElError ElDistPermutationHeight( ElConstDistPermutation P, ElInt* height )
{ EL_TRY( *height = CReflect(P)->Height() ) }

ElError ElPermutationWidth( ElConstPermutation P, ElInt* width )
{ EL_TRY( *width = CReflect(P)->Width() ) }
ElError ElDistPermutationWidth( ElConstDistPermutation P, ElInt* width )
{ EL_TRY( *width = CReflect(P)->Width() ) }

ElError ElPermutationParity( ElConstPermutation P, bool* parity )
{ EL_TRY( *parity = CReflect(P)->Parity() ) }
ElError ElDistPermutationParity( ElConstDistPermutation P, bool* parity )
{ EL_TRY( *parity = CReflect(P)->Parity() ) }

ElError ElPermutationIsSwapSequence
( ElConstPermutation P, bool* isSwap )
{ EL_TRY( *isSwap = CReflect(P)->IsSwapSequence() ) }
ElError ElDistPermutationIsSwapSequence
( ElConstDistPermutation P, bool* isSwap )
{ EL_TRY( *isSwap = CReflect(P)->IsSwapSequence() ) }

ElError ElPermutationIsImplicitSwapSequence
( ElConstPermutation P, bool* isImplicit )
{ EL_TRY( *isImplicit = CReflect(P)->IsImplicitSwapSequence() ) }
ElError ElDistPermutationIsImplicitSwapSequence
( ElConstDistPermutation P, bool* isImplicit )
{ EL_TRY( *isImplicit = CReflect(P)->IsImplicitSwapSequence() ) }

ElError ElPermutationPermuteRows_i
( ElConstPermutation P, ElMatrix_i A, ElInt offset )
{ EL_TRY( CReflect(P)->PermuteRows(*CReflect(A),offset) ) }
ElError ElPermutationPermuteRows_s
( ElConstPermutation P, ElMatrix_s A, ElInt offset )
{ EL_TRY( CReflect(P)->PermuteRows(*CReflect(A),offset) ) }
ElError ElPermutationPermuteRows_d
( ElConstPermutation P, ElMatrix_d A, ElInt offset )
{ EL_TRY( CReflect(P)->PermuteRows(*CReflect(A),offset) ) }
ElError ElPermutationPermuteRows_c
( ElConstPermutation P, ElMatrix_c A, ElInt offset )
{ EL_TRY( CReflect(P)->PermuteRows(*CReflect(A),offset) ) }
ElError ElPermutationPermuteRows_z
( ElConstPermutation P, ElMatrix_z A, ElInt offset )
{ EL_TRY( CReflect(P)->PermuteRows(*CReflect(A),offset) ) }

ElError ElDistPermutationPermuteRows_i
( ElConstDistPermutation P, ElDistMatrix_i A, ElInt offset )
{ EL_TRY( CReflect(P)->PermuteRows(*CReflect(A),offset) ) }
ElError ElDistPermutationPermuteRows_s
( ElConstDistPermutation P, ElDistMatrix_s A, ElInt offset )
{ EL_TRY( CReflect(P)->PermuteRows(*CReflect(A),offset) ) }
ElError ElDistPermutationPermuteRows_d
( ElConstDistPermutation P, ElDistMatrix_d A, ElInt offset )
{ EL_TRY( CReflect(P)->PermuteRows(*CReflect(A),offset) ) }
ElError ElDistPermutationPermuteRows_c
( ElConstDistPermutation P, ElDistMatrix_c A, ElInt offset )
{ EL_TRY( CReflect(P)->PermuteRows(*CReflect(A),offset) ) }
ElError ElDistPermutationPermuteRows_z
( ElConstDistPermutation P, ElDistMatrix_z A, ElInt offset )
{ EL_TRY( CReflect(P)->PermuteRows(*CReflect(A),offset) ) }

ElError ElPermutationPermuteCols_i
( ElConstPermutation P, ElMatrix_i A, ElInt offset )
{ EL_TRY( CReflect(P)->PermuteCols(*CReflect(A),offset) ) }
ElError ElPermutationPermuteCols_s
( ElConstPermutation P, ElMatrix_s A, ElInt offset )
{ EL_TRY( CReflect(P)->PermuteCols(*CReflect(A),offset) ) }
ElError ElPermutationPermuteCols_d
( ElConstPermutation P, ElMatrix_d A, ElInt offset )
{ EL_TRY( CReflect(P)->PermuteCols(*CReflect(A),offset) ) }
ElError ElPermutationPermuteCols_c
( ElConstPermutation P, ElMatrix_c A, ElInt offset )
{ EL_TRY( CReflect(P)->PermuteCols(*CReflect(A),offset) ) }
ElError ElPermutationPermuteCols_z
( ElConstPermutation P, ElMatrix_z A, ElInt offset )
{ EL_TRY( CReflect(P)->PermuteCols(*CReflect(A),offset) ) }

ElError ElDistPermutationPermuteCols_i
( ElConstDistPermutation P, ElDistMatrix_i A, ElInt offset )
{ EL_TRY( CReflect(P)->PermuteCols(*CReflect(A),offset) ) }
ElError ElDistPermutationPermuteCols_s
( ElConstDistPermutation P, ElDistMatrix_s A, ElInt offset )
{ EL_TRY( CReflect(P)->PermuteCols(*CReflect(A),offset) ) }
ElError ElDistPermutationPermuteCols_d
( ElConstDistPermutation P, ElDistMatrix_d A, ElInt offset )
{ EL_TRY( CReflect(P)->PermuteCols(*CReflect(A),offset) ) }
ElError ElDistPermutationPermuteCols_c
( ElConstDistPermutation P, ElDistMatrix_c A, ElInt offset )
{ EL_TRY( CReflect(P)->PermuteCols(*CReflect(A),offset) ) }
ElError ElDistPermutationPermuteCols_z
( ElConstDistPermutation P, ElDistMatrix_z A, ElInt offset )
{ EL_TRY( CReflect(P)->PermuteCols(*CReflect(A),offset) ) }

ElError ElPermutationPermuteSymmetrically_i
( ElConstPermutation P,
  ElUpperOrLower uplo, ElMatrix_i A, ElInt offset )
{ EL_TRY(
    CReflect(P)->PermuteSymmetrically
      (CReflect(uplo),*CReflect(A),false,offset)
  )
}
ElError ElPermutationPermuteSymmetrically_s
( ElConstPermutation P,
  ElUpperOrLower uplo, ElMatrix_s A, ElInt offset )
{ EL_TRY(
    CReflect(P)->PermuteSymmetrically
      (CReflect(uplo),*CReflect(A),false,offset)
  )
}
ElError ElPermutationPermuteSymmetrically_d
( ElConstPermutation P,
  ElUpperOrLower uplo, ElMatrix_d A, ElInt offset )
{ EL_TRY(
    CReflect(P)->PermuteSymmetrically
      (CReflect(uplo),*CReflect(A),false,offset)
  )
}
ElError ElPermutationPermuteSymmetrically_c
( ElConstPermutation P,
  ElUpperOrLower uplo, ElMatrix_c A, bool conjugate, ElInt offset )
{ EL_TRY(
    CReflect(P)->PermuteSymmetrically
      (CReflect(uplo),*CReflect(A),conjugate,offset)
  )
}
ElError ElPermutationPermuteSymmetrically_z
( ElConstPermutation P,
  ElUpperOrLower uplo, ElMatrix_z A, bool conjugate, ElInt offset )
{ EL_TRY(
    CReflect(P)->PermuteSymmetrically
      (CReflect(uplo),*CReflect(A),conjugate,offset)
  )
}

ElError ElDistPermutationPermuteSymmetrically_i
( ElConstDistPermutation P,
  ElUpperOrLower uplo, ElDistMatrix_i A, ElInt offset )
{ EL_TRY(
    CReflect(P)->PermuteSymmetrically
      (CReflect(uplo),*CReflect(A),false,offset)
  )
}
ElError ElDistPermutationPermuteSymmetrically_s
( ElConstDistPermutation P,
  ElUpperOrLower uplo, ElDistMatrix_s A, ElInt offset )
{ EL_TRY(
    CReflect(P)->PermuteSymmetrically
      (CReflect(uplo),*CReflect(A),false,offset)
  )
}
ElError ElDistPermutationPermuteSymmetrically_d
( ElConstDistPermutation P,
  ElUpperOrLower uplo, ElDistMatrix_d A, ElInt offset )
{ EL_TRY(
    CReflect(P)->PermuteSymmetrically
      (CReflect(uplo),*CReflect(A),false,offset)
  )
}
ElError ElDistPermutationPermuteSymmetrically_c
( ElConstDistPermutation P,
  ElUpperOrLower uplo, ElDistMatrix_c A, bool conjugate, ElInt offset )
{ EL_TRY(
    CReflect(P)->PermuteSymmetrically
      (CReflect(uplo),*CReflect(A),conjugate,offset)
  )
}
ElError ElDistPermutationPermuteSymmetrically_z
( ElConstDistPermutation P,
  ElUpperOrLower uplo, ElDistMatrix_z A, bool conjugate, ElInt offset )
{ EL_TRY(
    CReflect(P)->PermuteSymmetrically
      (CReflect(uplo),*CReflect(A),conjugate,offset)
  )
}

ElError ElPermutationExplicitVector
( ElConstPermutation P, ElMatrix_i pVec )
{ EL_TRY( CReflect(P)->ExplicitVector(*CReflect(pVec)) ) }
ElError ElDistPermutationExplicitVector
( ElConstDistPermutation P, ElDistMatrix_i pVec )
{ EL_TRY( CReflect(P)->ExplicitVector(*CReflect(pVec)) ) }

ElError ElPermutationExplicitMatrix
( ElConstPermutation P, ElMatrix_i pMat )
{ EL_TRY( CReflect(P)->ExplicitMatrix(*CReflect(pMat)) ) }
ElError ElDistPermutationExplicitMatrix
( ElConstDistPermutation P, ElDistMatrix_i pMat )
{ EL_TRY( CReflect(P)->ExplicitMatrix(*CReflect(pMat)) ) }

} // extern "C"
