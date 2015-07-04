/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./Trrk/Local.hpp"
#include "./Trrk/NN.hpp"
#include "./Trrk/NT.hpp"
#include "./Trrk/TN.hpp"
#include "./Trrk/TT.hpp"

namespace El {

template<typename Ring>
void Trrk
( UpperOrLower uplo, 
  Orientation orientA, Orientation orientB,
  Ring alpha, const Matrix<Ring>& A, const Matrix<Ring>& B,
  Ring beta,        Matrix<Ring>& C )
{
    DEBUG_ONLY(CSE cse("Trrk"))
    ScaleTrapezoid( beta, uplo, C );
    if( orientA==NORMAL && orientB==NORMAL )
        trrk::TrrkNN( uplo, alpha, A, B, C );
    else if( orientA==NORMAL )
        trrk::TrrkNT( uplo, orientB, alpha, A, B, C );
    else if( orientB==NORMAL )
        trrk::TrrkTN( uplo, orientA, alpha, A, B, C );
    else
        trrk::TrrkTT( uplo, orientA, orientB, alpha, A, B, C );
}

template<typename Ring>
void Trrk
( UpperOrLower uplo, Orientation orientA, Orientation orientB,
  Ring alpha, const AbstractDistMatrix<Ring>& A, 
              const AbstractDistMatrix<Ring>& B,
  Ring beta,        AbstractDistMatrix<Ring>& C )
{
    DEBUG_ONLY(CSE cse("Trrk"))
    ScaleTrapezoid( beta, uplo, C );
    if( orientA==NORMAL && orientB==NORMAL )
        trrk::TrrkNN( uplo, alpha, A, B, C );
    else if( orientA==NORMAL )
        trrk::TrrkNT( uplo, orientB, alpha, A, B, C );
    else if( orientB==NORMAL )
        trrk::TrrkTN( uplo, orientA, alpha, A, B, C );
    else
        trrk::TrrkTT( uplo, orientA, orientB, alpha, A, B, C );
}

#define PROTO(Ring) \
  template void Trrk \
  ( UpperOrLower uplo, \
    Orientation orientA, Orientation orientB, \
    Ring alpha, const Matrix<Ring>& A, \
                const Matrix<Ring>& B, \
    Ring beta,        Matrix<Ring>& C ); \
  template void Trrk \
  ( UpperOrLower uplo, \
    Orientation orientA, Orientation orientB, \
    Ring alpha, const AbstractDistMatrix<Ring>& A, \
                const AbstractDistMatrix<Ring>& B, \
    Ring beta,        AbstractDistMatrix<Ring>& C ); \
  template void LocalTrrk \
   ( UpperOrLower uplo, \
     Ring alpha, const DistMatrix<Ring,MC,  STAR>& A, \
                 const DistMatrix<Ring,STAR,MR  >& B, \
     Ring beta,        DistMatrix<Ring>& C ); \
  template void LocalTrrk \
  ( UpperOrLower uplo, Orientation orientB, \
    Ring alpha, const DistMatrix<Ring,MC,STAR>& A, \
                const DistMatrix<Ring,MR,STAR>& B, \
    Ring beta,        DistMatrix<Ring>& C ); \
  template void LocalTrrk \
  ( UpperOrLower uplo, Orientation orientA, \
    Ring alpha, const DistMatrix<Ring,STAR,MC>& A, \
                const DistMatrix<Ring,STAR,MR>& B, \
    Ring beta,        DistMatrix<Ring>& C ); \
  template void LocalTrrk \
  ( UpperOrLower uplo, \
    Orientation orientA, Orientation orientB, \
    Ring alpha, const DistMatrix<Ring,STAR,MC  >& A, \
                const DistMatrix<Ring,MR,  STAR>& B, \
    Ring beta,        DistMatrix<Ring>& C );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
