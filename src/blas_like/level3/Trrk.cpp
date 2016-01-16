/*
   Copyright (c) 2009-2016, Jack Poulson
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

template<typename T>
void Trrk
( UpperOrLower uplo, 
  Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
  T beta,        Matrix<T>& C )
{
    DEBUG_ONLY(CSE cse("Trrk"))
    ScaleTrapezoid( beta, uplo, C );
    if( orientationOfA==NORMAL && orientationOfB==NORMAL )
        trrk::TrrkNN( uplo, alpha, A, B, C );
    else if( orientationOfA==NORMAL )
        trrk::TrrkNT( uplo, orientationOfB, alpha, A, B, C );
    else if( orientationOfB==NORMAL )
        trrk::TrrkTN( uplo, orientationOfA, alpha, A, B, C );
    else
        trrk::TrrkTT( uplo, orientationOfA, orientationOfB, alpha, A, B, C );
}

template<typename T>
void Trrk
( UpperOrLower uplo, Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const ElementalMatrix<T>& A, const ElementalMatrix<T>& B,
  T beta,        ElementalMatrix<T>& C )
{
    DEBUG_ONLY(CSE cse("Trrk"))
    ScaleTrapezoid( beta, uplo, C );
    if( orientationOfA==NORMAL && orientationOfB==NORMAL )
        trrk::TrrkNN( uplo, alpha, A, B, C );
    else if( orientationOfA==NORMAL )
        trrk::TrrkNT( uplo, orientationOfB, alpha, A, B, C );
    else if( orientationOfB==NORMAL )
        trrk::TrrkTN( uplo, orientationOfA, alpha, A, B, C );
    else
        trrk::TrrkTT( uplo, orientationOfA, orientationOfB, alpha, A, B, C );
}

#define PROTO(T) \
  template void Trrk \
  ( UpperOrLower uplo, \
    Orientation orientationOfA, Orientation orientationOfB, \
    T alpha, const Matrix<T>& A, const Matrix<T>& B, \
    T beta,        Matrix<T>& C ); \
  template void Trrk \
  ( UpperOrLower uplo, \
    Orientation orientationOfA, Orientation orientationOfB, \
    T alpha, const ElementalMatrix<T>& A, const ElementalMatrix<T>& B, \
    T beta,        ElementalMatrix<T>& C ); \
  template void LocalTrrk \
   ( UpperOrLower uplo, \
     T alpha, const DistMatrix<T,MC,  STAR>& A, \
              const DistMatrix<T,STAR,MR  >& B, \
     T beta,        DistMatrix<T>& C ); \
  template void LocalTrrk \
  ( UpperOrLower uplo, Orientation orientationOfB, \
    T alpha, const DistMatrix<T,MC,STAR>& A, \
             const DistMatrix<T,MR,STAR>& B, \
    T beta,        DistMatrix<T>& C ); \
  template void LocalTrrk \
  ( UpperOrLower uplo, Orientation orientationOfA, \
    T alpha, const DistMatrix<T,STAR,MC>& A, \
             const DistMatrix<T,STAR,MR>& B, \
    T beta,        DistMatrix<T>& C ); \
  template void LocalTrrk \
  ( UpperOrLower uplo, \
    Orientation orientationOfA, Orientation orientationOfB, \
    T alpha, const DistMatrix<T,STAR,MC  >& A, \
             const DistMatrix<T,MR,  STAR>& B, \
    T beta,        DistMatrix<T>& C );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
