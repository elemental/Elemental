/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

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
    DEBUG_ONLY(CallStackEntry cse("Trrk"))
    if( orientationOfA==NORMAL && orientationOfB==NORMAL )
        trrk::TrrkNN( uplo, alpha, A, B, beta, C );
    else if( orientationOfA==NORMAL )
        trrk::TrrkNT( uplo, orientationOfB, alpha, A, B, beta, C );
    else if( orientationOfB==NORMAL )
        trrk::TrrkTN( uplo, orientationOfA, alpha, A, B, beta, C );
    else
        trrk::TrrkTT
        ( uplo, orientationOfA, orientationOfB, alpha, A, B, beta, C );
}

template<typename T>
void Trrk
( UpperOrLower uplo, Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
    DEBUG_ONLY(CallStackEntry cse("Trrk"))
    if( orientationOfA==NORMAL && orientationOfB==NORMAL )
        trrk::TrrkNN( uplo, alpha, A, B, beta, C );
    else if( orientationOfA==NORMAL )
        trrk::TrrkNT( uplo, orientationOfB, alpha, A, B, beta, C );
    else if( orientationOfB==NORMAL )
        trrk::TrrkTN( uplo, orientationOfA, alpha, A, B, beta, C );
    else
        trrk::TrrkTT
        ( uplo, orientationOfA, orientationOfB, alpha, A, B, beta, C );
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
    T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B, \
    T beta,        DistMatrix<T>& C ); \
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
    T beta,        DistMatrix<T>& C ); \
  template void trrk::TrrkNN \
  ( UpperOrLower uplo, \
    T alpha, const Matrix<T>& A, const Matrix<T>& B, \
    T beta,        Matrix<T>& C ); \
  template void trrk::TrrkNT \
  ( UpperOrLower uplo, Orientation orientationOfB, \
    T alpha, const Matrix<T>& A, const Matrix<T>& B, \
    T beta,        Matrix<T>& C ); \
  template void trrk::TrrkTN \
  ( UpperOrLower uplo, Orientation orientationOfA, \
    T alpha, const Matrix<T>& A, const Matrix<T>& B, \
    T beta,        Matrix<T>& C ); \
  template void trrk::TrrkTT \
  ( UpperOrLower uplo, \
    Orientation orientationOfA, Orientation orientationOfB, \
    T alpha, const Matrix<T>& A, const Matrix<T>& B, \
    T beta,        Matrix<T>& C );

#ifndef EL_DISABLE_FLOAT
PROTO(float);
#ifndef EL_DISABLE_COMPLEX
PROTO(Complex<float>);
#endif // ifndef EL_DISABLE_COMPLEX
#endif // ifndef EL_DISABLE_FLOAT

PROTO(double);
#ifndef EL_DISABLE_COMPLEX
PROTO(Complex<double>);
#endif // ifndef EL_DISABLE_COMPLEX

} // namespace El
