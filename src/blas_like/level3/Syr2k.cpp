/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like/level1.hpp>
#include <El/blas_like/level3.hpp>

#include "./Syr2k/LN.hpp"
#include "./Syr2k/LT.hpp"
#include "./Syr2k/UN.hpp"
#include "./Syr2k/UT.hpp"

namespace El {

template<typename T>
void Syr2k
( UpperOrLower uplo, Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
  T beta,        Matrix<T>& C,
  bool conjugate )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( orientation == NORMAL )
      {
          if( A.Height() != C.Height() || A.Height() != C.Width() ||
              B.Height() != C.Height() ||B.Height() != C.Width()    )
              LogicError("Nonconformal Syr2k");
      }
      else 
      {
          if( A.Width() != C.Height() || A.Width() != C.Width() ||
              B.Width() != C.Height() || B.Width() != C.Width()   )
              LogicError("Nonconformal Syr2k");
      }
    )
    const char uploChar = UpperOrLowerToChar( uplo );
    const char transChar = OrientationToChar( orientation );
    const Int k = ( orientation == NORMAL ? A.Width() : A.Height() );
    if( conjugate )
    {
        blas::Her2k
        ( uploChar, transChar, C.Height(), k,
          alpha,          A.LockedBuffer(), A.LDim(),
                          B.LockedBuffer(), B.LDim(),
          RealPart(beta), C.Buffer(),       C.LDim() );
    }
    else
    {
        blas::Syr2k
        ( uploChar, transChar, C.Height(), k,
          alpha, A.LockedBuffer(), A.LDim(),
                 B.LockedBuffer(), B.LDim(),
          beta,  C.Buffer(),       C.LDim() );
    }
}

template<typename T>
void Syr2k
( UpperOrLower uplo, Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& B,
                 Matrix<T>& C,
  bool conjugate )
{
    EL_DEBUG_CSE
    const Int n = ( orientation==NORMAL ? A.Height() : A.Width() );
    C.Resize( n, n );
    Zero( C );
    Syr2k( uplo, orientation, alpha, A, B, T(0), C, conjugate );
}

template<typename T>
void Syr2k
( UpperOrLower uplo, Orientation orientation,
  T alpha, const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B,
  T beta,        AbstractDistMatrix<T>& C, bool conjugate )
{
    EL_DEBUG_CSE
    ScaleTrapezoid( beta, uplo, C );
    if( uplo == LOWER && orientation == NORMAL )
        syr2k::LN( alpha, A, B, C, conjugate );
    else if( uplo == LOWER )
        syr2k::LT( alpha, A, B, C, conjugate );
    else if( orientation == NORMAL )
        syr2k::UN( alpha, A, B, C, conjugate );
    else
        syr2k::UT( alpha, A, B, C, conjugate );
}

template<typename T>
void Syr2k
( UpperOrLower uplo, Orientation orientation,
  T alpha, const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B,
                 AbstractDistMatrix<T>& C, bool conjugate )
{
    EL_DEBUG_CSE
    const Int n = ( orientation==NORMAL ? A.Height() : A.Width() );
    C.Resize( n, n );
    Zero( C );
    Syr2k( uplo, orientation, alpha, A, B, T(0), C, conjugate );
}

#define PROTO(T) \
  template void Syr2k \
  ( UpperOrLower uplo, Orientation orientation, \
    T alpha, const Matrix<T>& A, const Matrix<T>& B, \
    T beta,        Matrix<T>& C, bool conjugate ); \
  template void Syr2k \
  ( UpperOrLower uplo, Orientation orientation, \
    T alpha, const Matrix<T>& A, const Matrix<T>& B, \
                   Matrix<T>& C, bool conjugate ); \
  template void Syr2k \
  ( UpperOrLower uplo, Orientation orientation, \
    T alpha, const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B, \
    T beta,        AbstractDistMatrix<T>& C, bool conjugate ); \
  template void Syr2k \
  ( UpperOrLower uplo, Orientation orientation, \
    T alpha, const AbstractDistMatrix<T>& A, const AbstractDistMatrix<T>& B, \
                   AbstractDistMatrix<T>& C, bool conjugate );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
