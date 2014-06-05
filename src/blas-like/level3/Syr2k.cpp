/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

#include "./Syr2k/LN.hpp"
#include "./Syr2k/LT.hpp"
#include "./Syr2k/UN.hpp"
#include "./Syr2k/UT.hpp"

namespace El {

template<typename T>
void Syr2k
( UpperOrLower uplo, Orientation orientation,
  T alpha, const Matrix<T>& A, const Matrix<T>& B, T beta, Matrix<T>& C,
  bool conjugate )
{
    DEBUG_ONLY(
        CallStackEntry cse("Syr2k");
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
          alpha, A.LockedBuffer(), A.LDim(),
                 B.LockedBuffer(), B.LDim(),
          beta,  C.Buffer(),       C.LDim() );
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
  T alpha, const Matrix<T>& A, const Matrix<T>& B, Matrix<T>& C,
  bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("Syr2k"))
    const Int n = ( orientation==NORMAL ? A.Height() : A.Width() );
    Zeros( C, n, n );
    Syr2k( uplo, orientation, alpha, A, B, T(0), C, conjugate );
}

template<typename T>
void Syr2k
( UpperOrLower uplo, Orientation orientation,
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("Syr2k"))
    if( uplo == LOWER && orientation == NORMAL )
        syr2k::LN( alpha, A, B, beta, C, conjugate );
    else if( uplo == LOWER )
        syr2k::LT( alpha, A, B, beta, C, conjugate );
    else if( orientation == NORMAL )
        syr2k::UN( alpha, A, B, beta, C, conjugate );
    else
        syr2k::UT( alpha, A, B, beta, C, conjugate );
}

template<typename T>
void Syr2k
( UpperOrLower uplo, Orientation orientation,
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
                 DistMatrix<T>& C, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("Syr2k"))
    const Int n = ( orientation==NORMAL ? A.Height() : A.Width() );
    Zeros( C, n, n );
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
    T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B, \
    T beta,        DistMatrix<T>& C, bool conjugate ); \
  template void Syr2k \
  ( UpperOrLower uplo, Orientation orientation, \
    T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B, \
                   DistMatrix<T>& C, bool conjugate );

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
