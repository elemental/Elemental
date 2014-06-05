/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

#include "./Syrk/LN.hpp"
#include "./Syrk/LT.hpp"
#include "./Syrk/UN.hpp"
#include "./Syrk/UT.hpp"

namespace El {

template<typename T>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const Matrix<T>& A, T beta, Matrix<T>& C, bool conjugate )
{
    DEBUG_ONLY(
        CallStackEntry cse("Syrk");
        if( orientation == NORMAL )
        {
            if( A.Height() != C.Height() || A.Height() != C.Width() )
                LogicError("Nonconformal Syrk");
        }
        else
        {
            if( A.Width() != C.Height() || A.Width() != C.Width() )
                LogicError("Nonconformal Syrk");
        }
    )
    const char uploChar = UpperOrLowerToChar( uplo );
    const char transChar = OrientationToChar( orientation );
    const Int k = ( orientation == NORMAL ? A.Width() : A.Height() );
    if( conjugate )
    {
        blas::Herk
        ( uploChar, transChar, C.Height(), k,
          alpha, A.LockedBuffer(), A.LDim(),
          beta,  C.Buffer(),       C.LDim() );
    }
    else
    {
        blas::Syrk
        ( uploChar, transChar, C.Height(), k,
          alpha, A.LockedBuffer(), A.LDim(),
          beta,  C.Buffer(),       C.LDim() );
    }
}

template<typename T>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const Matrix<T>& A, Matrix<T>& C, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("Syrk"))
    const Int n = ( orientation==NORMAL ? A.Height() : A.Width() );
    Zeros( C, n, n );
    Syrk( uplo, orientation, alpha, A, T(0), C, conjugate );
}

template<typename T>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const DistMatrix<T>& A, T beta, DistMatrix<T>& C, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("Syrk"))
    if( uplo == LOWER && orientation == NORMAL )
        syrk::LN( alpha, A, beta, C, conjugate );
    else if( uplo == LOWER )
        syrk::LT( alpha, A, beta, C, conjugate );
    else if( orientation == NORMAL )
        syrk::UN( alpha, A, beta, C, conjugate );
    else
        syrk::UT( alpha, A, beta, C, conjugate );
}

template<typename T>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const DistMatrix<T>& A, DistMatrix<T>& C, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("Syrk"))
    const Int n = ( orientation==NORMAL ? A.Height() : A.Width() );
    Zeros( C, n, n );
    Syrk( uplo, orientation, alpha, A, T(0), C, conjugate );
}

#define PROTO(T) \
  template void Syrk \
  ( UpperOrLower uplo, Orientation orientation, \
    T alpha, const Matrix<T>& A, T beta, Matrix<T>& C, bool conjugate ); \
  template void Syrk \
  ( UpperOrLower uplo, Orientation orientation, \
    T alpha, const Matrix<T>& A, Matrix<T>& C, bool conjugate ); \
  template void Syrk \
  ( UpperOrLower uplo, Orientation orientation, \
    T alpha, const DistMatrix<T>& A, T beta, DistMatrix<T>& C, \
    bool conjugate ); \
  template void Syrk \
  ( UpperOrLower uplo, Orientation orientation, \
    T alpha, const DistMatrix<T>& A, DistMatrix<T>& C, bool conjugate );

// blas::Syrk is not yet supported for Int
//PROTO(Int)
PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
